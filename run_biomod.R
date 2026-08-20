# run_biomod.R  (shiny pipeline - PNG ONLY, direct draw, AUCroc fix + calibration-safe)
# UPDATED:
# - Adaptive modeling options based on predictor count (safe fallback to defaults)
# - Ensemble includes ONLY models passing BOTH AUCroc and TSS thresholds AND not "near-perfect"
# - BOYCE added to both BIOMOD_Modeling calls, ensemble evaluation and summary table
# - Ensemble algorithms: EMmean, EMca (committee averaging), EMwmean, EMcv (uncertainty)
# - Ensemble projections are now EXPORTED as GeoTIFF (continuous + MaxTSS binary)
# - MaxTSS cutoffs computed ONCE from current-day evaluation, held fixed across all scenarios
suppressPackageStartupMessages({
  library(jsonlite)
  library(terra)
  library(biomod2)
  library(dplyr)
})

# --- biomod2 compatibility: old getModelsBuilt() -> new get_built_models() ---
if (!exists("getModelsBuilt", mode = "function") && exists("get_built_models", mode = "function")) {
  getModelsBuilt <- get_built_models
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_biomod.R <config.json>")

# ---- keep nested lists nested (prevents atomic-vector future bug) ----
cfg <- jsonlite::read_json(args[1], simplifyVector = FALSE)

# ---- helpers ----
`%||%` <- function(a, b) {
  ok <- !is.null(a) && length(a) > 0 && !all(is.na(a))
  if (isTRUE(ok)) a else b
}

normalize_models <- function(x) {
  if (is.null(x)) return(character(0))
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x <- x[nzchar(x)]
  unique(x)
}

# ---- Normalize config types (robust for future blocks) ----
if (!is.null(cfg$future_n)) {
  if (is.character(cfg$future_n) && toupper(cfg$future_n) == "NA") cfg$future_n <- NULL
  if (is.character(cfg$future_n) && grepl("^[0-9]+$", cfg$future_n)) cfg$future_n <- as.integer(cfg$future_n)
  if (is.numeric(cfg$future_n)) cfg$future_n <- as.integer(cfg$future_n)
}

if (!is.null(cfg$future)) {
  if (!is.list(cfg$future)) cfg$future <- list(cfg$future)
  for (k in seq_along(cfg$future)) {
    fb <- cfg$future[[k]]
    if (!is.list(fb)) next
    if (!is.null(fb$scenarios)) {
      if (!is.list(fb$scenarios)) fb$scenarios <- list(fb$scenarios)
      cfg$future[[k]]$scenarios <- fb$scenarios
    }
  }
}

# ---- Normalize models (CRITICAL for biomod2) ----
cfg$models <- normalize_models(cfg$models)
if (length(cfg$models) == 0) stop("cfg$models is empty (no models selected).")

# ---- Normalize var_types (app sends named list: var -> type) ----
normalize_var_types <- function(cfg) {
  vt <- cfg$var_types
  if (is.null(vt)) return(NULL)
  if (!is.list(vt)) return(NULL)
  if (is.null(names(vt)) || any(!nzchar(names(vt)))) return(NULL)

  for (nm in names(vt)) {
    vv <- vt[[nm]]
    if (is.null(vv)) {
      vt[[nm]] <- "continuous"
    } else {
      vt[[nm]] <- tolower(as.character(vv)[1])
    }
  }
  vt
}
cfg$var_types <- normalize_var_types(cfg)


# ---------------------------------------------------------------------------
# KULLANICI CSV'SI: kodlama korumali okuma
# Windows'ta CSV'ler cogunlukla UTF-8 degil, Windows-1254 (Turkce) kaydedilir;
# bu durumda tur adlari bozulur ve tur suzgeci hicbir satir bulamaz.
# ---------------------------------------------------------------------------
read_user_csv <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("CSV not found: ", path)
  }
  for (enc in c("UTF-8", "", "windows-1254", "latin1")) {
    d <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = enc),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.data.frame(d) && nrow(d) > 0) {
      if (nzchar(enc) && !identical(enc, "UTF-8")) {
        message("CSV read using encoding: ", enc)
      }
      return(d)
    }
  }
  stop("Could not read the CSV (tried UTF-8, native, windows-1254, latin1): ", path)
}

# -------------------- STATUS --------------------
status_path <- file.path("runs", cfg$run_id, "status.txt")
dir.create(dirname(status_path), recursive = TRUE, showWarnings = FALSE)


# -------------------- ILERLEME GUNLUGU --------------------
# run.log kabuk yonlendirmesiyle (> dosya 2>&1) yazilir. Windows'ta bu yol
# iki sorun cikarir: (1) dosya is bitene kadar baska bir surec tarafindan
# okunamayabilir, (2) cikti tamponlandigi icin dosya uzun sure BOS gorunur.
# Bu yuzden kendi gunlugumuzu HER SATIRDA ACIP-KAPATARAK yaziyoruz; boylece
# arayuz ilerlemeyi aninda gorebilir. run.log yine tutulur (cokme teshisi icin).
log_path <- file.path("runs", cfg$run_id, "progress.log")
log_line <- function(...) {
  txt <- paste0(unlist(list(...)), collapse = "")
  try(cat(format(Sys.time(), "%H:%M:%S"), "  ", txt, "\n",
          sep = "", file = log_path, append = TRUE), silent = TRUE)
}
# message()/warning() ciktisi da gunluge dussun
message <- function(...) { log_line(...); base::message(...) }
warning <- function(...) {
  a <- list(...)
  if (length(a) > 0 && is.character(a[[1]])) log_line("WARNING: ", ...)
  do.call(base::warning, a)
}
log_line("=== run started: ", basename(commandArgs(trailingOnly = TRUE)[1]), " ===")


# -------------------- ORTAM KAYDI (yeniden uretilebilirlik) --------------------
# Hangi R ve paket surumleriyle uretildigini kaydeder. ODMAP ve cogu Q1 dergi
# yeniden uretilebilirlik icin bu bilgiyi ister.
write_session_info <- function(dir) {
  tryCatch({
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    si <- utils::sessionInfo()
    writeLines(utils::capture.output(print(si)), file.path(dir, "sessionInfo.txt"))

    ip <- utils::installed.packages()
    keep <- intersect(rownames(ip),
      c("biomod2","terra","dplyr","jsonlite","ggplot2","scales","R.utils",
        "randomForest","gbm","earth","mgcv","maxnet","xgboost","rpart",
        "nnet","mda","dismo","cito","shiny","shinyFiles"))
    if (length(keep) > 0) {
      utils::write.csv(
        data.frame(package = keep,
                   version = unname(ip[keep, "Version"]),
                   stringsAsFactors = FALSE),
        file.path(dir, "package_versions.csv"), row.names = FALSE)
    }
    log_line("Session info written (R ", getRversion(), ")")
  }, error = function(e) NULL)
}

set_status <- function(x) {
  try(writeLines(x, status_path), silent = TRUE)
  log_line("[STATUS] ", x)
}

# Ensure out_dir exists
if (is.null(cfg$out_dir) || !nzchar(cfg$out_dir)) cfg$out_dir <- file.path("runs", cfg$run_id, "out")
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
write_session_info(cfg$out_dir)

# --- Ensure R.utils exists (keep your behavior) ---
if (!requireNamespace("R.utils", quietly = TRUE)) {
  set_status("INSTALLING PACKAGE: R.utils")
  install.packages("R.utils", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(R.utils))

# -------------------- PAYLASILAN SECENEK MOTORU --------------------
# Literatur temelli parametre izgarasi run_biomod_cont.R ile ORTAKTIR.
opts_file <- file.path(getwd(), "biomod_opts.R")
if (!file.exists(opts_file)) {
  stop("biomod_opts.R bulunamadi. app.R, run_biomod.R ve biomod_opts.R ayni klasorde olmalidir.")
}
source(opts_file, encoding = "UTF-8")

# -------------------- CONFIG: Ensemble filtering thresholds --------------------
# - AUC veya TSS dusukse => ensemble disi
# - AUC veya TSS ~1 ise  => asiri uyum suphesi, ensemble disi
# Esikler artik arayuzden (config.json) okunur; yoksa varsayilan kullanilir.
as_num1 <- function(x, default) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE))[1])
  if (is.na(v)) default else v
}
ENS_MIN_AUC  <- as_num1(cfg$ens_min_auc,  0.70)
ENS_MIN_TSS  <- as_num1(cfg$ens_min_tss,  0.40)
ENS_NEAR_ONE <- as_num1(cfg$ens_near_one, 0.999)
message("Ensemble esikleri: AUCroc >= ", ENS_MIN_AUC,
        " | TSS >= ", ENS_MIN_TSS, " | near-one < ", ENS_NEAR_ONE)

# -------------------- HELPERS --------------------
# safe_png -> biomod_opts.R (paylasilan, platformdan bagimsiz cihaz zinciri)

# Apply variable types (continuous/categorical/nominal) to a SpatRaster
# categorical + nominal => factor
apply_var_types <- function(e, cfg) {
  vt <- cfg$var_types

  # Backward compatibility: if vt missing, fall back to cfg$cat_vars
  if (is.null(vt)) {
    cv <- cfg$cat_vars
    if (!is.null(cv)) cv <- cv[nzchar(cv)]
    if (!is.null(cv) && length(cv) > 0) {
      message("ENV LAYERS: ", paste(names(e), collapse = ", "))
      for (v in cv) {
        if (v %in% names(e)) {
          tv <- Sys.time()
          message("Converting to categorical (reads all cells): ", v, " ...")
          e[[v]] <- terra::as.factor(e[[v]])
          message("  done in ", round(as.numeric(difftime(Sys.time(), tv, units = "secs")), 1), " s")
        } else {
          warning("Categorical var '", v, "' not found among env layers. Skipping.")
        }
      }
    }
    return(e)
  }

  # var_types is named by env names AFTER make.names()
  for (v in intersect(names(e), names(vt))) {
    tp <- vt[[v]] %||% "continuous"
    tp <- tolower(as.character(tp)[1])
    if (tp %in% c("categorical", "nominal")) {
      tv <- Sys.time()
      message("Converting to categorical (reads all cells): ", v, " ...")
      e[[v]] <- terra::as.factor(e[[v]])
      message("  done in ", round(as.numeric(difftime(Sys.time(), tv, units = "secs")), 1), " s")
    }
  }
  e
}

# helper: pick sensible vars for response curves
pick_response_vars <- function(env, cfg, max_n = 50L, drop_categorical = TRUE) {
  vars <- names(env)
  if (length(vars) == 0) return(character(0))

  # allow user override via cfg$response_vars (optional)
  if (!is.null(cfg$response_vars) && length(cfg$response_vars) > 0) {
    rv <- unlist(cfg$response_vars, recursive = TRUE, use.names = FALSE)
    rv <- as.character(rv)
    rv <- rv[nzchar(rv)]
    rv <- rv[rv %in% vars]
    if (length(rv) > 0) vars <- rv
  }

  if (drop_categorical) {
    if (!is.null(cfg$var_types) && is.list(cfg$var_types) && !is.null(names(cfg$var_types))) {
      vt <- cfg$var_types
      tp <- tolower(vapply(vt, function(x) as.character(x)[1], ""))
      drop <- names(vt)[tp %in% c("categorical", "nominal")]
      vars <- setdiff(vars, drop)
    } else if (!is.null(cfg$cat_vars) && length(cfg$cat_vars) > 0) {
      cv <- as.character(unlist(cfg$cat_vars, recursive = TRUE, use.names = FALSE))
      cv <- cv[nzchar(cv)]
      vars <- setdiff(vars, cv)
    }
  }

  vars[seq_len(min(as.integer(max_n), length(vars)))]
}

# -------------------- RESPONSE CURVE DATA EXPORT (SAFE) --------------------
safe_write_csv <- function(df, path) {
  ok <- FALSE
  tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, path, row.names = FALSE)
    ok <- TRUE
  }, error = function(e) ok <<- FALSE)
  ok
}

export_response_curves_data <- function(bm_obj, obj_tag, env, vars, out_csv,
                                       fixed = "median", n_points = 200L) {
  rc_try <- try(
    bm_PlotResponseCurves(
      bm.out = bm_obj,
      models.chosen = "all",
      show.variables = vars,
      fixed.var = fixed,
      do.bivariate = FALSE,
      do.plot = FALSE
    ),
    silent = TRUE
  )

  if (!inherits(rc_try, "try-error") && !is.null(rc_try)) {
    if (is.data.frame(rc_try)) {
      rc_try$source <- "bm_PlotResponseCurves"
      rc_try$obj <- obj_tag
      return(safe_write_csv(rc_try, out_csv))
    }
    if (is.list(rc_try)) {
      dfs <- list()
      for (nm in names(rc_try)) {
        x <- rc_try[[nm]]
        if (is.data.frame(x)) {
          x$curve <- nm
          dfs[[length(dfs) + 1]] <- x
        }
      }
      if (length(dfs) > 0) {
        out <- do.call(rbind, dfs)
        out$source <- "bm_PlotResponseCurves"
        out$obj <- obj_tag
        return(safe_write_csv(out, out_csv))
      }
    }
  }

  fixed_fun <- switch(
    fixed,
    "median" = function(x) terra::global(x, "median", na.rm = TRUE)[1, ],
    "mean"   = function(x) terra::global(x, "mean",   na.rm = TRUE)[1, ],
    "min"    = function(x) terra::global(x, "min",    na.rm = TRUE)[1, ],
    "max"    = function(x) terra::global(x, "max",    na.rm = TRUE)[1, ],
    function(x) terra::global(x, "median", na.rm = TRUE)[1, ]
  )

  base_vals <- try(fixed_fun(env), silent = TRUE)
  if (inherits(base_vals, "try-error") || is.null(base_vals)) return(FALSE)

  base_df <- as.data.frame(base_vals)
  colnames(base_df) <- names(env)

  out_list <- list()

  for (v in vars) {
    if (!v %in% names(env)) next

    is_cat <- terra::is.factor(env[[v]])

    if (is_cat) {
      levs <- try(terra::levels(env[[v]])[[1]], silent = TRUE)
      if (!inherits(levs, "try-error") && !is.null(levs) && ncol(levs) >= 1) {
        ids <- levs[[1]]
        ids <- ids[!is.na(ids)]
        x_vals <- ids
      } else {
        x_vals <- unique(terra::values(env[[v]], mat = FALSE))
        x_vals <- x_vals[!is.na(x_vals)]
        if (length(x_vals) > 50) x_vals <- x_vals[seq_len(50)]
      }
      if (length(x_vals) < 2) next

      newdf <- base_df[rep(1, length(x_vals)), , drop = FALSE]
      newdf[[v]] <- x_vals
    } else {
      mm <- terra::minmax(env[[v]])
      x_min <- as.numeric(mm[1])
      x_max <- as.numeric(mm[2])
      if (!is.finite(x_min) || !is.finite(x_max) || x_min == x_max) next

      x_vals <- seq(x_min, x_max, length.out = as.integer(n_points))
      newdf <- base_df[rep(1, length(x_vals)), , drop = FALSE]
      newdf[[v]] <- x_vals
    }

    pred <- try(predict(bm_obj, newdata = newdf), silent = TRUE)
    if (inherits(pred, "try-error") || is.null(pred)) next

    if (is.data.frame(pred)) pred_vec <- suppressWarnings(as.numeric(pred[[1]]))
    else if (is.matrix(pred)) pred_vec <- suppressWarnings(as.numeric(pred[, 1]))
    else pred_vec <- suppressWarnings(as.numeric(pred))

    dfv <- data.frame(
      obj = obj_tag,
      variable = v,
      x = newdf[[v]],
      pred = pred_vec,
      source = "manual_predict",
      stringsAsFactors = FALSE
    )
    out_list[[length(out_list) + 1]] <- dfv
  }

  if (length(out_list) == 0) return(FALSE)

  out <- do.call(rbind, out_list)
  safe_write_csv(out, out_csv)
}

# -------------------- ADAPTIF MODELLEME SECENEKLERI --------------------
# Eski surumdeki make_adaptive_options() KALDIRILDI.
# Gerekce: biomod2 4.2+ imzasi bm_ModelingOptions(data.type=, models=, strategy=,
# user.val=, user.base=, bm.format=, calib.lines=) seklindedir; ANN=/GBM= gibi
# argumanlar taninmadigi icin eski cagri her zaman hata verip NULL donuyordu ve
# BIOMOD_Modeling varsayilan OPT.strategy="default" (ham paket varsayilanlari)
# ile calisiyordu. Ayrica models.options= diye bir arguman yoktur.
# Yeni yol: biomod_opts.R icindeki literatur temelli motor + OPT.user.val.

# -------------------- ENSEMBLE FILTER: choose models based on eval --------------------
choose_models_for_ensemble <- function(model_out, min_auc, min_tss, near_one) {

  ev <- try(get_evaluations(model_out), silent = TRUE)
  if (inherits(ev, "try-error") || is.null(ev)) return(character(0))

  df <- as.data.frame.table(ev, responseName = "score", stringsAsFactors = FALSE)

  # metric kolonu: icinde AUCroc/TSS gecen kolon hangisi?
  metric_col <- names(df)[sapply(df, function(x) any(as.character(x) %in% c("AUCroc","TSS")))]
  if (length(metric_col) == 0) return(character(0))
  metric_col <- metric_col[1]

  # model/algo kolonu: GLM/GBM/RF... iceren kolon
  algo_col <- names(df)[sapply(df, function(x) any(grepl("^(ANN|CTA|FDA|GAM|GBM|GLM|MARS|MAXENT|MAXNET|XGBOOST|SRE|RFd|RF)$", as.character(x))))]
  if (length(algo_col) == 0) return(character(0))
  algo_col <- algo_col[1]

  # calibration/validation skorlari bazen ayri kolon olur; bazen "score" tek kolondur.
  # Once validation varsa onu al, yoksa calibration, yoksa score.
  score_cols <- names(df)[sapply(df, is.numeric)]
  vcol <- score_cols[grepl("validation", score_cols, ignore.case = TRUE)]
  ccol <- score_cols[grepl("calibration", score_cols, ignore.case = TRUE)]

  if (length(vcol) > 0) use_score <- vcol[1]
  else if (length(ccol) > 0) use_score <- ccol[1]
  else use_score <- "score"

  sub <- df[df[[metric_col]] %in% c("AUCroc","TSS"), c(algo_col, metric_col, use_score)]
  names(sub) <- c("model","metric","value")
  sub$value <- as.numeric(sub$value)

  agg <- aggregate(value ~ model + metric, data = sub, FUN = mean, na.rm = TRUE)
  wide <- reshape(agg, idvar = "model", timevar = "metric", direction = "wide")

  auc_col <- "value.AUCroc"
  tss_col <- "value.TSS"
  if (!auc_col %in% names(wide)) wide[[auc_col]] <- NA_real_
  if (!tss_col %in% names(wide)) wide[[tss_col]] <- NA_real_

  keep <- wide$model[
    is.finite(wide[[auc_col]]) & is.finite(wide[[tss_col]]) &
      wide[[auc_col]] >= min_auc &
      wide[[tss_col]] >= min_tss &
      wide[[auc_col]] <  near_one &
      wide[[tss_col]] <  near_one
  ]

  as.character(keep)
}


# -------------------- LOAD DATA --------------------
set_status("LOADING DATA: listing raster files")

asc_files <- list.files(cfg$asc_dir, pattern = "\\.(asc|tif|tiff|grd)$",
                        full.names = TRUE, ignore.case = TRUE)
grd_f <- asc_files[grepl("\\.grd$", asc_files, ignore.case = TRUE)]
if (length(grd_f) > 0) {
  gri_f <- sub("\\.grd$", ".gri", grd_f, ignore.case = TRUE)
  asc_files <- unique(c(asc_files, gri_f[file.exists(gri_f)]))
}
stopifnot(length(asc_files) > 0)

set_status(paste0("LOADING DATA: opening ", length(asc_files), " raster file(s)"))
t0 <- Sys.time()
env <- rast(asc_files)
names(env) <- make.names(names(env), unique = TRUE)
message("Rasters opened in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " s | ", terra::nlyr(env), " layers | ",
        terra::ncol(env), " x ", terra::nrow(env), " cells")

# UYARI: terra::as.factor() katmanin TUM degerlerini okur. Buyuk bir raster,
# ozellikle .asc (metin) formatindaysa, bu adim cok uzun surebilir ve bellegi
# doldurabilir. Bu yuzden her kategorik katman ayri ayri bildirilir.
set_status("LOADING DATA: applying variable types")
env <- apply_var_types(env, cfg)

set_status("LOADING DATA: reading occurrence CSV")
occ <- read_user_csv(cfg$occ_file)
if ("Tur" %in% names(occ)) occ <- occ %>% filter(Tur == cfg$sp_name)

# Koordinat sutunlari "x"/"y" olmak ZORUNDA DEGIL. Yaygin adlar taninir
# (X/Y, lon/lat, boylam/enlem, POINT_X/POINT_Y, utm_x/utm_y, ...); arayuzden
# elle secilen ad (cfg$x_col / cfg$y_col) her zaman onceliklidir.
occ_xy <- get_xy_matrix(occ, cfg$x_col, cfg$y_col, what = "Occurrence CSV")
if (all(is.na(occ_xy))) {
  stop("Coordinate columns could not be read as numbers. Check the decimal ",
       "separator: R expects a dot (30.25), not a comma (30,25).")
}
ok_xy <- is.finite(occ_xy[, 1]) & is.finite(occ_xy[, 2])
if (any(!ok_xy)) {
  warning(sum(!ok_xy), " row(s) dropped: coordinates not numeric.")
  occ    <- occ[ok_xy, , drop = FALSE]
  occ_xy <- occ_xy[ok_xy, , drop = FALSE]
}
if (nrow(occ_xy) == 0) stop("No usable coordinates in the occurrence CSV.")

# -------------------- GERCEK YOKLUK (opsiyonel) --------------------
# cfg$pa_col verilmisse CSV'deki 1/0 sutunu GERCEK varlik-yokluk olarak okunur.
# Verilmemisse eski davranis: tum satirlar varlik, yokluklar PA ile uretilir.
pa_col <- as.character(cfg$pa_col %||% "")[1]
use_real_abs <- nzchar(pa_col) && pa_col %in% names(occ)

if (use_real_abs) {
  raw <- occ[[pa_col]]
  num <- suppressWarnings(as.numeric(as.character(raw)))
  resp <- ifelse(is.finite(num), as.integer(num > 0),
                 as.integer(!(trimws(tolower(as.character(raw))) %in%
                                c("", "0", "-", ".", "na", "no", "false", "yok", "absent"))))
  resp[is.na(resp)] <- 0L
  message("Real presence/absence column used: '", pa_col,
          "' (presences=", sum(resp == 1L), ", absences=", sum(resp == 0L), ")")
} else {
  if (nzchar(pa_col)) warning("pa_col not found in CSV, falling back to presence-only: ", pa_col)
  resp <- rep(1, nrow(occ_xy))
}

set_status("LOADING DATA: extracting predictor values at occurrences")
v1 <- terra::extract(env[[1]], occ_xy)
keep <- !is.na(v1[[1]])
occ_xy <- occ_xy[keep, , drop = FALSE]
resp <- resp[keep]
stopifnot(nrow(occ_xy) > 5)

n_pres <- sum(resp == 1L)
n_abs  <- sum(resp == 0L)

# PA gerekli mi?
#  - Gercek yokluk yoksa      -> PA zorunlu (eski davranis)
#  - Gercek yokluk CV katindan azsa -> PA ile takviye edilir; aksi halde
#    10 katli CV'de katlarin cogu yokluksuz kalir ve TSS/AUC hesaplanamaz.
CV_K <- 10L
need_pa <- (!use_real_abs) || (n_abs < CV_K)
if (use_real_abs) {
  message("Absence strategy: ",
          if (need_pa) paste0("real + pseudo-absence (only ", n_abs,
                              " real absences, fewer than CV.k=", CV_K, ")")
          else paste0("real absences only (", n_abs, ")"))
}

set_status(paste0("FORMATTING DATA (n=", nrow(occ_xy),
                  ", presences=", n_pres, ", absences=", n_abs, ")"))

# -------------------- BIOMOD formatting + PA --------------------
form_args <- list(
  resp.var = resp,
  expl.var = env,
  resp.xy = occ_xy,
  resp.name = cfg$sp_name,
  filter.raster = TRUE
)
if (need_pa) {
  form_args$PA.nb.rep     <- cfg$pa_rep
  form_args$PA.nb.absences <- cfg$pa_n
  form_args$PA.strategy   <- cfg$pa_strategy
  if (!is.null(cfg$pa_strategy) && cfg$pa_strategy == "disk") {
    form_args$PA.dist.min <- cfg$pa_dist_min
    form_args$PA.dist.max <- cfg$pa_dist_max
  }
}

biomod_data <- tryCatch(do.call(BIOMOD_FormatingData, form_args), error = function(e) e)
if (inherits(biomod_data, "error") && need_pa && use_real_abs) {
  # Bazi surumler gercek yokluk + PA karisimini reddeder -> yalniz varliklar + PA
  warning("Mixing real absences with PA failed (", conditionMessage(biomod_data),
          "); retrying as presence-only + PA.")
  fa <- form_args
  fa$resp.var <- rep(1, sum(resp == 1L))
  fa$resp.xy  <- occ_xy[resp == 1L, , drop = FALSE]
  biomod_data <- do.call(BIOMOD_FormatingData, fa)
}
if (inherits(biomod_data, "error")) stop(conditionMessage(biomod_data))

# -------------------- MODELING --------------------
set_status("MODELING (10-fold CV)")

models_to_run <- cfg$models
if (is.null(models_to_run) || length(models_to_run) == 0) stop("cfg$models is empty.")

req_pkgs <- list(
  ANN     = "nnet",
  CTA     = "rpart",
  DNN     = "cito",
  FDA     = "mda",
  GBM     = "gbm",
  MARS    = "earth",
  GAM     = "mgcv",
  MAXNET  = "maxnet",
  XGBOOST = "xgboost",
  RF      = "randomForest",
  RFd     = "randomForest"
)

for (m in intersect(names(req_pkgs), models_to_run)) {
  if (!requireNamespace(req_pkgs[[m]], quietly = TRUE)) {
    warning("Package '", req_pkgs[[m]], "' missing -> dropping model ", m)
    models_to_run <- setdiff(models_to_run, m)
  }
}

has_maxent <- any(grepl("^MAXENT$", models_to_run))
if (has_maxent) {
  if (is.null(cfg$maxent_jar) || !file.exists(cfg$maxent_jar)) {
    warning("MAXENT secili ancak maxent.jar bulunamadi -> MAXENT atiliyor.")
    models_to_run <- models_to_run[!grepl("^MAXENT$", models_to_run)]
  } else {
    # biomod2 java dosyasini simulasyon klasorunde (4.2) veya dismo/java (4.3) arar
    bmopt_stage_maxent(cfg$maxent_jar, verbose = TRUE)
  }
}

if (length(models_to_run) == 0) stop("No models left to run.")

cv_rep <- 1L
if (!is.null(cfg$cv$rep) && is.numeric(cfg$cv$rep) && cfg$cv$rep >= 1) {
  cv_rep <- as.integer(cfg$cv$rep)
} else if (!is.null(cfg$cv_rep) && is.numeric(cfg$cv_rep) && cfg$cv_rep >= 1) {
  cv_rep <- as.integer(cfg$cv_rep)
}

nb_cpu <- 1L
if (!is.null(cfg$cpu_n) && is.numeric(cfg$cpu_n) && cfg$cpu_n >= 1) {
  nb_cpu <- as.integer(cfg$cpu_n)
}

# ---- LITERATUR TEMELLI MODELLEME SECENEKLERI ----
p_env  <- terra::nlyr(env)
n_pres <- nrow(occ_xy)

set_status(paste0("BUILDING MODEL OPTIONS (p=", p_env, ", n=", n_pres, ")"))

built_opts <- bmopt_build(
  data_type  = "binary",
  models     = models_to_run,
  bm_format  = biomod_data,
  p          = p_env,
  n          = n_pres,                       # binary'de ETKIN orneklem = varlik sayisi
  maxent_jar = cfg$maxent_jar,
  strategy   = cfg$opt_strategy %||% "adaptive",
  transfer   = TRUE,                         # gelecek projeksiyonu -> daha guclu duzenlileme
  verbose    = TRUE
)

bmopt_report(built_opts, file.path(cfg$out_dir, "modeling_options_applied.csv"),
             data_type = "binary", p = p_env, n = n_pres)

mod_args <- list(
  bm.format   = biomod_data,
  modeling.id = "binary",
  models      = models_to_run,
  CV.strategy = "kfold",
  CV.k        = 10,
  CV.nb.rep   = cv_rep,
  metric.eval = c("AUCroc", "TSS", "BOYCE"),
  var.import  = 3,
  prevalence  = 0.5,        # Barbet-Massin et al. (2012): varlik/yokluk esit agirlik
  seed.val    = 42,
  nb.cpu      = nb_cpu,
  do.progress = TRUE
)
mod_args <- bmopt_inject(mod_args, built_opts)

model_out <- tryCatch(do.call(BIOMOD_Modeling, mod_args), error = function(e) e)

# Kademeli geri cekilme: user.defined -> bigboss -> default
if (inherits(model_out, "error")) {
  message("user.defined secenekleriyle modelleme basarisiz: ", conditionMessage(model_out))
  for (fb in c("bigboss", "default")) {
    message("OPT.strategy='", fb, "' ile tekrar deneniyor.")
    a2 <- mod_args
    a2$OPT.user.val <- NULL; a2$OPT.user.base <- NULL
    a2$OPT.strategy <- fb
    model_out <- tryCatch(do.call(BIOMOD_Modeling, a2), error = function(e) e)
    if (!inherits(model_out, "error")) break
  }
}
if (inherits(model_out, "error")) stop("BIOMOD_Modeling basarisiz: ", conditionMessage(model_out))

# -------------------- SAVE PLOTS (PNG ONLY, DIRECT DRAW) --------------------
set_status("SAVING DEFAULT BIOMOD PLOTS")

plots_dir <- file.path(cfg$out_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

dataset_plot <- "calibration"

ok_mean <- safe_png(
  file = file.path(plots_dir, "biomod_eval_mean.png"),
  expr = bm_PlotEvalMean(
    bm.out = model_out,
    metric.eval = c("AUCroc", "TSS"),
    dataset = dataset_plot,
    group.by = "algo",
    do.plot = TRUE
  )
)

ok_boyce <- safe_png(
  file = file.path(plots_dir, "biomod_eval_boyce.png"),
  expr = bm_PlotEvalMean(
    bm.out = model_out,
    metric.eval = c("BOYCE", "TSS"),
    dataset = dataset_plot,
    group.by = "algo",
    do.plot = TRUE
  )
)

ok_box <- safe_png(
  file = file.path(plots_dir, "biomod_eval_boxplot.png"),
  expr = bm_PlotEvalBoxplot(
    bm.out = model_out,
    dataset = dataset_plot,
    group.by = "algo",
    do.plot = TRUE
  )
)

ok_def <- safe_png(
  file = file.path(plots_dir, "biomod_default_plot.png"),
  expr = plot(model_out)
)

if (!ok_mean) set_status(paste0("SAVING DEFAULT BIOMOD PLOTS (FAILED: eval_mean; dataset=", dataset_plot, ")"))
if (!ok_box)  set_status(paste0("SAVING DEFAULT BIOMOD PLOTS (FAILED: eval_box; dataset=", dataset_plot, ")"))
if (!ok_def)  set_status("SAVING DEFAULT BIOMOD PLOTS (FAILED: default plot png)")

# -------------------- OPTIONAL: RESPONSE CURVES (single models) --------------------
set_status("RESPONSE CURVES (models)")

resp_vars_plot <- pick_response_vars(env, cfg, max_n = 50L, drop_categorical = FALSE)
resp_vars_csv  <- pick_response_vars(env, cfg, max_n = 50L, drop_categorical = FALSE)

if (length(resp_vars_plot) == 0 && length(resp_vars_csv) > 0) {
  resp_vars_plot <- resp_vars_csv[seq_len(min(50L, length(resp_vars_csv)))]
}

ok_rc_models <- safe_png(
  file = file.path(plots_dir, "models_response_curves.png"),
  expr = bm_PlotResponseCurves(
    bm.out = model_out,
    models.chosen = "all",
    show.variables = resp_vars_plot,
    fixed.var = "median",
    do.bivariate = FALSE,
    do.plot = TRUE
  )
)
if (!ok_rc_models) warning("Failed to save models response curves PNG.")

rc_dir <- file.path(cfg$out_dir, "response_curves")
dir.create(rc_dir, recursive = TRUE, showWarnings = FALSE)

ok_rc_models_csv <- export_response_curves_data(
  bm_obj = model_out,
  obj_tag = "models",
  env = env,
  vars = resp_vars_csv,
  out_csv = file.path(rc_dir, "models_response_curves_data.csv"),
  fixed = "median",
  n_points = 200L
)
if (!ok_rc_models_csv) warning("Failed to export models response-curve data CSV.")

# -------------------- EXPORT METRICS (MODEL) --------------------
set_status("EXPORTING METRICS (model)")

eval_mod <- get_evaluations(model_out)
eval_mod_df <- as.data.frame.table(eval_mod, responseName = "value", stringsAsFactors = FALSE)
write.csv(eval_mod_df, file.path(cfg$out_dir, "eval_model.csv"), row.names = FALSE)

# -------------------- ENSEMBLE (FILTERED) --------------------
set_status("ENSEMBLE")

keep_models <- choose_models_for_ensemble(
  model_out = model_out,
  min_auc   = ENS_MIN_AUC,
  min_tss   = ENS_MIN_TSS,
  near_one  = ENS_NEAR_ONE
)
# CLEAN (defensive)
keep_models <- keep_models[!tolower(keep_models) %in% c("model","none","na","null")]
# If filtering returns nothing, fall back safely (do not break pipeline)

# --- FIX: BIOMOD_EnsembleModeling needs FULL model names, not algo codes ---
built_full <- tryCatch(getModelsBuilt(model_out), error = function(e) character(0))
built_full <- unique(as.character(built_full))
built_full <- built_full[nzchar(built_full)]

keep_vec <- unique(as.character(keep_models))
keep_vec <- keep_vec[nzchar(keep_vec)]

keep_full <- character(0)

if (length(keep_vec) > 0) {
  # if keep_vec already contains RUN pattern, assume it's already full.name
  if (any(grepl("_RUN[0-9]+_", keep_vec))) {
    keep_full <- keep_vec
  } else {
    # keep_vec is algo list (ANN/GLM/...) -> map to full names by suffix
    algo_of_full <- sub(".*_", "", built_full)
    keep_full <- built_full[algo_of_full %in% keep_vec]
  }
}

models_chosen <- if (length(keep_full) > 0) keep_full else "all"


# prefer validation metrics for selection if available, else calibration
tmp_ev <- try(get_evaluations(model_out), silent = TRUE)
use_data <- "calibration"
if (!inherits(tmp_ev, "try-error") && !is.null(tmp_ev)) {
  df0 <- as.data.frame.table(tmp_ev, responseName = "value", stringsAsFactors = FALSE)
  # sutun sayisi biomod2 surumune gore degisebilir -> guvenli atama
  if (ncol(df0) == 6) colnames(df0) <- c("metric_set","metric","data","model","run","value")
  has_val <- any(vapply(df0, function(x) any(as.character(x) == "validation"), logical(1)))
  if (isTRUE(has_val)) use_data <- "validation"
}

ens <- BIOMOD_EnsembleModeling(
  bm.mod = model_out,
  models.chosen = models_chosen,
  em.by = "all",
  em.algo = c("EMmean", "EMca", "EMwmean", "EMcv"),
  metric.select = c("AUCroc", "TSS"),
  metric.select.thresh = c(ENS_MIN_AUC, ENS_MIN_TSS),
  metric.select.dataset = use_data,
  metric.eval = c("AUCroc", "TSS", "BOYCE")
)

# -------------------- ENSEMBLE: EXPORT SELECTED MODELS (CRITICAL) --------------------
# What YOU requested as candidate list (after your filtering)
write.csv(
  data.frame(model = if (length(keep_models) > 0) keep_models else character(0)),
  file = file.path(cfg$out_dir, "ensemble_candidates_after_thresholds.csv"),
  row.names = FALSE
)

# What BIOMOD2 ACTUALLY built into the ensemble
# -------------------- ENSEMBLE: DETECT ACTUALLY USED BASE MODELS --------------------
ens_built <- character(0)

# 1) Primary: BIOMOD native
ens_built_try <- try(getModelsBuilt(ens), silent = TRUE)
if (!inherits(ens_built_try, "try-error") &&
    !is.null(ens_built_try) &&
    length(ens_built_try) > 0) {
  ens_built <- as.character(ens_built_try)
}

# 2) Fallback: inspect ensemble predictions
if (length(ens_built) == 0) {

  ens_pred <- try(get_predictions(ens), silent = TRUE)

  if (!inherits(ens_pred, "try-error") && !is.null(ens_pred)) {

    # CASE A: data.frame output
    if (is.data.frame(ens_pred)) {

      if ("full.name" %in% names(ens_pred)) {
        ens_built <- unique(as.character(ens_pred$full.name))

      } else if ("model" %in% names(ens_pred)) {
        ens_built <- unique(as.character(ens_pred$model))
      }

    # CASE B: list output (names often carry model ids)
    } else if (is.list(ens_pred)) {

      nms <- names(ens_pred)
      if (!is.null(nms) && length(nms) > 0) {
        ens_built <- as.character(nms)
      }
    }
  }
}

# 3) Final cleanup: ensure character vector
ens_built <- unique(as.character(ens_built))
ens_built <- ens_built[nzchar(ens_built)]
# remove accidental header-like tokens
ens_built <- ens_built[!tolower(ens_built) %in% c("model", "none", "na", "null")]

# Optional: if only EM* names exist, record that fact explicitly
if (length(ens_built) > 0 && all(grepl("^EM", ens_built))) {
  write.csv(
    data.frame(model = as.character(models_chosen)),
    file = file.path(cfg$out_dir, "ensemble_models_chosen_argument.csv"),
    row.names = FALSE
  )
}

write.csv(
  data.frame(model = if (length(ens_built) > 0) ens_built else "NONE"),
  file = file.path(cfg$out_dir, "ensemble_selected_models.csv"),
  row.names = FALSE
)

# -------------------- VAR IMPORTANCE (ENSEMBLE-AVERAGED) [ROBUST] --------------------
set_status("EXPORTING VAR IMPORTANCE (ensemble-avg)")

# 1) get which base models were built (full.name)
built_full <- tryCatch(getModelsBuilt(model_out), error = function(e) character(0))
built_full <- unique(as.character(built_full))
built_full <- built_full[nzchar(built_full)]

# 2) decide which models are "selected for ensemble"
# keep_models senin choose_models_for_ensemble() ciktin -> genelde algo duzeyinde (ANN/RF/GBM...)
selected_algo <- if (length(keep_models) > 0) keep_models else character(0)
selected_algo <- unique(as.character(selected_algo))
selected_algo <- selected_algo[nzchar(selected_algo)]

# Eger keep_models algo duzeyindeyse, full.name icinden algoya gore sec
selected_full <- built_full
if (length(selected_algo) > 0) {
  selected_full <- built_full[sub(".*_", "", built_full) %in% selected_algo]
}

# fallback: hic eslesme yoksa tum built_full kullan
if (length(selected_full) == 0) selected_full <- built_full

# 3) get varimp table safely via bm_PlotVarImpBoxplot (returns a data.frame)
vip <- try(
  bm_PlotVarImpBoxplot(
    bm.out = model_out,
    group.by = c("full.name", "expl.var", "run"),
    do.plot = FALSE
  ),
  silent = TRUE
)

vip_df <- NULL
if (!inherits(vip, "try-error") && !is.null(vip)) {
  # biomod2 surumlerine gore isimler degisebiliyor: ilk eleman genelde data.frame
  if (is.data.frame(vip)) {
    vip_df <- vip
  } else if (is.list(vip)) {
    vip_df <- vip$tab %||% vip$data %||% vip[[1]]
    if (!is.data.frame(vip_df)) vip_df <- NULL
  }
}

if (is.null(vip_df) || nrow(vip_df) == 0) {
  warning("VarImp: bm_PlotVarImpBoxplot could not return a usable table (vip_df is empty).")
} else {

  # 4) detect columns robustly
  # variable column
  var_col <- if ("expl.var" %in% names(vip_df)) "expl.var" else NULL
  if (is.null(var_col)) {
    cand <- names(vip_df)[vapply(vip_df, function(x) any(as.character(x) %in% names(env)), logical(1))]
    if (length(cand) > 0) var_col <- cand[1]
  }

  # model column (full.name preferred)
  mod_col <- if ("full.name" %in% names(vip_df)) "full.name" else NULL
  if (is.null(mod_col) && "model" %in% names(vip_df)) mod_col <- "model"
  if (is.null(mod_col)) {
    cand <- names(vip_df)[vapply(vip_df, function(x) any(grepl("_(ANN|CTA|FDA|GAM|GBM|GLM|MARS|MAXENT|MAXNET|XGBOOST|SRE|RFd|RF)$", as.character(x))), logical(1))]
    if (length(cand) > 0) mod_col <- cand[1]
  }

  # importance numeric column
  num_cols <- names(vip_df)[vapply(vip_df, is.numeric, logical(1))]
  imp_col <- NULL
  if ("var.imp" %in% names(vip_df)) imp_col <- "var.imp"
  if (is.null(imp_col) && length(num_cols) > 0) {
    # choose the numeric col with mean in [0,1] if possible
    score <- sapply(num_cols, function(cc) {
      x <- vip_df[[cc]]
      x <- x[is.finite(x)]
      if (length(x) == 0) return(Inf)
      m <- mean(x, na.rm = TRUE)
      if (m >= 0 && m <= 1) abs(m - 0.5) else Inf
    })
    best <- names(sort(score))[1]
    if (!is.null(best) && is.finite(score[[best]])) imp_col <- best else imp_col <- num_cols[1]
  }

  if (is.null(var_col) || is.null(mod_col) || is.null(imp_col)) {
    warning("VarImp: could not detect required columns in vip_df (expl.var/full.name/var.imp).")
  } else {

    # 5) filter only selected models
    vip_sub <- vip_df
    vip_sub$..var <- as.character(vip_sub[[var_col]])
    vip_sub$..mod <- as.character(vip_sub[[mod_col]])
    vip_sub$..imp <- suppressWarnings(as.numeric(vip_sub[[imp_col]]))
    vip_sub <- vip_sub[is.finite(vip_sub$..imp), , drop = FALSE]

    # Keep selected full names; if vip uses algo instead, also allow algo match
    sel_algo2 <- unique(sub(".*_", "", selected_full))
    vip_keep <- (vip_sub$..mod %in% selected_full) | (vip_sub$..mod %in% sel_algo2)

    vip_sub <- vip_sub[vip_keep, , drop = FALSE]

    if (nrow(vip_sub) == 0) {
      warning("VarImp: after filtering, no rows matched selected ensemble models.")
    } else {

      # 6) ensemble-avg importance by variable
      ens_vi <- aggregate(..imp ~ ..var, data = vip_sub, FUN = mean, na.rm = TRUE)
      names(ens_vi) <- c("variable", "importance_mean")

      s <- sum(ens_vi$importance_mean, na.rm = TRUE)
      ens_vi$importance_percent <- if (is.finite(s) && s > 0) 100 * ens_vi$importance_mean / s else NA_real_

      ens_vi <- ens_vi[order(ens_vi$importance_mean, decreasing = TRUE), , drop = FALSE]

      write.csv(ens_vi, file.path(cfg$out_dir, "var_importance_ensemble_avg.csv"), row.names = FALSE)

      safe_png(
        file = file.path(plots_dir, "var_importance_ensemble_avg.png"),
        expr = {
          op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
          par(mar = c(6, 10, 3, 1))
          x <- ens_vi$importance_mean
          names(x) <- ens_vi$variable
          barplot(rev(x), horiz = TRUE, las = 1,
                  xlab = "Mean variable importance (selected models)")
          title("Ensemble-averaged variable importance")
        },
        width = 1600, height = 1200, res = 150
      )

      # debug raw table (optional but useful)
      write.csv(vip_df, file.path(cfg$out_dir, "var_importance_bmPlot_table_raw.csv"), row.names = FALSE)
      write.csv(vip_sub, file.path(cfg$out_dir, "var_importance_bmPlot_table_filtered.csv"), row.names = FALSE)
    }
  }
}


# -------------------- RESPONSE CURVES (ensemble) --------------------
set_status("RESPONSE CURVES (ensemble)")

ok_rc_ens <- safe_png(
  file = file.path(plots_dir, "ensemble_response_curves.png"),
  expr = bm_PlotResponseCurves(
    bm.out = ens,
    models.chosen = "all",
    show.variables = resp_vars_plot,
    fixed.var = "median",
    do.bivariate = FALSE,
    do.plot = TRUE
  )
)
if (!ok_rc_ens) warning("Failed to save ensemble response curves PNG.")

ok_rc_ens_csv <- export_response_curves_data(
  bm_obj = ens,
  obj_tag = "ensemble",
  env = env,
  vars = resp_vars_csv,
  out_csv = file.path(rc_dir, "ensemble_response_curves_data.csv"),
  fixed = "median",
  n_points = 200L
)
if (!ok_rc_ens_csv) warning("Failed to export ensemble response-curve data CSV.")

set_status("EXPORTING METRICS (ensemble)")

eval_ens <- get_evaluations(ens)
eval_ens_df <- as.data.frame.table(eval_ens, responseName = "value", stringsAsFactors = FALSE)
write.csv(eval_ens_df, file.path(cfg$out_dir, "eval_ensemble.csv"), row.names = FALSE)

# -------------------- EXPORT SUMMARY CSV (MEAN AUC/TSS + ENSEMBLE) --------------------
summarize_eval_long <- function(df_long) {
  if (is.null(df_long) || !is.data.frame(df_long) || nrow(df_long) == 0) return(NULL)

  # ---- helper: pick a numeric column that looks like a score in [0,1] (AUC/TSS) ----
  pick_01_col <- function(d) {
    num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    if (length(num_cols) == 0) return(NULL)

    # score each numeric col: prefer columns whose mean is within [0,1]
    score <- sapply(num_cols, function(cc) {
      x <- suppressWarnings(as.numeric(d[[cc]]))
      x <- x[is.finite(x)]
      if (length(x) == 0) return(Inf)

      m <- mean(x, na.rm = TRUE)
      # in [0,1] is what we want; closer to ~0.75 is typical for SDM metrics
      if (m >= 0 && m <= 1) return(abs(m - 0.75))
      Inf
    })

    best <- names(sort(score))[1]
    if (!is.null(best) && is.finite(score[[best]])) best else NULL
  }

  # --- pick value column (prefer 0..1 metric columns) ---
  val_col <- pick_01_col(df_long)

  # fallback: if no numeric column in [0,1], try classic names
  if (is.null(val_col)) {
    if ("value" %in% names(df_long)) val_col <- "value"
    if (is.null(val_col) && "Freq" %in% names(df_long)) val_col <- "Freq"

    # last fallback: first numeric column
    if (is.null(val_col)) {
      num_cols <- names(df_long)[vapply(df_long, is.numeric, logical(1))]
      if (length(num_cols) > 0) val_col <- num_cols[1]
    }
  }

  if (is.null(val_col)) return(NULL)

  # store numeric values in a consistent internal column
  val_col <- pick_01_col(df_long)
  if (is.null(val_col)) return(NULL)
  df_long$..value <- suppressWarnings(as.numeric(df_long[[val_col]]))



  # --- find metric, data, model columns robustly ---
  metric_col <- NULL
  if ("metric" %in% names(df_long)) metric_col <- "metric"
  if (is.null(metric_col) && "metric.eval" %in% names(df_long)) metric_col <- "metric.eval"
  if (is.null(metric_col) && "full.name" %in% names(df_long)) metric_col <- "full.name"
  if (is.null(metric_col)) {
    # try to detect: column containing AUCroc/TSS/ROC
    cand <- names(df_long)[vapply(df_long, function(x) any(as.character(x) %in% c("AUCroc","TSS","ROC","AUC")), logical(1))]
    if (length(cand) > 0) metric_col <- cand[1]
  }
  if (is.null(metric_col)) return(NULL)

  model_col <- NULL
  if ("model" %in% names(df_long)) model_col <- "model"
  if (is.null(model_col) && "algo" %in% names(df_long)) model_col <- "algo"
  if (is.null(model_col) && "model.name" %in% names(df_long)) model_col <- "model.name"
  if (is.null(model_col) && "full.name" %in% names(df_long)) model_col <- "full.name"
  if (is.null(model_col)) {
    # try: column containing ANN/CTA/GBM...
    cand <- names(df_long)[vapply(df_long, function(x) any(grepl("^(ANN|CTA|FDA|GAM|GBM|GLM|MARS|MAXENT|RFd|RF)$", as.character(x))), logical(1))]
    if (length(cand) > 0) model_col <- cand[1]
  }
  if (is.null(model_col)) return(NULL)

  data_col <- NULL
  if ("data" %in% names(df_long)) data_col <- "data"
  if (is.null(data_col) && "dataset" %in% names(df_long)) data_col <- "dataset"

  df_long$..metric <- as.character(df_long[[metric_col]])
  df_long$..model  <- as.character(df_long[[model_col]])

  # --- choose dataset: prefer validation if exists ---
  if (!is.null(data_col)) {
    df_long$..data <- as.character(df_long[[data_col]])
    use_data <- if ("validation" %in% unique(df_long$..data)) "validation" else "calibration"
    sub <- df_long[df_long$..data == use_data & df_long$..metric %in% c("AUCroc","ROC","AUC","TSS","BOYCE"), , drop = FALSE]
  } else {
    # no data column -> just use all
    use_data <- NA_character_
    sub <- df_long[df_long$..metric %in% c("AUCroc","ROC","AUC","TSS","BOYCE"), , drop = FALSE]
  }

  if (nrow(sub) == 0) return(NULL)

  agg <- aggregate(..value ~ ..model + ..metric, data = sub, FUN = mean, na.rm = TRUE)
  wide <- reshape(agg, idvar = "..model", timevar = "..metric", direction = "wide")

  # normalize names
  names(wide) <- gsub("^\\.\\.value\\.", "", names(wide))
  names(wide)[names(wide) == "..model"] <- "model"

  # unify AUC column name to AUCroc_mean later
  if ("ROC" %in% names(wide) && !"AUCroc" %in% names(wide)) wide$AUCroc <- wide$ROC
  if ("AUC" %in% names(wide) && !"AUCroc" %in% names(wide)) wide$AUCroc <- wide$AUC

  wide$dataset_used <- if (is.na(use_data)) "unknown" else use_data

  # rename to *_mean like your downstream expects
  if ("AUCroc" %in% names(wide)) names(wide)[names(wide) == "AUCroc"] <- "AUCroc_mean"
  if ("TSS" %in% names(wide))    names(wide)[names(wide) == "TSS"]    <- "TSS_mean"
  if ("BOYCE" %in% names(wide))  names(wide)[names(wide) == "BOYCE"]  <- "BOYCE_mean"

  wide
}


# Base model summary (from eval_model.csv in memory: eval_mod_df)
base_summary <- summarize_eval_long(eval_mod_df)
if (!is.null(base_summary)) {
  base_summary$type <- "base_model"

  # Mark whether each base model is in ensemble (needs ensemble_selected_models.csv from step #1)
  # keep_models icinde full.name gelirse algo'yu yakala (son parca genelde algo)
  # in_ensemble = threshold + near-one kuralina gore (AUC AND TSS)
  base_summary$in_ensemble <- with(base_summary,
  is.finite(AUCroc_mean) & is.finite(TSS_mean) &
  AUCroc_mean >= ENS_MIN_AUC &
  TSS_mean    >= ENS_MIN_TSS &
  AUCroc_mean <  ENS_NEAR_ONE &
  TSS_mean    <  ENS_NEAR_ONE
)

}

# Ensemble summary (from eval_ensemble.csv in memory: eval_ens_df)
ens_summary <- summarize_eval_long(eval_ens_df)
if (!is.null(ens_summary)) {
  ens_summary$type <- "ensemble"
  ens_summary$in_ensemble <- NA
}

# Bind + write
sum_out <- NULL

if (!is.null(base_summary) && !is.null(ens_summary)) {
  sum_out <- rbind(base_summary, ens_summary)
} else if (!is.null(base_summary)) {
  sum_out <- base_summary
} else if (!is.null(ens_summary)) {
  sum_out <- ens_summary
}

if (!is.null(sum_out)) {
  keep_cols <- intersect(c("type","model","AUCroc_mean","TSS_mean","BOYCE_mean","dataset_used","in_ensemble"), names(sum_out))
  sum_out <- sum_out[, keep_cols, drop = FALSE]
  sum_out <- sum_out[order(sum_out$type, sum_out$model), , drop = FALSE]

  write.csv(sum_out, file.path(cfg$out_dir, "auc_tss_summary.csv"), row.names = FALSE)
  cat("\nWrote summary: ", file.path(cfg$out_dir, "auc_tss_summary.csv"), "\n")
}


# -------------------- PROJECTION (current + optional future) --------------------
set_status("PROJECTION")

rasters_dir <- file.path(cfg$out_dir, "rasters")
dir.create(rasters_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------- CIKTI OLCEGI --------------------
# biomod2 0-1 olasiliklari 0-1000 TAM SAYI olcegine cevirerek saklar
# (on_0_1000 = TRUE varsayilani; 2 bayt/piksel -> ciddi disk tasarrufu).
# Dokumantasyonun kendi onerisi: 0-1'e donmek icin projeksiyonlari 1000'e bol.
#
# NOT: biomod2'ye on_0_1000 = FALSE GECIRMIYORUZ. Gerekce: bu yolda .grd/terra
# okumasinda bozulma bildirilmis (biomod2 issue #72). Bunun yerine biomod2 kendi
# denenmis tam sayi yolunda calisir, bolme YALNIZCA disa aktarimda yapilir.
# Ikili (MaxTSS) haritalar HER ZAMAN ham 0-1000 olcegi ve ham cutoff ile
# hesaplanir; boylece ondalik karsilastirma hatasi olusamaz.
out_scale <- tolower(as.character(cfg$out_scale %||% "0_1")[1])
SCALE_DIV <- if (identical(out_scale, "0_1")) 1000 else 1
message("Raster cikti olcegi: ", if (SCALE_DIV > 1) "0-1 (1000'e bolunerek)" else "0-1000 (biomod2 ham)")

# -------------------- MaxTSS CUTOFFS (ensemble) --------------------
# Esik DEGERI BIR KEZ gunumuz verisi uzerinden hesaplanir ve tum
# gelecek senaryolarinda DEGISTIRILMEDEN kullanilir.
get_em_cutoffs <- function(ens_obj, metric = "TSS") {
  ev <- try(get_evaluations(ens_obj), silent = TRUE)
  if (inherits(ev, "try-error") || is.null(ev)) return(NULL)

  df <- if (is.data.frame(ev)) ev else try(as.data.frame(ev), silent = TRUE)
  if (inherits(df, "try-error") || !is.data.frame(df)) return(NULL)

  nm_col <- if ("full.name" %in% names(df)) "full.name" else
            if ("model.name" %in% names(df)) "model.name" else
            if ("model" %in% names(df)) "model" else NULL
  mt_col <- if ("metric.eval" %in% names(df)) "metric.eval" else
            if ("metric" %in% names(df)) "metric" else NULL
  ct_col <- if ("cutoff" %in% names(df)) "cutoff" else NULL

  if (is.null(nm_col) || is.null(mt_col) || is.null(ct_col)) {
    warning("Cutoff tablosu okunamadi (sutun adlari eslesmedi).")
    return(NULL)
  }

  sub <- df[as.character(df[[mt_col]]) == metric, c(nm_col, ct_col), drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  names(sub) <- c("model", "cutoff")
  sub$model  <- as.character(sub$model)
  sub$cutoff <- suppressWarnings(as.numeric(sub$cutoff))
  sub <- sub[is.finite(sub$cutoff), , drop = FALSE]
  sub[!duplicated(sub$model), , drop = FALSE]
}

# -------------------- RASTER EXPORT (surekli + MaxTSS ikili) --------------------
# threshold_pattern: yalnizca EMwmean/EMmean esiklenir.
# EMca bir oy oranidir, EMcv belirsizlik katsayisidir -> esiklenmez.
write_proj_rasters <- function(proj_obj, tag, rasters_dir, cutoffs = NULL,
                               threshold_pattern = "EMwmean|EMmean",
                               scale_div = 1) {
  pred <- try(get_predictions(proj_obj), silent = TRUE)
  if (inherits(pred, "try-error") || is.null(pred)) {
    warning("get_predictions basarisiz: ", tag); return(invisible(NULL))
  }
  if (!inherits(pred, "SpatRaster")) {
    warning("Tahminler SpatRaster degil (atlaniyor): ", tag); return(invisible(NULL))
  }

  n_cont <- 0L; n_bin <- 0L

  for (i in seq_len(terra::nlyr(pred))) {
    lyr_name  <- names(pred)[i]
    safe_name <- make.names(lyr_name)

    ok_c <- try({
      r_out <- if (scale_div > 1) pred[[i]] / scale_div else pred[[i]]
      terra::writeRaster(
        r_out,
        file.path(rasters_dir, paste0(tag, "_", safe_name, "_continuous.tif")),
        overwrite = TRUE
      ); TRUE
    }, silent = TRUE)
    if (!inherits(ok_c, "try-error")) n_cont <- n_cont + 1L

    if (is.null(cutoffs) || nrow(cutoffs) == 0) next
    if (!grepl(threshold_pattern, lyr_name)) next

    idx <- which(vapply(cutoffs$model,
                        function(m) grepl(m, lyr_name, fixed = TRUE),
                        logical(1)))
    if (length(idx) == 0) {
      algo_lyr <- sub(".*(EM[a-zA-Z]+).*", "\\1", lyr_name)
      idx <- which(vapply(cutoffs$model,
                          function(m) sub(".*(EM[a-zA-Z]+).*", "\\1", m) == algo_lyr,
                          logical(1)))
    }
    if (length(idx) == 0) { warning("Cutoff eslesmedi: ", lyr_name); next }

    ct  <- cutoffs$cutoff[idx[1]]
    ok_b <- try({
      bin <- pred[[i]] >= ct
      terra::writeRaster(
        bin,
        file.path(rasters_dir, paste0(tag, "_", safe_name, "_binary_maxTSS.tif")),
        overwrite = TRUE, datatype = "INT1U"
      ); TRUE
    }, silent = TRUE)
    if (!inherits(ok_b, "try-error")) n_bin <- n_bin + 1L
  }

  message(tag, ": ", n_cont, " surekli, ", n_bin, " ikili raster yazildi.")
  invisible(NULL)
}

# Esikleri BIR KEZ hesapla (gunumuz degerlendirmesi uzerinden)
em_cutoffs <- get_em_cutoffs(ens, metric = "TSS")

if (!is.null(em_cutoffs) && nrow(em_cutoffs) > 0) {
  # Belirsizlik olmamasi icin her iki olcek de kaydedilir.
  em_cutoffs$cutoff_0_1000 <- em_cutoffs$cutoff
  em_cutoffs$cutoff_0_1    <- em_cutoffs$cutoff / 1000
  em_cutoffs$raster_scale  <- if (SCALE_DIV > 1) "0-1" else "0-1000"
  em_cutoffs$cutoff_for_raster <- em_cutoffs$cutoff / SCALE_DIV
  write.csv(em_cutoffs,
            file.path(cfg$out_dir, "ensemble_maxTSS_cutoffs.csv"),
            row.names = FALSE)
  message("MaxTSS cutoff araligi: ",
          paste(round(range(em_cutoffs$cutoff, na.rm = TRUE), 3), collapse = " - "))
} else {
  warning("MaxTSS cutoff tablosu bos -> ikili haritalar uretilmeyecek.")
}

# -------------------- ESIKLENDIRILMIS TEPKI EGRILERI (yayin formati) --------------------
# NOT: bu blok esik tablosu YAZILDIKTAN SONRA calismalidir; aksi halde esik
# bulunamaz ve egri medyanina geri cekilir.
set_status("RESPONSE CURVES (threshold-coloured)")
rct_file <- file.path(getwd(), "plot_response_threshold.R")
if (file.exists(rct_file)) {
  source(rct_file, encoding = "UTF-8")
  try(plot_response_threshold(
    out_dir    = cfg$out_dir,
    sp_name    = cfg$sp_name %||% "",
    var_types  = cfg$var_types,
    em_pick    = cfg$rc_em_pick %||% "EMwmean",
    var_labels = cfg$var_labels,
    y_lab      = "Habitat suitability (ensemble)"
  ), silent = TRUE)
} else {
  warning("plot_response_threshold.R bulunamadi -> esiklendirilmis egri atlandi.")
}


# helper: load env from folder and align names to reference (current env)
load_env_asc <- function(asc_dir, ref_names = NULL) {
  af <- list.files(
    asc_dir,
    pattern = "\\.(asc|tif|tiff|grd)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(af) == 0) stop("No raster files found in: ", asc_dir, " (expected .asc/.tif/.grd)")

  e <- rast(af)
  names(e) <- make.names(names(e), unique = TRUE)

  # Apply variable types to FUTURE env too
  e <- apply_var_types(e, cfg)

  # Ensure CRS matches current env (reference)
  ref_crs <- terra::crs(env)
  if (!is.na(ref_crs) && nzchar(ref_crs)) {
    if (is.na(terra::crs(e)) || !nzchar(terra::crs(e))) {
      terra::crs(e) <- ref_crs
    } else if (terra::crs(e) != ref_crs) {
      e <- terra::project(e, ref_crs)
    }
  }

  if (!is.null(ref_names)) {
    miss <- setdiff(ref_names, names(e))
    extra <- setdiff(names(e), ref_names)
    if (length(miss) > 0) stop("Future env is missing layers: ", paste(miss, collapse = ", "))
    e <- e[[ref_names]]
    if (length(extra) > 0) warning("Future env had extra layers dropped: ", paste(extra, collapse = ", "))
  }
  e
}

proj_cur <- BIOMOD_Projection(
  bm.mod  = model_out,
  new.env = env,
  proj.name = "current",
  models.chosen = "all",
  compress = "xz",
  build.clamping.mask = TRUE
)

ens_cur <- tryCatch({
  BIOMOD_EnsembleForecasting(
    bm.em = ens,
    bm.proj = proj_cur,
    proj.name = "current_ens"
  )
}, error = function(e) {
  warning("Ensemble forecasting (current) basarisiz: ", e$message)
  NULL
})

# ENSEMBLE rasterleri: surekli + MaxTSS ikili
if (!is.null(ens_cur)) {
  write_proj_rasters(ens_cur, "current_ensemble", rasters_dir, em_cutoffs,
                     scale_div = SCALE_DIV)
}

tryCatch({
  cur_pred <- get_predictions(proj_cur)
  if (!is.null(cur_pred)) {
    if (SCALE_DIV > 1) cur_pred <- cur_pred / SCALE_DIV
    terra::writeRaster(
      cur_pred,
      file.path(rasters_dir, "current_projection.tif"),
      overwrite = TRUE
    )
  }
}, error = function(e) {
  warning("Failed to export current projection TIF: ", e$message)
})

do_future <- !is.null(cfg$future_n) &&
  is.numeric(cfg$future_n) &&
  cfg$future_n >= 1 &&
  !is.null(cfg$future) &&
  length(cfg$future) >= 1

if (!do_future) {
  set_status("DONE (current only)")
} else {
  set_status("PROJECTION (future)")

  ref_names <- names(env)
  nf <- min(as.integer(cfg$future_n), length(cfg$future))

  for (i in seq_len(nf)) {
    fblock <- cfg$future[[i]]
    if (!is.list(fblock)) next
    if (is.null(fblock$scenarios) || length(fblock$scenarios) == 0) next

    period_label <- fblock$period_label
    if (is.null(period_label) || !nzchar(period_label)) period_label <- paste0("Future", i)

    scs <- fblock$scenarios
    if (!is.list(scs)) scs <- list(scs)

    for (j in seq_along(scs)) {
      sc <- scs[[j]]
      if (!is.list(sc)) next

      sc_label <- sc$label
      if (is.null(sc_label) || !nzchar(sc_label)) sc_label <- paste0("Scenario", j)

      asc_dir_future <- sc$asc_dir
      if (is.null(asc_dir_future) || !nzchar(asc_dir_future)) next

      set_status(paste0("FUTURE PROJECTION: ", period_label, " / ", sc_label))

      env_f <- load_env_asc(asc_dir_future, ref_names = ref_names)

      proj_name <- paste0(
        "future_", i, "_", j, "_",
        make.names(period_label), "_",
        make.names(sc_label)
      )

      proj_f <- BIOMOD_Projection(
        bm.mod  = model_out,
        new.env = env_f,
        proj.name = proj_name,
        models.chosen = "all",
        compress = "xz",
        build.clamping.mask = TRUE
      )

      ens_f <- tryCatch({
        BIOMOD_EnsembleForecasting(
          bm.em = ens,
          bm.proj = proj_f,
          proj.name = paste0(proj_name, "_ens")
        )
      }, error = function(e) {
        warning("Ensemble forecasting basarisiz (", proj_name, "): ", e$message)
        NULL
      })

      # ENSEMBLE rasterleri: surekli + MaxTSS ikili (esik gunumuzden, sabit)
      if (!is.null(ens_f)) {
        write_proj_rasters(ens_f, paste0(proj_name, "_ensemble"),
                           rasters_dir, em_cutoffs, scale_div = SCALE_DIV)
      }

      tryCatch({
        f_pred <- get_predictions(proj_f)
        if (!is.null(f_pred)) {
          if (SCALE_DIV > 1) f_pred <- f_pred / SCALE_DIV
          terra::writeRaster(
            f_pred,
            file.path(rasters_dir, paste0(proj_name, ".tif")),
            overwrite = TRUE
          )
        }
      }, error = function(e) {
        warning("Failed to export future projection TIF for ", proj_name, ": ", e$message)
      })
    }
  }

  set_status("DONE")
}




