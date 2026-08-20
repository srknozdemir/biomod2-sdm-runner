# run_biomod_cont.R
# =============================================================================
# biomod2 SUREKLI / NOMINAL pipeline
#   surekli : count | abundance | relative      -> regresyon
#   nominal : multiclass | ordinal              -> siniflandirma
#
# app.R bu betigi cfg$data_type != "binary" oldugunda calistirir.
# Var/Yok (binary) icin run_biomod.R kullanilmaya devam eder.
#
# Cikti duzeni run_biomod.R ile AYNIDIR (plots/, rasters/, eval_*.csv, status.txt)
# boylece Shiny izleme/gorsellestirme arayuzu degismeden calisir.
#
# biomod2 >= 4.3 gerekir (data.type destegi).
#
# TEMEL FARKLAR (binary'ye gore):
#  - Pseudo-absence YOK. Yanit degiskeni CSV'deki gercek olcum sutunudur.
#  - Model havuzu veri tipine gore daralir (MAXENT/MAXNET/SRE/ANN/RFd yok).
#  - Metrikler: RMSE/MAE/Rsquared/Rsquared_aj  (nicel)
#               Accuracy/F1/Recall/Precision   (nitel)
#  - RMSE/MAE/MSE/Max_error icin esik mantigi TERSTIR: "en iyi + tolerans".
#  - Ensemble algoritmalari: EMmean/EMmedian/EMwmean/EMcv/EMci (nicel)
#                            EMmode/EMfreq                    (nitel)
#  - MaxTSS ikili harita YOK. Yerine sabit kuantil sinif kesimleri (opsiyonel).
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(terra)
  library(biomod2)
})

# --- biomod2 uyumluluk ---
if (!exists("getModelsBuilt", mode = "function") && exists("get_built_models", mode = "function")) {
  getModelsBuilt <- get_built_models
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_biomod_cont.R <config.json>")

cfg <- jsonlite::read_json(args[1], simplifyVector = FALSE)

# -------------------- SURUM DENETIMI --------------------
# data.type (count/abundance/relative/multiclass/ordinal) destegi biomod2 4.3
# ile geldi. Daha eski surumde BIOMOD_FormatingData bu argumani tanimaz ve
# hata mesaji anlasilmaz olur; bu yuzden acik denetim yapiyoruz.
bm_ver <- tryCatch(utils::packageVersion("biomod2"), error = function(e) NULL)
if (is.null(bm_ver) || bm_ver < "4.3") {
  stop("Bu betik biomod2 >= 4.3 gerektirir (surekli/nominal veri tipi destegi). ",
       "Yuklu surum: ", if (is.null(bm_ver)) "yok" else as.character(bm_ver),
       ". Guncelleme: devtools::install_github('biomodhub/biomod2')")
}
message("biomod2 surumu: ", as.character(bm_ver))

# -------------------- PAYLASILAN SECENEK MOTORU --------------------
# Literatur temelli parametre izgarasi run_biomod.R ile ORTAKTIR.
opts_file <- file.path(getwd(), "biomod_opts.R")
if (!file.exists(opts_file)) {
  stop("biomod_opts.R bulunamadi. app.R, run_biomod_cont.R ve biomod_opts.R ayni klasorde olmalidir.")
}
source(opts_file, encoding = "UTF-8")

# -------------------- HELPERS --------------------
`%||%` <- function(a, b) {
  ok <- !is.null(a) && length(a) > 0 && !all(is.na(a))
  if (isTRUE(ok)) a else b
}

as_chr_vec <- function(x) {
  if (is.null(x)) return(character(0))
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  unique(x[nzchar(x)])
}

as_num1 <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE))[1])
  if (is.na(v)) default else v
}

as_int1 <- function(x, default = NA_integer_) {
  v <- as_num1(x, NA_real_)
  if (is.na(v)) default else as.integer(v)
}

first_present <- function(nms, cands) {
  hit <- cands[cands %in% nms]
  if (length(hit) > 0) hit[1] else NULL
}

# Bir fonksiyonu cagirir; hata verirse opsiyonel argumanlari sirayla dusurerek
# tekrar dener. biomod2 surum farkliliklarina karsi koruma saglar.
call_with_fallback <- function(fn, args_list, optional = character(0), label = "") {
  out <- tryCatch(do.call(fn, args_list), error = function(e) e)
  if (!inherits(out, "error")) return(out)
  if (length(optional) > 0) {
    for (k in seq_along(optional)) {
      drop_k <- optional[seq_len(k)]
      a2 <- args_list[setdiff(names(args_list), drop_k)]
      out2 <- tryCatch(do.call(fn, a2), error = function(e) e)
      if (!inherits(out2, "error")) {
        message("NOT (", label, "): su argumanlar dusuruldu -> ",
                paste(drop_k, collapse = ", "))
        return(out2)
      }
    }
  }
  stop("[", label, "] ", conditionMessage(out))
}

# safe_png -> biomod_opts.R (paylasilan, platformdan bagimsiz cihaz zinciri)

safe_write_csv <- function(df, path) {
  ok <- FALSE
  tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, path, row.names = FALSE)
    ok <- TRUE
  }, error = function(e) ok <<- FALSE)
  ok
}

# -------------------- VERI TIPI TABLOLARI --------------------
DT_QUANT <- c("count", "abundance", "relative")   # nicel (regresyon)
DT_QUAL  <- c("multiclass", "ordinal")            # nitel (siniflandirma)

# biomod2 4.3 model uygunluk matrisi (binary disi)
ALLOWED_MODELS <- list(
  count      = c("CTA", "DNN", "GAM", "GBM", "GLM", "MARS", "RF", "XGBOOST"),
  abundance  = c("CTA", "DNN", "GAM", "GBM", "GLM", "MARS", "RF", "XGBOOST"),
  relative   = c("CTA", "DNN", "GAM", "GBM", "GLM", "MARS", "RF", "XGBOOST"),
  ordinal    = c("CTA", "DNN", "FDA", "GAM", "GLM", "MARS", "RF", "XGBOOST"),
  multiclass = c("CTA", "DNN", "FDA", "MARS", "RF", "XGBOOST")
)

REQ_PKGS <- list(
  CTA = "rpart", DNN = "cito", FDA = "mda", GAM = "mgcv", GBM = "gbm",
  MARS = "earth", RF = "randomForest", XGBOOST = "xgboost"
)

METRICS_QUANT <- c("Rsquared", "Rsquared_aj", "RMSE", "MAE")
METRICS_QUAL  <- c("Accuracy", "F1", "Recall", "Precision")
LOWER_BETTER  <- c("RMSE", "MSE", "MAE", "Max_error")

EM_ALGO_QUANT <- c("EMmean", "EMmedian", "EMwmean", "EMcv", "EMci")
EM_ALGO_QUAL  <- c("EMmode", "EMfreq")

# -------------------- CONFIG NORMALIZASYONU --------------------
data_type <- tolower(as.character(cfg$data_type %||% "abundance")[1])
if (!data_type %in% c(DT_QUANT, DT_QUAL)) {
  stop("Gecersiz cfg$data_type: '", data_type,
       "'. Beklenen: count/abundance/relative/multiclass/ordinal.")
}
is_quant <- data_type %in% DT_QUANT
is_qual  <- data_type %in% DT_QUAL

cfg$models <- as_chr_vec(cfg$models)
if (length(cfg$models) == 0) stop("cfg$models bos.")

if (!is.null(cfg$future_n)) {
  if (is.character(cfg$future_n) && toupper(cfg$future_n) == "NA") cfg$future_n <- NULL
  if (!is.null(cfg$future_n)) cfg$future_n <- as_int1(cfg$future_n, NA_integer_)
  if (!is.null(cfg$future_n) && is.na(cfg$future_n)) cfg$future_n <- NULL
}
if (!is.null(cfg$future)) {
  if (!is.list(cfg$future)) cfg$future <- list(cfg$future)
  for (k in seq_along(cfg$future)) {
    fb <- cfg$future[[k]]
    if (!is.list(fb)) next
    if (!is.null(fb$scenarios) && !is.list(fb$scenarios)) {
      cfg$future[[k]]$scenarios <- list(fb$scenarios)
    }
  }
}

normalize_var_types <- function(cfg) {
  vt <- cfg$var_types
  if (is.null(vt) || !is.list(vt)) return(NULL)
  if (is.null(names(vt)) || any(!nzchar(names(vt)))) return(NULL)
  for (nm in names(vt)) {
    vv <- vt[[nm]]
    vt[[nm]] <- if (is.null(vv)) "continuous" else tolower(as.character(vv)[1])
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

if (is.null(cfg$out_dir) || !nzchar(cfg$out_dir)) {
  cfg$out_dir <- file.path("runs", cfg$run_id, "out")
}
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
write_session_info(cfg$out_dir)
plots_dir   <- file.path(cfg$out_dir, "plots");   dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
rasters_dir <- file.path(cfg$out_dir, "rasters"); dir.create(rasters_dir, recursive = TRUE, showWarnings = FALSE)
rc_dir      <- file.path(cfg$out_dir, "response_curves"); dir.create(rc_dir, recursive = TRUE, showWarnings = FALSE)

set_status(paste0("STARTED (data.type=", data_type, ")"))

# -------------------- ESIKLER --------------------
# Nicel: Rsquared alt sinir + RMSE toleransi (en iyi RMSE + tol)
# Nitel: Accuracy / F1 alt sinirlari
ENS_MIN_R2   <- as_num1(cfg$ens_min_r2,   0.30)
ENS_RMSE_TOL <- as_num1(cfg$ens_rmse_tol, NA_real_)   # NA -> veriden turetilir
ENS_MIN_ACC  <- as_num1(cfg$ens_min_acc,  0.50)
ENS_MIN_F1   <- as_num1(cfg$ens_min_f1,   0.40)
ENS_NEAR_ONE <- as_num1(cfg$ens_near_one, 0.999)      # asiri uyum supheli

# -------------------- ENV --------------------
set_status("LOADING DATA")

apply_var_types <- function(e, cfg) {
  vt <- cfg$var_types
  if (is.null(vt)) {
    cv <- as_chr_vec(cfg$cat_vars)
    for (v in intersect(cv, names(e))) e[[v]] <- terra::as.factor(e[[v]])
    return(e)
  }
  for (v in intersect(names(e), names(vt))) {
    tp <- tolower(as.character(vt[[v]] %||% "continuous")[1])
    if (tp %in% c("categorical", "nominal")) {
      tv <- Sys.time()
      message("Converting to categorical (reads all cells): ", v, " ...")
      e[[v]] <- terra::as.factor(e[[v]])
      message("  done in ", round(as.numeric(difftime(Sys.time(), tv, units = "secs")), 1), " s")
    }
  }
  e
}

read_env_dir <- function(dir_path) {
  if (is.null(dir_path) || !nzchar(dir_path)) stop("Raster klasor yolu bos.")
  if (!dir.exists(dir_path)) stop("Raster klasoru bulunamadi: ", dir_path)
  files <- list.files(dir_path, pattern = "\\.(asc|tif|tiff|grd)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("Klasorde raster yok (.asc/.tif/.grd): ", dir_path)
  grd <- files[grepl("\\.grd$", files, ignore.case = TRUE)]
  if (length(grd) > 0) {
    gri <- sub("\\.grd$", ".gri", grd, ignore.case = TRUE)
    files <- unique(c(files, gri[file.exists(gri)]))
  }
  e <- terra::rast(files)
  names(e) <- make.names(names(e), unique = TRUE)
  e
}

env <- read_env_dir(cfg$asc_dir)
env <- apply_var_types(env, cfg)
message("ENV katmanlari: ", paste(names(env), collapse = ", "))

# -------------------- YANIT DEGISKENI --------------------
occ <- read_user_csv(cfg$occ_file)

sp_col <- first_present(names(occ), c("Tur", "tur", "species", "Species", "sp", "taxon"))
if (!is.null(sp_col) && !is.null(cfg$sp_name) && nzchar(cfg$sp_name)) {
  sel <- as.character(occ[[sp_col]]) == cfg$sp_name
  if (any(sel, na.rm = TRUE)) occ <- occ[which(sel), , drop = FALSE]
}

resp_col <- as.character(cfg$resp_col %||% "")[1]
if (!nzchar(resp_col) || !(resp_col %in% names(occ))) {
  stop("Yanit sutunu bulunamadi. cfg$resp_col = '", resp_col,
       "'. CSV sutunlari: ", paste(names(occ), collapse = ", "))
}

resp_raw <- occ[[resp_col]]
# Koordinat sutunlari "x"/"y" olmak zorunda degil (bkz. biomod_opts.R)
occ_xy   <- get_xy_matrix(occ, cfg$x_col, cfg$y_col, what = "Occurrence CSV")

# Cevresel degiskenlerde NA olan noktalari at
ex <- terra::extract(env, occ_xy)
if ("ID" %in% names(ex)) ex$ID <- NULL
keep <- stats::complete.cases(ex) & !is.na(resp_raw) &
  is.finite(occ_xy[, 1]) & is.finite(occ_xy[, 2])
occ_xy   <- occ_xy[keep, , drop = FALSE]
resp_raw <- resp_raw[keep]
if (nrow(occ_xy) < 20) stop("Gecerli gozlem sayisi cok dusuk (n=", nrow(occ_xy), ").")

# Veri tipine gore yanit vektorunu kur + dogrula
build_response <- function(v, dt, ord_levels = NULL) {
  if (dt %in% DT_QUANT) {
    x <- suppressWarnings(as.numeric(v))
    if (all(is.na(x))) stop("Yanit sutunu sayisala cevrilemedi (dt=", dt, ").")
    if (dt == "count") {
      if (any(x < 0, na.rm = TRUE)) stop("count verisi negatif deger iceremez.")
      if (any(abs(x - round(x)) > 1e-8, na.rm = TRUE)) {
        warning("count verisinde ondalik degerler var -> yuvarlaniyor.")
        x <- round(x)
      }
      x <- as.integer(x)
    }
    if (dt == "relative" && (any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE))) {
      stop("relative verisi 0-1 araliginda olmalidir (yuzde ise 100'e bolun).")
    }
    if (dt == "abundance" && any(x < 0, na.rm = TRUE)) {
      stop("abundance verisi negatif deger iceremez.")
    }
    return(x)
  }
  # nitel
  x <- as.character(v)
  if (dt == "ordinal") {
    if (!is.null(ord_levels) && length(ord_levels) > 1) {
      miss <- setdiff(unique(x), ord_levels)
      if (length(miss) > 0) stop("Ordinal siralamada tanimsiz sinif(lar): ", paste(miss, collapse = ", "))
      return(factor(x, levels = ord_levels, ordered = TRUE))
    }
    num_try <- suppressWarnings(as.numeric(x))
    lev <- if (!any(is.na(num_try))) as.character(sort(unique(num_try))) else sort(unique(x))
    return(factor(x, levels = lev, ordered = TRUE))
  }
  factor(x, levels = sort(unique(x)))
}

resp <- build_response(resp_raw, data_type, as_chr_vec(cfg$ordinal_levels))

n_obs <- length(resp)
if (is_qual) {
  cls_tab <- table(resp)
  if (length(cls_tab) < 2) stop("Nitel yanit tek sinif iceriyor.")
  message("Sinif dagilimi: ", paste(names(cls_tab), cls_tab, sep = "=", collapse = ", "))
} else {
  message("Yanit ozeti: n=", n_obs,
          " min=", round(min(resp, na.rm = TRUE), 4),
          " ortalama=", round(mean(resp, na.rm = TRUE), 4),
          " max=", round(max(resp, na.rm = TRUE), 4),
          " sifir orani=", round(mean(resp == 0, na.rm = TRUE), 3))
}

set_status(paste0("FORMATTING DATA (n=", n_obs, ", data.type=", data_type, ")"))

# -------------------- BIOMOD_FormatingData (PA YOK) --------------------
fmt_args <- list(
  resp.name    = cfg$sp_name,
  resp.var     = resp,
  resp.xy      = occ_xy,
  expl.var     = env,
  data.type    = data_type,
  na.rm        = TRUE,
  filter.raster = isTRUE(cfg$filter_raster)
)

biomod_data <- call_with_fallback(
  BIOMOD_FormatingData, fmt_args,
  optional = c("filter.raster", "na.rm"),
  label = "BIOMOD_FormatingData"
)

# -------------------- MODEL HAVUZU --------------------
set_status("MODELING")

models_to_run <- cfg$models
allowed <- ALLOWED_MODELS[[data_type]]

dropped_unsupported <- setdiff(models_to_run, allowed)
if (length(dropped_unsupported) > 0) {
  warning("'", data_type, "' veri tipinde desteklenmeyen modeller atildi: ",
          paste(dropped_unsupported, collapse = ", "))
}
models_to_run <- intersect(models_to_run, allowed)

for (m in intersect(names(REQ_PKGS), models_to_run)) {
  if (!requireNamespace(REQ_PKGS[[m]], quietly = TRUE)) {
    warning("Paket yok ('", REQ_PKGS[[m]], "') -> model atildi: ", m)
    models_to_run <- setdiff(models_to_run, m)
  }
}
if (length(models_to_run) == 0) {
  stop("Calistirilacak model kalmadi. '", data_type, "' icin uygun olanlar: ",
       paste(allowed, collapse = ", "))
}
message("Calistirilacak modeller: ", paste(models_to_run, collapse = ", "))

# -------------------- CV --------------------
cv_k   <- as_int1(cfg$cv$k, 10L)
cv_rep <- as_int1(cfg$cv$rep, 1L)
if (is.na(cv_k) || cv_k < 2) cv_k <- 10L
if (is.na(cv_rep) || cv_rep < 1) cv_rep <- 1L

# Nitel veride en kucuk sinif k'dan kucukse fold sayisini dusur
if (is_qual) {
  min_cls <- min(as.integer(table(resp)))
  if (min_cls < cv_k) {
    cv_k <- max(2L, min_cls)
    warning("En kucuk sinif ", min_cls, " gozlem iceriyor -> CV.k = ", cv_k)
  }
}

nb_cpu <- as_int1(cfg$cpu_n, 1L); if (is.na(nb_cpu) || nb_cpu < 1) nb_cpu <- 1L

metric_eval <- if (is_quant) METRICS_QUANT else METRICS_QUAL

# -------------------- LITERATUR TEMELLI MODELLEME SECENEKLERI --------------------
# Parametre izgarasi biomod_opts.R'dedir ve run_biomod.R ile ORTAKTIR.
p_env <- terra::nlyr(env)

set_status(paste0("BUILDING MODEL OPTIONS (p=", p_env, ", n=", n_obs, ")"))

built_opts <- bmopt_build(
  data_type  = data_type,
  models     = models_to_run,
  bm_format  = biomod_data,
  p          = p_env,
  n          = n_obs,
  maxent_jar = NULL,                 # MAXENT binary disinda kullanilamaz
  strategy   = cfg$opt_strategy %||% "adaptive",
  transfer   = TRUE,
  verbose    = TRUE
)

bmopt_report(built_opts, file.path(cfg$out_dir, "modeling_options_applied.csv"),
             data_type = data_type, p = p_env, n = n_obs)

# -------------------- BIOMOD_Modeling --------------------
mod_args <- list(
  bm.format   = biomod_data,
  modeling.id = paste0("ab_", data_type),
  models      = models_to_run,
  CV.strategy = "kfold",
  CV.k        = cv_k,
  CV.nb.rep   = cv_rep,
  metric.eval = metric_eval,
  var.import  = 3,
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

# -------------------- DEGERLENDIRME TABLOSU --------------------
set_status("EXPORTING METRICS (model)")

# biomod2 >= 4.2 get_evaluations() bir data.frame dondurur.
eval_to_long <- function(ev) {
  if (is.null(ev)) return(NULL)
  df <- if (is.data.frame(ev)) as.data.frame(ev, stringsAsFactors = FALSE)
        else tryCatch(as.data.frame(ev, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  nm <- names(df)
  mod_col  <- first_present(nm, c("full.name", "model.name", "model"))
  algo_col <- first_present(nm, c("algo", "model"))
  met_col  <- first_present(nm, c("metric.eval", "metric"))
  run_col  <- first_present(nm, c("run", "RUN"))
  if (is.null(mod_col) || is.null(met_col)) return(NULL)

  val_cols <- intersect(c("calibration", "validation", "evaluation"), nm)
  if (length(val_cols) == 0) {
    val_cols <- nm[vapply(df, is.numeric, logical(1))]
    val_cols <- setdiff(val_cols, c("cutoff", "sensitivity", "specificity"))
    if (length(val_cols) == 0) return(NULL)
  }

  out <- do.call(rbind, lapply(val_cols, function(cc) {
    data.frame(
      full.name = as.character(df[[mod_col]]),
      algo      = if (is.null(algo_col)) sub(".*_", "", as.character(df[[mod_col]]))
                  else as.character(df[[algo_col]]),
      run       = if (is.null(run_col)) NA_character_ else as.character(df[[run_col]]),
      metric    = as.character(df[[met_col]]),
      dataset   = cc,
      value     = suppressWarnings(as.numeric(df[[cc]])),
      stringsAsFactors = FALSE
    )
  }))
  out[is.finite(out$value), , drop = FALSE]
}

eval_mod_long <- eval_to_long(get_evaluations(model_out))
if (is.null(eval_mod_long)) warning("Model degerlendirme tablosu okunamadi.")
if (!is.null(eval_mod_long)) {
  safe_write_csv(eval_mod_long, file.path(cfg$out_dir, "eval_model.csv"))
}

# Hangi veri seti kullanilacak: validation varsa o, yoksa calibration
pick_dataset <- function(el) {
  if (is.null(el)) return("calibration")
  if (any(el$dataset == "validation" & is.finite(el$value))) "validation" else "calibration"
}
use_data <- pick_dataset(eval_mod_long)
message("Model secimi icin kullanilan veri seti: ", use_data)

# algo x metric ortalama tablosu
eval_wide <- function(el, ds) {
  if (is.null(el)) return(NULL)
  sub <- el[el$dataset == ds, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  agg <- stats::aggregate(value ~ algo + metric, data = sub, FUN = mean, na.rm = TRUE)
  w <- stats::reshape(agg, idvar = "algo", timevar = "metric", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w
}
mod_wide <- eval_wide(eval_mod_long, use_data)

# -------------------- ENSEMBLE ICIN MODEL SECIMI --------------------
# NICEL : Rsquared >= esik  VE  Rsquared < near_one  VE  RMSE <= en_iyi_RMSE + tolerans
# NITEL : Accuracy >= esik  VE  Accuracy < near_one  VE  F1 >= esik
set_status("ENSEMBLE")

rmse_tol <- ENS_RMSE_TOL
if (is_quant && (is.na(rmse_tol) || rmse_tol <= 0)) {
  # Veriden turetilen makul tolerans: yanit std sapmasinin %25'i
  sdv <- stats::sd(as.numeric(resp), na.rm = TRUE)
  rmse_tol <- if (is.finite(sdv) && sdv > 0) 0.25 * sdv else 1
  message("RMSE toleransi otomatik: ", round(rmse_tol, 4))
}

choose_models_nb <- function(w) {
  if (is.null(w) || nrow(w) == 0) return(character(0))
  g <- function(cc) if (cc %in% names(w)) suppressWarnings(as.numeric(w[[cc]])) else rep(NA_real_, nrow(w))

  if (is_quant) {
    r2   <- g("Rsquared")
    rmse <- g("RMSE")
    ok <- rep(TRUE, nrow(w))
    if (any(is.finite(r2))) ok <- ok & is.finite(r2) & r2 >= ENS_MIN_R2 & r2 < ENS_NEAR_ONE
    if (any(is.finite(rmse))) {
      best <- min(rmse[is.finite(rmse)], na.rm = TRUE)
      ok <- ok & is.finite(rmse) & rmse <= best + rmse_tol
    }
  } else {
    acc <- g("Accuracy")
    f1  <- g("F1")
    ok <- rep(TRUE, nrow(w))
    if (any(is.finite(acc))) ok <- ok & is.finite(acc) & acc >= ENS_MIN_ACC & acc < ENS_NEAR_ONE
    if (any(is.finite(f1)))  ok <- ok & is.finite(f1)  & f1  >= ENS_MIN_F1
  }
  as.character(w$algo[which(ok)])
}

keep_algo <- choose_models_nb(mod_wide)
keep_algo <- keep_algo[!tolower(keep_algo) %in% c("model", "none", "na", "null")]

built_full <- tryCatch(as.character(getModelsBuilt(model_out)), error = function(e) character(0))
built_full <- unique(built_full[nzchar(built_full)])

keep_full <- character(0)
if (length(keep_algo) > 0 && length(built_full) > 0) {
  keep_full <- built_full[sub(".*_", "", built_full) %in% keep_algo]
}
models_chosen <- if (length(keep_full) > 0) keep_full else "all"

safe_write_csv(
  data.frame(model = if (length(keep_algo) > 0) keep_algo else character(0)),
  file.path(cfg$out_dir, "ensemble_candidates_after_thresholds.csv")
)

# -------------------- BIOMOD_EnsembleModeling --------------------
em_algo <- as_chr_vec(cfg$em_algo)
em_algo <- intersect(em_algo, if (is_quant) EM_ALGO_QUANT else EM_ALGO_QUAL)
if (length(em_algo) == 0) em_algo <- if (is_quant) c("EMmean", "EMmedian", "EMwmean", "EMcv") else EM_ALGO_QUAL

if (is_quant) {
  metric_select       <- c("Rsquared", "RMSE")
  metric_select_thresh <- c(ENS_MIN_R2, rmse_tol)
} else {
  metric_select       <- c("Accuracy", "F1")
  metric_select_thresh <- c(ENS_MIN_ACC, ENS_MIN_F1)
}
metric_select        <- intersect(metric_select, metric_eval)
metric_select_thresh <- metric_select_thresh[seq_along(metric_select)]

ens_args <- list(
  bm.mod               = model_out,
  models.chosen        = models_chosen,
  em.by                = "all",
  em.algo              = em_algo,
  metric.select        = metric_select,
  metric.select.thresh = metric_select_thresh,
  metric.select.dataset = use_data,
  metric.eval          = metric_eval,
  var.import           = 3,
  EMci.alpha           = 0.05
)

ens <- tryCatch(
  call_with_fallback(BIOMOD_EnsembleModeling, ens_args,
                     optional = c("EMci.alpha", "var.import", "metric.select.dataset"),
                     label = "BIOMOD_EnsembleModeling"),
  error = function(e) { warning("Ensemble kurulamadi: ", conditionMessage(e)); NULL }
)

ens_built <- character(0)
if (!is.null(ens)) {
  ens_built <- tryCatch(as.character(getModelsBuilt(ens)), error = function(e) character(0))
  ens_built <- unique(ens_built[nzchar(ens_built)])
}
safe_write_csv(
  data.frame(model = if (length(ens_built) > 0) ens_built else "NONE"),
  file.path(cfg$out_dir, "ensemble_selected_models.csv")
)

if (!is.null(ens)) {
  eval_ens_long <- eval_to_long(get_evaluations(ens))
  if (!is.null(eval_ens_long)) {
    safe_write_csv(eval_ens_long, file.path(cfg$out_dir, "eval_ensemble.csv"))
  }
} else {
  eval_ens_long <- NULL
}

# -------------------- OZET TABLO --------------------
build_summary <- function(el, ds, type_lab) {
  w <- eval_wide(el, ds)
  if (is.null(w)) return(NULL)
  names(w)[names(w) == "algo"] <- "model"
  keep_metrics <- intersect(metric_eval, names(w))
  out <- w[, c("model", keep_metrics), drop = FALSE]
  names(out)[-1] <- paste0(names(out)[-1], "_mean")
  out$dataset_used <- ds
  out$type <- type_lab
  out
}

base_summary <- build_summary(eval_mod_long, use_data, "base_model")
if (!is.null(base_summary)) {
  base_summary$in_ensemble <- base_summary$model %in% keep_algo
}
ens_summary <- build_summary(eval_ens_long, use_data, "ensemble")
if (!is.null(ens_summary)) ens_summary$in_ensemble <- NA

sum_out <- NULL
if (!is.null(base_summary) && !is.null(ens_summary)) {
  common <- union(names(base_summary), names(ens_summary))
  for (cc in setdiff(common, names(base_summary))) base_summary[[cc]] <- NA
  for (cc in setdiff(common, names(ens_summary)))  ens_summary[[cc]]  <- NA
  sum_out <- rbind(base_summary[, common, drop = FALSE], ens_summary[, common, drop = FALSE])
} else {
  sum_out <- base_summary %||% ens_summary
}
if (!is.null(sum_out)) {
  first_cols <- intersect(c("type", "model"), names(sum_out))
  sum_out <- sum_out[, c(first_cols, setdiff(names(sum_out), first_cols)), drop = FALSE]
  sum_out <- sum_out[order(sum_out$type, sum_out$model), , drop = FALSE]
  safe_write_csv(sum_out, file.path(cfg$out_dir, "auc_tss_summary.csv"))  # ayni dosya adi: arayuz uyumu
}

# -------------------- PLOTLAR --------------------
set_status("SAVING DEFAULT BIOMOD PLOTS")

m1 <- metric_eval[seq_len(min(2L, length(metric_eval)))]
m2 <- if (length(metric_eval) >= 4) metric_eval[3:4] else m1

safe_png(file.path(plots_dir, "biomod_eval_mean.png"),
         bm_PlotEvalMean(bm.out = model_out, metric.eval = m1,
                         dataset = use_data, group.by = "algo", do.plot = TRUE))

safe_png(file.path(plots_dir, "biomod_eval_boyce.png"),
         bm_PlotEvalMean(bm.out = model_out, metric.eval = m2,
                         dataset = use_data, group.by = "algo", do.plot = TRUE))

safe_png(file.path(plots_dir, "biomod_eval_boxplot.png"),
         bm_PlotEvalBoxplot(bm.out = model_out, dataset = use_data,
                            group.by = c("algo", "algo"), do.plot = TRUE))

safe_png(file.path(plots_dir, "biomod_default_plot.png"), plot(model_out))

# Artik/uyum analizi (yalnizca nicel veri) - biomod2 >= 4.3
if (is_quant) {
  safe_png(file.path(plots_dir, "model_analysis_residuals.png"),
           bm_ModelAnalysis(bm.mod = model_out, models.chosen = built_full))
}

# -------------------- DEGISKEN ONEMI --------------------
set_status("EXPORTING VAR IMPORTANCE")

vip <- tryCatch(
  bm_PlotVarImpBoxplot(bm.out = model_out,
                       group.by = c("expl.var", "algo", "run"), do.plot = FALSE),
  error = function(e) NULL
)
vip_df <- NULL
if (!is.null(vip)) {
  vip_df <- if (is.data.frame(vip)) vip else (vip$tab %||% vip$data %||% vip[[1]])
  if (!is.data.frame(vip_df)) vip_df <- NULL
}
if (is.null(vip_df)) {
  vi <- tryCatch(get_variables_importance(model_out), error = function(e) NULL)
  if (is.data.frame(vi)) vip_df <- vi
}

if (!is.null(vip_df) && nrow(vip_df) > 0) {
  nmv <- names(vip_df)
  var_col <- first_present(nmv, c("expl.var", "variable"))
  mod_col <- first_present(nmv, c("full.name", "algo", "model"))
  imp_col <- first_present(nmv, c("var.imp", "importance", "value"))
  if (is.null(imp_col)) {
    numc <- nmv[vapply(vip_df, is.numeric, logical(1))]
    if (length(numc) > 0) imp_col <- numc[1]
  }

  if (!is.null(var_col) && !is.null(mod_col) && !is.null(imp_col)) {
    d <- data.frame(
      var = as.character(vip_df[[var_col]]),
      mod = as.character(vip_df[[mod_col]]),
      imp = suppressWarnings(as.numeric(vip_df[[imp_col]])),
      stringsAsFactors = FALSE
    )
    d <- d[is.finite(d$imp), , drop = FALSE]

    if (length(keep_algo) > 0) {
      sel <- (d$mod %in% keep_algo) | (sub(".*_", "", d$mod) %in% keep_algo)
      if (any(sel)) d <- d[sel, , drop = FALSE]
    }

    if (nrow(d) > 0) {
      ens_vi <- stats::aggregate(imp ~ var, data = d, FUN = mean, na.rm = TRUE)
      names(ens_vi) <- c("variable", "importance_mean")
      s <- sum(ens_vi$importance_mean, na.rm = TRUE)
      ens_vi$importance_percent <- if (is.finite(s) && s > 0) 100 * ens_vi$importance_mean / s else NA_real_
      ens_vi <- ens_vi[order(ens_vi$importance_mean, decreasing = TRUE), , drop = FALSE]

      safe_write_csv(ens_vi, file.path(cfg$out_dir, "var_importance_ensemble_avg.csv"))
      safe_write_csv(vip_df, file.path(cfg$out_dir, "var_importance_bmPlot_table_raw.csv"))

      safe_png(file.path(plots_dir, "var_importance_ensemble_avg.png"), {
        op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
        par(mar = c(6, 10, 3, 1))
        x <- ens_vi$importance_mean; names(x) <- ens_vi$variable
        barplot(rev(x), horiz = TRUE, las = 1,
                xlab = "Ortalama degisken onemi (secili modeller)")
        title(paste0("Ensemble-ortalamali degisken onemi (", data_type, ")"))
      })
    }
  }
}

# -------------------- TEPKI EGRILERI --------------------
set_status("RESPONSE CURVES (models)")

pick_response_vars <- function(env, cfg, max_n = 50L) {
  vars <- names(env)
  ov <- as_chr_vec(cfg$response_vars)
  if (length(ov) > 0) {
    ov <- ov[ov %in% vars]
    if (length(ov) > 0) vars <- ov
  }
  vars[seq_len(min(as.integer(max_n), length(vars)))]
}
resp_vars_plot <- pick_response_vars(env, cfg, 50L)

safe_png(file.path(plots_dir, "models_response_curves.png"),
         bm_PlotResponseCurves(bm.out = model_out, models.chosen = "all",
                               show.variables = resp_vars_plot, fixed.var = "median",
                               do.bivariate = FALSE, do.plot = TRUE))

rc_dat <- tryCatch(
  bm_PlotResponseCurves(bm.out = model_out, models.chosen = "all",
                        show.variables = resp_vars_plot, fixed.var = "median",
                        do.bivariate = FALSE, do.plot = FALSE),
  error = function(e) NULL
)
if (!is.null(rc_dat)) {
  d <- if (is.data.frame(rc_dat)) rc_dat else (rc_dat$tab %||% rc_dat[[1]])
  if (is.data.frame(d)) safe_write_csv(d, file.path(rc_dir, "models_response_curves_data.csv"))
}

if (!is.null(ens)) {
  set_status("RESPONSE CURVES (ensemble)")
  safe_png(file.path(plots_dir, "ensemble_response_curves.png"),
           bm_PlotResponseCurves(bm.out = ens, models.chosen = "all",
                                 show.variables = resp_vars_plot, fixed.var = "median",
                                 do.bivariate = FALSE, do.plot = TRUE))
  rc_e <- tryCatch(
    bm_PlotResponseCurves(bm.out = ens, models.chosen = "all",
                          show.variables = resp_vars_plot, fixed.var = "median",
                          do.bivariate = FALSE, do.plot = FALSE),
    error = function(e) NULL
  )
  if (!is.null(rc_e)) {
    d <- if (is.data.frame(rc_e)) rc_e else (rc_e$tab %||% rc_e[[1]])
    if (is.data.frame(d)) safe_write_csv(d, file.path(rc_dir, "ensemble_response_curves_data.csv"))
  }
}


# -------------------- PROJEKSIYON --------------------
set_status("PROJECTION")

proj_digits <- as_int1(cfg$proj_digits, if (data_type == "count") 0L else 3L)

# ---- CIKTI OLCEGI ----
# Yalnizca 'relative' (0-1 oran) biomod2 tarafindan 0-1000'e cevrilir.
# count / abundance kendi biriminde kalir -> asla olceklenmez.
out_scale <- tolower(as.character(cfg$out_scale %||% "0_1")[1])
SCALE_DIV <- if (identical(data_type, "relative") && identical(out_scale, "0_1")) 1000 else 1
if (SCALE_DIV > 1) message("Raster cikti olcegi: 0-1 (1000'e bolunerek)")

do_project <- function(bm_mod, new_env, proj_name) {
  a <- list(
    bm.mod  = bm_mod,
    new.env = new_env,
    proj.name = proj_name,
    models.chosen = "all",
    build.clamping.mask = TRUE,
    digits = proj_digits
  )
  # NOT: on_0_1000 = FALSE GECIRILMEZ (biomod2 issue #72). biomod2 kendi tam sayi
  # yolunda calisir; 0-1'e indirgeme yalnizca disa aktarimda yapilir.
  call_with_fallback(BIOMOD_Projection, a,
                     optional = c("digits", "build.clamping.mask"),
                     label = paste0("BIOMOD_Projection/", proj_name))
}

# Sabit sinif kesimleri: GUNUMUZ ensemble ortalamasindan BIR KEZ hesaplanir,
# tum gelecek senaryolarinda DEGISTIRILMEDEN kullanilir (binary MaxTSS mantiginin
# nicel karsiligi).
get_quantile_breaks <- function(r, probs) {
  if (is.null(r) || !inherits(r, "SpatRaster")) return(NULL)
  v <- terra::values(r[[1]], mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) < 10) return(NULL)
  q <- unname(stats::quantile(v, probs = probs, na.rm = TRUE))
  q <- unique(round(q, 6))
  if (length(q) < 1) return(NULL)
  data.frame(prob = probs[seq_along(q)], break_value = q)
}

class_probs <- suppressWarnings(as.numeric(unlist(cfg$class_probs %||% list(0.25, 0.5, 0.75))))
class_probs <- class_probs[is.finite(class_probs) & class_probs > 0 & class_probs < 1]
do_breaks <- isTRUE(cfg$export_classes) && is_quant && length(class_probs) > 0

write_proj_rasters <- function(proj_obj, tag, breaks = NULL, scale_div = 1) {
  pred <- tryCatch(get_predictions(proj_obj), error = function(e) NULL)
  if (is.null(pred) || !inherits(pred, "SpatRaster")) {
    warning("Tahmin rasteri alinamadi: ", tag); return(invisible(NULL))
  }
  n_c <- 0L; n_b <- 0L
  for (i in seq_len(terra::nlyr(pred))) {
    lyr <- names(pred)[i]
    safe_name <- make.names(lyr)
    ok <- tryCatch({
      r_out <- if (scale_div > 1) pred[[i]] / scale_div else pred[[i]]
      terra::writeRaster(r_out,
                         file.path(rasters_dir, paste0(tag, "_", safe_name, "_continuous.tif")),
                         overwrite = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (isTRUE(ok)) n_c <- n_c + 1L

    if (!is.null(breaks) && nrow(breaks) > 0 && grepl("EMmean|EMmedian|EMwmean", lyr)) {
      ok2 <- tryCatch({
        b <- c(-Inf, breaks$break_value, Inf)
        # from | to | becomes  -> sinif kodlari 1..n olarak acikca tanimlanir
        rcl <- cbind(b[-length(b)], b[-1], seq_len(length(b) - 1L))
        cls <- terra::classify(pred[[i]], rcl = rcl,
                               include.lowest = TRUE, right = TRUE)
        terra::writeRaster(cls,
                           file.path(rasters_dir, paste0(tag, "_", safe_name, "_classes.tif")),
                           overwrite = TRUE, datatype = "INT1U")
        TRUE
      }, error = function(e) FALSE)
      if (isTRUE(ok2)) n_b <- n_b + 1L
    }
  }
  message(tag, ": ", n_c, " surekli, ", n_b, " sinif rasteri yazildi.")
  invisible(NULL)
}

proj_cur <- do_project(model_out, env, "current")

tryCatch({
  cp <- get_predictions(proj_cur)
  if (inherits(cp, "SpatRaster")) {
    if (SCALE_DIV > 1) cp <- cp / SCALE_DIV
    terra::writeRaster(cp, file.path(rasters_dir, "current_projection.tif"), overwrite = TRUE)
  }
}, error = function(e) warning("Gunumuz tekil model rasteri yazilamadi: ", conditionMessage(e)))

ens_cur <- NULL
if (!is.null(ens)) {
  ens_cur <- tryCatch(
    BIOMOD_EnsembleForecasting(bm.em = ens, bm.proj = proj_cur, proj.name = "current_ens"),
    error = function(e) { warning("Ensemble forecasting (gunumuz) basarisiz: ", conditionMessage(e)); NULL }
  )
}

em_breaks <- NULL
if (do_breaks && !is.null(ens_cur)) {
  pc <- tryCatch(get_predictions(ens_cur), error = function(e) NULL)
  if (inherits(pc, "SpatRaster")) {
    idx <- grep("EMmean|EMmedian|EMwmean", names(pc))
    if (length(idx) > 0) {
      em_breaks <- get_quantile_breaks(pc[[idx[1]]], class_probs)
      if (!is.null(em_breaks)) {
        em_breaks$layer <- names(pc)[idx[1]]
        em_breaks$raster_scale <- if (SCALE_DIV > 1) "0-1" else "native"
        em_breaks$break_for_raster <- em_breaks$break_value / SCALE_DIV
        safe_write_csv(em_breaks, file.path(cfg$out_dir, "ensemble_class_breaks.csv"))
        message("Sabit sinif kesimleri: ",
                paste(round(em_breaks$break_value, 4), collapse = " | "))
      }
    }
  }
}

if (!is.null(ens_cur)) write_proj_rasters(ens_cur, "current_ensemble", em_breaks,
                                         scale_div = SCALE_DIV)

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
    y_lab      = if (is_quant) "Predicted response (ensemble)" else "Class probability (ensemble)"
  ), silent = TRUE)
} else {
  warning("plot_response_threshold.R bulunamadi -> esiklendirilmis egri atlandi.")
}


# Nitel veride sinif kodlari icin sozluk
if (is_qual) {
  safe_write_csv(
    data.frame(code = seq_along(levels(resp)), class = levels(resp)),
    file.path(cfg$out_dir, "class_levels.csv")
  )
}

# -------------------- GELECEK SENARYOLARI --------------------
load_env_future <- function(dir_path, ref_names) {
  e <- read_env_dir(dir_path)
  e <- apply_var_types(e, cfg)
  ref_crs <- terra::crs(env)
  if (!is.na(ref_crs) && nzchar(ref_crs)) {
    ce <- terra::crs(e)
    if (is.na(ce) || !nzchar(ce)) terra::crs(e) <- ref_crs
    else if (ce != ref_crs) e <- terra::project(e, ref_crs)
  }
  miss <- setdiff(ref_names, names(e))
  if (length(miss) > 0) stop("Gelecek katmanlarda eksik degiskenler: ", paste(miss, collapse = ", "))
  extra <- setdiff(names(e), ref_names)
  if (length(extra) > 0) warning("Fazla katmanlar atildi: ", paste(extra, collapse = ", "))
  e[[ref_names]]
}

do_future <- !is.null(cfg$future_n) && !is.na(cfg$future_n) && cfg$future_n >= 1 &&
  !is.null(cfg$future) && length(cfg$future) >= 1

if (!do_future) {
  set_status("DONE (current only)")
} else {
  set_status("PROJECTION (future)")
  ref_names <- names(env)
  nf <- min(as.integer(cfg$future_n), length(cfg$future))

  for (i in seq_len(nf)) {
    fb <- cfg$future[[i]]
    if (!is.list(fb) || is.null(fb$scenarios) || length(fb$scenarios) == 0) next
    period_label <- fb$period_label %||% paste0("Future", i)
    scs <- if (is.list(fb$scenarios)) fb$scenarios else list(fb$scenarios)

    for (j in seq_along(scs)) {
      sc <- scs[[j]]
      if (!is.list(sc)) next
      sc_label <- sc$label %||% paste0("Scenario", j)
      dirf <- sc$asc_dir %||% ""
      if (!nzchar(dirf)) next

      set_status(paste0("FUTURE PROJECTION: ", period_label, " / ", sc_label))

      env_f <- tryCatch(load_env_future(dirf, ref_names), error = function(e) e)
      if (inherits(env_f, "error")) {
        warning("Gelecek katmanlari okunamadi (", period_label, "/", sc_label, "): ",
                conditionMessage(env_f))
        next
      }

      proj_name <- paste0("future_", i, "_", j, "_",
                          make.names(period_label), "_", make.names(sc_label))

      proj_f <- tryCatch(do_project(model_out, env_f, proj_name), error = function(e) e)
      if (inherits(proj_f, "error")) {
        warning("Projeksiyon basarisiz (", proj_name, "): ", conditionMessage(proj_f)); next
      }

      tryCatch({
        fp <- get_predictions(proj_f)
        if (inherits(fp, "SpatRaster")) {
          if (SCALE_DIV > 1) fp <- fp / SCALE_DIV
          terra::writeRaster(fp, file.path(rasters_dir, paste0(proj_name, ".tif")),
                             overwrite = TRUE)
        }
      }, error = function(e) warning("Tekil model rasteri yazilamadi (", proj_name, "): ",
                                     conditionMessage(e)))

      if (!is.null(ens)) {
        ens_f <- tryCatch(
          BIOMOD_EnsembleForecasting(bm.em = ens, bm.proj = proj_f,
                                     proj.name = paste0(proj_name, "_ens")),
          error = function(e) { warning("Ensemble forecasting basarisiz (", proj_name, "): ",
                                        conditionMessage(e)); NULL }
        )
        if (!is.null(ens_f)) {
          write_proj_rasters(ens_f, paste0(proj_name, "_ensemble"), em_breaks,
                             scale_div = SCALE_DIV)
        }
      }
    }
  }
  set_status("DONE")
}
