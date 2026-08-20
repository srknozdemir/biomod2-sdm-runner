# run_ssdm.R
# =============================================================================
# STACKED SPECIES DISTRIBUTION MODELLING (S-SDM)
# Potential species richness by stacking per-species binary (MaxTSS) maps.
#
# Called by app.R when cfg$ssdm_enable is TRUE.
#
# Method
#   1. Each species is modelled separately with the selected algorithms.
#      The 0 cells of the species x plot matrix are TRUE ABSENCES, so no
#      pseudo-absences are generated (biomod2 does not need them here).
#   2. Models passing the AUCroc / TSS thresholds enter an ensemble.
#   3. The ensemble map is binarised at its own MaxTSS cut-off.
#   4. Binary maps are summed  ->  bSSDM richness (binary stacking).
#      Continuous maps are summed -> pSSDM richness (probability stacking),
#      written as a secondary layer for comparison.
#
# Outputs (all file names and column names are ENGLISH by design)
#   rasters/richness_current_bSSDM.tif
#   rasters/richness_current_pSSDM.tif
#   rasters/richness_<scenario>_bSSDM.tif
#   rasters/richness_change_<scenario>_bSSDM.tif
#   rasters/species/<species>_current_binary.tif  (per species)
#   species_summary.csv , skipped_species.csv , ssdm_richness_stats.csv
#
# Reference for the stacking approach:
#   Guisan A, Rahbek C (2011) J Biogeogr 38:1433-1444
#   Calabrese JM et al. (2014) Global Ecol Biogeogr 23:99-112
#   D'Amen M et al. (2015) Biol Rev 90:1248-1263
#   Schmitt S et al. (2017) Methods Ecol Evol 8:1795-1803  (SSDM package)
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(terra)
  library(biomod2)
})

if (!exists("getModelsBuilt", mode = "function") && exists("get_built_models", mode = "function")) {
  getModelsBuilt <- get_built_models
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_ssdm.R <config.json>")
cfg <- jsonlite::read_json(args[1], simplifyVector = FALSE)

opts_file <- file.path(getwd(), "biomod_opts.R")
if (!file.exists(opts_file)) {
  stop("biomod_opts.R not found. app.R, run_ssdm.R and biomod_opts.R must share a folder.")
}
source(opts_file, encoding = "UTF-8")

# -------------------- helpers --------------------
`%||%` <- function(a, b) {
  ok <- !is.null(a) && length(a) > 0 && !all(is.na(a))
  if (isTRUE(ok)) a else b
}
as_chr_vec <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(unlist(x, recursive = TRUE, use.names = FALSE))
  unique(x[nzchar(x)])
}
as_num1 <- function(x, d = NA_real_) {
  if (is.null(x)) return(d)
  v <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE))[1])
  if (is.na(v)) d else v
}
as_int1 <- function(x, d = NA_integer_) {
  v <- as_num1(x, NA_real_); if (is.na(v)) d else as.integer(v)
}
first_present <- function(nms, cands) { h <- cands[cands %in% nms]; if (length(h)) h[1] else NULL }


read_user_csv <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("CSV not found: ", path)
  for (enc in c("UTF-8", "", "windows-1254", "latin1")) {
    d <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = enc),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.data.frame(d) && nrow(d) > 0) {
      if (nzchar(enc) && !identical(enc, "UTF-8")) message("CSV read using encoding: ", enc)
      return(d)
    }
  }
  stop("Could not read the CSV (tried UTF-8, native, windows-1254, latin1): ", path)
}

safe_write_csv <- function(df, path) {
  tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, path, row.names = FALSE); TRUE
  }, error = function(e) FALSE)
}

# safe_png -> biomod_opts.R (paylasilan, platformdan bagimsiz cihaz zinciri)

# -------------------- status --------------------
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

if (is.null(cfg$out_dir) || !nzchar(cfg$out_dir)) cfg$out_dir <- file.path("runs", cfg$run_id, "out")
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
write_session_info(cfg$out_dir)
plots_dir   <- file.path(cfg$out_dir, "plots");           dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
rasters_dir <- file.path(cfg$out_dir, "rasters");         dir.create(rasters_dir, recursive = TRUE, showWarnings = FALSE)
sp_ras_dir  <- file.path(rasters_dir, "species");         dir.create(sp_ras_dir, recursive = TRUE, showWarnings = FALSE)
sp_out_dir  <- file.path(cfg$out_dir, "species");         dir.create(sp_out_dir, recursive = TRUE, showWarnings = FALSE)

set_status("STARTED (S-SDM)")

ENS_MIN_AUC  <- as_num1(cfg$ens_min_auc,  0.70)
ENS_MIN_TSS  <- as_num1(cfg$ens_min_tss,  0.40)
ENS_NEAR_ONE <- as_num1(cfg$ens_near_one, 0.999)
# Tek olcut: frekans yuzdesi. 0 = hicbir tur cikarilmaz.
MIN_FREQ_PCT <- as_num1(cfg$ssdm_min_freq, 0)
if (!is.finite(MIN_FREQ_PCT) || MIN_FREQ_PCT < 0) MIN_FREQ_PCT <- 0
CV_K <- 10L
DO_PSSDM     <- isTRUE(cfg$ssdm_pssdm)
DO_CHANGE    <- isTRUE(cfg$ssdm_change)

out_scale <- tolower(as.character(cfg$out_scale %||% "0_1")[1])
SCALE_DIV <- if (identical(out_scale, "0_1")) 1000 else 1

# -------------------- environment --------------------
set_status("LOADING ENVIRONMENT")

apply_var_types <- function(e, cfg) {
  vt <- cfg$var_types
  if (is.null(vt)) {
    for (v in intersect(as_chr_vec(cfg$cat_vars), names(e))) e[[v]] <- terra::as.factor(e[[v]])
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
  if (is.null(dir_path) || !nzchar(dir_path)) stop("Raster folder path is empty.")
  if (!dir.exists(dir_path)) stop("Raster folder not found: ", dir_path)
  files <- list.files(dir_path, pattern = "\\.(asc|tif|tiff|grd)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No rasters (.asc/.tif/.grd) in: ", dir_path)
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
p_env <- terra::nlyr(env)
message("Predictors: ", paste(names(env), collapse = ", "))

# -------------------- long-format occurrence table --------------------
# app.R writes the converted table; this script re-reads it so that the run is
# fully reproducible from config.json + the CSV alone.
long_csv <- cfg$ssdm_long_file %||% file.path(cfg$out_dir, "ssdm_long_format.csv")
if (!file.exists(long_csv)) stop("Converted S-SDM table not found: ", long_csv)

occ <- read_user_csv(long_csv)
miss <- setdiff(c("species", "presence"), names(occ))
if (length(miss) > 0) {
  stop("Converted table is missing column(s): ", paste(miss, collapse = ", "))
}
xy0 <- get_xy_matrix(occ, cfg$x_col, cfg$y_col, what = "Converted S-SDM table")

occ$species  <- as.character(occ$species)
occ$x        <- xy0[, 1]
occ$y        <- xy0[, 2]
occ$presence <- suppressWarnings(as.integer(occ$presence))
occ <- occ[is.finite(occ$x) & is.finite(occ$y) & occ$presence %in% c(0L, 1L), , drop = FALSE]
if (nrow(occ) == 0) stop("No usable rows in the converted table.")

# Drop plots whose environment is NA
ex <- terra::extract(env, as.matrix(occ[, c("x", "y")]))
if ("ID" %in% names(ex)) ex$ID <- NULL
keep <- stats::complete.cases(ex)
if (sum(keep) < nrow(occ)) {
  message("Rows dropped for NA predictors: ", nrow(occ) - sum(keep))
}
occ <- occ[keep, , drop = FALSE]

sp_all <- sort(unique(occ$species))
message("Species in table: ", length(sp_all), " | rows: ", nrow(occ))

# -------------------- species screening --------------------
# The app already applies this filter before writing the table; it is repeated
# here so that running the script directly on an unfiltered CSV is still safe.
n_plots_total <- length(unique(occ$plot %||% occ$x))
if ("plot" %in% names(occ)) n_plots_total <- length(unique(occ$plot))

occ_tab <- do.call(rbind, lapply(sp_all, function(s) {
  d <- occ[occ$species == s, , drop = FALSE]
  data.frame(species = s,
             n_presence = sum(d$presence == 1L),
             n_absence  = sum(d$presence == 0L),
             stringsAsFactors = FALSE)
}))
occ_tab$frequency_pct <- 100 * occ_tab$n_presence / n_plots_total

if (MIN_FREQ_PCT <= 0) {
  message("Species filter: OFF (0%) - every species enters modelling.")
  ok_sp <- occ_tab$species[occ_tab$n_presence >= 1]
} else {
  freq_n <- ceiling(MIN_FREQ_PCT / 100 * n_plots_total)
  message("Species filter: frequency >= ", MIN_FREQ_PCT, "% of ", n_plots_total,
          " plots (= ", freq_n, " plots)")
  ok_sp <- occ_tab$species[occ_tab$frequency_pct >= MIN_FREQ_PCT & occ_tab$n_presence >= 1]
}
skip_sp <- setdiff(sp_all, ok_sp)

if (length(skip_sp) > 0) {
  sk <- occ_tab[occ_tab$species %in% skip_sp, , drop = FALSE]
  sk$reason <- ifelse(sk$n_presence < 1, "no presence records",
                      sprintf("frequency < %g%%", MIN_FREQ_PCT))
  safe_write_csv(sk, file.path(cfg$out_dir, "skipped_species.csv"))
  message("Skipped species: ", length(skip_sp))
}

# Yokluk stratejisi ozeti (asagida tur bazinda uygulanir)
message("Absence strategy: real absences are used when a species has at least ",
        CV_K, " of them; otherwise pseudo-absences are generated.")
if (length(ok_sp) < 2) {
  stop("At least 2 modellable species are required for stacking (found ", length(ok_sp), ").")
}

# -------------------- model pool --------------------
models_to_run <- as_chr_vec(cfg$models)
models_to_run <- intersect(models_to_run, ALLOWED_MODELS_BINARY)
req_pkgs <- list(ANN = "nnet", CTA = "rpart", DNN = "cito", FDA = "mda", GBM = "gbm",
                 MARS = "earth", GAM = "mgcv", MAXNET = "maxnet", XGBOOST = "xgboost",
                 RF = "randomForest", RFd = "randomForest")
for (m in intersect(names(req_pkgs), models_to_run)) {
  if (!requireNamespace(req_pkgs[[m]], quietly = TRUE)) {
    warning("Package missing ('", req_pkgs[[m]], "') -> dropping model: ", m)
    models_to_run <- setdiff(models_to_run, m)
  }
}
if (any(grepl("^MAXENT$", models_to_run))) {
  if (is.null(cfg$maxent_jar) || !file.exists(cfg$maxent_jar)) {
    warning("MAXENT selected but maxent.jar not found -> dropping MAXENT.")
    models_to_run <- setdiff(models_to_run, "MAXENT")
  } else {
    bmopt_stage_maxent(cfg$maxent_jar, verbose = TRUE)
  }
}
if (length(models_to_run) == 0) stop("No models left to run.")
message("Algorithms: ", paste(models_to_run, collapse = ", "))

cv_rep <- as_int1(cfg$cv$rep, 1L); if (is.na(cv_rep) || cv_rep < 1) cv_rep <- 1L
nb_cpu <- as_int1(cfg$cpu_n, 1L);  if (is.na(nb_cpu) || nb_cpu < 1) nb_cpu <- 1L

# -------------------- future scenarios --------------------
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
  if (length(miss) > 0) stop("Missing predictors in future layers: ", paste(miss, collapse = ", "))
  e[[ref_names]]
}

scenarios <- list()
if (!is.null(cfg$future_n) && !is.na(as_int1(cfg$future_n, NA_integer_)) &&
    as_int1(cfg$future_n, 0L) >= 1 && !is.null(cfg$future)) {
  fl <- if (is.list(cfg$future)) cfg$future else list(cfg$future)
  nf <- min(as_int1(cfg$future_n, 0L), length(fl))
  for (i in seq_len(nf)) {
    fb <- fl[[i]]; if (!is.list(fb) || is.null(fb$scenarios)) next
    scs <- if (is.list(fb$scenarios)) fb$scenarios else list(fb$scenarios)
    for (j in seq_along(scs)) {
      sc <- scs[[j]]; if (!is.list(sc)) next
      d <- sc$asc_dir %||% ""; if (!nzchar(d)) next
      scenarios[[length(scenarios) + 1]] <- list(
        tag = make.names(paste0(fb$period_label %||% paste0("Future", i), "_",
                                sc$label %||% paste0("Scenario", j))),
        dir = d
      )
    }
  }
}
message("Future scenarios: ", length(scenarios))

# -------------------- per-species modelling --------------------
# Accumulators: binary sum (bSSDM) and continuous sum (pSSDM)
acc <- new.env(parent = emptyenv())
acc$bin <- list()      # tag -> SpatRaster (running sum)
acc$prob <- list()
acc$n_sp <- 0L

add_layer <- function(store, tag, r) {
  if (is.null(acc[[store]][[tag]])) acc[[store]][[tag]] <- r
  else acc[[store]][[tag]] <- acc[[store]][[tag]] + r
  invisible(NULL)
}

pick_em_layer <- function(pred, pattern = "EMwmean|EMmean") {
  idx <- grep(pattern, names(pred))
  if (length(idx) == 0) idx <- 1L
  idx[1]
}

get_maxtss_cutoff <- function(ens, layer_name) {
  ev <- tryCatch(get_evaluations(ens), error = function(e) NULL)
  if (is.null(ev) || !is.data.frame(ev)) return(NA_real_)
  nm <- names(ev)
  mcol <- first_present(nm, c("full.name", "model.name", "model"))
  ecol <- first_present(nm, c("metric.eval", "metric"))
  if (is.null(mcol) || is.null(ecol) || !("cutoff" %in% nm)) return(NA_real_)
  sub <- ev[ev[[ecol]] == "TSS", , drop = FALSE]
  if (nrow(sub) == 0) return(NA_real_)
  hit <- which(as.character(sub[[mcol]]) == layer_name)
  if (length(hit) == 0) hit <- grep("EMwmean|EMmean", as.character(sub[[mcol]]))
  if (length(hit) == 0) hit <- seq_len(nrow(sub))
  v <- suppressWarnings(as.numeric(sub$cutoff[hit]))
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else stats::median(v)
}

summary_rows <- list()
n_ok <- 0L

# biomod2 resp.name'i KLASOR ADI olarak kullanir. make.names() tek basina
# yetmez: "Pinus nigra" ile "Pinus-nigra" ayni ada duser ve iki tur ayni
# klasore yazip birbirini ezer. make.unique() ile benzersizlestiriyoruz.
sp_ids <- make.unique(make.names(ok_sp))
if (any(make.names(ok_sp) != sp_ids)) {
  clash <- ok_sp[make.names(ok_sp) != sp_ids]
  warning("Species names collapsing to the same folder name were disambiguated: ",
          paste(clash, collapse = ", "))
}

for (si in seq_along(ok_sp)) {
  sp <- ok_sp[si]
  sp_id <- sp_ids[si]
  prefix <- paste0("SPECIES ", si, "/", length(ok_sp), ": ", sp)
  set_status(paste0(prefix, " - FORMATTING"))
  message("\n===== ", prefix, " =====")

  d <- occ[occ$species == sp, , drop = FALSE]
  resp <- as.integer(d$presence)
  xy   <- as.matrix(d[, c("x", "y")]); storage.mode(xy) <- "double"

  n_pres_sp <- sum(resp == 1L); n_abs_sp <- sum(resp == 0L)
  # Gercek yokluk CV kat sayisindan azsa 10 katli CV'de katlarin cogu yokluksuz
  # kalir ve TSS/AUC hesaplanamaz -> sozde-yokluk ile takviye edilir.
  need_pa_sp <- n_abs_sp < CV_K
  abs_strategy <- if (!need_pa_sp) "real"
                  else if (n_abs_sp > 0) "real+PA" else "PA"

  row <- data.frame(species = sp, model_id = sp_id,
                    n_presence = n_pres_sp, n_absence = n_abs_sp,
                    frequency_pct = round(100 * n_pres_sp / n_plots_total, 3),
                    absence_strategy = abs_strategy,
                    status = "ok", models_in_ensemble = NA_character_,
                    ensemble_AUCroc = NA_real_, ensemble_TSS = NA_real_,
                    maxTSS_cutoff_0_1000 = NA_real_, stacked = FALSE,
                    stringsAsFactors = FALSE)

  res <- try({
    fa <- list(resp.name = sp_id, resp.var = resp, resp.xy = xy,
               expl.var = env, data.type = "binary", na.rm = TRUE)
    if (need_pa_sp) {
      fa$PA.nb.rep      <- as_int1(cfg$pa_rep, 3L)
      fa$PA.nb.absences <- as_int1(cfg$pa_n, 10000L)
      fa$PA.strategy    <- as.character(cfg$pa_strategy %||% "random")[1]
      if (identical(fa$PA.strategy, "disk")) {
        fa$PA.dist.min <- as_num1(cfg$pa_dist_min, 0)
        fa$PA.dist.max <- as_num1(cfg$pa_dist_max, 0)
      }
      message("  absence strategy: ", abs_strategy,
              " (", n_abs_sp, " real absences < CV.k=", CV_K, ")")
    }
    bmf <- tryCatch(do.call(BIOMOD_FormatingData, fa), error = function(e) e)
    if (inherits(bmf, "error") && need_pa_sp && n_abs_sp > 0) {
      # Bazi surumler gercek yokluk + PA karisimini reddeder -> varlik + PA
      fa$resp.var <- rep(1, n_pres_sp)
      fa$resp.xy  <- xy[resp == 1L, , drop = FALSE]
      bmf <- do.call(BIOMOD_FormatingData, fa)
      abs_strategy <- "PA"
    }
    if (inherits(bmf, "error")) stop(conditionMessage(bmf))

    set_status(paste0(prefix, " - MODELING"))
    built <- bmopt_build(
      data_type = "binary", models = models_to_run, bm_format = bmf,
      p = p_env, n = sum(resp == 1L), maxent_jar = cfg$maxent_jar,
      strategy = cfg$opt_strategy %||% "adaptive", transfer = TRUE, verbose = FALSE
    )

    margs <- list(
      bm.format = bmf, modeling.id = sp_id, models = models_to_run,
      CV.strategy = "kfold", CV.k = 10, CV.nb.rep = cv_rep,
      metric.eval = c("AUCroc", "TSS"), var.import = 3,
      prevalence = 0.5, seed.val = 42, nb.cpu = nb_cpu, do.progress = FALSE
    )
    margs <- bmopt_inject(margs, built)

    mo <- tryCatch(do.call(BIOMOD_Modeling, margs), error = function(e) e)
    if (inherits(mo, "error")) {
      margs$OPT.user.val <- NULL; margs$OPT.user.base <- NULL; margs$OPT.strategy <- "bigboss"
      mo <- do.call(BIOMOD_Modeling, margs)
    }

    # ---- model selection on AUC + TSS ----
    ev <- get_evaluations(mo)
    ds <- if ("validation" %in% names(ev) && any(is.finite(ev$validation))) "validation" else "calibration"
    ag <- stats::aggregate(ev[[ds]], by = list(algo = ev$algo, metric = ev$metric.eval),
                           FUN = mean, na.rm = TRUE)
    names(ag)[3] <- "value"
    auc <- ag$value[ag$metric == "AUCroc"]; names(auc) <- ag$algo[ag$metric == "AUCroc"]
    tss <- ag$value[ag$metric == "TSS"];    names(tss) <- ag$algo[ag$metric == "TSS"]
    cand <- intersect(names(auc), names(tss))
    keep_algo <- cand[is.finite(auc[cand]) & is.finite(tss[cand]) &
                        auc[cand] >= ENS_MIN_AUC & tss[cand] >= ENS_MIN_TSS &
                        auc[cand] < ENS_NEAR_ONE & tss[cand] < ENS_NEAR_ONE]
    if (length(keep_algo) == 0) stop("no algorithm passed the AUC/TSS thresholds")

    bf <- as.character(getModelsBuilt(mo))
    keep_full <- bf[sub(".*_", "", bf) %in% keep_algo]
    row$models_in_ensemble <- paste(keep_algo, collapse = "|")

    set_status(paste0(prefix, " - ENSEMBLE"))
    ens <- BIOMOD_EnsembleModeling(
      bm.mod = mo, models.chosen = if (length(keep_full)) keep_full else "all",
      em.by = "all", em.algo = c("EMmean", "EMca", "EMwmean", "EMcv"),
      metric.select = c("AUCroc", "TSS"),
      metric.select.thresh = c(ENS_MIN_AUC, ENS_MIN_TSS),
      metric.eval = c("AUCroc", "TSS"), var.import = 0
    )

    eev <- get_evaluations(ens)
    if (is.data.frame(eev)) {
      dsx <- if ("validation" %in% names(eev) && any(is.finite(eev$validation))) "validation" else "calibration"
      row$ensemble_AUCroc <- suppressWarnings(mean(eev[[dsx]][eev$metric.eval == "AUCroc"], na.rm = TRUE))
      row$ensemble_TSS    <- suppressWarnings(mean(eev[[dsx]][eev$metric.eval == "TSS"],    na.rm = TRUE))
    }
    if (!is.null(eev)) safe_write_csv(eev, file.path(sp_out_dir, paste0(sp_id, "_eval_ensemble.csv")))

    # ---- current projection ----
    set_status(paste0(prefix, " - PROJECTION (current)"))
    pj <- BIOMOD_Projection(bm.mod = mo, new.env = env, proj.name = paste0(sp_id, "_cur"),
                            models.chosen = "all", build.clamping.mask = FALSE)
    ec <- BIOMOD_EnsembleForecasting(bm.em = ens, bm.proj = pj,
                                     proj.name = paste0(sp_id, "_cur_ens"))
    pc <- get_predictions(ec)
    if (!inherits(pc, "SpatRaster")) stop("current ensemble raster unavailable")

    li <- pick_em_layer(pc)
    lyr_name <- names(pc)[li]
    ct <- get_maxtss_cutoff(ens, lyr_name)
    if (!is.finite(ct)) stop("MaxTSS cut-off unavailable")
    row$maxTSS_cutoff_0_1000 <- ct

    bin_cur <- terra::classify(pc[[li]] >= ct, cbind(NA, 0))
    terra::writeRaster(bin_cur,
                       file.path(sp_ras_dir, paste0(sp_id, "_current_binary.tif")),
                       overwrite = TRUE, datatype = "INT1U")

    add_layer("bin", "current", bin_cur)
    if (DO_PSSDM) add_layer("prob", "current", pc[[li]] / SCALE_DIV)

    # ---- future projections ----
    for (sc in scenarios) {
      set_status(paste0(prefix, " - PROJECTION (", sc$tag, ")"))
      ef <- try({
        envf <- load_env_future(sc$dir, names(env))
        pjf <- BIOMOD_Projection(bm.mod = mo, new.env = envf,
                                 proj.name = paste0(sp_id, "_", sc$tag),
                                 models.chosen = "all", build.clamping.mask = FALSE)
        BIOMOD_EnsembleForecasting(bm.em = ens, bm.proj = pjf,
                                   proj.name = paste0(sp_id, "_", sc$tag, "_ens"))
      }, silent = TRUE)
      if (inherits(ef, "try-error")) {
        warning("Future projection failed for ", sp, " / ", sc$tag); next
      }
      pf <- get_predictions(ef)
      if (!inherits(pf, "SpatRaster")) next
      lif <- pick_em_layer(pf)
      # The SAME MaxTSS cut-off is reused across scenarios so that richness
      # differences reflect habitat change, not a moving threshold.
      bin_f <- terra::classify(pf[[lif]] >= ct, cbind(NA, 0))
      add_layer("bin", sc$tag, bin_f)
      if (DO_PSSDM) add_layer("prob", sc$tag, pf[[lif]] / SCALE_DIV)
    }

    row$stacked <- TRUE
    TRUE
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    row$status <- paste0("failed: ", trimws(conditionMessage(attr(res, "condition"))))
    warning("Species failed: ", sp, " - ", row$status)
  } else {
    n_ok <- n_ok + 1L
  }
  summary_rows[[length(summary_rows) + 1]] <- row
  safe_write_csv(do.call(rbind, summary_rows), file.path(cfg$out_dir, "species_summary.csv"))
}

if (n_ok < 2) stop("Fewer than 2 species were modelled successfully; cannot stack.")
message("\nSpecies successfully stacked: ", n_ok, " / ", length(ok_sp))

# -------------------- stacking --------------------
set_status("STACKING RICHNESS MAPS")

stats_rows <- list()
richness_stats <- function(r, label) {
  v <- terra::values(r, mat = FALSE); v <- v[is.finite(v)]
  if (length(v) == 0) return(NULL)
  data.frame(map = label, n_cells = length(v),
             mean_richness = mean(v), sd_richness = stats::sd(v),
             min_richness = min(v), max_richness = max(v),
             stringsAsFactors = FALSE)
}

cur_bin <- acc$bin[["current"]]
if (is.null(cur_bin)) stop("No current binary stack was produced.")

terra::writeRaster(cur_bin, file.path(rasters_dir, "richness_current_bSSDM.tif"),
                   overwrite = TRUE, datatype = "INT2U")
stats_rows[[length(stats_rows) + 1]] <- richness_stats(cur_bin, "richness_current_bSSDM")

if (DO_PSSDM && !is.null(acc$prob[["current"]])) {
  terra::writeRaster(acc$prob[["current"]],
                     file.path(rasters_dir, "richness_current_pSSDM.tif"), overwrite = TRUE)
  stats_rows[[length(stats_rows) + 1]] <- richness_stats(acc$prob[["current"]], "richness_current_pSSDM")
}

for (sc in scenarios) {
  rb <- acc$bin[[sc$tag]]
  if (is.null(rb)) next
  terra::writeRaster(rb, file.path(rasters_dir, paste0("richness_", sc$tag, "_bSSDM.tif")),
                     overwrite = TRUE, datatype = "INT2U")
  stats_rows[[length(stats_rows) + 1]] <- richness_stats(rb, paste0("richness_", sc$tag, "_bSSDM"))

  if (DO_CHANGE) {
    ch <- rb - cur_bin
    terra::writeRaster(ch, file.path(rasters_dir, paste0("richness_change_", sc$tag, "_bSSDM.tif")),
                       overwrite = TRUE, datatype = "INT2S")
    stats_rows[[length(stats_rows) + 1]] <- richness_stats(ch, paste0("richness_change_", sc$tag, "_bSSDM"))
  }
  if (DO_PSSDM && !is.null(acc$prob[[sc$tag]])) {
    terra::writeRaster(acc$prob[[sc$tag]],
                       file.path(rasters_dir, paste0("richness_", sc$tag, "_pSSDM.tif")),
                       overwrite = TRUE)
  }
}

if (length(stats_rows) > 0) {
  safe_write_csv(do.call(rbind, stats_rows), file.path(cfg$out_dir, "ssdm_richness_stats.csv"))
}

# -------------------- plots (English labels by design) --------------------
safe_png(file.path(plots_dir, "richness_current_bSSDM.png"), {
  terra::plot(cur_bin, main = "Potential species richness (current, bSSDM)")
})
for (sc in scenarios) {
  rb <- acc$bin[[sc$tag]]; if (is.null(rb)) next
  safe_png(file.path(plots_dir, paste0("richness_", sc$tag, "_bSSDM.png")), {
    terra::plot(rb, main = paste0("Potential species richness (", sc$tag, ", bSSDM)"))
  })
  if (DO_CHANGE) {
    safe_png(file.path(plots_dir, paste0("richness_change_", sc$tag, "_bSSDM.png")), {
      terra::plot(rb - cur_bin,
                  main = paste0("Richness change (", sc$tag, " - current)"))
    })
  }
}

set_status("DONE")
message("S-SDM finished. Richness maps written to: ", rasters_dir)
