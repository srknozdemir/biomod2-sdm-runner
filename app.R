# app.R  (PATH PICKERS + CHECK BUTTONS; variable types UI) -- UPDATED
# YENI (bu surum): BAGIMLI DEGISKEN TIPI SECIMI
# - En ustte hedef tur verisi tipi secilir: Var/Yok | Surekli | Nominal
# - Secime gore TUM alt konfigurasyon degisir:
#     * Var/Yok  -> data.type = "binary"      -> run_biomod.R      (PA, AUC/TSS/BOYCE, MaxTSS)
#     * Surekli  -> count/abundance/relative  -> run_biomod_cont.R   (regresyon; RMSE/Rsquared)
#     * Nominal  -> multiclass/ordinal        -> run_biomod_cont.R   (siniflandirma; Accuracy/F1)
# - Surekli/nominal icin CSV'den YANIT SUTUNU secilir (PA uretimi yoktur)
# - Model listesi, metrikler, ensemble algoritmalari ve esikler veri tipine gore filtrelenir
# - Optimizasyon: bm_ModelingOptions(data.type=..., strategy=...) + OPT.user
#
# ONCEKI YENILIKLER:
# - shinyFiles ile tikla-sec dosya/klasor secicileri (textInput'lar korundu, elle yapistirma da calisir)
# - CPU cekirdek sayisi (cfg$cpu_n) -> biomod2 paralel calisir (onceden hep 1 idi)
# - BOYCE metrigi metrics tablosunda ve grafik listesinde
# - MaxTSS cutoff tablosu (0-1000 olcek kontrolu icin)
# - Windows/Linux uyumlu arka plan calistirma
# - RUN oncesi on kosul denetimi
# - MAXNET / XGBOOST model listesine eklendi
# - Manual path inputs (paste)
# - Check buttons for CSV / current / future / maxent.jar
# - Auto species detection from CSV
# - Current raster check reads stack + lists variable names
# - Variable types UI after current check (continuous/categorical/nominal)
# - Apply variable types saves mapping into cfg$var_types
# - RUN writes robust config.json (models forced to character; var_types forced list)
# - Keeps run_biomod.R launch behavior (background Rscript)
# - FIX: Shiny image serving + cache busting; dynamic Image Viewer choices
# - FIX (NEW): Raster viewer refresh-safe: dropdown doesn't reset; selection persists; auto-updates when new tif arrives

library(shiny)
library(jsonlite)
library(terra)

# Dosya/klasor secici (tikla-sec). Yoksa kurulum uyarisi verir.
if (!requireNamespace("shinyFiles", quietly = TRUE)) {
  stop("shinyFiles paketi gerekli. Kurulum: install.packages('shinyFiles')")
}
library(shinyFiles)

# FIX: vector-safe null-coalesce
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (is.atomic(a)) {
    if (!any(!is.na(a))) return(b)
  }
  a
}

# =========================================================
# BASLANGIC SISTEM DENETIMI
# Arka plandaki Rscript sessizce olurse tanisi zor olur; bu yuzden
# gerekli dosya/paketleri UYGULAMA ACILISINDA kontrol ediyoruz.
# =========================================================
REQUIRED_FILES <- c("run_biomod.R", "run_biomod_cont.R", "run_ssdm.R",
                    "biomod_opts.R", "plot_response_threshold.R",
                    "i18n.R", "ssdm_convert.R")

for (.f in c("i18n.R", "ssdm_convert.R", "biomod_opts.R")) {
  if (file.exists(file.path(getwd(), .f))) source(file.path(getwd(), .f), encoding = "UTF-8")
}
if (!exists("tr_")) tr_ <- function(key, lang = "en", ...) key

MODEL_PKGS <- c(
  ANN = "nnet", CTA = "rpart", DNN = "cito", FDA = "mda", GAM = "mgcv",
  GBM = "gbm", MARS = "earth", MAXNET = "maxnet", RF = "randomForest",
  RFd = "randomForest", XGBOOST = "xgboost"
)

setup_status <- function() {
  files_missing <- REQUIRED_FILES[!file.exists(file.path(getwd(), REQUIRED_FILES))]

  core <- c("shiny", "jsonlite", "terra", "shinyFiles", "biomod2")
  core_missing <- core[!vapply(core, function(p) requireNamespace(p, quietly = TRUE), logical(1))]

  bm_ver <- tryCatch(as.character(utils::packageVersion("biomod2")), error = function(e) NA_character_)
  bm_ok43 <- !is.na(bm_ver) && utils::compareVersion(bm_ver, "4.3") >= 0

  mp <- unique(MODEL_PKGS)
  mp_missing <- mp[!vapply(mp, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
  models_blocked <- names(MODEL_PKGS)[MODEL_PKGS %in% mp_missing]

  java_ok <- nzchar(Sys.which("java"))

  # Arka plan islerini baslatan Rscript bulunabiliyor mu?
  # Windows'ta R'in bin klasoru cogu kurulumda PATH'e eklenmez; bu durumda
  # is sessizce baslamaz. R.home() her zaman dogru sonuc verir.
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  rs_cand <- c(file.path(R.home("bin"), exe),
               file.path(R.home(), "bin", exe),
               file.path(R.home(), "bin", "x64", exe),
               unname(Sys.which("Rscript")))
  rs_cand <- rs_cand[nzchar(rs_cand)]
  rs_hit  <- rs_cand[file.exists(rs_cand)]
  rscript <- if (length(rs_hit) > 0) rs_hit[1] else NA_character_

  list(
    wd = getwd(),
    files_missing = files_missing,
    core_missing = core_missing,
    bm_ver = bm_ver,
    bm_ok43 = bm_ok43,
    mp_missing = mp_missing,
    models_blocked = models_blocked,
    java_ok = java_ok,
    rscript = rscript
  )
}

# =========================================================
# DAYANIKLI DOSYA OKUMA  (Windows dosya kilidi icin)
#
# Windows'ta arka plandaki Rscript run.log'a yazarken dosyayi kilitleyebilir;
# bu sirada readLines() "cannot open the connection" hatasi verir ve Shiny
# ciktisi coker. Linux/macOS eszamanli okumaya izin verdigi icin bu hata
# yalnizca Windows'ta gorulur.
#
# Cozum: once dogrudan oku; kilitliyse gecici bir kopya alip kopyayi oku;
# o da olmazsa bos don. Hicbir kosulda hata FIRLATMAZ.
# =========================================================
safe_read_lines <- function(path, n = Inf) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(character(0))

  out <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"),
                  error = function(e) NULL, warning = function(w) NULL)

  if (is.null(out)) {
    out <- tryCatch({
      tmp <- tempfile(fileext = ".txt")
      on.exit(unlink(tmp), add = TRUE)
      if (isTRUE(file.copy(path, tmp, overwrite = TRUE))) {
        readLines(tmp, warn = FALSE, encoding = "UTF-8")
      } else NULL
    }, error = function(e) NULL, warning = function(w) NULL)
  }

  if (is.null(out)) return(character(0))
  if (is.finite(n) && length(out) > n) out <- utils::tail(out, n)
  out
}

# CSV okuma: kilit + kodlama korumali.
# Kullanici dosyalari Windows'ta genellikle UTF-8 DEGIL, Windows-1254
# (Turkce) olarak kaydedilir; bu durumda tur adlari bozuk gorunur.
safe_read_csv <- function(path, ...) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)

  attempt <- function(p, enc) {
    tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, fileEncoding = enc, ...),
             error = function(e) NULL, warning = function(w) NULL)
  }
  encs <- c("UTF-8", "", "windows-1254", "latin1")

  for (e in encs) {
    d <- attempt(path, e)
    if (is.data.frame(d) && nrow(d) > 0) return(d)
  }
  # Kilitliyse kopya uzerinden dene
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  if (isTRUE(tryCatch(file.copy(path, tmp, overwrite = TRUE), error = function(e) FALSE))) {
    for (e in encs) {
      d <- attempt(tmp, e)
      if (is.data.frame(d) && nrow(d) > 0) return(d)
    }
  }
  NULL
}

detect_species_name <- function(csv_path) {
  df <- safe_read_csv(csv_path)
  if (is.null(df) || nrow(df) == 0) return("")

  cand <- c(
    "Tur","tur","species","Species","sp","Sp","taxon","Taxon",
    "scientificName","ScientificName","name","Name"
  )
  col <- cand[cand %in% names(df)]
  if (length(col) == 0) return("")

  v <- as.character(df[[col[1]]])
  v <- v[!is.na(v)]
  v <- v[nzchar(v)]
  if (length(v) == 0) return("")

  tab <- sort(table(v), decreasing = TRUE)
  nm <- names(tab)[1]
  if (is.null(nm)) "" else nm
}

# =========================================================
# BAGIMLI DEGISKEN (data.type) YARDIMCILARI
# =========================================================

# biomod2 >= 4.3 model uygunluk matrisi
ALLOWED_MODELS <- list(
  binary     = c("ANN","CTA","DNN","FDA","GAM","GBM","GLM","MARS",
                 "MAXENT","MAXNET","RF","RFd","SRE","XGBOOST"),
  count      = c("CTA","DNN","GAM","GBM","GLM","MARS","RF","XGBOOST"),
  abundance  = c("CTA","DNN","GAM","GBM","GLM","MARS","RF","XGBOOST"),
  relative   = c("CTA","DNN","GAM","GBM","GLM","MARS","RF","XGBOOST"),
  ordinal    = c("CTA","DNN","FDA","GAM","GLM","MARS","RF","XGBOOST"),
  multiclass = c("CTA","DNN","FDA","MARS","RF","XGBOOST")
)

DEFAULT_MODELS <- list(
  binary     = c("GLM","GAM","GBM","RF","MAXNET"),
  count      = c("GLM","GAM","GBM","RF","XGBOOST"),
  abundance  = c("GLM","GAM","GBM","RF","XGBOOST"),
  relative   = c("GLM","GAM","GBM","RF","XGBOOST"),
  ordinal    = c("CTA","GLM","MARS","RF","XGBOOST"),
  multiclass = c("CTA","FDA","MARS","RF","XGBOOST")
)

EM_ALGOS <- list(
  binary = c("EMmean","EMca","EMwmean","EMcv"),
  quant  = c("EMmean","EMmedian","EMwmean","EMcv","EMci"),
  qual   = c("EMmode","EMfreq")
)

DT_QUANT <- c("count","abundance","relative")
DT_QUAL  <- c("multiclass","ordinal")

dt_family <- function(dt) {
  if (identical(dt, "binary")) "binary" else if (dt %in% DT_QUANT) "quant" else "qual"
}

metrics_for <- function(dt) {
  switch(dt_family(dt),
         binary = c("AUCroc","TSS","BOYCE"),
         quant  = c("Rsquared","Rsquared_aj","RMSE","MAE"),
         qual   = c("Accuracy","F1","Recall","Precision"))
}

# Surekli alt tipi tahmini: 0-1 -> relative, negatif olmayan tam sayi -> count, diger -> abundance
guess_continuous_type <- function(v) {
  x <- suppressWarnings(as.numeric(v))
  x <- x[is.finite(x)]
  if (length(x) == 0) return("abundance")
  if (all(x >= 0) && all(x <= 1) && length(unique(x)) > 2) return("relative")
  if (all(x >= 0) && all(abs(x - round(x)) < 1e-8)) return("count")
  "abundance"
}

# Nominal alt tipi tahmini: sayisal/siralanabilir -> ordinal, metinsel -> multiclass
guess_nominal_type <- function(v) {
  x <- v[!is.na(v)]
  if (length(x) == 0) return("multiclass")
  num <- suppressWarnings(as.numeric(as.character(x)))
  if (!any(is.na(num))) return("ordinal")
  "multiclass"
}

# CSV sutunlarini ve aday yanit sutunlarini okur
read_csv_head <- function(csv_path) {
  df <- safe_read_csv(csv_path)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

resp_candidates <- function(df) {
  if (is.null(df)) return(character(0))
  drop <- c("x","y","X","Y","lon","lat","longitude","latitude",
            "Tur","tur","species","Species","sp","Sp","taxon","Taxon",
            "scientificName","ScientificName","name","Name","id","ID")
  setdiff(names(df), drop)
}

# Yanit sutunu ozeti (kullaniciya gosterilir)
resp_summary_text <- function(v) {
  x <- v[!is.na(v)]
  if (length(x) == 0) return("Sutun bos veya tamami NA.")
  num <- suppressWarnings(as.numeric(as.character(x)))
  if (!any(is.na(num))) {
    u <- unique(num)
    paste0(
      "n = ", length(num),
      " | min = ", signif(min(num), 4),
      " | ortalama = ", signif(mean(num), 4),
      " | max = ", signif(max(num), 4),
      " | benzersiz = ", length(u),
      " | sifir orani = ", round(mean(num == 0), 3),
      " | tam sayi: ", if (all(abs(num - round(num)) < 1e-8)) "evet" else "hayir"
    )
  } else {
    tb <- sort(table(as.character(x)), decreasing = TRUE)
    paste0("Metinsel/kategorik | n = ", length(x), " | sinif sayisi = ", length(tb),
           " | siniflar: ", paste(utils::head(names(tb), 12), collapse = ", "))
  }
}

# Read env rasters from folder; supports asc/tif/grd.
# IMPORTANT: if .grd exists, also include matching .gri (more robust for terra).
read_env_stack <- function(dir_path) {
  if (is.null(dir_path) || !nzchar(dir_path)) stop("Empty folder path.")
  if (!dir.exists(dir_path)) stop("Folder not found: ", dir_path)

  exts <- c("asc", "tif", "tiff", "grd")
  files <- list.files(
    dir_path,
    pattern = paste0("\\.(", paste(exts, collapse = "|"), ")$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    stop("No raster files found in folder. Expected: .asc/.tif/.grd")
  }

  # If any .grd, append its .gri (terra often expects both)
  grd <- files[grepl("\\.grd$", files, ignore.case = TRUE)]
  if (length(grd) > 0) {
    gri <- sub("\\.grd$", ".gri", grd, ignore.case = TRUE)
    gri <- gri[file.exists(gri)]
    files <- unique(c(files, gri))
  }

  r <- try(terra::rast(files), silent = TRUE)
  if (inherits(r, "try-error")) stop("terra::rast() failed to read rasters in folder.")
  names(r) <- make.names(names(r), unique = TRUE)
  r
}

ui <- fluidPage(
  tags$head(tags$style(
    HTML(
      "
      .progress { height: 22px; }
      .progress-bar { font-weight: 700; }
      .smallnote { font-size: 12px; color: #555; }
      .boxpad { padding: 10px; border: 1px solid #eee; border-radius: 8px; margin-bottom: 10px; }
      .ok { color: #0a7a0a; font-weight: 700; }
      .bad { color: #b00020; font-weight: 700; }
      code { font-size: 12px; }
      .vartype-row { display:flex; gap:10px; align-items:center; margin:6px 0; }
      .vartype-name { min-width: 160px; font-family: ui-monospace, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; font-size: 12px; }
      .vartype-select { flex: 1; }
      .dtbox { border: 2px solid #2c6fbb; background: #f4f8fd; }
      .dtbadge { display:inline-block; padding:2px 8px; border-radius:10px;
                 background:#2c6fbb; color:#fff; font-size:11px; font-weight:700; }
      .langbox { background:#fafafa; }
      .ssdmbox { border: 2px solid #6a3d9a; background: #f7f3fb; }
      .previewtbl table { font-size: 11px; }
      .previewtbl { max-height: 260px; overflow: auto; margin-bottom: 6px; }
    "
    )
  )),
  titlePanel("biomod2 SDM Runner (Manual Paths + Check Buttons + Variable Types)"),

  sidebarLayout(
    sidebarPanel(

      tags$div(
        class = "boxpad langbox",
        selectInput("lang", "Language / Dil",
                    choices = c("English" = "en", "T\u00fcrk\u00e7e" = "tr"),
                    selected = "en", width = "100%")
      ),

      tags$div(class = "boxpad", uiOutput("setup_ui")),

      # =====================================================
      # S-SDM (opsiyonel) - tek tur akisini devre disi birakir
      # =====================================================
      tags$div(
        class = "boxpad ssdmbox",
        uiOutput("ssdm_hdr"),
        checkboxInput("ssdm_enable", "Enable stacked species distribution modelling",
                      value = FALSE),
        uiOutput("ssdm_note_ui"),

        conditionalPanel(
          condition = "input.ssdm_enable == true",

          uiOutput("ssdm_mat_hdr_ui"),
          textInput("ssdm_mat_file", "Species matrix CSV (full path)", value = ""),
          shinyFilesButton("ssdm_mat_btn", "CSV", "Select species matrix", multiple = FALSE),

          uiOutput("ssdm_xy_hdr_ui"),
          textInput("ssdm_xy_file", "Coordinate CSV (full path)", value = ""),
          shinyFilesButton("ssdm_xy_btn", "CSV", "Select coordinate file", multiple = FALSE),

          tags$br(), tags$br(),
          actionButton("ssdm_check", "Check both files"),
          uiOutput("ssdm_check_ui"),

          uiOutput("ssdm_xcol_ui"),
          uiOutput("ssdm_ycol_ui"),
          checkboxInput("ssdm_latlon",
                        "Coordinates are geographic degrees (EPSG:4326)", value = FALSE),

          tags$br(),
          actionButton("ssdm_convert", "Convert to long format", class = "btn-primary"),
          uiOutput("ssdm_convert_ui"),

          uiOutput("ssdm_filter_hdr_ui"),
          numericInput("ssdm_min_freq", "Minimum frequency (% of plots)",
                       value = 0, min = 0, max = 100, step = 0.5),
          uiOutput("ssdm_filter_note_ui"),
          uiOutput("ssdm_filter_ui"),
          tags$hr(),
          checkboxInput("ssdm_pssdm", "Also produce pSSDM (sum of suitabilities)", value = TRUE),
          checkboxInput("ssdm_change", "Produce richness change maps (future - current)",
                        value = TRUE)
        )
      ),

      # =====================================================
      # 0) HEDEF TUR / BAGIMLI DEGISKEN TIPI  (ILK ADIM)
      # S-SDM acikken gizlenir: orada yanit matristen gelir.
      # =====================================================
      conditionalPanel(
        condition = "input.ssdm_enable != true",
      tags$div(
        class = "boxpad dtbox",
        uiOutput("hdr_resp_box"),
        selectInput("resp_kind", "Response variable type",
                    choices = c("Presence / Absence" = "binary",
                                "Continuous"         = "continuous",
                                "Nominal"            = "nominal"),
                    selected = "binary"),

        conditionalPanel(
          condition = "input.resp_kind == 'continuous'",
          selectInput("dt_cont", "Continuous sub-type",
                      choices = c("Automatic" = "auto", "count" = "count",
                                  "abundance" = "abundance", "relative" = "relative"),
                      selected = "auto")
        ),

        conditionalPanel(
          condition = "input.resp_kind == 'nominal'",
          selectInput("dt_nom", "Nominal sub-type",
                      choices = c("Automatic" = "auto", "multiclass" = "multiclass",
                                  "ordinal" = "ordinal"),
                      selected = "auto"),
          textInput(
            "ord_levels",
            "Ordinal class order", value = "", placeholder = ""
          )
        ),

        uiOutput("dt_note")
      )
      ),

      uiOutput("inputs_hdr"),

      tags$div(
        class = "boxpad",
        textInput(
          "asc_dir",
          "Raster folder (current)",
          value = ""
        ),
        uiOutput("btn_asc_dir"),
        actionButton("check_current", "\u2705 Check current raster folder"),
        uiOutput("current_check_ui"),
        tags$div(
          class = "smallnote",
          "Check sonucu: klasor okunur, stack olusturulur, degisken isimleri listelenir."
        ),
        uiOutput("current_vars_ui")
      ),

      tags$div(
        class = "boxpad",
        h4("Predictor types (from current ASC)"),
        tags$p(
          class = "smallnote",
          tr_("vartypes_hint", "en")
        ),
        uiOutput("ui_var_types"),
        actionButton("apply_var_types", "\u2705 Apply variable types"),
        uiOutput("var_types_msg")
      ),

      tags$div(
        class = "boxpad",
        textInput("occ_file", "Occurrence CSV (full path)", value = ""),
        uiOutput("btn_occ_file"),
        actionButton("check_csv", "\u2705 Check CSV & detect species"),
        uiOutput("csv_check_ui"),
        textInput("sp_name", "Species name (auto-filled from CSV)", value = ""),

        conditionalPanel(
          condition = "input.ssdm_enable != true",
          tags$hr(),
          uiOutput("xy_hdr_ui"),
          uiOutput("xy_status_ui"),
          uiOutput("ui_x_col"),
          uiOutput("ui_y_col")
        ),

        conditionalPanel(
          condition = "input.resp_kind == 'binary' && input.ssdm_enable != true",
          tags$hr(),
          uiOutput("ui_pa_col"),
          uiOutput("pa_col_note_ui")
        ),

        conditionalPanel(
          condition = "input.resp_kind != 'binary'",
          tags$hr(),
          uiOutput("hdr_resp_col"),
          uiOutput("note_resp_col"),
          uiOutput("ui_resp_col"),
          uiOutput("resp_summary_ui"),
          checkboxInput("filter_raster", "filter.raster", value = FALSE)
        )
      ),

      conditionalPanel(
        condition = "input.resp_kind == 'binary'",
        tags$div(
          class = "boxpad",
          h4("MaxEnt"),
          textInput("maxent_jar", "MaxEnt jar (maxent.jar)", value = ""),
          uiOutput("btn_maxent_jar"),
          actionButton("check_jar", "\u2705 Check maxent.jar"),
          uiOutput("jar_check_ui")
        ),

        tags$div(
          class = "boxpad",
          h4("Pseudo-absence"),
          uiOutput("note_pa"),
          selectInput(
            "pa_strategy",
            "PA strategy",
            choices = c("disk", "random", "sre"),
            selected = "disk"
          ),
          numericInput("pa_rep", "PA.nb.rep", value = 3, min = 1, step = 1),
          numericInput("pa_n", "PA.nb.absences", value = 10000, min = 1000, step = 1000),
          numericInput("pa_dist_min", "PA.dist.min (m)", value = 2000, min = 0, step = 500),
          numericInput("pa_dist_max", "PA.dist.max (m)", value = 200000, min = 0, step = 5000)
        )
      ),

      tags$div(
        class = "boxpad",
        h4("Cross-validation"),
        helpText("k = 10 fixed. Choose repetitions (repeat full 10-fold)."),
        numericInput("kfold_rep", "k-fold repetitions", value = 1, min = 1, step = 1)
      ),

      tags$div(
        class = "boxpad",
        h4("Performance"),
        # Windows'ta bu alan GIZLENIR ve deger 1'e sabitlenir.
        # Gerekce: biomod2 paralel calisirken PSOCK kumesi kurar (fork yok).
        # Her isci ayri bir R sureci acip rasterleri yeniden yukler; bir iscide
        # paket yuklemesi veya java (MAXENT) cagrisi takilirsa is HATA VERMEDEN
        # asili kalir. Linux/macOS'ta fork kullanildigi icin bu risk yoktur.
        if (.Platform$OS.type == "windows") {
          tagList(
            numericInput("cpu_n", NULL, value = 1, min = 1, max = 1),
            tags$script(HTML(
              "$(document).ready(function(){$('#cpu_n').closest('.form-group').hide();});")),
            tags$p(class = "smallnote",
                   "Parallel processing is disabled on Windows (nb.cpu = 1) ",
                   "because biomod2 workers can hang without reporting an error.")
          )
        } else {
          nc <- suppressWarnings(parallel::detectCores(logical = FALSE))
          if (!is.finite(nc)) nc <- suppressWarnings(parallel::detectCores())
          if (!is.finite(nc)) nc <- 1L
          numericInput("cpu_n", "CPU cores (nb.cpu)",
                              value = max(1L, as.integer(nc) - 1L), min = 1, step = 1)
        }
      ),

      tags$div(
        class = "boxpad",
        h4("Models"),
        uiOutput("ui_models"),
        uiOutput("models_note")
      ),

      tags$div(
        class = "boxpad",
        uiOutput("hdr_thresh"),
        uiOutput("ui_thresholds"),
        uiOutput("ui_em_algo")
      ),

      tags$div(
        class = "boxpad",
        uiOutput("hdr_scale"),
        selectInput("out_scale", "Output scale",
                    choices = c("0 - 1" = "0_1", "0 - 1000" = "0_1000"),
                    selected = "0_1"),
        uiOutput("note_scale")
      ),

      tags$div(
        class = "boxpad",
        uiOutput("hdr_rc"),
        selectInput("rc_em_pick", "Ensemble curve",
                    choices = c("EMwmean", "EMmean", "EMmedian", "EMca"),
                    selected = "EMwmean"),
        textInput("var_labels", "Variable name translation",
                  value = "", placeholder = "anakaya=bedrock, yukselti=elevation"),
        uiOutput("note_rc")
      ),

      tags$div(
        class = "boxpad",
        h4("Optimizasyon (modeling options)"),
        selectInput(
          "opt_strategy", "Strateji",
          choices = c("adaptive" = "adaptive", "bigboss" = "bigboss",
                      "tuned" = "tuned", "default" = "default"),
          selected = "adaptive"
        ),
        uiOutput("note_opt")
      ),

      conditionalPanel(
        condition = "input.resp_kind != 'binary'",
        tags$div(
          class = "boxpad",
          uiOutput("hdr_proj"),
          numericInput("proj_digits", "digits", value = 3, min = 0, step = 1),
          conditionalPanel(
            condition = "input.resp_kind == 'continuous'",
            checkboxInput("export_classes", "Quantile class maps", value = TRUE),
            textInput("class_probs", "Quantile breaks", value = "0.25, 0.5, 0.75")
          )
        )
      ),

      tags$div(
        class = "boxpad",
        h4("Future (optional)"),
        selectInput(
          "future_n",
          "How many future periods?",
          choices = c("NA", "1", "2", "3", "4", "5"),
          selected = "NA"
        ),
        uiOutput("ui_future_blocks")
      ),

      hr(),
      actionButton("run", "\u25b6 Run biomod2", class = "btn-primary"),
      actionButton("refresh", "\u21bb Refresh now"),
      tags$p(class = "smallnote", "Auto-refresh every 5 seconds after you start a run."),
      hr(),

      uiOutput("progress_ui"),
      verbatimTextOutput("status")
    ),

    mainPanel(tabsetPanel(
      tabPanel("Maps & Plots", fluidRow(
        column(
          6,
          h4("Saved PNG plots"),
          tags$p(class = "smallnote", "Saved to: results/<run_id>/plots/ (PNG only)."),
          uiOutput("ui_plot_images"),

          hr(),
          h4("Image viewer"),
          uiOutput("ui_img_choice"),
          uiOutput("img_viewer_ui")
        ),
        column(
          6,
          h4("Raster viewer"),
          tags$p(class = "smallnote", "Uses TIF exports: results/<run_id>/rasters/*.tif"),
          uiOutput("ui_raster_selectors"),
          plotOutput("raster_plot", height = "520px"),
          verbatimTextOutput("raster_info")
        )
      )),
      tabPanel(
        "Metrics",
        uiOutput("metrics_header"),
        h4("Model metrics (mean across CV folds/runs)"),
        tableOutput("metrics_model"),
        hr(),
        h4("Ensemble metrics (if produced)"),
        tableOutput("metrics_ens"),

        hr(),
        h4("Ensemble membership (which models were used)"),
        tableOutput("ens_candidates_tbl"),
        tableOutput("ens_selected_tbl"),

        hr(),
        uiOutput("hdr_summary"),
        tableOutput("auc_tss_summary_tbl"),

        hr(),
        uiOutput("hdr_cutoffs"),
        uiOutput("cutoffs_note"),
        tableOutput("cutoffs_tbl")
      ),
      tabPanel("Log", verbatimTextOutput("logtail"))
    ))
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    run_id = NULL,
    run_dir = NULL,
    res_dir = NULL,
    log_file = NULL,
    csv_ok = FALSE,
    csv_msg = "Not checked yet.",
    current_ok = FALSE,
    current_msg = "Not checked yet.",
    current_vars = character(0),
    jar_ok = FALSE,
    jar_msg = "Not checked yet.",
    future_ok = list(),
    future_msg = list(),
    future_vars = list(),
    var_types = NULL,
    selected_raster = NULL,
    csv_df = NULL,
    resp_cols = character(0),
    dt_guess_cont = NULL,
    dt_guess_nom = NULL,
    ssdm_mat = NULL,
    ssdm_xy = NULL,
    ssdm_long = NULL
  )

  make_run_id <- function()
    format(Sys.time(), "%Y%m%d_%H%M%S")

  # =========================================================
  # DIL / LANGUAGE
  # Arayuz metinleri cevrilir; grafik ve CSV ciktilarinin dili
  # HER ZAMAN Ingilizcedir (analiz zinciri tutarliligi icin).
  # Deger kaybi olmamasi icin girdi ETIKETLERI guncellenir,
  # girdiler yeniden olusturulmaz.
  # =========================================================
  lang <- reactive(input$lang %||% "en")

  LABELS <- list(
    list("asc_dir","text","asc_dir"),            list("occ_file","text","occ_file"),
    list("sp_name","text","sp_name"),            list("maxent_jar","text","maxent_jar"),
    list("var_labels","text","var_labels"),      list("class_probs","text","class_probs"),
    list("ord_levels","text","ord_levels"),
    list("ssdm_mat_file","text","ssdm_mat_file"),list("ssdm_xy_file","text","ssdm_xy_file"),
    list("resp_kind","select","resp_kind"),      list("dt_cont","select","dt_cont"),
    list("dt_nom","select","dt_nom"),            list("out_scale","select","out_scale"),
    list("rc_em_pick","select","rc_em_pick"),    list("opt_strategy","select","opt_strategy"),
    list("pa_strategy","select","pa_strategy"),  list("raster_pick","select","raster_pick"),
    list("pa_rep","numeric","pa_rep"),           list("pa_n","numeric","pa_n"),
    list("pa_dist_min","numeric","pa_dist_min"), list("pa_dist_max","numeric","pa_dist_max"),
    list("kfold_rep","numeric","kfold_rep"),     list("cpu_n","numeric","cpu_n"),
    list("proj_digits","numeric","proj_digits"), list("ssdm_min_occ","numeric","ssdm_min_occ"),
    list("ssdm_min_freq","numeric","ssdm_min_freq"),
    list("ssdm_min_abs","numeric","ssdm_min_abs"),
    list("ens_min_auc","numeric","min_auc"),     list("ens_min_tss","numeric","min_tss"),
    list("ens_min_r2","numeric","min_r2"),       list("ens_rmse_tol","numeric","rmse_tol"),
    list("ens_min_acc","numeric","min_acc"),     list("ens_min_f1","numeric","min_f1"),
    list("filter_raster","check","filter_raster"),
    list("export_classes","check","export_classes"),
    list("ssdm_enable","check","ssdm_enable"),   list("ssdm_latlon","check","ssdm_latlon"),
    list("ssdm_pssdm","check","ssdm_pssdm"),     list("ssdm_change","check","ssdm_change"),
    list("check_current","action","check_current"),
    list("apply_var_types","action","apply_vartypes"),
    list("check_csv","action","check_csv"),      list("check_jar","action","check_jar"),
    list("ssdm_check","action","ssdm_check"),    list("ssdm_convert","action","ssdm_convert"),
    list("run","action","run_btn"),              list("refresh","action","refresh_btn")
  )

  observeEvent(lang(), {
    lg <- lang()
    for (e in LABELS) {
      id <- e[[1]]; ty <- e[[2]]; key <- e[[3]]
      lb <- tr_(key, lg)
      switch(ty,
        text    = updateTextInput(session, id, label = lb),
        numeric = updateNumericInput(session, id, label = lb),
        check   = updateCheckboxInput(session, id, label = lb),
        action  = updateActionButton(session, id, label = lb),
        select  = updateSelectInput(session, id, label = lb))
    }
    # Secenek ETIKETLERI de cevrilir, DEGERLER korunur
    updateSelectInput(session, "resp_kind",
      choices = tr_choices(c("resp_binary","resp_continuous","resp_nominal"),
                           c("binary","continuous","nominal"), lg),
      selected = isolate(input$resp_kind) %||% "binary")
    updateSelectInput(session, "dt_cont",
      choices = tr_choices(c("dt_cont_auto","dt_cont_count","dt_cont_abund","dt_cont_rel"),
                           c("auto","count","abundance","relative"), lg),
      selected = isolate(input$dt_cont) %||% "auto")
    updateSelectInput(session, "dt_nom",
      choices = tr_choices(c("dt_nom_auto","dt_nom_multi","dt_nom_ord"),
                           c("auto","multiclass","ordinal"), lg),
      selected = isolate(input$dt_nom) %||% "auto")
    updateSelectInput(session, "out_scale",
      choices = tr_choices(c("scale_01","scale_01000"), c("0_1","0_1000"), lg),
      selected = isolate(input$out_scale) %||% "0_1")
    updateSelectInput(session, "opt_strategy",
      choices = tr_choices(c("opt_adaptive","opt_bigboss","opt_tuned","opt_default"),
                           c("adaptive","bigboss","tuned","default"), lg),
      selected = isolate(input$opt_strategy) %||% "adaptive")
  }, ignoreInit = FALSE)

  hdr <- function(key) renderUI(h4(tr_(key, lang())))
  note <- function(key) renderUI(tags$p(class = "smallnote", tr_(key, lang())))


  # ---- Dile bagli statik metinler (uiOutput; dil degisince yenilenir) ----
  output$hdr_resp_box <- hdr("resp_box")
  output$hdr_thresh   <- hdr("thresh_box")
  output$hdr_scale    <- hdr("scale_box")
  output$hdr_rc       <- hdr("rc_box")
  output$hdr_proj     <- hdr("proj_box")
  output$hdr_summary  <- hdr("summary_tbl")
  output$hdr_cutoffs  <- hdr("cutoffs_tbl")
  output$hdr_resp_col <- renderUI(tags$strong(tr_("resp_col_hdr", lang())))
  output$note_resp_col <- note("resp_col_note")
  output$note_pa       <- note("pa_note")
  output$note_scale    <- note("scale_note")
  output$note_rc       <- note("rc_note")
  output$note_opt      <- note("opt_note")

  # shinyFiles dugmeleri: etiket ve baslik cevrilir
  output$btn_asc_dir <- renderUI(
    shinyDirButton("asc_dir_btn", tr_("pick_folder", lang()),
                   tr_("pick_folder_ttl", lang())))
  output$btn_occ_file <- renderUI(
    shinyFilesButton("occ_file_btn", tr_("pick_csv", lang()),
                     tr_("pick_csv_ttl", lang()), multiple = FALSE))
  output$btn_maxent_jar <- renderUI(
    shinyFilesButton("maxent_jar_btn", tr_("pick_jar", lang()),
                     tr_("pick_jar_ttl", lang()), multiple = FALSE))

  output$inputs_hdr        <- hdr("inputs_box")
  output$ssdm_hdr          <- hdr("ssdm_box")
  output$ssdm_mat_hdr_ui   <- renderUI(tagList(tags$hr(), tags$strong(tr_("ssdm_mat_hdr", lang())),
                                               tags$p(class="smallnote", tr_("ssdm_mat_note", lang()))))
  output$ssdm_xy_hdr_ui    <- renderUI(tagList(tags$hr(), tags$strong(tr_("ssdm_xy_hdr", lang())),
                                               tags$p(class="smallnote", tr_("ssdm_xy_note", lang()))))
  output$ssdm_note_ui      <- note("ssdm_note")
  # ---- Koordinat sutunlari: otomatik sapta, gerekirse elle sectir ----
  xy_detect <- reactive({
    df <- rv$csv_df
    if (is.null(df)) return(NULL)
    find_xy_cols(df, input$x_col, input$y_col)
  })

  output$xy_hdr_ui <- renderUI(tags$strong(tr_("xy_hdr", lang())))

  output$xy_status_ui <- renderUI({
    f <- xy_detect(); lg <- lang()
    if (is.null(f)) return(tags$div(class = "smallnote", tr_("ssdm_not_checked", lg)))
    if (isTRUE(f$ok)) tags$div(class = "ok", tr_("xy_auto", lg, f$x, f$y))
    else tags$div(class = "bad", tr_("xy_missing", lg))
  })

  output$ui_x_col <- renderUI({
    df <- rv$csv_df; if (is.null(df)) return(NULL)
    f <- find_xy_cols(df, NULL, NULL)
    selectInput("x_col", tr_("x_col", lang()),
                choices = c("(auto)" = "", names(df)),
                selected = isolate(input$x_col) %||% "")
  })
  output$ui_y_col <- renderUI({
    df <- rv$csv_df; if (is.null(df)) return(NULL)
    selectInput("y_col", tr_("y_col", lang()),
                choices = c("(auto)" = "", names(df)),
                selected = isolate(input$y_col) %||% "")
  })

  # Tek-tur binary yolunda GERCEK yokluk sutunu (opsiyonel)
  output$ui_pa_col <- renderUI({
    cols <- rv$resp_cols %||% character(0)
    selectInput("pa_col", tr_("pa_col", lang()),
                choices = c("(none - presence only)" = "", cols),
                selected = isolate(input$pa_col) %||% "")
  })
  output$pa_col_note_ui <- note("pa_col_note")

  output$ssdm_filter_hdr_ui   <- renderUI(tagList(tags$hr(),
                                    tags$strong(tr_("ssdm_filter_hdr", lang()))))
  output$ssdm_filter_note_ui  <- note("ssdm_filter_note")

  # Suzgec REAKTIF: yuzdeyi degistirince yeniden donusturmeye gerek yok.
  ssdm_filtered <- reactive({
    lf <- rv$ssdm_long
    if (is.null(lf) || !isTRUE(lf$ok)) return(NULL)
    ssdm_filter_species(lf, min_freq_pct = input$ssdm_min_freq %||% 0)
  })

  output$ssdm_filter_ui <- renderUI({
    f <- ssdm_filtered(); lg <- lang()
    if (is.null(f)) return(NULL)
    show <- function(d, cols) {
      d <- d[, cols, drop = FALSE]
      d$frequency_pct <- round(d$frequency_pct, 2)
      utils::head(d, 20)
    }
    tagList(
      tags$div(class = "smallnote",
               if (f$min_freq_pct <= 0) tr_("ssdm_filter_off", lg)
               else tr_("ssdm_filter_eff", lg, f$min_freq_pct, f$n_plots, f$freq_n)),
      if (f$n_pa > 0) tags$div(class = "smallnote", tr_("ssdm_pa_note", lg, f$n_pa)),
      tags$div(class = if (f$n_kept >= 2) "ok" else "bad",
               tr_("ssdm_filter_res", lg, f$n_kept, f$n_total, f$n_dropped)),
      if (f$n_dropped > 0) tagList(
        tags$div(class = "smallnote", tr_("ssdm_filter_drop", lg)),
        tags$div(class = "previewtbl",
                 HTML(renderTable(show(f$dropped,
                   c("species","n_presence","n_absence","frequency_pct","reason")),
                   rownames = FALSE)()))
      ),
      tags$div(class = "smallnote", tr_("ssdm_filter_kept", lg)),
      tags$div(class = "previewtbl",
               HTML(renderTable(show(f$kept,
                 c("species","n_presence","n_absence","frequency_pct","absence_strategy")),
                 rownames = FALSE)()))
    )
  })

  # =========================================================
  # S-SDM: dosya kontrolu, onizleme, donusum
  # =========================================================
  shinyFileChoose(input, "ssdm_mat_btn", roots = volumes, session = session,
                  filetypes = c("csv", "txt"))
  shinyFileChoose(input, "ssdm_xy_btn", roots = volumes, session = session,
                  filetypes = c("csv", "txt"))
  observeEvent(input$ssdm_mat_btn, {
    f <- parseFilePaths(volumes, input$ssdm_mat_btn)
    if (nrow(f) > 0) updateTextInput(session, "ssdm_mat_file", value = as.character(f$datapath[1]))
  })
  observeEvent(input$ssdm_xy_btn, {
    f <- parseFilePaths(volumes, input$ssdm_xy_btn)
    if (nrow(f) > 0) updateTextInput(session, "ssdm_xy_file", value = as.character(f$datapath[1]))
  })

  preview_tbl <- function(d, n = 10) {
    d <- utils::head(as.data.frame(d), n)
    if (ncol(d) > 12) d <- d[, seq_len(12), drop = FALSE]
    d
  }

  # DIKKAT: renderTable(...)() bir HTML METNI dondurur. Bunu dogrudan bir tag
  # icine koyarsan Shiny metni kacirir ve ekranda ham <table> etiketleri gorunur.
  # HTML() ile sarmalamak sart.
  preview_html <- function(d, n = 10) {
    tags$div(class = "previewtbl",
             HTML(renderTable(preview_tbl(d, n), rownames = FALSE)()))
  }

  observeEvent(input$ssdm_check, {
    lg <- lang()
    m  <- ssdm_read_matrix(input$ssdm_mat_file %||% "")
    cc <- ssdm_read_coords(input$ssdm_xy_file %||% "")
    rv$ssdm_mat <- m; rv$ssdm_xy <- cc; rv$ssdm_long <- NULL

    output$ssdm_check_ui <- renderUI({
      parts <- list()
      if (isTRUE(m$ok)) {
        parts <- c(parts, list(
          tags$div(class = "ok", paste0("\u2714 ", tr_("ssdm_preview_mat", lg), " \u2014 ", m$msg)),
          preview_html(m$data)
        ))
      } else parts <- c(parts, list(tags$div(class = "bad", paste0("\u2716 ", m$msg))))

      if (isTRUE(cc$ok)) {
        parts <- c(parts, list(
          tags$div(class = "ok", paste0("\u2714 ", tr_("ssdm_preview_xy", lg), " \u2014 ", cc$msg)),
          preview_html(cc$data)
        ))
        if (length(cc$duplicate_plots) > 0) {
          parts <- c(parts, list(tags$div(class = "bad",
            paste0(tr_("dup_plots_warn", lg),
                   paste(cc$duplicate_plots, collapse = ", ")))))
        }
      } else parts <- c(parts, list(tags$div(class = "bad", paste0("\u2716 ", cc$msg))))
      do.call(tagList, parts)
    })
  }, ignoreInit = TRUE)

  output$ssdm_xcol_ui <- renderUI({
    cc <- rv$ssdm_xy
    if (is.null(cc) || !isTRUE(cc$ok)) return(NULL)
    g <- ssdm_guess_xy(cc$coord_cols, cc$data)
    selectInput("ssdm_xcol", tr_("ssdm_xcol", lang()), choices = cc$coord_cols, selected = g$x)
  })
  output$ssdm_ycol_ui <- renderUI({
    cc <- rv$ssdm_xy
    if (is.null(cc) || !isTRUE(cc$ok)) return(NULL)
    g <- ssdm_guess_xy(cc$coord_cols, cc$data)
    selectInput("ssdm_ycol", tr_("ssdm_ycol", lang()), choices = cc$coord_cols, selected = g$y)
  })

  observeEvent(input$ssdm_convert, {
    lg <- lang()
    if (is.null(rv$ssdm_mat) || is.null(rv$ssdm_xy)) {
      output$ssdm_convert_ui <- renderUI(tags$div(class = "bad", tr_("ssdm_not_checked", lg)))
      return(NULL)
    }
    r <- ssdm_convert(input$ssdm_mat_file %||% "", input$ssdm_xy_file %||% "",
                      x_col = input$ssdm_xcol, y_col = input$ssdm_ycol)
    if (!isTRUE(r$ok)) {
      rv$ssdm_long <- NULL
      output$ssdm_convert_ui <- renderUI(tags$div(class = "bad", paste0("\u2716 ", r$msg)))
      return(NULL)
    }
    rv$ssdm_long <- r

    output$ssdm_convert_ui <- renderUI({
      warn <- list()
      if (length(r$missing_xy) > 0)
        warn <- c(warn, list(tags$div(class = "smallnote",
          paste0(tr_("dropped_no_xy", lg),
                 paste(utils::head(r$missing_xy, 8), collapse = ", "),
                 if (length(r$missing_xy) > 8) " ..." else ""))))
      if (length(r$missing_mat) > 0)
        warn <- c(warn, list(tags$div(class = "smallnote",
          paste0(tr_("dropped_no_mat", lg),
                 paste(utils::head(r$missing_mat, 8), collapse = ", "),
                 if (length(r$missing_mat) > 8) " ..." else ""))))
      if (length(r$duplicates) > 0)
        warn <- c(warn, list(tags$div(class = "smallnote",
          paste0(tr_("dup_merged", lg), paste(r$duplicates, collapse = ", ")))))
      if (length(r$duplicate_plots) > 0)
        warn <- c(warn, list(tags$div(class = "bad",
          paste0(tr_("dup_plots_warn", lg),
                 paste(r$duplicate_plots, collapse = ", ")))))

      do.call(tagList, c(list(
        tags$div(class = "ok",
                 tr_("ssdm_conv_ok", lg, r$n_species, r$n_plots, r$n_rows)),
        tags$div(class = "smallnote",
                 paste0(tr_("ssdm_preview_lng", lg), " \u2014 ", r$msg,
                        "  |  X = ", r$x_col, ", Y = ", r$y_col)),
        preview_html(r$data),
        preview_html(r$per_species, 15)
      ), warn))
    })
  }, ignoreInit = TRUE)


  # =========================================================
  # SISTEM DENETIMI PANELI
  # =========================================================
  output$setup_ui <- renderUI({
    lg <- lang()
    s <- setup_status()
    rows <- list()
    add <- function(ok, txt) rows[[length(rows) + 1]] <<-
      tags$div(class = if (ok) "ok" else "bad", paste0(if (ok) "OK " else "X ", txt))

    add(length(s$files_missing) == 0,
        if (length(s$files_missing) == 0) tr_("files_ok", lg, length(REQUIRED_FILES), length(REQUIRED_FILES))
        else tr_("files_missing", lg, paste(s$files_missing, collapse = ", ")))

    add(length(s$core_missing) == 0,
        if (length(s$core_missing) == 0) tr_("core_ok", lg)
        else tr_("core_missing", lg, paste(s$core_missing, collapse = ", ")))

    if (!is.na(s$bm_ver)) {
      add(s$bm_ok43,
          if (s$bm_ok43) tr_("bm_ok", lg, s$bm_ver)
          else tr_("bm_old", lg, s$bm_ver))
    }

    if (length(s$models_blocked) > 0) {
      add(FALSE, tr_("models_missing", lg, paste(s$models_blocked, collapse = ", "), paste(s$mp_missing, collapse = ", ")))
    } else {
      add(TRUE, tr_("models_ok", lg))
    }

    if (!s$java_ok) add(FALSE, tr_("java_missing", lg))
    if (is.na(s$rscript)) add(FALSE, tr_("rscript_missing", lg, R.home("bin")))
    else add(TRUE, tr_("rscript_ok", lg))

    if (length(s$core_missing) > 0 || length(s$mp_missing) > 0) {
      rows[[length(rows) + 1]] <- tags$div(
        class = "smallnote",
        tr_("install_hint", lg))
    }

    tagList(
      tags$strong(tr_("setup_status", lg)),
      tags$div(class = "smallnote", paste0(tr_("working_dir", lg), ": ", s$wd)),
      rows
    )
  })

  # =========================================================
  # ETKIN VERI TIPI (data.type)
  # =========================================================
  data_type <- reactive({
    k <- input$resp_kind %||% "binary"
    if (identical(k, "binary")) return("binary")

    if (identical(k, "continuous")) {
      s <- input$dt_cont %||% "auto"
      if (!identical(s, "auto")) return(s)
      return(rv$dt_guess_cont %||% "abundance")
    }

    s <- input$dt_nom %||% "auto"
    if (!identical(s, "auto")) return(s)
    rv$dt_guess_nom %||% "multiclass"
  })

  output$dt_note <- renderUI({
    dt  <- data_type()
    fam <- dt_family(dt)
    mets <- paste(metrics_for(dt), collapse = ", ")

    desc <- switch(
      fam,
      binary = tr_("dt_note_binary", lg),
      quant  = tr_("dt_note_quant", lg),
      qual   = tr_("dt_note_qual", lg)
    )

    tagList(
      tags$p(tags$span(class = "dtbadge", paste0("data.type = ", dt))),
      tags$div(class = "smallnote", desc),
      tags$div(class = "smallnote", tags$strong("Metrikler: "), mets)
    )
  })

  # =========================================================
  # DOSYA / KLASOR SECICILER (shinyFiles)
  #
  # WINDOWS NOTU: shinyFiles::getVolumes() surucu listesini
  # "wmic logicaldisk get Caption, VolumeName" ciktisini ayristirarak kurar.
  # Guncel Windows surumlerinde WMIC kaldirilmaktadir; ayrica surucu etiketinde
  # bosluk varsa ayristirma bozulur ve su uyari cikar:
  #   "data length [7] is not a sub-multiple or multiple of the number of rows [4]"
  # Bu uyari YALNIZCA secici dugmelerini etkiler, modelleme isini etkilemez.
  # Yine de burada guvenli bir yedek surucu listesi kuruyoruz.
  # =========================================================
  safe_volumes <- function() {
    base_v <- c("Working directory" = getwd(), "Home" = path.expand("~"))

    sys_v <- tryCatch(
      withCallingHandlers(
        shinyFiles::getVolumes()(),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    sys_v <- sys_v[nzchar(names(sys_v) %||% "") & nzchar(sys_v)]
    sys_v <- sys_v[vapply(sys_v, function(p) isTRUE(dir.exists(p)), logical(1))]

    # Yedek: getVolumes() bos donerse platforma gore kokleri dogrudan sina
    if (length(sys_v) == 0) {
      if (.Platform$OS.type == "windows") {
        cand <- paste0(LETTERS, ":/")
        nm   <- sub("/$", "", cand)
      } else {
        # macOS: /Volumes altinda takili diskler; Linux: /media, /mnt
        mounts <- unlist(lapply(c("/Volumes", "/media", "/mnt"), function(d)
          if (dir.exists(d)) list.dirs(d, recursive = FALSE, full.names = TRUE) else character(0)))
        cand <- c("/", mounts)
        nm   <- c("Root", basename(mounts))
      }
      keep <- vapply(cand, function(p) isTRUE(dir.exists(p)), logical(1))
      if (any(keep)) sys_v <- setNames(cand[keep], nm[keep])
    }

    v <- c(base_v, sys_v)
    v[!duplicated(v)]
  }

  volumes <- safe_volumes()

  # --- Guncel raster klasoru ---
  shinyDirChoose(input, "asc_dir_btn", roots = volumes, session = session)
  observeEvent(input$asc_dir_btn, {
    pth <- try(parseDirPath(volumes, input$asc_dir_btn), silent = TRUE)
    if (!inherits(pth, "try-error") && length(pth) > 0 && nzchar(pth)) {
      updateTextInput(session, "asc_dir", value = as.character(pth))
    }
  }, ignoreInit = TRUE)

  # --- Presence CSV ---
  shinyFileChoose(input, "occ_file_btn", roots = volumes, session = session,
                  filetypes = c("csv", "CSV", "txt"))
  observeEvent(input$occ_file_btn, {
    f <- try(parseFilePaths(volumes, input$occ_file_btn), silent = TRUE)
    if (!inherits(f, "try-error") && nrow(f) > 0) {
      updateTextInput(session, "occ_file", value = as.character(f$datapath[1]))
    }
  }, ignoreInit = TRUE)

  # --- maxent.jar ---
  shinyFileChoose(input, "maxent_jar_btn", roots = volumes, session = session,
                  filetypes = c("jar", "JAR"))
  observeEvent(input$maxent_jar_btn, {
    f <- try(parseFilePaths(volumes, input$maxent_jar_btn), silent = TRUE)
    if (!inherits(f, "try-error") && nrow(f) > 0) {
      updateTextInput(session, "maxent_jar", value = as.character(f$datapath[1]))
    }
  }, ignoreInit = TRUE)

  # --- Gelecek senaryo klasorleri (5 donem x 4 senaryo) ---
  # Handler'lar renderUI'dan bagimsiz olarak ONCEDEN kurulur.
  for (i in 1:5) {
    for (j in 1:4) {
      local({
        ii <- i; jj <- j
        btn <- paste0("future_dir_btn_", ii, "_", jj)
        txt <- paste0("asc_future_dir_", ii, "_", jj)

        shinyDirChoose(input, btn, roots = volumes, session = session)

        observeEvent(input[[btn]], {
          pth <- try(parseDirPath(volumes, input[[btn]]), silent = TRUE)
          if (!inherits(pth, "try-error") && length(pth) > 0 && nzchar(pth)) {
            updateTextInput(session, txt, value = as.character(pth))
          }
        }, ignoreInit = TRUE)
      })
    }
  }

  # ---- Future UI blocks (manual paths + check) ----
  output$ui_future_blocks <- renderUI({
    if (input$future_n == "NA") {
      return(tags$div(tags$em("Future is NA \u2192 only current maps will be produced.")))
    }
    nf <- as.integer(input$future_n)
    tagList(lapply(seq_len(nf), function(i) {
      tags$div(
        class = "boxpad",
        tags$h5(paste0("Future Period ", i)),
        textInput(paste0("period_label_", i), "Period label (optional)", value = paste0("Future", i)),
        selectInput(paste0("future_scen_n_", i), "How many scenarios for this period?",
                    choices = c("1", "2", "3", "4"), selected = "1"),
        uiOutput(paste0("ui_scenarios_", i))
      )
    }))
  })

  for (i in 1:5) {
    local({
      ii <- i
      output[[paste0("ui_scenarios_", ii)]] <- renderUI({
        req(input$future_n)
        if (input$future_n == "NA") return(NULL)
        if (as.integer(input$future_n) < ii) return(NULL)

        nsc <- as.integer(input[[paste0("future_scen_n_", ii)]] %||% "1")

        tagList(lapply(seq_len(nsc), function(j) {
          key <- paste0("P", ii, "_S", j)
          tags$div(
            class = "boxpad",
            tags$strong(paste0("Scenario ", j)),
            br(),
            selectInput(
              paste0("ssp_", ii, "_", j),
              "SSP label",
              choices = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"),
              selected = c("SSP2-4.5", "SSP5-8.5", "SSP1-2.6", "SSP3-7.0")[(j - 1) %% 4 + 1]
            ),
            textInput(paste0("sc_label_", ii, "_", j),
                      "Label override (optional; leave blank to use SSP label)",
                      value = ""),
            textInput(paste0("asc_future_dir_", ii, "_", j),
                      "Raster folder (future scenario)",
                      value = ""),
            shinyDirButton(paste0("future_dir_btn_", ii, "_", j),
                           "Browse",
                           tr_("scen_pick_ttl", lang(), j)),
            actionButton(paste0("check_future_", ii, "_", j), "\u2705 Check this future folder"),
            uiOutput(paste0("future_check_ui_", ii, "_", j)),
            uiOutput(paste0("future_vars_ui_", ii, "_", j)),
            tags$div(class = "smallnote", paste0("Key: ", key))
          )
        }))
      })
    })
  }

  # ---- Check: CSV ----
  observeEvent(input$check_csv, {
    p <- input$occ_file %||% ""
    if (!nzchar(p)) {
      rv$csv_ok <- FALSE
      rv$csv_msg <- "CSV path is empty."
      rv$csv_df <- NULL; rv$resp_cols <- character(0)
    } else if (!file.exists(p)) {
      rv$csv_ok <- FALSE
      rv$csv_msg <- paste0("CSV not found: ", p)
      rv$csv_df <- NULL; rv$resp_cols <- character(0)
    } else {
      sp <- detect_species_name(p)
      df <- read_csv_head(p)

      if (is.null(df)) {
        rv$csv_ok <- FALSE
        rv$csv_msg <- "CSV okunamadi veya bos."
        rv$csv_df <- NULL; rv$resp_cols <- character(0)
      } else {
        rv$csv_ok <- TRUE
        rv$csv_df <- df
        rv$resp_cols <- resp_candidates(df)

        xy_ok <- all(c("x", "y") %in% names(df))
        rv$csv_msg <- paste0(
          "OK. CSV okundu (", nrow(df), " satir, ", ncol(df), " sutun). ",
          "Tur: ", if (nzchar(sp)) sp else "(saptanamadi)",
          if (xy_ok) " | x,y bulundu." else " | UYARI: x,y sutunlari bulunamadi!"
        )
        if (nzchar(sp)) updateTextInput(session, "sp_name", value = sp)
      }
    }
  }, ignoreInit = TRUE)

  output$csv_check_ui <- renderUI({
    cls <- if (isTRUE(rv$csv_ok)) "ok" else "bad"
    tags$div(class = cls, rv$csv_msg)
  })

  # ---- Yanit sutunu secimi (surekli / nominal) ----
  output$ui_resp_col <- renderUI({
    cols <- rv$resp_cols %||% character(0)
    if (!isTRUE(rv$csv_ok) || length(cols) == 0) {
      return(tags$div(class = "bad", "Once 'Check CSV & detect species' calistirin."))
    }
    sel <- input$resp_col %||% cols[1]
    if (!(sel %in% cols)) sel <- cols[1]
    selectInput("resp_col", tr_("resp_col", lang()), choices = cols, selected = sel)
  })

  # Secilen yanit sutununa gore alt tip tahmini
  observeEvent(list(input$resp_col, rv$csv_df, input$resp_kind), {
    df <- rv$csv_df
    cc <- input$resp_col %||% ""
    if (is.null(df) || !nzchar(cc) || !(cc %in% names(df))) return(NULL)
    v <- df[[cc]]
    rv$dt_guess_cont <- guess_continuous_type(v)
    rv$dt_guess_nom  <- guess_nominal_type(v)
  }, ignoreInit = TRUE)

  output$resp_summary_ui <- renderUI({
    df <- rv$csv_df
    cc <- input$resp_col %||% ""
    if (is.null(df) || !nzchar(cc) || !(cc %in% names(df))) return(NULL)

    txt <- resp_summary_text(df[[cc]])
    gs <- if (identical(input$resp_kind, "nominal")) rv$dt_guess_nom else rv$dt_guess_cont

    tagList(
      tags$div(class = "smallnote", tags$strong(tr_("col_summary_lbl", lang())), txt),
      tags$div(class = "smallnote",
               tags$strong(tr_("guess_lbl", lang())), gs %||% "-",
               tr_("guess_hint", lang()))
    )
  })

  # =========================================================
  # VERI TIPINE BAGLI DINAMIK UI
  # =========================================================
  # NOT: isolate() sart. Aksi halde renderUI kendi urettigi input'a bagimli olur
  # ve sonsuz yeniden cizim dongusu olusur.
  output$ui_models <- renderUI({
    dt <- data_type()
    ch <- ALLOWED_MODELS[[dt]]
    prev <- isolate(input$models) %||% character(0)
    sel <- intersect(prev, ch)
    if (length(sel) == 0) sel <- intersect(DEFAULT_MODELS[[dt]], ch)

    checkboxGroupInput("models", "Select models", choices = ch, selected = sel)
  })

  output$models_note <- renderUI({
    dt <- data_type()
    dropped <- setdiff(ALLOWED_MODELS$binary, ALLOWED_MODELS[[dt]])
    if (length(dropped) == 0) {
      return(tags$p(class = "smallnote", tr_("models_pkg_note", lang())))
    }
    tags$p(class = "smallnote",
           tr_("models_drop_note", lang(), dt, paste(dropped, collapse = ", ")))
  })

  output$ui_thresholds <- renderUI({
    fam <- dt_family(data_type())
    lg <- lang()
    near_one <- tagList(
      numericInput("ens_near_one", tr_("ens_near_one", lg),
                   value = isolate(input$ens_near_one) %||% 0.999,
                   min = 0.5, max = 1, step = 0.001),
      tags$p(class = "smallnote", tr_("near_one_note", lg))
    )

    if (identical(fam, "binary")) {
      return(tagList(
        numericInput("ens_min_auc", "Minimum AUCroc", value = 0.70, min = 0, max = 1, step = 0.01),
        numericInput("ens_min_tss", "Minimum TSS", value = 0.40, min = 0, max = 1, step = 0.01),
        near_one,
        tags$p(class = "smallnote", tr_("thresh_bin_note", lg))
      ))
    }

    if (identical(fam, "quant")) {
      return(tagList(
        numericInput("ens_min_r2", "Minimum R\u00b2", value = 0.30, min = 0, max = 1, step = 0.05),
        numericInput("ens_rmse_tol", tr_("rmse_tol", lg), value = 0, min = 0, step = 0.1),
        near_one,
        tags$p(class = "smallnote", tr_("thresh_q_note", lg))
      ))
    }

    tagList(
      numericInput("ens_min_acc", "Minimum Accuracy", value = 0.50, min = 0, max = 1, step = 0.05),
      numericInput("ens_min_f1", "Minimum F1", value = 0.40, min = 0, max = 1, step = 0.05),
      near_one,
      tags$p(class = "smallnote", tr_("thresh_bin_note", lg))
    )
  })

  output$ui_em_algo <- renderUI({
    fam <- dt_family(data_type())
    ch <- EM_ALGOS[[fam]]
    prev <- isolate(input$em_algo) %||% character(0)
    sel <- intersect(prev, ch)
    if (length(sel) == 0) sel <- setdiff(ch, "EMci")

    tagList(
      checkboxGroupInput("em_algo", tr_("em_algo", lang()), choices = ch, selected = sel),
      tags$p(class = "smallnote",
             switch(fam,
                    binary = tr_("em_note_bin", lang()),
                    quant  = tr_("em_note_quant", lang()),
                    qual   = tr_("em_note_qual", lang())))
    )
  })

  output$metrics_header <- renderUI({
    dt <- data_type()
    tags$p(tags$span(class = "dtbadge", paste0("data.type = ", dt)),
           tags$span(class = "smallnote",
                     paste0("  \u2022  metrikler: ", paste(metrics_for(dt), collapse = ", "))))
  })

  output$cutoffs_note <- renderUI({
    if (identical(dt_family(data_type()), "binary")) {
      tags$p(class = "smallnote",
             tr_("cutoff_note_bin", lang()))
    } else {
      tags$p(class = "smallnote",
             tr_("cutoff_note_cont", lang()))
    }
  })

  # ---- Check: MaxEnt jar ----
  observeEvent(input$check_jar, {
    p <- input$maxent_jar %||% ""
    if (!nzchar(p)) {
      rv$jar_ok <- FALSE
      rv$jar_msg <- "Jar path is empty."
    } else if (!file.exists(p)) {
      rv$jar_ok <- FALSE
      rv$jar_msg <- paste0("Jar not found: ", p)
    } else if (!grepl("\\.jar$", p, ignore.case = TRUE)) {
      rv$jar_ok <- TRUE
      rv$jar_msg <- paste0("Found file (not .jar extension, but exists): ", p)
    } else {
      rv$jar_ok <- TRUE
      rv$jar_msg <- paste0("OK. maxent.jar found: ", p)
    }
  }, ignoreInit = TRUE)

  output$jar_check_ui <- renderUI({
    cls <- if (isTRUE(rv$jar_ok)) "ok" else "bad"
    tags$div(class = cls, rv$jar_msg)
  })

  # ---- Check: Current raster folder ----
  observeEvent(input$check_current, {
    p <- input$asc_dir %||% ""
    if (!nzchar(p)) {
      rv$current_ok <- FALSE
      rv$current_msg <- "Current raster folder path is empty."
      rv$current_vars <- character(0)
      rv$var_types <- NULL
    } else if (!dir.exists(p)) {
      rv$current_ok <- FALSE
      rv$current_msg <- paste0("Folder not found: ", p)
      rv$current_vars <- character(0)
      rv$var_types <- NULL
    } else {
      r <- try(read_env_stack(p), silent = TRUE)
      if (inherits(r, "try-error")) {
        rv$current_ok <- FALSE
        rv$current_msg <- paste0("Failed to read rasters: ", as.character(r))
        rv$current_vars <- character(0)
        rv$var_types <- NULL
      } else {
        rv$current_ok <- TRUE
        rv$current_msg <- paste0("OK. Read ", terra::nlyr(r), " layers.")
        rv$current_vars <- names(r)

        vars <- rv$current_vars
        vt <- setNames(as.list(rep("continuous", length(vars))), vars)
        for (v in vars) {
          if (tolower(v) %in% c("anakaya", "geology", "lithology", "landuse", "soil")) {
            vt[[v]] <- "categorical"
          }
        }
        rv$var_types <- vt
      }
    }
  }, ignoreInit = TRUE)

  output$current_check_ui <- renderUI({
    cls <- if (isTRUE(rv$current_ok)) "ok" else "bad"
    tags$div(class = cls, rv$current_msg)
  })

  output$current_vars_ui <- renderUI({
    if (!isTRUE(rv$current_ok)) return(NULL)
    tags$div(tags$strong("Detected current variables:"), tags$pre(paste(rv$current_vars, collapse = "\n")))
  })

  # ---- Variable types UI ----
  output$ui_var_types <- renderUI({
    vars <- rv$current_vars %||% character(0)
    if (!isTRUE(rv$current_ok) || length(vars) == 0) {
      return(tags$em("No variables yet. Run 'Check current ASC folder' first."))
    }

    tagList(lapply(vars, function(v) {
      sel <- "continuous"
      if (!is.null(rv$var_types) && !is.null(rv$var_types[[v]])) sel <- rv$var_types[[v]]

      tags$div(
        class = "vartype-row",
        tags$div(class = "vartype-name", v),
        tags$div(
          class = "vartype-select",
          selectInput(
            inputId = paste0("vartype__", v),
            label = NULL,
            choices = c("continuous", "categorical", "nominal"),
            selected = sel
          )
        )
      )
    }))
  })

  observeEvent(input$apply_var_types, {
    vars <- rv$current_vars %||% character(0)
    if (!isTRUE(rv$current_ok) || length(vars) == 0) return(NULL)

    vt <- setNames(as.list(rep("continuous", length(vars))), vars)
    for (v in vars) vt[[v]] <- input[[paste0("vartype__", v)]] %||% "continuous"
    rv$var_types <- vt
  }, ignoreInit = TRUE)

  output$var_types_msg <- renderUI({
    if (!isTRUE(rv$current_ok) || length(rv$current_vars %||% character(0)) == 0) {
      return(tags$div(class = "smallnote", "Check current ASC folder first."))
    }
    if (is.null(rv$var_types)) return(tags$div(class = "bad", "Variable types not set yet."))
    tags$div(class = "ok", "Variable types saved. They will be applied to BOTH current and future.")
  })

  # ---- Check: Future scenario folders ----
  for (i in 1:5) {
    for (j in 1:4) {
      local({
        ii <- i
        jj <- j
        btn_id <- paste0("check_future_", ii, "_", jj)
        ui_id1 <- paste0("future_check_ui_", ii, "_", jj)
        ui_id2 <- paste0("future_vars_ui_", ii, "_", jj)
        key <- paste0("P", ii, "_S", jj)

        observeEvent(input[[btn_id]], {
          p <- input[[paste0("asc_future_dir_", ii, "_", jj)]] %||% ""

          if (!nzchar(p)) {
            rv$future_ok[[key]] <- FALSE
            rv$future_msg[[key]] <- "Future folder path is empty."
            rv$future_vars[[key]] <- character(0)
          } else if (!dir.exists(p)) {
            rv$future_ok[[key]] <- FALSE
            rv$future_msg[[key]] <- paste0("Folder not found: ", p)
            rv$future_vars[[key]] <- character(0)
          } else {
            r <- try(read_env_stack(p), silent = TRUE)
            if (inherits(r, "try-error")) {
              rv$future_ok[[key]] <- FALSE
              rv$future_msg[[key]] <- paste0("Failed to read rasters: ", as.character(r))
              rv$future_vars[[key]] <- character(0)
            } else {
              rv$future_ok[[key]] <- TRUE
              rv$future_msg[[key]] <- paste0("OK. Read ", terra::nlyr(r), " layers.")
              rv$future_vars[[key]] <- names(r)

              if (isTRUE(rv$current_ok) && length(rv$current_vars) > 0) {
                miss <- setdiff(rv$current_vars, rv$future_vars[[key]])
                extra <- setdiff(rv$future_vars[[key]], rv$current_vars)
                if (length(miss) > 0 || length(extra) > 0) {
                  rv$future_msg[[key]] <- paste0(
                    rv$future_msg[[key]],
                    " (WARNING: var mismatch vs current; missing: ",
                    if (length(miss) > 0) paste(miss, collapse = ", ") else "none",
                    "; extra: ",
                    if (length(extra) > 0) paste(extra, collapse = ", ") else "none",
                    ")"
                  )
                }
              }
            }
          }
        }, ignoreInit = TRUE)

        output[[ui_id1]] <- renderUI({
          ok <- isTRUE(rv$future_ok[[key]])
          msg <- rv$future_msg[[key]] %||% "Not checked yet."
          cls <- if (ok) "ok" else "bad"
          tags$div(class = cls, msg)
        })

        output[[ui_id2]] <- renderUI({
          ok <- isTRUE(rv$future_ok[[key]])
          vars <- rv$future_vars[[key]] %||% character(0)
          if (!ok || length(vars) == 0) return(NULL)
          tags$div(tags$strong("Detected variables:"), tags$pre(paste(vars, collapse = "\n")))
        })
      })
    }
  }

  # ---- Run job (background Rscript) ----
  # Arka planda calistir. Windows'ta '&' calismadigi icin platforma gore ayrilir.
  # ---------------------------------------------------------------------
  # Rscript'i PATH'e GUVENMEDEN bul.
  # Windows'ta R'in bin klasoru cogu kurulumda PATH'e eklenmez; bu durumda
  # system2("Rscript", ...) sessizce basarisiz olur, log bos kalir ve is
  # "STARTED"ta asili kalmis gibi gorunur. R.home() her zaman dogrudur.
  # ---------------------------------------------------------------------
  find_rscript <- function() {
    exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    cands <- c(
      file.path(R.home("bin"), exe),
      file.path(R.home(), "bin", exe),
      file.path(R.home(), "bin", "x64", exe),
      unname(Sys.which("Rscript"))
    )
    cands <- cands[nzchar(cands)]
    hit <- cands[file.exists(cands)]
    if (length(hit) > 0) normalizePath(hit[1], winslash = "/") else NA_character_
  }

  start_job <- function(cfg_path, log_path, script_name = "run_biomod.R") {
    run_script <- file.path(getwd(), script_name)
    if (!file.exists(run_script)) {
      showNotification(paste0(script_name, " bulunamadi: ", run_script),
                       type = "error", duration = NULL)
      return(invisible(FALSE))
    }
    run_script <- normalizePath(run_script)
    cfg_p <- normalizePath(cfg_path)
    log_p <- normalizePath(log_path, mustWork = FALSE)

    rscript <- find_rscript()
    if (is.na(rscript)) {
      showNotification(
        paste0("Rscript not found. Add R's bin folder to PATH, or reinstall R. ",
               "Looked in: ", R.home("bin")),
        type = "error", duration = NULL)
      try(writeLines("FAILED: Rscript executable not found", 
                     file.path(dirname(cfg_p), "status.txt")), silent = TRUE)
      return(invisible(FALSE))
    }

    # stdout ve stderr AYNI dosyaya verilir; system2 bunu "2>&1" olarak kurar,
    # boylece Windows'ta iki ayri yonlendirmenin dosyayi kilitlemesi olusmaz.
    ok <- tryCatch({
      system2(rscript, c(shQuote(run_script), shQuote(cfg_p)),
              stdout = log_p, stderr = log_p, wait = FALSE)
      TRUE
    }, error = function(e) {
      showNotification(paste0("Could not start the job: ", conditionMessage(e)),
                       type = "error", duration = NULL)
      FALSE
    })
    if (!isTRUE(ok)) return(invisible(FALSE))

    # Is gercekten basladi mi? 6 sn icinde log dosyasi olusmazsa uyar.
    later_check <- function() {
      if (!file.exists(log_p) || file.info(log_p)$size == 0) {
        showNotification(
          paste0("No output after 6 s. Check runs/", basename(dirname(cfg_p)),
                 "/run.log and confirm that Rscript can start."),
          type = "warning", duration = 15)
      }
    }
    try(shiny::observe({
      invalidateLater(6000, session)
      isolate(later_check())
    }), silent = TRUE)

    invisible(TRUE)
  }

  progress_map <- c(
    "STARTED" = 5,
    "LOADING DATA" = 10,
    "FORMATTING DATA" = 20,
    "MODELING" = 50,
    "SAVING DEFAULT BIOMOD PLOTS" = 70,
    "RESPONSE CURVES (models)" = 75,
    "EXPORTING METRICS (model)" = 80,
    "ENSEMBLE" = 85,
    "EXPORTING VAR IMPORTANCE" = 86,
    "RESPONSE CURVES (ensemble)" = 88,
    "EXPORTING METRICS (ensemble)" = 90,
    "PROJECTION" = 95,
    "DONE" = 100
  )

  ALL_METRICS <- c("AUCroc","ROC","AUC","TSS","BOYCE",
                   "RMSE","MSE","MAE","Max_error","Rsquared","Rsquared_aj",
                   "Accuracy","F1","Recall","Precision")

  summarize_eval <- function(csv_path) {
    if (!file.exists(csv_path)) return(NULL)
    df <- safe_read_csv(csv_path)
    if (is.null(df) || nrow(df) == 0) return(NULL)

    # --- YENI: run_biomod_cont.R duzenli (tidy) format yazar ---
    if (all(c("algo", "metric", "dataset", "value") %in% names(df))) {
      ds <- if (any(df$dataset == "validation" & is.finite(df$value))) "validation" else "calibration"
      sub <- df[df$dataset == ds & is.finite(df$value), , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
      agg <- aggregate(value ~ algo + metric, data = sub, FUN = mean, na.rm = TRUE)
      wide <- reshape(agg, idvar = "algo", timevar = "metric", direction = "wide")
      names(wide) <- sub("^value\\.", "", names(wide))
      names(wide)[names(wide) == "algo"] <- "model"
      wide$dataset <- ds
      return(wide[order(wide$model), , drop = FALSE])
    }

    if (!("value" %in% names(df))) {
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      if (length(num_cols) == 0) return(data.frame(message = "eval CSV found but no numeric value column."))
      names(df)[names(df) == num_cols[1]] <- "value"
    }

    metric_col <- NULL
    cand <- names(df)[vapply(df, function(x) any(as.character(x) %in% ALL_METRICS), logical(1))]
    if (length(cand) > 0) metric_col <- cand[1]

    model_col <- NULL
    algo_pat <- "^(ANN|CTA|DNN|FDA|GAM|GBM|GLM|MARS|MAXENT|SRE|RFd|RF|MAXNET|XGBOOST)$"
    cand2 <- names(df)[vapply(df, function(x) any(grepl(algo_pat, as.character(x))), logical(1))]
    if (length(cand2) > 0) model_col <- cand2[1]

    if (is.null(metric_col) || is.null(model_col)) {
      return(data.frame(message = "eval CSV found but columns not recognized."))
    }

    df$..metric <- as.character(df[[metric_col]])
    df$..model  <- as.character(df[[model_col]])
    df$..value  <- suppressWarnings(as.numeric(df[["value"]]))

    df <- df[df$..metric %in% ALL_METRICS, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    agg <- aggregate(..value ~ ..model + ..metric, data = df, FUN = mean, na.rm = TRUE)
    wide <- reshape(agg, idvar = "..model", timevar = "..metric", direction = "wide")
    names(wide) <- gsub("^\\.\\.value\\.", "", names(wide))
    names(wide)[names(wide) == "..model"] <- "model"
    wide[order(wide$model), , drop = FALSE]
  }

  # ---- candidates (fixed names remain; Image viewer uses dynamic list too) ----
  plot_candidates <- function(res_dir) {
    pdir <- file.path(res_dir, "plots")
    c(
      eval_mean = file.path(pdir, "biomod_eval_mean.png"),
      eval_boyce = file.path(pdir, "biomod_eval_boyce.png"),
      eval_box  = file.path(pdir, "biomod_eval_boxplot.png"),
      rc_models = file.path(pdir, "models_response_curves.png"),
      rc_ens    = file.path(pdir, "ensemble_response_curves.png"),
      varimp    = file.path(pdir, "var_importance_ensemble_avg.png"),
      rc_thr    = file.path(pdir, "ensemble_response_curves_threshold.png")
    )
  }

  raster_candidates <- function(res_dir) {
    rdir <- file.path(res_dir, "rasters")
    if (!dir.exists(rdir)) return(character(0))
    list.files(rdir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  }

  # ---- IMPORTANT: serve results/ as static ----
  ensure_results_resource <- function() {
    res_root <- normalizePath("results", winslash = "/", mustWork = FALSE)
    if (!dir.exists(res_root)) dir.create(res_root, recursive = TRUE, showWarnings = FALSE)
    addResourcePath("results", res_root)
  }
  ensure_results_resource()

  # ---- helper: file mtime-based cache-bust ----
  cache_ver <- function(path) {
    if (!file.exists(path)) return(as.integer(Sys.time()))
    as.integer(as.POSIXct(file.info(path)$mtime))
  }

  # =========================================================
  # \u2705 RASTER VIEWER FIX (refresh-safe)
  # =========================================================

  rasters_poll <- reactivePoll(
    intervalMillis = 1500,
    session = session,
    checkFunc = function() {
      if (is.null(rv$res_dir)) return(0)
      rdir <- file.path(rv$res_dir, "rasters")
      if (!dir.exists(rdir)) return(0)
      files <- list.files(rdir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
      if (length(files) == 0) return(0)
      max(file.info(files)$mtime, na.rm = TRUE)
    },
    valueFunc = function() {
      if (is.null(rv$res_dir)) return(character(0))
      raster_candidates(rv$res_dir)
    }
  )

  output$ui_raster_selectors <- renderUI({
    selectInput("raster_pick", "Select raster (.tif)", choices = character(0))
  })

  observe({
    ras <- rasters_poll()

    if (length(ras) == 0) {
      updateSelectInput(session, "raster_pick", choices = character(0), selected = character(0))
      return()
    }

    labs <- basename(ras)
    choices <- setNames(ras, labs)

    sel <- rv$selected_raster %||% input$raster_pick %||% ras[1]
    if (!nzchar(sel) || !(sel %in% ras)) sel <- ras[1]

    updateSelectInput(session, "raster_pick", choices = choices, selected = sel)
  })

  observeEvent(input$raster_pick, {
    rv$selected_raster <- input$raster_pick
  }, ignoreInit = TRUE)

  output$raster_plot <- renderPlot({
    fp <- rv$selected_raster %||% input$raster_pick
    if (is.null(fp) || !nzchar(fp)) return(NULL)

    if (!file.exists(fp)) {
      plot.new(); text(0.5, 0.5, "Raster file not found.")
      return()
    }

    r <- try(terra::rast(fp), silent = TRUE)
    if (inherits(r, "try-error")) {
      plot.new(); text(0.5, 0.5, "Failed to read raster.")
      return()
    }

    if (terra::nlyr(r) > 1) r <- r[[1]]
    terra::plot(r, main = basename(fp))
  })

  output$raster_info <- renderText({
    fp <- rv$selected_raster %||% input$raster_pick
    if (is.null(fp) || !nzchar(fp)) return("No raster selected yet.")
    if (!file.exists(fp)) return(paste0("Raster not found: ", fp))

    r <- try(terra::rast(fp), silent = TRUE)
    if (inherits(r, "try-error")) return("Could not read raster.")

    crs_txt <- try(as.character(terra::crs(r))[1], silent = TRUE)
    if (inherits(crs_txt, "try-error") || is.na(crs_txt) || !nzchar(crs_txt)) crs_txt <- "NA"

    ex_txt <- paste(as.vector(terra::ext(r)), collapse = ", ")

    paste(
      paste0("File: ", fp),
      paste0("Layers: ", terra::nlyr(r)),
      paste0("Extent: ", ex_txt),
      paste0("CRS: ", crs_txt),
      sep = "\n"
    )
  })

  # =========================================================
  # do_refresh (NO raster selector re-render here!)
  # =========================================================
  do_refresh <- function() {
    req(rv$run_dir)

    status_path <- file.path(rv$run_dir, "status.txt")
    st_lines <- safe_read_lines(status_path)
    st <- if (length(st_lines) > 0) paste(st_lines, collapse = "\n") else "-"
    last_line <- if (length(st_lines) > 0) st_lines[length(st_lines)] else "-"

    output$status <- renderText(paste0("run_id: ", rv$run_id, "\nstatus: ", st))

    output$logtail <- renderText({
      # Once KENDI gunlugumuz (progress.log): her satirda acilip kapandigi icin
      # Windows'ta da aninda okunabilir. run.log kabuk yonlendirmesiyle yazilir;
      # tamponlama ve dosya kilidi yuzunden is bitene kadar bos gorunebilir.
      prog <- if (!is.null(rv$run_dir)) file.path(rv$run_dir, "progress.log") else NULL
      ll <- safe_read_lines(prog, 180)
      src <- "progress.log"

      if (length(ll) == 0) {
        ll <- safe_read_lines(rv$log_file, 180)
        src <- "run.log"
      }
      if (length(ll) == 0) {
        st <- safe_read_lines(if (!is.null(rv$run_dir))
                                file.path(rv$run_dir, "status.txt") else NULL)
        if (length(st) > 0) {
          paste0("Status: ", paste(st, collapse = " "),
                 paste0("\n\n", tr_("job_running", lang()), "\n"),
                 "On Windows run.log stays buffered until the job ends; progress.log\n",
                 "appears as soon as the first message is produced.")
        } else tr_("no_log", lang())
      } else {
        paste0("[", src, "]\n", paste(ll, collapse = "\n"))
      }
    })

    model_csv <- file.path(rv$res_dir, "eval_model.csv")
    ens_csv   <- file.path(rv$res_dir, "eval_ensemble.csv")

    cand_csv <- file.path(rv$res_dir, "ensemble_candidates_after_thresholds.csv")
    sel_csv  <- file.path(rv$res_dir, "ensemble_selected_models.csv")
    sum_csv  <- file.path(rv$res_dir, "auc_tss_summary.csv")
    cut_csv  <- file.path(rv$res_dir, "ensemble_maxTSS_cutoffs.csv")
    brk_csv  <- file.path(rv$res_dir, "ensemble_class_breaks.csv")

    output$ens_candidates_tbl <- renderTable({
      if (!file.exists(cand_csv)) return(data.frame(message = "No ensemble_candidates_after_thresholds.csv yet."))
      df <- safe_read_csv(cand_csv)
      if ("model" %in% names(df)) df <- df[nzchar(df$model) & tolower(df$model) != "model", , drop = FALSE]
      df
    }, rownames = FALSE)

    output$ens_selected_tbl <- renderTable({
      if (!file.exists(sel_csv)) return(data.frame(message = "No ensemble_selected_models.csv yet."))
      df <- safe_read_csv(sel_csv)
      if ("model" %in% names(df)) df <- df[nzchar(df$model) & tolower(df$model) != "model", , drop = FALSE]
      df
    }, rownames = FALSE)

    output$auc_tss_summary_tbl <- renderTable({
      if (!file.exists(sum_csv)) return(data.frame(message = "No auc_tss_summary.csv yet."))
      safe_read_csv(sum_csv)
    }, rownames = FALSE)

    output$cutoffs_tbl <- renderTable({
      if (file.exists(cut_csv)) return(safe_read_csv(cut_csv))
      if (file.exists(brk_csv)) return(safe_read_csv(brk_csv))
      data.frame(message = tr_("no_thresh_tbl", lang()))
    }, rownames = FALSE)

    output$metrics_model <- renderTable({
      x <- summarize_eval(model_csv)
      if (is.null(x)) return(data.frame(message = "No eval_model.csv yet."))
      x
    }, rownames = FALSE)

    output$metrics_ens <- renderTable({
      x <- summarize_eval(ens_csv)
      if (is.null(x)) return(data.frame(message = "No eval_ensemble.csv yet (ensemble may still be running)."))
      x
    }, rownames = FALSE)

    output$progress_ui <- renderUI({
      if (!file.exists(status_path)) return(NULL)
      pct <- 0
      for (k in names(progress_map)) {
        if (grepl(k, last_line, fixed = TRUE)) pct <- max(pct, progress_map[[k]])
      }
      tags$div(
        tags$div(class = "smallnote", paste0("Stage: ", last_line)),
        tags$div(
          class = "progress",
          tags$div(
            class = "progress-bar progress-bar-striped progress-bar-animated",
            role = "progressbar",
            style = paste0("width:", pct, "%"),
            paste0(pct, "%")
          )
        )
      )
    })

    # ----- Saved PNG plots list (cache bust added) -----
    output$ui_plot_images <- renderUI({
      req(rv$run_id)
      files <- plot_candidates(rv$res_dir)

      to_src <- function(abs_path) {
        rel <- file.path("results", rv$run_id, "plots", basename(abs_path))
        gsub("\\\\", "/", rel)
      }

      img_block <- function(path, title) {
        if (!file.exists(path)) return(NULL)
        ver <- cache_ver(path)
        tags$div(
          class = "boxpad",
          tags$strong(title),
          br(),
          tags$img(src = paste0(to_src(path), "?v=", ver), style = "width: 100%; height: auto;")
        )
      }

      tagList(
        img_block(files["eval_mean"], "AUC/TSS Mean (Calibration)"),
        img_block(files["eval_boyce"], "Boyce / TSS Mean (Calibration)"),
        img_block(files["eval_box"], "AUC/TSS Boxplot (Calibration)"),
        img_block(files["rc_models"], "Response curves (models)"),
        img_block(files["rc_ens"], "Response curves (ensemble)"),
        img_block(files["rc_thr"], tr_("rc_thr_caption", lang())),
        img_block(files["varimp"], "Variable Importance Plot (ensemble-avg)")
      )
    })

    # ----- Image viewer: choices auto from plots dir (dynamic) -----
    output$ui_img_choice <- renderUI({
      req(rv$run_id)
      plots_dir <- file.path(rv$res_dir, "plots")
      if (!dir.exists(plots_dir)) return(tags$div(class = "smallnote", "No plots folder yet. Run the model first."))

      pngs <- list.files(plots_dir, pattern = "\\.png$", ignore.case = TRUE, full.names = FALSE)
      if (length(pngs) == 0) {
        return(tags$div(class = "smallnote", "No PNG exported yet. Wait for the run to reach plotting stage."))
      }

      sel <- input$img_choice %||% pngs[1]
      if (!(sel %in% pngs)) sel <- pngs[1]

      selectInput("img_choice", label = NULL, choices = pngs, selected = sel)
    })

    output$img_viewer_ui <- renderUI({
      req(rv$run_id)
      plots_dir <- file.path(rv$res_dir, "plots")
      img_name <- input$img_choice %||% "var_importance_ensemble_avg.png"
      img_path <- file.path(plots_dir, img_name)

      if (!dir.exists(plots_dir)) return(tags$div(class = "smallnote", "No plots folder yet. Run the model first."))
      if (!file.exists(img_path)) {
        return(tags$div(class = "smallnote",
                        paste0("Image not found yet: ", img_name, " (will appear after run_biomod.R exports it).")))
      }

      rel <- file.path("results", rv$run_id, "plots", basename(img_path))
      rel <- gsub("\\\\", "/", rel)
      ver <- cache_ver(img_path)

      tags$div(
        class = "boxpad",
        tags$strong(img_name),
        br(),
        tags$img(src = paste0(rel, "?v=", ver), style = "width: 100%; height: auto;")
      )
    })
  }

  # ---- RUN ----
  observeEvent(input$run, {

    # ---- On kosul denetimi: Rscript baslatilmadan once ----
    problems <- character(0)
    dt  <- data_type()
    fam <- dt_family(dt)

    if (!isTRUE(rv$current_ok)) {
      problems <- c(problems, "Guncel raster klasoru dogrulanmadi ('Check current raster folder').")
    }
    if (!isTRUE(rv$csv_ok) && !isTRUE(input$ssdm_enable)) {
      problems <- c(problems, tr_("err_no_csv", lang()))
    }
    if (!nzchar(input$sp_name %||% "") && !isTRUE(input$ssdm_enable)) {
      problems <- c(problems, tr_("err_no_sp", lang()))
    }
    if (length(input$models %||% character(0)) == 0) {
      problems <- c(problems, "Hicbir model secilmedi.")
    }
    if (!file.exists(file.path(getwd(), "biomod_opts.R"))) {
      problems <- c(problems, tr_("err_no_opts", lang()))
    }

    ssdm_on <- isTRUE(input$ssdm_enable)

    if (ssdm_on) {
      # S-SDM modunda tek-tur denetimleri gecersizdir
      problems <- setdiff(problems, tr_("err_no_csv", lang()))
      problems <- problems[!grepl("CSV dogrulanmadi|CSV not verified|Tur adi bos|Species name is empty",
                                  problems)]
      if (!nzchar(input$ssdm_mat_file %||% "") || !nzchar(input$ssdm_xy_file %||% "")) {
        problems <- c(problems, tr_("err_ssdm_files", lang()))
      }
      if (is.null(rv$ssdm_long) || !isTRUE(rv$ssdm_long$ok)) {
        problems <- c(problems, tr_("err_ssdm_convert", lang()))
      } else {
        fk <- ssdm_filtered()
        if (is.null(fk) || fk$n_kept < 2) {
          problems <- c(problems, sprintf(
            "Species filter leaves fewer than 2 species (%s kept). Lower the minimum frequency.",
            if (is.null(fk)) "0" else as.character(fk$n_kept)))
        }
      }
      if (any(grepl("^MAXENT$", input$models %||% character(0))) && !isTRUE(rv$jar_ok)) {
        problems <- c(problems, tr_("err_no_jar", lang()))
      }
    } else

    # --- Bagimli degisken tipine ozgu denetimler ---
    if (TRUE) {
      # Koordinat sutunlari her iki tek-tur yolunda da zorunlu
      fxy <- xy_detect()
      if (!is.null(rv$csv_df) && (is.null(fxy) || !isTRUE(fxy$ok))) {
        problems <- c(problems, tr_("err_no_xy_pick", lang(),
                                    paste(names(rv$csv_df), collapse = ", ")))
      }
    }
    if (identical(fam, "binary")) {
      if (any(grepl("^MAXENT$", input$models %||% character(0))) && !isTRUE(rv$jar_ok)) {
        problems <- c(problems, tr_("err_no_jar", lang()))
      }
    } else {
      rc <- input$resp_col %||% ""
      if (!nzchar(rc)) {
        problems <- c(problems, "Yanit sutunu secilmedi (surekli/nominal modelleme icin zorunlu).")
      } else if (!is.null(rv$csv_df) && !(rc %in% names(rv$csv_df))) {
        problems <- c(problems, paste0("Yanit sutunu CSV'de bulunamadi: ", rc))
      }

      bad_models <- setdiff(input$models %||% character(0), ALLOWED_MODELS[[dt]])
      if (length(bad_models) > 0) {
        problems <- c(problems, paste0("'", dt, "' veri tipinde desteklenmeyen model(ler): ",
                                       paste(bad_models, collapse = ", ")))
      }
      if (identical(dt, "ordinal") && !is.null(rv$csv_df) && nzchar(input$resp_col %||% "")) {
        lv <- trimws(unlist(strsplit(input$ord_levels %||% "", ",")))
        lv <- lv[nzchar(lv)]
        if (length(lv) > 0) {
          obs <- unique(as.character(rv$csv_df[[input$resp_col]]))
          obs <- obs[!is.na(obs)]
          miss <- setdiff(obs, lv)
          if (length(miss) > 0) {
            problems <- c(problems, paste0("Ordinal siralamada tanimsiz sinif(lar): ",
                                           paste(miss, collapse = ", ")))
          }
        }
      }
    }
    if (!is.null(input$future_n) && input$future_n != "NA") {
      nf <- as.integer(input$future_n)
      for (i in seq_len(nf)) {
        nsc <- as.integer(input[[paste0("future_scen_n_", i)]] %||% "1")
        for (j in seq_len(nsc)) {
          d <- input[[paste0("asc_future_dir_", i, "_", j)]] %||% ""
          if (!nzchar(d)) {
            problems <- c(problems, paste0("Gelecek P", i, "_S", j, ": klasor yolu bos."))
          } else if (!dir.exists(d)) {
            problems <- c(problems, paste0("Gelecek P", i, "_S", j, ": klasor bulunamadi."))
          }
        }
      }
    }

    if (length(problems) > 0) {
      showModal(modalDialog(
        title = "Calistirma baslatilamadi",
        tags$ul(lapply(problems, tags$li)),
        easyClose = TRUE, footer = modalButton("Tamam")
      ))
      return(invisible(NULL))
    }

    run_id  <- make_run_id()
    run_dir <- file.path("runs", run_id)
    res_dir <- file.path("results", run_id)

    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(res_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(res_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)

    var_types <- rv$var_types
    if (!is.null(var_types) && !is.list(var_types)) var_types <- as.list(var_types)
    if (!is.null(var_types) && length(var_types) == 0) var_types <- NULL

    if (is.null(input$future_n) || input$future_n == "NA") {
      future_n_num <- NULL
      future_list  <- NULL
    } else {
      future_n_num <- as.integer(input$future_n)
      future_list <- vector("list", future_n_num)

      for (i in seq_len(future_n_num)) {
        period_label <- input[[paste0("period_label_", i)]] %||% paste0("Future", i)
        scen_n <- as.integer(input[[paste0("future_scen_n_", i)]] %||% "1")

        scenarios <- vector("list", scen_n)
        for (j in seq_len(scen_n)) {
          ssp_lab   <- input[[paste0("ssp_", i, "_", j)]] %||% paste0("Scenario", j)
          lab_over  <- input[[paste0("sc_label_", i, "_", j)]] %||% ""
          sc_label  <- if (nzchar(lab_over)) lab_over else ssp_lab
          asc_dir_f <- input[[paste0("asc_future_dir_", i, "_", j)]] %||% ""
          scenarios[[j]] <- list(label = sc_label, asc_dir = asc_dir_f)
        }
        future_list[[i]] <- list(period_label = period_label, scenarios = scenarios)
      }
    }

    # "anakaya=bedrock, yukselti=elevation" -> list(anakaya="bedrock", ...)
    var_labels_list <- NULL
    vl_raw <- trimws(unlist(strsplit(input$var_labels %||% "", ",")))
    vl_raw <- vl_raw[nzchar(vl_raw) & grepl("=", vl_raw)]
    if (length(vl_raw) > 0) {
      kv <- strsplit(vl_raw, "=", fixed = TRUE)
      keys <- trimws(vapply(kv, function(z) z[1], ""))
      vals <- trimws(vapply(kv, function(z) paste(z[-1], collapse = "="), ""))
      keep <- nzchar(keys) & nzchar(vals)
      if (any(keep)) var_labels_list <- as.list(setNames(vals[keep], keys[keep]))
    }

    ord_levels <- trimws(unlist(strsplit(input$ord_levels %||% "", ",")))
    ord_levels <- ord_levels[nzchar(ord_levels)]

    class_probs <- suppressWarnings(as.numeric(trimws(unlist(
      strsplit(input$class_probs %||% "0.25,0.5,0.75", ",")))))
    class_probs <- class_probs[is.finite(class_probs) & class_probs > 0 & class_probs < 1]
    if (length(class_probs) == 0) class_probs <- c(0.25, 0.5, 0.75)

    cfg <- list(
      run_id   = run_id,
      asc_dir  = input$asc_dir,
      occ_file = input$occ_file,
      sp_name  = input$sp_name,
      out_dir  = res_dir,

      # ---- BAGIMLI DEGISKEN ----
      resp_kind      = input$resp_kind %||% "binary",
      data_type      = dt,
      resp_col       = if (identical(fam, "binary")) NULL else (input$resp_col %||% NULL),
      pa_col         = if (identical(fam, "binary") && nzchar(input$pa_col %||% ""))
                          input$pa_col else NULL,
      ordinal_levels = if (length(ord_levels) > 0) as.list(ord_levels) else NULL,
      filter_raster  = if (identical(fam, "binary")) TRUE else isTRUE(input$filter_raster),

      var_types = var_types,

      pa_strategy = input$pa_strategy,
      pa_rep      = as.integer(input$pa_rep %||% 3L),
      pa_n        = as.integer(input$pa_n %||% 10000L),
      pa_dist_min = as.integer(input$pa_dist_min %||% 0L),
      pa_dist_max = as.integer(input$pa_dist_max %||% 0L),

      cv = list(strategy = "kfold", k = 10, rep = as.integer(input$kfold_rep), perc = 70),

      cpu_n = as.integer(input$cpu_n %||% 1L),

      models  = as.character(input$models %||% character(0)),
      em_algo = as.character(input$em_algo %||% character(0)),

      # ---- ENSEMBLE ESIKLERI ----
      opt_strategy = input$opt_strategy %||% "adaptive",
      ens_min_auc  = as.numeric(input$ens_min_auc  %||% 0.70),
      ens_min_tss  = as.numeric(input$ens_min_tss  %||% 0.40),
      ens_min_r2   = as.numeric(input$ens_min_r2   %||% 0.30),
      ens_rmse_tol = as.numeric(input$ens_rmse_tol %||% 0),
      ens_min_acc  = as.numeric(input$ens_min_acc  %||% 0.50),
      ens_min_f1   = as.numeric(input$ens_min_f1   %||% 0.40),
      ens_near_one = as.numeric(input$ens_near_one %||% 0.999),
      x_col        = if (nzchar(input$x_col %||% "")) input$x_col else NULL,
      y_col        = if (nzchar(input$y_col %||% "")) input$y_col else NULL,

      # ---- PROJEKSIYON ----
      proj_digits    = as.integer(input$proj_digits %||% 3L),
      export_classes = isTRUE(input$export_classes),
      class_probs    = as.list(class_probs),

      out_scale  = input$out_scale %||% "0_1",
      rc_em_pick = input$rc_em_pick %||% "EMwmean",
      var_labels = var_labels_list,

      ssdm_enable   = isTRUE(input$ssdm_enable),
      ssdm_min_freq = as.numeric(input$ssdm_min_freq %||% 0),
      ssdm_pssdm    = isTRUE(input$ssdm_pssdm),
      ssdm_change   = isTRUE(input$ssdm_change),
      ssdm_long_file = if (isTRUE(input$ssdm_enable))
        file.path(res_dir, "ssdm_long_format.csv") else NULL,

      maxent_jar = input$maxent_jar,

      future_n = future_n_num,
      future   = future_list
    )

    if (length(cfg$models) == 0) cfg$models <- DEFAULT_MODELS[[dt]]

    cfg_path <- file.path(run_dir, "config.json")
    log_path <- file.path(run_dir, "run.log")
    write_json(cfg, cfg_path, auto_unbox = TRUE, pretty = TRUE)

    rv$run_id   <- run_id
    rv$run_dir  <- run_dir
    rv$res_dir  <- res_dir
    rv$log_file <- log_path

    writeLines("STARTED", file.path(run_dir, "status.txt"))

    if (isTRUE(input$ssdm_enable)) {
      # Donusturulmus tabloyu diske yaz: betik onu okur, calistirma yeniden uretilebilir olur
      lf <- rv$ssdm_long
      fk <- ssdm_filtered()
      d <- if (!is.null(fk)) fk$data else lf$data
      if (isTRUE(input$ssdm_latlon)) {
        rc <- tryCatch(terra::crs(terra::rast(list.files(input$asc_dir, full.names = TRUE,
                        pattern = "\\.(asc|tif|tiff|grd)$", ignore.case = TRUE)[1])),
                       error = function(e) NULL)
        if (!is.null(rc) && nzchar(rc)) d <- ssdm_reproject(d, 4326, rc)
      }
      utils::write.csv(d, file.path(res_dir, "ssdm_long_format.csv"), row.names = FALSE)
      if (!is.null(fk)) {
        utils::write.csv(fk$per_species, file.path(res_dir, "ssdm_species_counts.csv"),
                         row.names = FALSE)
        if (fk$n_dropped > 0) {
          utils::write.csv(fk$dropped, file.path(res_dir, "skipped_species.csv"),
                           row.names = FALSE)
        }
      } else {
        utils::write.csv(lf$per_species, file.path(res_dir, "ssdm_species_counts.csv"),
                         row.names = FALSE)
      }
      script_name <- "run_ssdm.R"
      dt <- "binary (S-SDM)"
    } else {
      script_name <- if (identical(fam, "binary")) "run_biomod.R" else "run_biomod_cont.R"
    }
    showNotification(tr_("msg_running", lang(), script_name, dt),
                     type = "message", duration = 8)
    start_job(cfg_path, log_path, script_name)

    do_refresh()
  }, ignoreInit = TRUE)

  observeEvent(input$refresh, {
    if (!is.null(rv$run_dir)) do_refresh()
  }, ignoreInit = TRUE)

  auto <- reactiveTimer(5000, session = session)
  observe({
    auto()
    if (!is.null(rv$run_dir)) do_refresh()
  })

  output$status <- renderText("Idle.")
  output$logtail <- renderText("-")
  output$progress_ui <- renderUI(NULL)
  output$metrics_model <- renderTable(data.frame(message = "Run a model, then metrics will appear."), rownames = FALSE)
  output$metrics_ens <- renderTable(data.frame(message = "Run a model, then metrics will appear."), rownames = FALSE)
}

shinyApp(ui, server)