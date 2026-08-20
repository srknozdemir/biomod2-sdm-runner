# biomod_opts.R
# =============================================================================
# biomod2 (>= 4.2) icin LITERATUR TEMELLI modelleme secenekleri motoru.
# run_biomod.R (binary) ve run_biomod_cont.R (surekli/nominal) tarafindan
# source() edilir. Tek kaynak => iki hat arasinda ayar farki olusamaz.
#
# ---------------------------------------------------------------------------
# PARAMETRE GEREKCELERI (kaynaklar)
#
#  GBM / BRT   Elith J, Leathwick JR, Hastie T (2008) A working guide to boosted
#              regression trees. Journal of Animal Ecology 77:802-813.
#              -> yavas ogrenme orani (lr) + >1000 agac; agac karmasikligi (tc)
#                 orneklem buyuklugu ile olceklenir; bag.fraction 0.5-0.75.
#
#  RF          Breiman L (2001) Random forests. Machine Learning 45:5-32.
#              Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G (2021)
#              Modelling species presence-only data with random forests.
#              Ecography 44:1731-1742.
#              -> mtry = sqrt(p) siniflandirmada, p/3 regresyonda;
#                 dengesiz PA verisinde asagi-orneklemeli RF (biomod2: RFd).
#
#  MAXENT /    Phillips SJ, Dudik M (2008) Modeling of species distributions with
#  MAXNET      Maxent: new extensions and a comprehensive evaluation.
#              Ecography 31:161-175.  -> "auto features" esikleri.
#              Merow C, Smith MJ, Silander JA (2013) A practical guide to MaxEnt.
#              Ecography 36:1058-1069.
#              Radosavljevic A, Anderson RP (2014) Making better Maxent models of
#              species distributions: complexity, overfitting and evaluation.
#              Journal of Biogeography 41:629-643.
#              -> beta (regularization) > 1, aktarilabilirlik (transferability)
#                 gerektiren gelecek projeksiyonlarinda daha duz tepki egrileri.
#
#  GAM         Wood SN (2011) Fast stable restricted maximum likelihood ...
#              JRSS-B 73:3-36.  -> REML, GCV'ye gore daha az az-duzlestirme.
#              Guisan A, Edwards TC, Hastie T (2002) Ecological Modelling 157:89-100.
#              Merow C et al. (2014) Ecography 37:1267-1281.
#              -> dusuk k (temel fonksiyon sayisi) + select=TRUE + gamma=1.4
#                 asiri uyumu sinirlar.
#
#  GLM         Austin M (2007) Ecological Modelling 200:1-19.
#              -> tur-cevre iliskileri tipik olarak tek tepeli; kuadratik terim.
#
#  MARS        Leathwick JR, Elith J, Hastie T (2006) Ecological Modelling
#              199:188-196.  -> degree <= 2, GCV ile budama.
#
#  PA/agirlik  Barbet-Massin M et al. (2012) Selecting pseudo-absences for species
#              distribution models. Methods in Ecology and Evolution 3:327-338.
#              -> prevalence = 0.5 (varlik/yokluk esit agirlik).
#
#  XGBOOST     Chen T, Guestrin C (2016) XGBoost: A scalable tree boosting
#              system. KDD '16.  -> dusuk eta + cok tur; satir/sutun altornekleme.
# =============================================================================

BMOPT_VERSION <- "1.1"

# Binary veri tipinde biomod2'nin destekledigi algoritmalar.
# run_ssdm.R bu listeyi model havuzunu suzmek icin kullanir.
ALLOWED_MODELS_BINARY <- c("ANN","CTA","DNN","FDA","GAM","GBM","GLM","MARS",
                           "MAXENT","MAXNET","RF","RFd","SRE","XGBOOST")

# ---------------------------------------------------------------------------
# GRAFIK CIHAZI (Windows / macOS / Linux)
#
# png() varsayilan tipi platforma gore degisir: Windows'ta "windows" (GDI),
# macOS'ta "quartz", Linux'ta "cairo" veya "Xlib". Bassiz (headless) bir Linux
# sunucusunda X11 yoksa ve R cairo'suz derlenmisse varsayilan png() COKER ve
# hicbir grafik uretilmez. Bu yuzden once platform varsayilani denenir,
# calismazsa sirayla alternatiflere gecilir.
# ---------------------------------------------------------------------------
bmopt_open_png <- function(file, width = 1600, height = 1200, res = 150) {
  attempts <- list(
    function() grDevices::png(file, width = width, height = height, res = res)
  )
  if (isTRUE(unname(capabilities("cairo")))) {
    attempts[[length(attempts) + 1]] <-
      function() grDevices::png(file, width = width, height = height, res = res,
                                type = "cairo")
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    attempts[[length(attempts) + 1]] <-
      function() ragg::agg_png(file, width = width, height = height, res = res)
  }
  for (f in attempts) {
    ok <- tryCatch({ f(); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(TRUE)
  }
  warning("No usable graphics device found. On a headless Linux server install ",
          "cairo support or the 'ragg' package.")
  FALSE
}

# ---------------------------------------------------------------------------
# GUVENLI PNG YAZIMI
#
# ONEMLI: cihaz, dosya boyutu denetlenmeden ONCE kapatilmalidir. on.exit ile
# kapatilirsa dev.off() fonksiyon DONERKEN calisir; boyut denetimi ise daha
# once yapilir ve dosya heniz diske yazilmadigi icin grafik sessizce atilir.
# ---------------------------------------------------------------------------
safe_png <- function(file, expr, width = 1600, height = 1200, res = 150) {
  tmp <- paste0(file, ".tmp.png")
  if (!bmopt_open_png(tmp, width = width, height = height, res = res)) return(FALSE)
  dev_id <- grDevices::dev.cur()

  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    message("Plot failed (", basename(file), "): ", conditionMessage(e)); FALSE
  })

  try(grDevices::dev.off(which = dev_id), silent = TRUE)   # ONCE kapat

  if (isTRUE(ok) && file.exists(tmp) && file.info(tmp)$size > 2000) {
    file.rename(tmp, file)
    return(TRUE)
  }
  if (file.exists(tmp)) unlink(tmp)
  FALSE
}

# ---------------------------------------------------------------------------
# KOORDINAT SUTUNU TESPITI
# CSV'lerde koordinat sutunu her zaman "x"/"y" adinda olmaz. Sabit isim
# beklemek yerine yaygin adlari taniyoruz; kullanici arayuzden yine de
# elle secebilir (cfg$x_col / cfg$y_col).
# ---------------------------------------------------------------------------
XCOL_PAT <- paste0("^(x|x_?coord|coord_?x|lon|lng|long|longitude|boylam|",
                   "easting|east|utm_?x|point_?x|decimal_?longitude|dec_?long)$")
YCOL_PAT <- paste0("^(y|y_?coord|coord_?y|lat|latitude|enlem|",
                   "northing|north|utm_?y|point_?y|decimal_?latitude|dec_?lat)$")

find_xy_cols <- function(df, x_col = NULL, y_col = NULL) {
  nm <- names(df)
  lo <- tolower(trimws(nm))

  use <- function(given, pat, exact) {
    # 1) Arayuzden elle secilen ad her zaman oncelikli
    if (!is.null(given) && nzchar(given) && given %in% nm) return(given)
    hit <- which(grepl(pat, lo))
    if (length(hit) == 0) return(NULL)
    # 2) Birden fazla aday varsa tam "x"/"y" tercih edilir
    if (length(hit) > 1) {
      ex <- which(lo == exact)
      if (length(ex) > 0) hit <- ex
    }
    nm[hit[1]]
  }
  x <- use(x_col, XCOL_PAT, "x")
  y <- use(y_col, YCOL_PAT, "y")

  list(x = x, y = y, ok = !is.null(x) && !is.null(y), columns = nm)
}

# Koordinatlari cikarir; bulunamazsa ACIKLAYICI hata verir.
get_xy_matrix <- function(df, x_col = NULL, y_col = NULL, what = "CSV") {
  f <- find_xy_cols(df, x_col, y_col)
  if (!isTRUE(f$ok)) {
    stop(what, ": coordinate columns not found. Detected columns: ",
         paste(f$columns, collapse = ", "),
         ". Rename them to x / y, or pick them in the interface.")
  }
  m <- as.matrix(data.frame(x = suppressWarnings(as.numeric(df[[f$x]])),
                            y = suppressWarnings(as.numeric(df[[f$y]]))))
  storage.mode(m) <- "double"
  message("Coordinate columns: x = '", f$x, "', y = '", f$y, "'")
  attr(m, "x_col") <- f$x; attr(m, "y_col") <- f$y
  m
}

# --------------------------------------------------------------------------
# Orneklem buyuklugu sinifi
# binary  : varlik sayisi (n_pres)
# nonbinary: toplam gozlem sayisi
# --------------------------------------------------------------------------
bmopt_nclass <- function(n) {
  n <- as.integer(n)
  if (is.na(n)) return("moderate")
  if (n <  50)  return("small")
  if (n < 250)  return("moderate")
  if (n < 1000) return("large")
  "verylarge"
}

# --------------------------------------------------------------------------
# MAXENT "auto features" (Phillips & Dudik 2008)
#   n <  10 : yalniz dogrusal
#   10-14   : dogrusal + kuadratik
#   15-79   : + hinge
#   n >= 80 : + product + threshold
# --------------------------------------------------------------------------
bmopt_maxent_features <- function(n_pres) {
  n <- as.integer(n_pres)
  list(
    linear    = TRUE,
    quadratic = n >= 10,
    hinge     = n >= 15,
    product   = n >= 80,
    threshold = n >= 80
  )
}

# --------------------------------------------------------------------------
# Regularization carpani (beta).
# Varsayilan 1 asiri uyum egilimlidir; gelecek iklim projeksiyonu gibi
# AKTARILABILIRLIK gerektiren kullanimlarda literatur beta > 1 onerir
# (Radosavljevic & Anderson 2014; Merow et al. 2013).
# --------------------------------------------------------------------------
bmopt_beta <- function(n_pres, transfer = TRUE) {
  n <- as.integer(n_pres)
  b <- if (n < 30) 3 else if (n < 100) 2 else 1.5
  if (!isTRUE(transfer)) b <- max(1, b - 0.5)
  b
}

# --------------------------------------------------------------------------
# Formul uretimi (GLM kuadratik, GAM s_smoother + sinirli k)
# bm_MakeFormula yoksa NULL doner; cagiran taraf sessizce atlar.
# --------------------------------------------------------------------------
bmopt_formula <- function(bm_format, type = "quadratic", interaction.level = 0, k = NULL) {
  if (is.null(bm_format)) return(NULL)
  if (!exists("bm_MakeFormula", mode = "function")) return(NULL)

  sp <- tryCatch(bm_format@sp.name, error = function(e) NULL)
  ev <- tryCatch(utils::head(bm_format@data.env.var), error = function(e) NULL)
  if (is.null(sp) || is.null(ev)) return(NULL)

  a <- list(resp.name = sp, expl.var = ev, type = type,
            interaction.level = as.integer(interaction.level))
  if (!is.null(k)) a$k <- as.integer(k)

  out <- tryCatch(do.call(bm_MakeFormula, a), error = function(e) NULL)
  if (is.null(out) && !is.null(k)) {
    a$k <- NULL
    out <- tryCatch(do.call(bm_MakeFormula, a), error = function(e) NULL)
  }
  out
}

# --------------------------------------------------------------------------
# ANA PARAMETRE IZGARASI
#   data_type : binary | count | abundance | relative | multiclass | ordinal
#   p         : ongorucu (katman) sayisi
#   n         : etkin orneklem (binary -> varlik sayisi)
# --------------------------------------------------------------------------
bmopt_params <- function(data_type, p, n, bm_format = NULL,
                         maxent_jar = NULL, transfer = TRUE) {

  p <- max(1L, as.integer(p))
  n <- max(10L, as.integer(n))
  nc <- bmopt_nclass(n)

  binary <- identical(data_type, "binary")
  quant  <- data_type %in% c("count", "abundance", "relative")
  qual   <- data_type %in% c("multiclass", "ordinal")

  out <- list()

  # ---- GBM (Elith et al. 2008) --------------------------------------------
  # lr, n.trees > 1000 olacak sekilde secilir; tc orneklemle olceklenir.
  gbm_lr <- switch(nc, small = 0.001, moderate = 0.005, large = 0.01, verylarge = 0.01)
  gbm_nt <- switch(nc, small = 5000L, moderate = 5000L, large = 3000L, verylarge = 3000L)
  gbm_tc <- switch(nc, small = 2L,    moderate = 3L,    large = 4L,    verylarge = 5L)
  out$GBM <- list(
    n.trees           = gbm_nt,
    interaction.depth = gbm_tc,
    shrinkage         = gbm_lr,
    bag.fraction      = 0.5,
    n.minobsinnode    = max(2L, min(10L, as.integer(floor(n / 20)))),
    cv.folds          = 0L,          # capraz gecerleme biomod2 tarafinda
    keep.data         = FALSE,
    verbose           = FALSE
  )

  # ---- RF / RFd (Breiman 2001; Valavi et al. 2021) ------------------------
  out$RF <- list(
    ntree    = 1000L,
    mtry     = if (quant) max(1L, as.integer(floor(p / 3))) else max(1L, as.integer(floor(sqrt(p)))),
    nodesize = if (quant) 5L else 1L
  )
  out$RFd <- out$RF

  # ---- CTA (rpart) --------------------------------------------------------
  out$CTA <- list(
    control = list(
      xval      = 10L,               # rpart kendi CV budamasi
      minbucket = max(2L, as.integer(floor(n / 40))),
      minsplit  = max(5L, as.integer(floor(n / 20))),
      cp        = 0.001,
      maxdepth  = 25L
    )
  )

  # ---- GLM (Austin 2007) --------------------------------------------------
  # Tek tepeli tur-cevre iliskisi icin kuadratik terim; etkilesim yalnizca
  # bol orneklem + az degisken durumunda.
  glm_int <- if (nc %in% c("large", "verylarge") && p <= 12) 1L else 0L
  glm_form <- bmopt_formula(bm_format, type = "quadratic", interaction.level = glm_int)
  out$GLM <- list(control = stats::glm.control(maxit = 100))
  if (!is.null(glm_form)) out$GLM$formula <- glm_form

  # ---- GAM (Wood 2011; Merow et al. 2014) --------------------------------
  # Dusuk k + select=TRUE + gamma=1.4 => aktarilabilir, az-uyumlu duzlestirme.
  gam_k <- switch(nc, small = 3L, moderate = 4L, large = 5L, verylarge = 5L)
  gam_form <- bmopt_formula(bm_format, type = "s_smoother", interaction.level = 0, k = gam_k)
  out$GAM <- list(method = "REML", select = TRUE, gamma = 1.4)
  if (!is.null(gam_form)) out$GAM$formula <- gam_form

  # ---- MARS (Leathwick et al. 2006) --------------------------------------
  out$MARS <- list(
    degree  = if (nc == "small") 1L else 2L,
    nk      = min(200L, max(21L, as.integer(2 * p + 1))),
    penalty = if (nc == "small") 3 else 2,
    pmethod = "backward",
    nprune  = max(4L, min(25L, as.integer(floor(n / 8))))
  )

  # ---- ANN (nnet) --------------------------------------------------------
  if (binary) {
    out$ANN <- list(
      size     = max(2L, min(8L, as.integer(round(p / 2)))),
      decay    = if (nc %in% c("small", "moderate")) 0.1 else 0.01,
      rang     = 0.1,
      maxit    = if (nc == "small") 200L else 500L,
      MaxNWts  = max(10000L, as.integer(2000L * p)),
      trace    = FALSE
    )
  }

  # ---- MAXNET (Phillips & Dudik 2008; Radosavljevic & Anderson 2014) -----
  if (binary) {
    out$MAXNET <- list(regmult = bmopt_beta(n, transfer))
  }

  # ---- MAXENT (java) -----------------------------------------------------
  if (binary) {
    ft <- bmopt_maxent_features(n)
    out$MAXENT <- list(
      linear            = ft$linear,
      quadratic         = ft$quadratic,
      hinge             = ft$hinge,
      product           = ft$product,
      threshold         = ft$threshold,
      betamultiplier    = bmopt_beta(n, transfer),
      defaultprevalence = 0.5
    )
    if (!is.null(maxent_jar) && nzchar(maxent_jar) && file.exists(maxent_jar)) {
      # 4.2 serisinde gecerli; 4.3'te reddedilirse otomatik atlanir.
      out$MAXENT$path_to_maxent.jar <- dirname(normalizePath(maxent_jar))
    }
  }

  # ---- SRE ---------------------------------------------------------------
  if (binary) out$SRE <- list(quant = 0.025)   # %95 zarf

  # ---- FDA (mda) ---------------------------------------------------------
  if (!quant) out$FDA <- list(method = "mars")

  # ---- XGBOOST (Chen & Guestrin 2016) ------------------------------------
  out$XGBOOST <- list(
    nrounds          = if (nc == "small") 300L else 1000L,
    max_depth        = switch(nc, small = 2L, moderate = 3L, large = 4L, verylarge = 6L),
    eta              = if (nc %in% c("small", "moderate")) 0.01 else 0.05,
    subsample        = 0.75,
    colsample_bytree = 0.8,
    min_child_weight = max(1L, as.integer(floor(n / 100))),
    gamma            = 0,
    verbose          = 0L
  )

  out
}

# --------------------------------------------------------------------------
# SECENEK NESNESI KURULUMU
#
# Onemli: kurulan user.val, BIOMOD_Modeling'e OPT.user DEGIL,
# OPT.strategy="user.defined" + OPT.user.val + OPT.user.base ile verilir.
# Boylece veri seti anahtarlarini (_PA1_RUN1, _allData_allRun ...) capraz
# gecerleme kurgusuna gore biomod2 kendisi acar; surum farkindan kaynaklanan
# "names(OPT.user@options[[...]]@args.values) must be ..." hatasi olusmaz.
#
# Doner: list(mode, user.val, accepted, rejected, dropped_args)
# --------------------------------------------------------------------------
bmopt_build <- function(data_type, models, bm_format, p, n,
                        maxent_jar = NULL, strategy = "adaptive",
                        transfer = TRUE, verbose = TRUE) {

  say <- function(...) if (isTRUE(verbose)) message(...)

  if (!exists("bm_ModelingOptions", mode = "function")) {
    say("bm_ModelingOptions bulunamadi -> OPT.strategy='bigboss'.")
    return(list(mode = "bigboss", user.val = NULL,
                accepted = character(0), rejected = character(0),
                dropped_args = character(0)))
  }

  strategy <- tolower(as.character(strategy)[1])
  if (strategy %in% c("bigboss", "default", "tuned")) {
    say("Modelleme secenekleri: OPT.strategy='", strategy, "' (adaptif ayar uygulanmayacak).")
    return(list(mode = strategy, user.val = NULL,
                accepted = character(0), rejected = character(0),
                dropped_args = character(0)))
  }

  # 1) Secenek ADLARINI surumden ogren (RF.binary.randomForest.randomForest gibi)
  probe <- tryCatch(
    bm_ModelingOptions(data.type = data_type, models = models,
                       strategy = "bigboss", bm.format = bm_format),
    error = function(e) NULL
  )
  if (is.null(probe)) {
    probe <- tryCatch(
      bm_ModelingOptions(data.type = data_type, models = models, strategy = "bigboss"),
      error = function(e) NULL
    )
  }
  if (is.null(probe)) {
    say("bm_ModelingOptions calistirilamadi -> OPT.strategy='bigboss'.")
    return(list(mode = "bigboss", user.val = NULL,
                accepted = character(0), rejected = character(0),
                dropped_args = character(0)))
  }

  opt_names <- tryCatch(names(probe@options), error = function(e) character(0))
  if (length(opt_names) == 0) {
    return(list(mode = "bigboss", user.val = NULL,
                accepted = character(0), rejected = character(0),
                dropped_args = character(0)))
  }

  # Bir modelin veri seti anahtarlari (yedek yol)
  ds_keys <- function(nm) {
    k <- tryCatch(names(probe@options[[nm]]@args.values), error = function(e) NULL)
    if (is.null(k) || length(k) == 0) "_allData_allRun" else k
  }

  # user.val adayini dogrula
  try_opt <- function(uv) {
    tryCatch(
      bm_ModelingOptions(data.type = data_type, models = models,
                         strategy = "user.defined", user.val = uv,
                         user.base = "bigboss", bm.format = bm_format),
      error = function(e) NULL
    )
  }

  # Bir algoritma blogunu iki anahtarlama bicimiyle dener
  wrap_val <- function(nm, v, style) {
    if (identical(style, "all")) list('_for_all_datasets' = v)
    else stats::setNames(rep(list(v), length(ds_keys(nm))), ds_keys(nm))
  }

  grid <- bmopt_params(data_type, p = p, n = n, bm_format = bm_format,
                       maxent_jar = maxent_jar, transfer = transfer)

  user.val <- list()
  accepted <- character(0)
  rejected <- character(0)
  dropped_args <- character(0)

  for (algo in names(grid)) {
    v <- grid[[algo]]
    if (length(v) == 0) next

    nm <- opt_names[startsWith(opt_names, paste0(algo, "."))]
    if (length(nm) == 0) next          # bu veri tipinde model yok
    nm <- nm[1]

    placed <- FALSE
    for (style in c("all", "keys")) {
      trial <- user.val
      trial[[nm]] <- wrap_val(nm, v, style)
      if (!is.null(try_opt(trial))) {
        user.val <- trial; accepted <- c(accepted, nm); placed <- TRUE
        say("OPT [", algo, "] uygulandi (", length(v), " parametre).")
        break
      }
    }

    # Tum blok reddedildiyse: argumanlari TEK TEK ele, kabul edilenleri tut.
    if (!placed) {
      keep <- list()
      for (a in names(v)) {
        cand <- c(keep, v[a])
        ok <- FALSE
        for (style in c("all", "keys")) {
          trial <- user.val
          trial[[nm]] <- wrap_val(nm, cand, style)
          if (!is.null(try_opt(trial))) { ok <- TRUE; break }
        }
        if (ok) keep <- cand else dropped_args <- c(dropped_args, paste0(algo, "$", a))
      }
      if (length(keep) > 0) {
        for (style in c("all", "keys")) {
          trial <- user.val
          trial[[nm]] <- wrap_val(nm, keep, style)
          if (!is.null(try_opt(trial))) {
            user.val <- trial; accepted <- c(accepted, nm); placed <- TRUE
            say("OPT [", algo, "] kismen uygulandi (", length(keep), "/", length(v),
                " parametre).")
            break
          }
        }
      }
      if (!placed) {
        rejected <- c(rejected, nm)
        say("OPT [", algo, "] reddedildi -> bigboss degerlerinde birakildi.")
      }
    }
  }

  if (length(user.val) == 0) {
    return(list(mode = "bigboss", user.val = NULL,
                accepted = character(0), rejected = rejected,
                dropped_args = dropped_args))
  }

  list(mode = "user.defined", user.val = user.val,
       accepted = accepted, rejected = rejected, dropped_args = dropped_args)
}

# --------------------------------------------------------------------------
# BIOMOD_Modeling argumanlarina secenekleri yerlestirir.
# --------------------------------------------------------------------------
bmopt_inject <- function(mod_args, built) {
  mod_args$OPT.user <- NULL
  mod_args$OPT.user.val <- NULL
  mod_args$models.options <- NULL

  if (identical(built$mode, "user.defined")) {
    mod_args$OPT.strategy  <- "user.defined"
    mod_args$OPT.user.val  <- built$user.val
    mod_args$OPT.user.base <- "bigboss"
  } else {
    mod_args$OPT.strategy <- built$mode
  }
  mod_args
}

# --------------------------------------------------------------------------
# Uygulanan ayarlarin denetim kaydi (CSV)
# --------------------------------------------------------------------------
bmopt_report <- function(built, path, data_type, p, n) {
  flat <- list()
  if (!is.null(built$user.val)) {
    for (nm in names(built$user.val)) {
      block <- built$user.val[[nm]]
      v <- if (length(block) > 0) block[[1]] else list()
      for (a in names(v)) {
        val <- v[[a]]
        txt <- tryCatch({
          if (inherits(val, "formula")) paste(deparse(val), collapse = " ")
          else if (is.list(val)) paste(names(val), unlist(val), sep = "=", collapse = "; ")
          else paste(as.character(val), collapse = ", ")
        }, error = function(e) "<?>")
        flat[[length(flat) + 1]] <- data.frame(
          option_set = nm, parameter = a, value = txt, stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- if (length(flat) > 0) do.call(rbind, flat) else
    data.frame(option_set = NA_character_, parameter = NA_character_,
               value = NA_character_, stringsAsFactors = FALSE)

  df$mode          <- built$mode
  df$data_type     <- data_type
  df$n_predictors  <- p
  df$n_effective   <- n
  df$size_class    <- bmopt_nclass(n)
  df$engine        <- BMOPT_VERSION

  if (length(built$rejected) > 0) {
    df <- rbind(df, data.frame(
      option_set = built$rejected, parameter = "<REJECTED>", value = "bigboss",
      mode = built$mode, data_type = data_type, n_predictors = p,
      n_effective = n, size_class = bmopt_nclass(n), engine = BMOPT_VERSION,
      stringsAsFactors = FALSE))
  }
  if (length(built$dropped_args) > 0) {
    df <- rbind(df, data.frame(
      option_set = "<DROPPED_ARGS>", parameter = built$dropped_args, value = "bigboss",
      mode = built$mode, data_type = data_type, n_predictors = p,
      n_effective = n, size_class = bmopt_nclass(n), engine = BMOPT_VERSION,
      stringsAsFactors = FALSE))
  }

  tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, path, row.names = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

# --------------------------------------------------------------------------
# MAXENT java dosyasini biomod2'nin bulabilecegi yerlere yerlestirir.
# 4.2 : simulasyon klasoru / path_to_maxent.jar
# 4.3 : ENMevaluate -> dismo paketinin java klasoru
# --------------------------------------------------------------------------
bmopt_stage_maxent <- function(maxent_jar, verbose = TRUE) {
  say <- function(...) if (isTRUE(verbose)) message(...)
  if (is.null(maxent_jar) || !nzchar(maxent_jar) || !file.exists(maxent_jar)) {
    say("maxent.jar bulunamadi; MAXENT atlanacak.")
    return(invisible(FALSE))
  }
  src <- normalizePath(maxent_jar)
  ok <- FALSE

  dst1 <- file.path(getwd(), "maxent.jar")
  if (!file.exists(dst1) || !identical(normalizePath(dst1), src)) {
    if (isTRUE(try(file.copy(src, dst1, overwrite = TRUE), silent = TRUE))) {
      say("maxent.jar calisma dizinine kopyalandi: ", dst1); ok <- TRUE
    }
  } else ok <- TRUE

  if (requireNamespace("dismo", quietly = TRUE)) {
    jdir <- system.file("java", package = "dismo")
    if (nzchar(jdir) && dir.exists(jdir)) {
      dst2 <- file.path(jdir, "maxent.jar")
      if (!file.exists(dst2)) {
        if (isTRUE(try(file.copy(src, dst2), silent = TRUE))) {
          say("maxent.jar dismo/java klasorune kopyalandi: ", dst2); ok <- TRUE
        } else {
          say("dismo/java klasorune yazilamadi (yazma izni yok): ", jdir)
        }
      } else ok <- TRUE
    }
  }
  invisible(ok)
}
