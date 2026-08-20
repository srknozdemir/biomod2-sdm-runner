# check_code.R -- STATIK KOD DENETIMI
# "nesnesi bulunamadi" (no visible binding) hatalarini, kod CALISTIRILMADAN
# once yakalar. Bu hata sinifi parse denetiminden gecer, cunku R degiskenleri
# ancak calisma aninda arar.
# ONEMLI: yan etkili kod CALISTIRILMAZ; yalnizca "ad <- function/c/list"
# bicimindeki tanimlar degerlendirilir (install.packages() asla tetiklenmez).
if (!requireNamespace("codetools", quietly = TRUE)) stop("codetools gerekli")

FILES <- c("i18n.R","ssdm_convert.R","biomod_opts.R","plot_response_threshold.R",
           "run_biomod.R","run_biomod_cont.R","run_ssdm.R","app.R","install_packages.R")

KNOWN <- c("fluidPage","titlePanel","sidebarLayout","sidebarPanel","mainPanel","tabsetPanel",
  "tabPanel","fluidRow","column","selectInput","textInput","numericInput","checkboxInput",
  "checkboxGroupInput","actionButton","uiOutput","tableOutput","plotOutput","renderTable",
  "renderPlot","renderText","renderUI","verbatimTextOutput","textOutput","conditionalPanel",
  "tagList","HTML","hr","h3","h4","br","helpText","reactiveValues","reactive","observe",
  "observeEvent","isolate","invalidateLater","showNotification","updateTextInput",
  "updateNumericInput","updateSelectInput","updateCheckboxInput","updateActionButton",
  "addResourcePath","shinyApp","req","icon","wellPanel","downloadButton","fileInput",
  "shinyFilesButton","shinyDirButton","shinyFileChoose","shinyDirChoose",
  "parseFilePaths","parseDirPath","getVolumes",
  "BIOMOD_FormatingData","BIOMOD_Modeling","BIOMOD_EnsembleModeling","BIOMOD_Projection",
  "BIOMOD_EnsembleForecasting","bm_ModelingOptions","bm_MakeFormula","bm_PlotEvalMean",
  "bm_PlotEvalBoxplot","bm_PlotVarImpBoxplot","bm_PlotResponseCurves","bm_ModelAnalysis",
  "get_evaluations","get_predictions","get_variables_importance","get_built_models",
  "getModelsBuilt","rast","ggplot","aes","filter")

IGNORE <- c("input","output","session","rv","tags","volumes","lang",".data","cfg","env","occ")

# Dosyada HERHANGI bir yerde atanan TUM adlari topla (degerlendirmeden).
# Betik duzeyindeki degiskenler (status_path, rasters_dir ...) calisma aninda
# global ortamda bulunur; codetools bunlari goremedigi icin yanlis alarm verir.
# YALNIZCA BETIK DUZEYI (top-level) atamalari toplanir.
# Fonksiyon govdelerine GIRILMEZ: bir fonksiyonda atanan ad, baska bir
# fonksiyondan gorunmez. Kapsami genis tutmak, "rscript" tipi gercek
# hatalari eleyip denetimi ise yaramaz hale getirir.
collect_names <- function(expr, acc) {
  if (!is.call(expr)) return(acc)
  h <- as.character(expr[[1]])[1]

  if (h == "function") return(acc)          # fonksiyon govdesine GIRME

  if (h %in% c("<-", "=", "<<-") && length(expr) == 3 && is.name(expr[[2]])) {
    acc <- c(acc, as.character(expr[[2]]))
    rhs <- expr[[3]]
    if (is.call(rhs) && identical(as.character(rhs[[1]])[1], "function")) return(acc)
    return(collect_names(rhs, acc))
  }
  if (h == "for" && length(expr) >= 3 && is.name(expr[[2]]))
    acc <- c(acc, as.character(expr[[2]]))

  # yalnizca kontrol akisi bloklarina in
  if (h %in% c("{", "if", "for", "while", "repeat", "(")) {
    for (i in seq_along(expr))
      acc <- tryCatch({ s <- expr[[i]]; if (is.call(s)) collect_names(s, acc) else acc },
                      error = function(e) acc)
  }
  acc
}

# Bir ifade YAN ETKISIZ mi? Yalnizca sabit deger (sayi/metin/mantiksal) ve
# bunlardan olusan c()/list() cagrilari guvenlidir. Icinde baska bir fonksiyon
# cagrisi varsa DEGERLENDIRILMEZ.
# Gerekce: install_packages.R icinde "failed <- c(failed, install_if_missing(...))"
# gibi satirlar var; RHS bir c() cagrisi oldugu icin "sabit" sanilip
# calistiriliyor ve gercekten paket kurmaya girisiyordu.
is_pure <- function(x) {
  if (is.name(x)) return(TRUE)
  if (!is.call(x)) return(is.atomic(x) || is.null(x))
  h <- as.character(x[[1]])[1]
  if (!h %in% c("c", "list")) return(FALSE)
  for (i in seq_along(x)[-1]) {
    ok <- tryCatch(is_pure(x[[i]]), error = function(e) FALSE)
    if (!isTRUE(ok)) return(FALSE)
  }
  TRUE
}

# YALNIZCA UST DUZEY fonksiyonlari degerlendir.
# Ic ice tanimli fonksiyonlar (closure) dis kapsamdaki degiskenleri yakalar;
# onlari yalitilmis denetlemek yanlis alarm uretir. Ust duzey fonksiyonu
# denetlerken codetools ic kapsamlari zaten dogru cozer.
collect_defs <- function(expr, env) {
  if (!is.call(expr)) return(invisible(NULL))
  h <- as.character(expr[[1]])[1]
  if (h == "function") return(invisible(NULL))

  if (h %in% c("<-","=") && length(expr) == 3 && is.name(expr[[2]])) {
    rhs <- expr[[3]]
    if (is.call(rhs)) {
      h2 <- as.character(rhs[[1]])[1]
      if (identical(h2, "function")) {
        try(eval(expr, envir = env), silent = TRUE)
      } else if (h2 %in% c("c", "list") && isTRUE(is_pure(rhs))) {
        try(eval(expr, envir = env), silent = TRUE)
      }
    }
    return(invisible(NULL))
  }
  if (h %in% c("{","if","for","while","repeat","(")) {
    for (i in seq_along(expr))
      tryCatch({ s <- expr[[i]]; if (is.call(s)) collect_defs(s, env) }, error = function(e) NULL)
  }
  invisible(NULL)
}

total <- 0L; checked <- 0L
for (f in FILES) {
  if (!file.exists(f)) {
    cat("EKSIK DOSYA: ", f, "  <- bu dosyayi klasore kopyalayin\n", sep = "")
    total <- total + 1L; next
  }
  # keep.source = TRUE -> mesajlarda dosya:satir referansi cikar
  ex <- tryCatch(parse(f, encoding = "UTF-8", keep.source = TRUE),
                 error = function(e) e)
  if (inherits(ex, "error")) {
    cat("PARSE HATASI:", f, "-", conditionMessage(ex), "\n"); total <- total + 1L; next }

  # Dosyada atanan tum adlar -> yanlis alarmlari eler
  assigned <- character(0)
  for (e in ex) assigned <- collect_names(e, assigned)
  assigned <- unique(assigned)

  env <- new.env(parent = globalenv())
  for (k in KNOWN)  assign(k, function(...) NULL, envir = env)
  for (k in IGNORE) assign(k, list(), envir = env)
  for (e in ex) collect_defs(e, env)

  msgs <- character(0)
  for (nm in ls(env)) {
    obj <- get(nm, envir = env)
    if (!is.function(obj) || nm %in% c(KNOWN, IGNORE)) next
    checked <- checked + 1L
    tryCatch(codetools::checkUsage(obj, name = nm,
      report = function(x) msgs <<- c(msgs, x)), error = function(e) NULL)
  }
  bad <- grep("no visible binding for global variable", msgs, value = TRUE)

  # codetools degisken adini tirnak icinde verir. R'in "fancy quotes" ayari
  # ve yerel ayarlara gore bu tirnak DUZ ('x') ya da EGIK (\u2018x\u2019)
  # olabilir. Mesaja gore metin eslestirmek yerine ADI cikarip karsilastiriyoruz;
  # aksi halde egik tirnak kullanan sistemlerde eleme calismaz.
  # DIKKAT: codetools mesajlari SATIR SONU ile biter. Once bosluk/satir sonu
  # temizlenmeli, sonra tirnaklar soyulmalidir; aksi halde kapanis tirnagi
  # ada yapisir ve eleme calismaz.
  # codetools mesaj bicimi ortama gore DEGISIR:
  #   'x'                          (duz tirnak, satir referansi yok)
  #   \u2018x\u2019 (dosya.R:111)          (egik tirnak + satir referansi)
  # Satir referansi keep.source ayarina bagli olarak gelir; RStudio'da gelir,
  # Rscript'te varsayilan olarak gelmez. Bu yuzden kapanis tirnagindan
  # SONRASINI da atiyoruz.
  # codetools mesaj bicimi surume/yerel ayara gore DEGISIR:
  #   - tirnak duz ('x') veya egik (\u2018x\u2019) olabilir
  #   - satir sonu ile bitebilir
  #   - yeni surumler sona "(dosya:satir)" ekler
  # Bu yuzden tirnak soymak yerine, on ekten SONRAKI ILK R TANIMLAYICISINI
  # aliyoruz. Bicim degisse de calisir.
  var_of <- function(msg) {
    m <- sub(".*no visible binding for global variable", "", msg)
    hit <- regmatches(m, regexpr("[A-Za-z._][A-Za-z0-9._]*", m))
    if (length(hit) == 0) "" else hit
  }
  keep <- vapply(bad, function(b) !(var_of(b) %in% c(KNOWN, IGNORE, assigned)),
                 logical(1))
  bad <- unique(bad[keep])
  if (length(bad) > 0) {
    cat("\n== ", f, " ==\n", sep = ""); cat(paste0("  ", bad, collapse = "\n"), "\n")
    total <- total + length(bad)
  }
}
cat("\n", checked, " fonksiyon denetlendi. ",
    if (total == 0) "Tanimsiz degisken referansi YOK." else
      paste0(total, " potansiyel sorun."), "\n", sep = "")
