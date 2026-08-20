# install_packages.R
# =============================================================================
# Kurulum betigi / Setup script
#
# KULLANIM
#   RStudio'da :  source("install_packages.R")
#   Terminalde :  Rscript install_packages.R
#
# Bu betik kurulu olanlari ATLAR, yalnizca eksikleri kurar. Tekrar tekrar
# calistirmak guvenlidir. Sonunda hangi modellerin kullanilabilir oldugunu
# gosteren bir ozet tablo basar.
#
# NOT: Once kurulum, sonra app.R. Uygulama acildiginda sol ustteki sistem
# paneli ayni denetimi tekrar yapar.
# =============================================================================

# ---- AYARLAR ---------------------------------------------------------------
INSTALL_DNN      <- FALSE   # TRUE yaparsan 'cito' + 'torch' kurulur (BUYUK indirme,
                            # ~2 GB; yalnizca DNN algoritmasi icin gerekli)
INSTALL_OPTIONAL <- TRUE    # dismo (MAXENT jar yerlesimi icin yardimci)
CRAN_MIRROR      <- "https://cloud.r-project.org"

BIOMOD_MIN <- "4.3"         # surekli/nominal veri tipi destegi bu surumle geldi

# ---- PAKET LISTELERI -------------------------------------------------------
CORE <- c(
  shiny       = "user interface",
  shinyFiles  = "file/folder pickers",
  jsonlite    = "config.json read/write",
  terra       = "raster handling",
  biomod2     = "modelling engine",
  dplyr       = "data manipulation",
  ggplot2     = "plots",
  scales      = "plot axis breaks",
  R.utils     = "timeouts / helpers"
)

# Algoritma -> paket. GLM, SRE ve MAXENT ek paket gerektirmez
# (MAXENT java + maxent.jar ister).
MODEL_PKGS <- c(
  ANN     = "nnet",
  CTA     = "rpart",
  FDA     = "mda",
  GAM     = "mgcv",
  GBM     = "gbm",
  MARS    = "earth",
  MAXNET  = "maxnet",
  RF      = "randomForest",
  RFd     = "randomForest",
  XGBOOST = "xgboost",
  DNN     = "cito"
)

OPTIONAL <- c(dismo = "helps biomod2 locate maxent.jar",
              ragg  = "fallback graphics device (headless Linux)")

# ---- YARDIMCILAR -----------------------------------------------------------
has_pkg <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
ver_of  <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                error = function(e) NA_character_)

hr <- function(ch = "-") cat(strrep(ch, 74), "\n")
say <- function(...) cat(..., "\n", sep = "")

install_if_missing <- function(pkgs) {
  pkgs <- unique(pkgs[nzchar(pkgs)])
  todo <- pkgs[!vapply(pkgs, has_pkg, logical(1))]
  if (length(todo) == 0) {
    say("  All present, nothing to install.")
    return(invisible(character(0)))
  }
  say("  Installing: ", paste(todo, collapse = ", "))
  failed <- character(0)
  for (p in todo) {
    ok <- tryCatch({
      utils::install.packages(p, repos = CRAN_MIRROR, quiet = FALSE)
      has_pkg(p)
    }, error = function(e) {
      message("    ERROR for '", p, "': ", conditionMessage(e)); FALSE
    }, warning = function(w) {
      message("    WARNING for '", p, "': ", conditionMessage(w)); has_pkg(p)
    })
    if (!isTRUE(ok)) failed <- c(failed, p)
  }
  invisible(failed)
}

# ---- BASLA -----------------------------------------------------------------
hr("=")
say("biomod2 SDM Runner - dependency setup")
hr("=")
# NOT: ust duzeyde 'else' yeni satira YAZILAMAZ (R ifadeyi bitmis sayar).
os <- {
  if (.Platform$OS.type == "windows") "Windows"
  else if (grepl("darwin", R.version$os, ignore.case = TRUE)) "macOS"
  else "Linux"
}
say("Platform  : ", os)
say("R version : ", R.version.string)
say("Platform  : ", R.version$platform)
say("Library   : ", .libPaths()[1])
say("CRAN      : ", CRAN_MIRROR)
cat("\n")

if (getRversion() < "4.1.0") {
  warning("R 4.1 or newer is strongly recommended (native pipe / terra requirements).")
}

# Kutuphane klasoru yazilabilir mi? (Windows'ta sik karsilasilan sorun)
libdir <- .libPaths()[1]
if (file.access(libdir, 2) != 0) {
  say("!! The library folder is not writable: ", libdir)
  say("   Run R as administrator, or set a personal library:")
  say('   dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)')
  say('   .libPaths(Sys.getenv("R_LIBS_USER"))')
  cat("\n")
}

# Platforma ozgu on bilgi
if (os == "Linux") {
  say("Linux note: packages are compiled from source. 'terra' needs the system")
  say("  libraries GDAL, PROJ and GEOS. If terra fails to build, install them first:")
  say("    Debian/Ubuntu: sudo apt install libgdal-dev libproj-dev libgeos-dev libudunits2-dev")
  say("    Fedora/RHEL  : sudo dnf install gdal-devel proj-devel geos-devel udunits2-devel")
  if (!isTRUE(unname(capabilities("cairo")))) {
    say("  WARNING: this R has no cairo support. Plots may fail on a headless server.")
    say("           Install the 'ragg' package as a fallback: install.packages(\"ragg\")")
  }
  cat("\n")
} else if (os == "macOS") {
  say("macOS note: CRAN provides binaries, no compiler needed for the packages above.")
  say("  If a package tries to build from source, install Xcode command line tools:")
  say("    xcode-select --install")
  cat("\n")
} else {
  say("Windows note: CRAN provides binaries, no Rtools needed for these packages.")
  cat("\n")
}

failed <- character(0)

hr()
say("1/4  Core packages")
hr()
failed <- c(failed, install_if_missing(names(CORE)))

hr()
say("2/4  Model algorithm packages")
hr()
model_wanted <- unique(unname(MODEL_PKGS[names(MODEL_PKGS) != "DNN"]))
failed <- c(failed, install_if_missing(model_wanted))

if (isTRUE(INSTALL_DNN)) {
  say("  DNN requested -> installing 'cito' (this also pulls in 'torch', ~2 GB)")
  failed <- c(failed, install_if_missing("cito"))
  if (has_pkg("torch")) {
    say("  Note: torch may need a one-off backend download on first use:")
    say("        torch::install_torch()")
  }
} else {
  say("  Skipping 'cito' (DNN). Set INSTALL_DNN <- TRUE at the top to include it.")
}

hr()
say("3/4  Optional packages")
hr()
if (isTRUE(INSTALL_OPTIONAL)) {
  failed <- c(failed, install_if_missing(names(OPTIONAL)))
} else {
  say("  Skipped.")
}

# ---- biomod2 surum denetimi ------------------------------------------------
hr()
say("4/4  Version checks")
hr()

bm <- ver_of("biomod2")
if (is.na(bm)) {
  say("  biomod2 : NOT INSTALLED")
} else if (utils::compareVersion(bm, BIOMOD_MIN) < 0) {
  say("  biomod2 : ", bm, "  -> presence/absence works, but continuous/nominal")
  say("            data types need >= ", BIOMOD_MIN, ".")
  say("            To upgrade from the development source:")
  say('              install.packages("remotes")')
  say('              remotes::install_github("biomodhub/biomod2")')
} else {
  say("  biomod2 : ", bm, "  (continuous / nominal supported)")
}

tv <- ver_of("terra")
say("  terra   : ", if (is.na(tv)) "NOT INSTALLED" else tv)

# ---- java / MAXENT ---------------------------------------------------------
# Rscript bulunabiliyor mu? Arayuz arka plan islerini bununla baslatir.
exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
rs_c <- c(file.path(R.home("bin"), exe), file.path(R.home(), "bin", exe),
          file.path(R.home(), "bin", "x64", exe), unname(Sys.which("Rscript")))
rs_c <- rs_c[nzchar(rs_c)]
rs_hit <- rs_c[file.exists(rs_c)]
say("  Rscript : ", if (length(rs_hit) > 0) rs_hit[1] else "NOT FOUND - jobs cannot start")

java_ok <- nzchar(Sys.which("java"))
say("  java    : ", if (java_ok) "found" else "NOT FOUND - MAXENT cannot run")
if (!java_ok) {
  say("            MAXNET is a pure-R alternative and needs no java.")
}

# ---- OZET ------------------------------------------------------------------
cat("\n")
hr("=")
say("SUMMARY")
hr("=")

core_state <- vapply(names(CORE), has_pkg, logical(1))
say("Core packages:")
for (p in names(CORE)) {
  cat(sprintf("  %-12s %-8s %s\n", p,
              if (core_state[[p]]) "OK" else "MISSING",
              if (core_state[[p]]) ver_of(p) else CORE[[p]]))
}

cat("\nModel availability:\n")
algos <- names(MODEL_PKGS)
extra <- c("GLM", "SRE", "MAXENT")
rows <- c(algos, extra)
for (a in rows) {
  if (a %in% algos) {
    pk <- MODEL_PKGS[[a]]
    st <- if (has_pkg(pk)) "available" else paste0("needs '", pk, "'")
  } else if (a == "MAXENT") {
    st <- if (java_ok) "available (also needs maxent.jar)" else "needs java + maxent.jar"
  } else {
    st <- "available (no extra package)"
  }
  cat(sprintf("  %-8s %s\n", a, st))
}

failed <- unique(failed[nzchar(failed)])
cat("\n")
if (length(failed) > 0) {
  say("!! Could not install: ", paste(failed, collapse = ", "))
  say("   Common causes on Windows:")
  say("     - No internet / proxy blocking CRAN")
  say("     - Library folder not writable (run R as administrator)")
  say("     - A package is loaded in another R session - close it and retry")
  say("   Try one at a time, e.g.:  install.packages(\"", failed[1], "\")")
} else if (all(core_state)) {
  say("All core packages are ready. You can now start the app:")
  say('  shiny::runApp("', normalizePath(getwd(), winslash = "/"), '")')
} else {
  say("Some core packages are still missing - see the table above.")
}
hr("=")
