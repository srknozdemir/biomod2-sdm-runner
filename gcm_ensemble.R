# gcm_ensemble.R
# =============================================================================
# GCM ENSEMBLE + CHANGE MAPS (post-processing of future projections)
#
# Called by run_biomod.R (presence/absence) and run_biomod_cont.R (continuous)
# AFTER all future scenarios have been projected. Works only on the GeoTIFFs the
# pipelines have already written to rasters/, so it never refits or reprojects.
#
# Why average PREDICTIONS (not climate layers): each GCM is projected separately
# and the predictions are combined afterwards, so the spread among GCMs is kept
# as an explicit uncertainty layer (Araujo & New 2007; Thuiller et al. 2019).
#
# Grouping: every future projection carries (period, ssp, gcm). Projections that
# share period + ssp are one group; the GCMs inside a group are averaged.
# A group with a single member (e.g. manual mode, one folder per scenario)
# still gets change maps; only the GCM spread layers are skipped.
#
# Outputs (all in rasters/ and plots/, so the app viewers list them):
#   CONTINUOUS response or habitat suitability
#     gcm_<P>_<SSP>_mean.tif          mean of the GCM predictions
#     gcm_<P>_<SSP>_sd.tif            SD among GCMs (n >= 2)
#     change_<P>_<SSP>_abs.tif        mean future - current
#     change_<P>_<SSP>_pct.tif        100 * (future - current) / |current|
#     change_<P>_<SSP>_agreement.tif  share of GCMs whose change has the sign of the mean change
#     change_<P>_<SSP>_class.tif      1 robust decrease | 2 decrease, low agreement |
#                                     3 no change (|change| < tolerance) |
#                                     4 increase, low agreement | 5 robust increase
#     change_<P>_<SSP>_<GCM>_abs.tif  per-GCM change
#   PRESENCE/ABSENCE (MaxTSS binary maps), in addition to the above on suitability
#     range_<P>_<SSP>_<GCM>.tif       -2 loss | -1 stable presence | 0 stable absence | 1 gain
#                                     (coding of biomod2::BIOMOD_RangeSize)
#     range_<P>_<SSP>_presence_agreement.tif  share of GCMs predicting presence
#     range_<P>_<SSP>_consensus.tif   range change of the majority-vote future map
#   TABLES (out_dir/)
#     change_summary.csv              area increase / decrease / no change, mean change
#     range_change_summary.csv        Loss, Stable0, Stable1, Gain, PercLoss, PercGain,
#                                     SpeciesRangeChange (binary only)
#     gcm_scenario_registry.csv       which folder / GCM / SSP / period went where
#   FIGURES (plots/)
#     change_<P>_<SSP>.png, range_<P>_<SSP>.png, change_overview.png
#
# Settings read from cfg (all optional):
#   cfg$change_em          ensemble layer to use: "EMwmean" (default), "EMmean", "EMmedian"
#   cfg$change_agree       agreement threshold for "robust" change (default 0.8)
#   cfg$change_tol_frac    no-change tolerance as a fraction of the SD of the current map (default 0.05)
#   cfg$change_tol_abs     absolute no-change tolerance (overrides the fraction)
#   cfg$gcm_exclude        GCM names for an additional "excluding" ensemble (sensitivity)
#   cfg$change_mask        polygon (gpkg/shp) or raster; summaries are ALSO given inside it
#   cfg$binary_consensus   share of GCMs needed for presence in the consensus map (default 0.5)
# =============================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    ok <- !is.null(a) && length(a) > 0 && !all(is.na(a))
    if (isTRUE(ok)) a else b
  }
}

# ---- labels -----------------------------------------------------------------
# "SSP2-4.5", "ssp245", "SSP 245" -> "ssp245"
gcm_norm_ssp <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_character_)
  m <- regmatches(x, regexpr("(?i)ssp\\s*[-_]?\\s*([0-9])[-_.]?([0-9])[.]?([0-9])", x, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return(NA_character_)
  paste0("ssp", gsub("[^0-9]", "", m))
}

# Folder names written by prepare_climate_for_biomod2.R:
#   future_<period>_<ssp>_<GCM>     e.g. future_2041-2060_ssp245_MPI-ESM1-2-HR
gcm_parse_folder <- function(path, pattern = "^future_(.+?)_(ssp[0-9]{3})_(.+)$") {
  b <- basename(normalizePath(path, winslash = "/", mustWork = FALSE))
  m <- regmatches(b, regexec(pattern, b, ignore.case = TRUE))[[1]]
  if (length(m) == 4) return(list(period = m[2], ssp = tolower(m[3]), gcm = m[4]))
  list(period = NA_character_, ssp = gcm_norm_ssp(b), gcm = NA_character_)
}

gcm_key <- function(...) {
  k <- paste(unlist(list(...)), collapse = "_")
  gsub("[^A-Za-z0-9._-]+", "-", k)
}

# ---- locate the ensemble GeoTIFFs written by the pipelines -------------------
gcm_find_layer <- function(rasters_dir, tag, em_pick = "EMwmean",
                           suffix = "_continuous.tif") {
  f <- list.files(rasters_dir, pattern = "\\.tif$", full.names = FALSE)
  f <- f[startsWith(f, paste0(tag, "_")) & endsWith(f, suffix)]
  if (length(f) == 0) return(NA_character_)
  # prefer the R2-filtered ensemble when several selection metrics exist
  f <- f[order(!grepl("ByRsquared", f, fixed = TRUE), f)]
  for (p in unique(c(em_pick, "EMwmean", "EMmean", "EMmedian"))) {
    hit <- f[grepl(p, f, fixed = TRUE)]
    if (length(hit) > 0) return(file.path(rasters_dir, hit[1]))
  }
  NA_character_
}

# ---- small raster helpers -----------------------------------------------------
gcm_area_km2 <- function(r) {
  a <- tryCatch(terra::cellSize(r[[1]], unit = "km", mask = FALSE), error = function(e) NULL)
  if (is.null(a)) { a <- r[[1]]; terra::values(a) <- prod(terra::res(r)) / 1e6 }
  a
}

gcm_read_mask <- function(path, ref) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  m <- tryCatch({
    if (grepl("\\.(gpkg|shp|geojson|kml)$", path, ignore.case = TRUE)) {
      v <- terra::vect(path)
      if (!terra::same.crs(v, ref)) v <- terra::project(v, terra::crs(ref))
      terra::rasterize(v, ref[[1]], field = 1, background = NA)
    } else {
      r <- terra::rast(path)[[1]]
      if (!terra::compareGeom(r, ref, stopOnError = FALSE)) r <- terra::resample(r, ref, method = "near")
      terra::ifel(r > 0, 1, NA)
    }
  }, error = function(e) { message("Change mask could not be read (ignored): ", conditionMessage(e)); NULL })
  m
}

gcm_png <- function(path, expr, width = 2000, height = 1600, res = 200) {
  if (exists("safe_png", mode = "function"))
    return(safe_png(path, expr, width = width, height = height, res = res))
  ok <- tryCatch({
    grDevices::png(path, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr); TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}

# area-weighted summaries of a change raster (inside a region mask or everywhere)
gcm_change_stats <- function(cur, fut, delta, cls, area, region = NULL) {
  m <- !is.na(cur) & !is.na(fut)
  if (!is.null(region)) m <- m & !is.na(region)
  v_cur <- terra::values(cur, mat = FALSE); v_fut <- terra::values(fut, mat = FALSE)
  v_d   <- terra::values(delta, mat = FALSE); v_a <- terra::values(area, mat = FALSE)
  v_c   <- if (!is.null(cls)) terra::values(cls, mat = FALSE) else rep(NA, length(v_d))
  k <- which(terra::values(m, mat = FALSE) == 1)
  if (length(k) == 0) return(NULL)
  A <- sum(v_a[k], na.rm = TRUE)
  wm <- function(x) sum(x[k] * v_a[k], na.rm = TRUE) / A
  ar <- function(sel) sum(v_a[k][sel], na.rm = TRUE)
  out <- data.frame(
    area_km2          = A,
    mean_current      = wm(v_cur),
    mean_future       = wm(v_fut),
    mean_change       = wm(v_d),
    pct_change_of_mean = 100 * (wm(v_fut) - wm(v_cur)) / abs(wm(v_cur)),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(v_c[k]))) {
    cc <- v_c[k]
    out$area_decrease_km2        <- ar(cc %in% c(1, 2))
    out$area_nochange_km2        <- ar(cc == 3)
    out$area_increase_km2        <- ar(cc %in% c(4, 5))
    out$area_robust_decrease_km2 <- ar(cc == 1)
    out$area_robust_increase_km2 <- ar(cc == 5)
    out$pct_area_decrease        <- 100 * out$area_decrease_km2 / A
    out$pct_area_nochange        <- 100 * out$area_nochange_km2 / A
    out$pct_area_increase        <- 100 * out$area_increase_km2 / A
    out$pct_area_robust_decrease <- 100 * out$area_robust_decrease_km2 / A
    out$pct_area_robust_increase <- 100 * out$area_robust_increase_km2 / A
  }
  out
}

gcm_range_stats <- function(code, area, region = NULL) {
  v <- terra::values(code, mat = FALSE); a <- terra::values(area, mat = FALSE)
  k <- !is.na(v)
  if (!is.null(region)) k <- k & !is.na(terra::values(region, mat = FALSE))
  n  <- function(x) sum(k & v == x)
  ak <- function(x) sum(a[k & v == x], na.rm = TRUE)
  loss <- n(-2); s1 <- n(-1); s0 <- n(0); gain <- n(1)
  cur <- loss + s1
  data.frame(Loss = loss, Stable0 = s0, Stable1 = s1, Gain = gain,
             PercLoss = if (cur > 0) 100 * loss / cur else NA_real_,
             PercGain = if (cur > 0) 100 * gain / cur else NA_real_,
             SpeciesRangeChange = if (cur > 0) 100 * (gain - loss) / cur else NA_real_,
             CurrentRangeSize = cur, FutureRangeSize = s1 + gain,
             Loss_km2 = ak(-2), Stable1_km2 = ak(-1), Gain_km2 = ak(1))
}

# ---- plotting -------------------------------------------------------------------
GCM_CLASS_COLS <- c("#B2182B", "#F4A582", "#E8E8E8", "#92C5DE", "#2166AC")
GCM_CLASS_LABS <- c("Robust decrease", "Decrease (low agreement)", "No change",
                    "Increase (low agreement)", "Robust increase")
GCM_RANGE_COLS <- c("-2" = "#D7301F", "-1" = "#1A9850", "0" = "#EFEFEF", "1" = "#2C7BB6")
GCM_RANGE_LABS <- c("-2" = "Loss", "-1" = "Stable presence", "0" = "Stable absence", "1" = "Gain")

gcm_plot_change <- function(path, cur, fut_mean, delta, cls, title, resp_lab, n = 2) {
  rng <- range(c(terra::values(cur, mat = FALSE), terra::values(fut_mean, mat = FALSE)), na.rm = TRUE)
  dmax <- max(abs(terra::values(delta, mat = FALSE)), na.rm = TRUE)
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  seq_cols <- grDevices::hcl.colors(100, "viridis")
  div_cols <- grDevices::hcl.colors(101, "Blue-Red 3", rev = TRUE)
  gcm_png(path, {
    graphics::par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0))
    terra::plot(cur, range = rng, col = seq_cols, main = paste0("Current (", resp_lab, ")"), axes = FALSE)
    terra::plot(fut_mean, range = rng, col = seq_cols, main = if (n > 1) "Future (GCM mean)" else "Future", axes = FALSE)
    terra::plot(delta, range = c(-dmax, dmax), col = div_cols, main = "Change (future - current)", axes = FALSE)
    present <- sort(unique(stats::na.omit(terra::values(cls, mat = FALSE))))
    labs <- if (n > 1) GCM_CLASS_LABS else c("Decrease", "Decrease", "No change", "Increase", "Increase")
    terra::plot(cls, col = GCM_CLASS_COLS[present], type = "classes",
                levels = labs[present], main = "Change class", axes = FALSE,
                mar = c(1.5, 1, 2.5, 11), plg = list(cex = 0.8))
    graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  })
}

gcm_plot_range <- function(path, prop, code, title) {
  gcm_png(path, {
    graphics::par(mfrow = c(1, if (is.null(prop)) 1 else 2), oma = c(0, 0, 2.5, 0))
    if (!is.null(prop)) terra::plot(prop, range = c(0, 1), type = "continuous", col = grDevices::hcl.colors(100, "YlGnBu", rev = TRUE),
                main = "Share of GCMs predicting presence", axes = FALSE)
    present <- as.character(sort(unique(stats::na.omit(terra::values(code, mat = FALSE)))))
    terra::plot(code, col = unname(GCM_RANGE_COLS[present]), type = "classes",
                levels = unname(GCM_RANGE_LABS[present]),
                main = if (is.null(prop)) "Range change" else "Range change (GCM consensus)",
                axes = FALSE, mar = c(1.5, 1, 2.5, 9), plg = list(cex = 0.8))
    graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  }, height = 1000)
}

gcm_plot_overview <- function(path, summ, resp_lab) {
  s <- summ[summ$region == "full", , drop = FALSE]
  ens <- s[(s$member == "GCM ensemble" & s$set == "all GCMs") | s$set == "single projection", , drop = FALSE]
  if (nrow(ens) == 0) return(invisible(NULL))
  exc <- s[s$member == "GCM ensemble" & s$set != "all GCMs", , drop = FALSE]
  ind <- s[s$member != "GCM ensemble" & s$set == "all GCMs", , drop = FALSE]
  ens <- ens[order(ens$period, ens$ssp), , drop = FALSE]
  x <- seq_len(nrow(ens))
  yl <- range(c(0, ens$mean_change, ind$mean_change, exc$mean_change), na.rm = TRUE)
  gcm_png(path, {
    graphics::par(mar = c(7, 5, 3, 1))
    graphics::plot(x, ens$mean_change, ylim = yl, xlim = c(0.5, length(x) + 0.5), pch = 19, cex = 1.6,
                   xaxt = "n", xlab = "", ylab = paste0("Mean change in ", resp_lab),
                   main = "Projected change (area-weighted mean)")
    graphics::abline(h = 0, lty = 2, col = "grey40")
    graphics::axis(1, at = x, labels = paste(ens$period, ens$ssp, sep = "\n"), las = 2, cex.axis = 0.75)
    for (i in x) {
      d <- ind[ind$period == ens$period[i] & ind$ssp == ens$ssp[i], , drop = FALSE]
      if (nrow(d) > 0) graphics::points(rep(i, nrow(d)) + seq(-0.15, 0.15, length.out = nrow(d)),
                                        d$mean_change, pch = 1, col = "grey30")
      e <- exc[exc$period == ens$period[i] & exc$ssp == ens$ssp[i], , drop = FALSE]
      if (nrow(e) > 0) graphics::points(rep(i + 0.25, nrow(e)), e$mean_change, pch = 17, col = "#D95F02", cex = 1.3)
    }
    lg <- c("GCM ensemble", "Individual GCMs")
    lp <- c(19, 1); lc <- c("black", "grey30")
    if (nrow(exc) > 0) { lg <- c(lg, paste0("Ensemble, ", unique(exc$set)[1])); lp <- c(lp, 17); lc <- c(lc, "#D95F02") }
    graphics::legend("bottomleft", legend = lg, pch = lp, col = lc, bty = "n", cex = 0.85)
  }, height = 1300)
}

gcm_rbind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols, drop = FALSE]
  }))
}

# ---- main ---------------------------------------------------------------------
# registry: data.frame(period, ssp, gcm, tag, folder)  -- one row per future projection
# mode    : "quant" (continuous response) or "binary" (presence/absence)
gcm_postprocess <- function(cfg, registry, rasters_dir, plots_dir, out_dir, mode = "quant",
                            resp_lab = "response") {
  if (is.null(registry) || nrow(registry) == 0) return(invisible(NULL))
  em_pick  <- as.character(cfg$change_em %||% "EMwmean")[1]
  agree    <- as.numeric(cfg$change_agree %||% 0.8)[1]
  tol_frac <- as.numeric(cfg$change_tol_frac %||% 0.05)[1]
  tol_abs  <- suppressWarnings(as.numeric(cfg$change_tol_abs %||% NA)[1])
  cons_thr <- as.numeric(cfg$binary_consensus %||% 0.5)[1]
  excl     <- unique(as.character(unlist(cfg$gcm_exclude %||% character(0))))

  # fill missing labels from folder names
  # (manual mode: SSP from the folder name or the scenario label; the label then
  #  identifies the member, so two folders with the same SSP and period are averaged)
  for (i in seq_len(nrow(registry))) {
    pf <- gcm_parse_folder(registry$folder[i])
    if (is.na(registry$ssp[i]) || !nzchar(registry$ssp[i]))
      registry$ssp[i] <- pf$ssp %||% gcm_norm_ssp(registry$label[i]) %||% registry$label[i]
    if (is.na(registry$gcm[i]) || !nzchar(registry$gcm[i]))
      registry$gcm[i] <- pf$gcm %||% registry$label[i] %||% paste0("member", i)
  }

  registry$cont_file <- vapply(registry$tag, function(t) gcm_find_layer(rasters_dir, t, em_pick), "")
  registry$bin_file  <- if (mode == "binary")
    vapply(registry$tag, function(t) gcm_find_layer(rasters_dir, t, em_pick, "_binary_maxTSS.tif"), "")
    else NA_character_
  utils::write.csv(registry, file.path(out_dir, "gcm_scenario_registry.csv"), row.names = FALSE)

  cur_file <- gcm_find_layer(rasters_dir, "current_ensemble", em_pick)
  if (is.na(cur_file)) { message("GCM ensemble: current ensemble raster not found -> skipped."); return(invisible(NULL)) }
  cur  <- terra::rast(cur_file)[[1]]
  area <- gcm_area_km2(cur)
  region <- gcm_read_mask(cfg$change_mask %||% "", cur)
  tol <- if (is.finite(tol_abs) && tol_abs >= 0) tol_abs else
    tol_frac * stats::sd(terra::values(cur, mat = FALSE), na.rm = TRUE)
  message("GCM ensemble: layer '", basename(cur_file), "', no-change tolerance = ", signif(tol, 3),
          ", agreement threshold = ", agree)

  cur_bin <- NULL
  if (mode == "binary") {
    cb <- gcm_find_layer(rasters_dir, "current_ensemble", em_pick, "_binary_maxTSS.tif")
    if (!is.na(cb)) cur_bin <- terra::rast(cb)[[1]]
  }

  summ <- list(); rsumm <- list()
  add_stats <- function(lst, df, ...) {
    if (is.null(df)) return(lst)
    meta <- data.frame(..., stringsAsFactors = FALSE)
    lst[[length(lst) + 1]] <- cbind(meta[rep(1, nrow(df)), , drop = FALSE], df)
    lst
  }

  groups <- unique(registry[, c("period", "ssp")])
  for (g in seq_len(nrow(groups))) {
    P <- groups$period[g]; S <- groups$ssp[g]
    mem_all <- registry[registry$period == P & registry$ssp == S & !is.na(registry$cont_file), , drop = FALSE]
    if (nrow(mem_all) == 0) next
    sets <- list(all = mem_all)
    if (length(excl) > 0 && any(mem_all$gcm %in% excl) && sum(!mem_all$gcm %in% excl) >= 1)
      sets[[paste0("excl_", paste(excl[excl %in% mem_all$gcm], collapse = "+"))]] <- mem_all[!mem_all$gcm %in% excl, , drop = FALSE]

    for (sn in names(sets)) {
      mem <- sets[[sn]]
      key <- gcm_key(P, S, if (sn == "all") NULL else sn)
      single  <- nrow(mem) == 1
      set_lab <- if (single) "single projection" else if (sn == "all") "all GCMs" else sn
      ens_lab <- if (single) mem$gcm[1] else "GCM ensemble"
      fut <- terra::rast(mem$cont_file)
      names(fut) <- make.unique(mem$gcm)
      n <- terra::nlyr(fut)
      fmean <- if (n > 1) terra::app(fut, mean, na.rm = TRUE) else fut[[1]]
      delta <- fmean - cur
      dl    <- fut - cur
      sgn   <- sign(delta)
      agr   <- if (n > 1) terra::app(sign(dl) == sgn, mean, na.rm = TRUE) else terra::ifel(is.na(delta), NA, 1)
      cls   <- terra::ifel(abs(delta) < tol, 3,
               terra::ifel(delta < 0, terra::ifel(agr >= agree, 1, 2),
                                      terra::ifel(agr >= agree, 5, 4)))
      pct   <- terra::ifel(abs(cur) > 0, 100 * delta / abs(cur), NA)

      wr <- function(r, nm, dt = "FLT4S") terra::writeRaster(r, file.path(rasters_dir, nm), overwrite = TRUE, datatype = dt)
      wr(fmean, paste0("gcm_", key, "_mean.tif"))
      if (n > 1) wr(terra::app(fut, stats::sd, na.rm = TRUE), paste0("gcm_", key, "_sd.tif"))
      wr(delta, paste0("change_", key, "_abs.tif"))
      wr(pct,   paste0("change_", key, "_pct.tif"))
      if (n > 1) wr(agr, paste0("change_", key, "_agreement.tif"))
      wr(cls,   paste0("change_", key, "_class.tif"), "INT1U")
      if (sn == "all" && n > 1) for (k in seq_len(n)) wr(dl[[k]], paste0("change_", gcm_key(P, S, mem$gcm[k]), "_abs.tif"))

      gcm_plot_change(file.path(plots_dir, paste0("change_", key, ".png")), cur, fmean, delta, cls,
                      paste0(P, " | ", S, " | ", set_lab, " (n = ", n, ")"), resp_lab, n)

      for (reg in c("full", if (!is.null(region)) "mask")) {
        rr <- if (reg == "mask") region else NULL
        summ <- add_stats(summ, gcm_change_stats(cur, fmean, delta, cls, area, rr),
                          period = P, ssp = S, set = set_lab, member = ens_lab, n_gcm = n, region = reg)
        if (sn == "all" && n > 1) for (k in seq_len(n)) {
          cls_k <- terra::ifel(abs(dl[[k]]) < tol, 3, terra::ifel(dl[[k]] < 0, 2, 4))
          summ <- add_stats(summ, gcm_change_stats(cur, fut[[k]], dl[[k]], cls_k, area, rr),
                            period = P, ssp = S, set = set_lab, member = mem$gcm[k], n_gcm = 1, region = reg)
        }
      }

      # ---- presence/absence: range change ----
      if (mode == "binary" && !is.null(cur_bin)) {
        bm <- mem[!is.na(mem$bin_file), , drop = FALSE]
        if (nrow(bm) > 0) {
          fb <- terra::rast(bm$bin_file)
          prop <- if (terra::nlyr(fb) > 1) terra::app(fb, mean, na.rm = TRUE) else fb[[1]]
          cons <- prop >= cons_thr
          rcode <- function(fbin) terra::ifel(cur_bin == 1, terra::ifel(fbin == 1, -1, -2), terra::ifel(fbin == 1, 1, 0))
          code_c <- rcode(cons)
          if (!single) wr(prop, paste0("range_", key, "_presence_agreement.tif"))
          wr(code_c, paste0("range_", key, if (single) "" else "_consensus", ".tif"), "INT1S")
          gcm_plot_range(file.path(plots_dir, paste0("range_", key, ".png")), if (single) NULL else prop, code_c,
                         paste0(P, " | ", S, " | ", set_lab, " (n = ", terra::nlyr(fb), ")"))
          for (reg in c("full", if (!is.null(region)) "mask")) {
            rr <- if (reg == "mask") region else NULL
            rsumm <- add_stats(rsumm, gcm_range_stats(code_c, area, rr), period = P, ssp = S, set = set_lab,
                               member = if (single) mem$gcm[1] else "GCM consensus", region = reg)
            if (sn == "all" && !single) for (k in seq_len(terra::nlyr(fb))) {
              ck <- rcode(fb[[k]])
              if (reg == "full") wr(ck, paste0("range_", gcm_key(P, S, bm$gcm[k]), ".tif"), "INT1S")
              rsumm <- add_stats(rsumm, gcm_range_stats(ck, area, rr), period = P, ssp = S, set = set_lab,
                                 member = bm$gcm[k], region = reg)
            }
          }
        }
      }
    }
  }

  if (length(summ) > 0) {
    s <- gcm_rbind_fill(summ)
    utils::write.csv(s, file.path(out_dir, "change_summary.csv"), row.names = FALSE)
    gcm_plot_overview(file.path(plots_dir, "change_overview.png"), s, resp_lab)
    message("GCM ensemble: change_summary.csv written (", nrow(s), " rows).")
  }
  if (length(rsumm) > 0) {
    utils::write.csv(gcm_rbind_fill(rsumm), file.path(out_dir, "range_change_summary.csv"), row.names = FALSE)
    message("GCM ensemble: range_change_summary.csv written.")
  }
  invisible(TRUE)
}
