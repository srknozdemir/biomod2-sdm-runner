# plot_response_threshold.R
# =============================================================================
# ESIKLENDIRILMIS TEPKI EGRILERI (yayin formati)
#
# Kullanicinin daha once Jenks kesimiyle urettigi grafigin AYNI gorsel dili,
# fakat esik degeri artik MaxTSS'ten (ensemble degerlendirmesinden) gelir.
#
# Girdi : results/<run_id>/response_curves/ensemble_response_curves_data.csv
#         results/<run_id>/ensemble_maxTSS_cutoffs.csv        (binary)
#         results/<run_id>/ensemble_class_breaks.csv          (surekli)
# Cikti : plots/ensemble_response_curves_threshold.png   (arayuz izleyicisi)
#         plots/ensemble_response_curves_threshold.tiff  (300 dpi, yayin)
#         response_curves/response_curves_threshold_data.csv (yeniden cizim icin)
#
# ---------------------------------------------------------------------------
# OLCEK UYARISI (onemli)
# biomod2 tahminleri genellikle 0-1000 tam sayi olceginde saklar; get_evaluations()
# cutoff degerleri de ayni olcektedir. Ancak bm_PlotResponseCurves()$tab ciktisi
# surume ve nesneye gore 0-1 veya 0-1000 gelebilir (bkz. biomod2 issue #493).
# Bu yuzden egri ve esik BAGIMSIZ olarak olcek denetiminden gecirilip 0-1'e
# indirgenir. Aksi halde esik cizgisi grafigin disina duser.
#
# ---------------------------------------------------------------------------
# YORUM NOTU
# biomod2 egrileri "evaluation strip" yontemidir (Elith et al. 2005 Ecological
# Modelling 186:280-289): diger degiskenler medyanda SABIT tutulur, tek degisken
# gradyan boyunca degistirilir. Bu, her degisken icin AYRI tek-degiskenli model
# kurmakla ayni sey DEGILDIR. Sekil altyazisi buna gore yazilmalidir.
# =============================================================================

# --------------------------------------------------------------------------
# 1) VERI HAZIRLIGI  (ggplot2 gerektirmez -> ayrica test edilebilir)
# --------------------------------------------------------------------------
rc_read <- function(rc_csv) {
  if (!file.exists(rc_csv)) return(NULL)
  d <- tryCatch(utils::read.csv(rc_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)

  nm <- names(d)
  pick <- function(cands) { h <- cands[cands %in% nm]; if (length(h)) h[1] else NULL }

  # biomod2 yolu: expl.name / expl.val / pred.val / pred.name
  # yedek yol   : variable  / x        / pred     / obj
  v_col <- pick(c("expl.name", "variable"))
  x_col <- pick(c("expl.val", "x"))
  y_col <- pick(c("pred.val", "pred"))
  m_col <- pick(c("pred.name", "obj", "curve"))

  if (is.null(v_col) || is.null(x_col) || is.null(y_col)) return(NULL)

  out <- data.frame(
    variable = as.character(d[[v_col]]),
    env_raw  = d[[x_col]],
    suit_raw = suppressWarnings(as.numeric(d[[y_col]])),
    model    = if (is.null(m_col)) NA_character_ else as.character(d[[m_col]]),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$suit_raw), , drop = FALSE]
}

# 0-1000 -> 0-1 indirgeme (deger araligindan cikarim)
rc_to01 <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  mx <- suppressWarnings(max(v, na.rm = TRUE))
  if (!is.finite(mx)) return(v)
  if (mx > 1.001) v / 1000 else v
}

# Ensemble uyesini sec (EMwmean varsayilan). Bulunamazsa uyari + ortalama.
rc_pick_model <- function(d, em_pick = "EMwmean") {
  if (all(is.na(d$model))) return(list(d = d, used = "(tek egri)"))
  mods <- unique(d$model)
  hit <- mods[grepl(em_pick, mods, ignore.case = TRUE)]
  if (length(hit) > 0) {
    return(list(d = d[d$model %in% hit[1], , drop = FALSE], used = hit[1]))
  }
  warning("'", em_pick, "' egrisi bulunamadi; mevcut egriler ortalanacak: ",
          paste(utils::head(mods, 5), collapse = ", "))
  agg <- stats::aggregate(suit_raw ~ variable + env_raw, data = d, FUN = mean, na.rm = TRUE)
  agg$model <- "mean_of_all"
  list(d = agg, used = "mean_of_all")
}

# Esik degerini cutoff tablosundan al
rc_threshold <- function(cutoff_csv, em_pick = "EMwmean", fallback = NA_real_) {
  if (!file.exists(cutoff_csv)) return(fallback)
  d <- tryCatch(utils::read.csv(cutoff_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(fallback)

  # binary: model / cutoff   |  surekli: break_value
  if (all(c("model", "cutoff") %in% names(d))) {
    hit <- grep(em_pick, d$model, ignore.case = TRUE)
    val <- if (length(hit) > 0) d$cutoff[hit[1]] else stats::median(d$cutoff, na.rm = TRUE)
    return(suppressWarnings(as.numeric(val)))
  }
  if ("break_value" %in% names(d)) {
    return(suppressWarnings(as.numeric(stats::median(d$break_value, na.rm = TRUE))))
  }
  fallback
}

# Kategorik degiskenleri sayisal kodlara cevir.
# Gerekce: facet_wrap tek bir x estetigi kullanir; kategorik panelde faktor
# kullanilirsa surekli panellerle catisir. Bu yuzden kod + sozluk uretilir.
rc_encode_categorical <- function(d, cat_vars) {
  lookup <- list()
  for (v in intersect(cat_vars, unique(d$variable))) {
    idx <- d$variable == v
    raw <- d$env_raw[idx]
    num <- suppressWarnings(as.numeric(as.character(raw)))
    if (!any(is.na(num))) {
      d$env[idx] <- num
      lv <- sort(unique(num))
      lookup[[v]] <- data.frame(variable = v, code = lv, label = as.character(lv),
                                stringsAsFactors = FALSE)
    } else {
      lv <- sort(unique(as.character(raw)))
      d$env[idx] <- match(as.character(raw), lv)
      lookup[[v]] <- data.frame(variable = v, code = seq_along(lv), label = lv,
                                stringsAsFactors = FALSE)
    }
  }
  cont <- setdiff(unique(d$variable), cat_vars)
  for (v in cont) {
    idx <- d$variable == v
    d$env[idx] <- suppressWarnings(as.numeric(as.character(d$env_raw[idx])))
  }
  list(d = d, lookup = if (length(lookup)) do.call(rbind, lookup) else NULL)
}

# Esik USTU ardisik segmentler (siyah egri uzerine kirmizi bindirme)
rc_segments <- function(d, threshold) {
  out <- list()
  for (v in unique(d$variable)) {
    s <- d[d$variable == v, , drop = FALSE]
    s <- s[order(s$env), , drop = FALSE]
    if (nrow(s) < 2) next
    n <- nrow(s)
    seg <- data.frame(
      variable  = v,
      env       = s$env[-n],
      env_next  = s$env[-1],
      suit      = s$suit[-n],
      suit_next = s$suit[-1],
      stringsAsFactors = FALSE
    )
    seg <- seg[seg$suit >= threshold & seg$suit_next >= threshold, , drop = FALSE]
    if (nrow(seg) > 0) out[[length(out) + 1]] <- seg
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# --------------------------------------------------------------------------
# Tum veri hazirligini tek adimda yapar (test edilebilir cekirdek)
# --------------------------------------------------------------------------
rc_prepare <- function(rc_csv, cutoff_csv, var_types = NULL,
                       em_pick = "EMwmean", var_labels = NULL,
                       threshold_override = NA_real_) {

  d <- rc_read(rc_csv)
  if (is.null(d)) return(NULL)

  sel <- rc_pick_model(d, em_pick)
  d <- sel$d

  # Olcek: egri ve esik BAGIMSIZ denetlenir
  d$suit <- rc_to01(d$suit_raw)

  thr <- if (is.finite(threshold_override)) threshold_override else
    rc_threshold(cutoff_csv, em_pick, fallback = NA_real_)
  if (is.finite(thr)) thr <- rc_to01(thr)

  # Esik hala yoksa: egrinin medyani (acikca uyarilir)
  thr_source <- "MaxTSS"
  if (!is.finite(thr)) {
    thr <- stats::median(d$suit, na.rm = TRUE)
    thr_source <- "medyan (MaxTSS bulunamadi)"
    warning("Esik tablosu okunamadi; gecici olarak egri medyani kullanildi.")
  }

  # Kategorik degiskenler
  cat_vars <- character(0)
  if (!is.null(var_types) && is.list(var_types) && !is.null(names(var_types))) {
    tp <- tolower(vapply(var_types, function(x) as.character(x)[1], ""))
    cat_vars <- names(var_types)[tp %in% c("categorical", "nominal")]
  }
  cat_vars <- intersect(cat_vars, unique(d$variable))

  d$env <- NA_real_
  enc <- rc_encode_categorical(d, cat_vars)
  d <- enc$d
  d <- d[is.finite(d$env) & is.finite(d$suit), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  # Etiket cevirisi (ornek: anakaya -> bedrock)
  d$variable_lab <- d$variable
  if (!is.null(var_labels) && length(var_labels) > 0) {
    hit <- d$variable %in% names(var_labels)
    d$variable_lab[hit] <- unlist(var_labels[d$variable[hit]], use.names = FALSE)
  }

  d$above <- d$suit >= thr

  d_cont <- d[!(d$variable %in% cat_vars), , drop = FALSE]
  segs <- if (nrow(d_cont) > 0) rc_segments(d_cont, thr) else NULL
  if (!is.null(segs)) {
    segs$variable_lab <- d$variable_lab[match(segs$variable, d$variable)]
  }

  list(
    data       = d,
    segments   = segs,
    threshold  = thr,
    thr_source = thr_source,
    cat_vars   = cat_vars,
    lookup     = enc$lookup,
    em_used    = sel$used
  )
}

# --------------------------------------------------------------------------
# 2) CIZIM
# --------------------------------------------------------------------------
rc_plot <- function(prep, sp_name = "", y_lab = "Suitability", subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 yok -> esiklendirilmis tepki egrisi cizilemedi.")
    return(NULL)
  }
  gg <- ggplot2::ggplot
  aes <- ggplot2::aes

  d    <- prep$data
  thr  <- prep$threshold
  segs <- prep$segments
  catv <- prep$cat_vars

  d_cat  <- d[d$variable %in% catv, , drop = FALSE]
  d_cont <- d[!(d$variable %in% catv), , drop = FALSE]

  lab_lo <- paste0("< ",  format(round(thr, 3), nsmall = 3))
  lab_hi <- paste0(">= ", format(round(thr, 3), nsmall = 3))

  p <- gg() +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed")

  if (nrow(d_cat) > 0) {
    p <- p + ggplot2::geom_col(
      data = d_cat,
      aes(x = .data$env, y = .data$suit, fill = .data$above),
      width = 0.7
    )
  }

  if (nrow(d_cont) > 0) {
    p <- p + ggplot2::geom_line(
      data = d_cont,
      aes(x = .data$env, y = .data$suit, group = .data$variable_lab),
      color = "black", linewidth = 0.7
    )
  }

  if (!is.null(segs) && nrow(segs) > 0) {
    p <- p + ggplot2::geom_segment(
      data = segs,
      aes(x = .data$env, xend = .data$env_next,
          y = .data$suit, yend = .data$suit_next,
          group = .data$variable_lab),
      color = "red", linewidth = 0.9
    )
  }

  p <- p +
    ggplot2::facet_wrap(~ variable_lab, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "grey40", `TRUE` = "red"),
      name   = "Threshold",
      labels = c(`FALSE` = lab_lo, `TRUE` = lab_hi)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text    = ggplot2::element_text(size = 11),
      axis.title   = ggplot2::element_text(size = 12),
      strip.text   = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 11),
      legend.text  = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(
      x = "Environmental gradient",
      y = y_lab,
      title = if (nzchar(sp_name))
        paste0(sp_name, " - ensemble response curves") else "Ensemble response curves",
      subtitle = subtitle
    )
  p
}

# --------------------------------------------------------------------------
# 3) TEK CAGRILIK SARMALAYICI (run betikleri bunu kullanir)
# --------------------------------------------------------------------------
plot_response_threshold <- function(out_dir, sp_name = "", var_types = NULL,
                                    em_pick = "EMwmean", var_labels = NULL,
                                    y_lab = "Suitability (ensemble)",
                                    threshold_override = NA_real_,
                                    width_cm = 20, height_cm = 12, dpi = 300) {

  rc_csv <- file.path(out_dir, "response_curves", "ensemble_response_curves_data.csv")
  if (!file.exists(rc_csv)) {
    rc_csv <- file.path(out_dir, "response_curves", "models_response_curves_data.csv")
  }
  cut_csv <- file.path(out_dir, "ensemble_maxTSS_cutoffs.csv")
  if (!file.exists(cut_csv)) cut_csv <- file.path(out_dir, "ensemble_class_breaks.csv")

  prep <- tryCatch(
    rc_prepare(rc_csv, cut_csv, var_types = var_types, em_pick = em_pick,
               var_labels = var_labels, threshold_override = threshold_override),
    error = function(e) { warning("Tepki egrisi hazirligi basarisiz: ",
                                  conditionMessage(e)); NULL }
  )
  if (is.null(prep)) return(invisible(FALSE))

  message("Esiklendirilmis tepki egrileri | esik = ", round(prep$threshold, 4),
          " (", prep$thr_source, ") | egri = ", prep$em_used)

  # Yeniden cizim icin duzenli veri
  tryCatch({
    utils::write.csv(prep$data[, c("variable", "variable_lab", "env", "suit", "above")],
                     file.path(out_dir, "response_curves",
                               "response_curves_threshold_data.csv"), row.names = FALSE)
    if (!is.null(prep$lookup)) {
      utils::write.csv(prep$lookup,
                       file.path(out_dir, "response_curves",
                                 "response_curves_class_labels.csv"), row.names = FALSE)
    }
  }, error = function(e) NULL)

  sub <- paste0("Threshold = ", round(prep$threshold, 3),
                " (", prep$thr_source, "); other predictors fixed at median")
  p <- rc_plot(prep, sp_name = sp_name, y_lab = y_lab, subtitle = sub)
  if (is.null(p)) return(invisible(FALSE))

  plots_dir <- file.path(out_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Bassiz Linux'ta varsayilan png cihazi olmayabilir -> sirayla dene
  save_try <- function(path, ...) {
    args <- list(filename = path, plot = p, width = width_cm, height = height_cm,
                 units = "cm", dpi = dpi, ...)
    tryCatch({ do.call(ggplot2::ggsave, args); TRUE },
             error = function(e) FALSE)
  }
  ok <- save_try(file.path(plots_dir, "ensemble_response_curves_threshold.png"))
  if (!ok && isTRUE(unname(capabilities("cairo")))) {
    ok <- save_try(file.path(plots_dir, "ensemble_response_curves_threshold.png"),
                   type = "cairo")
  }
  if (!ok && requireNamespace("ragg", quietly = TRUE)) {
    ok <- save_try(file.path(plots_dir, "ensemble_response_curves_threshold.png"),
                   device = ragg::agg_png)
  }
  if (!ok) warning("PNG could not be saved (no usable graphics device).")

  tryCatch({
    ggplot2::ggsave(file.path(plots_dir, "ensemble_response_curves_threshold.tiff"),
                    plot = p, width = width_cm, height = height_cm,
                    units = "cm", dpi = dpi, compression = "lzw")
  }, error = function(e) warning("TIFF kaydedilemedi: ", conditionMessage(e)))

  invisible(ok)
}
