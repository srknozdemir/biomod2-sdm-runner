# ssdm_convert.R
# =============================================================================
# Species x plot matrix  +  plot coordinates   ->   long format table
#
#   matrix.csv                        coords.csv
#   Species , P1 , P2 , P3            Plot , lat  , lon
#   Pinus nigra , 1 , 0 , 1           P1   , 37.5 , 30.2
#   Cedrus libani , 0 , 1 , 1         P2   , 37.6 , 30.4
#
#   ->  species | plot | x | y | presence
#       Pinus nigra   | P1 | 30.2 | 37.5 | 1
#       Pinus nigra   | P2 | 30.4 | 37.6 | 0
#       ...
#
# Non-binary cell values (cover %, Braun-Blanquet codes, counts) are converted
# to 1/0: anything that means "recorded" becomes 1.
# =============================================================================

# --- Braun-Blanquet and similar codes that mean "present" -------------------
BB_PRESENT <- c("+", "r", "1", "2", "3", "4", "5",
                "2a", "2b", "2m", "i", "p", "x", "y", "yes", "true", "var", "v")
BB_ABSENT  <- c("", "0", "-", ".", "na", "n/a", "no", "false", "yok", "a")

ssdm_read_matrix <- function(path) {
  if (!file.exists(path)) return(list(ok = FALSE, msg = paste0("File not found: ", path)))
  d <- try(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           silent = TRUE)
  if (inherits(d, "try-error") || is.null(d)) {
    return(list(ok = FALSE, msg = "Could not read the CSV."))
  }
  if (ncol(d) < 2) return(list(ok = FALSE, msg = "Matrix must have at least 2 columns."))
  if (nrow(d) < 2) return(list(ok = FALSE, msg = "Matrix must have at least 2 rows (species)."))

  sp <- trimws(as.character(d[[1]]))
  bad <- which(!nzchar(sp) | is.na(sp))
  if (length(bad) > 0) d <- d[-bad, , drop = FALSE]
  sp <- trimws(as.character(d[[1]]))

  dup <- unique(sp[duplicated(sp)])
  list(ok = TRUE, data = d, species = sp,
       plots = trimws(names(d)[-1]),
       n_species = length(unique(sp)), n_plots = ncol(d) - 1L,
       duplicates = dup,
       msg = sprintf("%d rows x %d columns", nrow(d), ncol(d)))
}

ssdm_read_coords <- function(path) {
  if (!file.exists(path)) return(list(ok = FALSE, msg = paste0("File not found: ", path)))
  d <- try(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           silent = TRUE)
  if (inherits(d, "try-error") || is.null(d)) {
    return(list(ok = FALSE, msg = "Could not read the CSV."))
  }
  if (ncol(d) < 3) {
    return(list(ok = FALSE, msg = "Coordinate file needs at least 3 columns (plot, and two coordinates)."))
  }
  num_cols <- names(d)[vapply(d, function(z) {
    z <- suppressWarnings(as.numeric(as.character(z))); sum(is.finite(z)) > 0
  }, logical(1))]
  num_cols <- setdiff(num_cols, names(d)[1])

  pl <- trimws(as.character(d[[1]]))
  list(ok = TRUE, data = d, plots = pl,
       coord_cols = num_cols,
       duplicate_plots = unique(pl[duplicated(pl)]),
       msg = sprintf("%d rows x %d columns", nrow(d), ncol(d)))
}

# Guess which numeric column is longitude (x) and which is latitude (y)
ssdm_guess_xy <- function(coord_cols, d = NULL) {
  lo <- tolower(coord_cols)
  x <- coord_cols[grepl("^(x|lon|long|longitude|boylam|easting|utm_?x)$", lo)]
  y <- coord_cols[grepl("^(y|lat|latitude|enlem|northing|utm_?y)$", lo)]

  if (length(x) == 0 || length(y) == 0) {
    # Fall back to value ranges: latitude is bounded by +/-90
    if (!is.null(d) && length(coord_cols) >= 2) {
      rng <- vapply(coord_cols, function(cc) {
        v <- suppressWarnings(as.numeric(as.character(d[[cc]])))
        max(abs(v[is.finite(v)]), na.rm = TRUE)
      }, numeric(1))
      if (length(x) == 0 && length(y) == 0 && all(is.finite(rng)) && length(rng) >= 2) {
        o <- order(rng)
        y <- coord_cols[o[1]]   # smaller range -> latitude
        x <- coord_cols[o[2]]
      }
    }
  }
  if (length(x) == 0) x <- coord_cols[min(2L, length(coord_cols))]
  if (length(y) == 0) y <- coord_cols[1]
  list(x = x[1], y = y[1])
}

# Cell value -> 1/0
ssdm_to_binary <- function(v) {
  ch <- trimws(tolower(as.character(v)))
  num <- suppressWarnings(as.numeric(ch))

  out <- integer(length(ch))
  is_num <- is.finite(num)
  out[is_num] <- as.integer(num[is_num] > 0)

  nn <- !is_num
  if (any(nn)) {
    out[nn] <- as.integer(!(ch[nn] %in% BB_ABSENT) & nzchar(ch[nn]) & !is.na(ch[nn]))
  }
  out[is.na(out)] <- 0L
  out
}

# --- main conversion -------------------------------------------------------
ssdm_convert <- function(mat_path, xy_path, x_col = NULL, y_col = NULL) {
  m <- ssdm_read_matrix(mat_path)
  if (!isTRUE(m$ok)) return(list(ok = FALSE, msg = m$msg))
  cc <- ssdm_read_coords(xy_path)
  if (!isTRUE(cc$ok)) return(list(ok = FALSE, msg = cc$msg))

  if (is.null(x_col) || is.null(y_col) || !nzchar(x_col) || !nzchar(y_col)) {
    g <- ssdm_guess_xy(cc$coord_cols, cc$data)
    x_col <- g$x; y_col <- g$y
  }
  if (!(x_col %in% names(cc$data)) || !(y_col %in% names(cc$data))) {
    return(list(ok = FALSE, msg = "Selected coordinate columns not found."))
  }

  mat_plots <- m$plots
  xy_plots  <- cc$plots
  common <- intersect(mat_plots, xy_plots)
  if (length(common) == 0) {
    return(list(ok = FALSE, msg = paste0(
      "No plot names match between the two files. Matrix header example: '",
      paste(utils::head(mat_plots, 3), collapse = ", "),
      "' | coordinate file example: '",
      paste(utils::head(xy_plots, 3), collapse = ", "), "'")))
  }

  missing_xy  <- setdiff(mat_plots, xy_plots)
  missing_mat <- setdiff(xy_plots, mat_plots)

  # Koordinat dosyasinda AYNI ornek alan birden fazla kez gecerse, isimle
  # erisimde R yalnizca ILK satiri dondurur ve digerleri sessizce yok sayilir.
  # Bu, farkina varilmadan yanlis koordinat kullanilmasina yol acabilir.
  dup_plots <- unique(xy_plots[duplicated(xy_plots)])

  xmap <- suppressWarnings(as.numeric(as.character(cc$data[[x_col]])))
  ymap <- suppressWarnings(as.numeric(as.character(cc$data[[y_col]])))
  names(xmap) <- xy_plots; names(ymap) <- xy_plots

  sp_vec <- m$species
  n_sp <- length(sp_vec); n_pl <- length(common)

  species <- rep(sp_vec, each = n_pl)
  plot_id <- rep(common, times = n_sp)

  vals <- integer(n_sp * n_pl)
  for (i in seq_len(n_sp)) {
    rowvals <- ssdm_to_binary(unlist(m$data[i, common, drop = TRUE], use.names = FALSE))
    vals[((i - 1L) * n_pl + 1L):(i * n_pl)] <- rowvals
  }

  out <- data.frame(
    species  = species,
    plot     = plot_id,
    x        = unname(xmap[plot_id]),
    y        = unname(ymap[plot_id]),
    presence = vals,
    stringsAsFactors = FALSE
  )

  bad_xy <- !is.finite(out$x) | !is.finite(out$y)
  n_bad <- sum(bad_xy)
  if (n_bad > 0) out <- out[!bad_xy, , drop = FALSE]

  # Collapse duplicated species rows (same species listed twice): presence wins
  if (length(m$duplicates) > 0) {
    agg <- stats::aggregate(presence ~ species + plot + x + y, data = out, FUN = max)
    out <- agg[order(agg$species, agg$plot), c("species", "plot", "x", "y", "presence"),
               drop = FALSE]
  }

  per_sp <- stats::aggregate(presence ~ species, data = out, FUN = sum)
  names(per_sp)[2] <- "n_presence"
  per_sp$n_absence <- n_pl - per_sp$n_presence

  list(
    ok = TRUE, data = out, per_species = per_sp,
    n_species = length(unique(out$species)),
    n_plots = length(unique(out$plot)),
    n_rows = nrow(out),
    x_col = x_col, y_col = y_col,
    missing_xy = missing_xy, missing_mat = missing_mat,
    duplicate_plots = dup_plots,
    n_bad_xy = n_bad,
    duplicates = m$duplicates,
    msg = sprintf("%d rows x %d columns", nrow(out), ncol(out))
  )
}

# Reproject geographic degrees to the raster CRS if requested
ssdm_reproject <- function(df, from_epsg = 4326, to_crs = NULL) {
  if (is.null(to_crs) || !nzchar(to_crs)) return(df)
  if (!requireNamespace("terra", quietly = TRUE)) return(df)
  v <- terra::vect(df, geom = c("x", "y"), crs = paste0("EPSG:", from_epsg))
  v <- terra::project(v, to_crs)
  xy <- terra::crds(v)
  df$x <- xy[, 1]; df$y <- xy[, 2]
  df
}

# =============================================================================
# FREKANS TABANLI TUR SUZGECI  (yalnizca S-SDM icin)
#
# Frekans (fitososyolojide "constancy"): turun goruldugu ornek alan sayisinin
# toplam ornek alan sayisina orani. Mutlak sayi yerine yuzde kullanmak, farkli
# buyuklukteki veri setleri arasinda karsilastirilabilir bir olcut verir.
#
# min_freq_pct = 0  ->  HICBIR TUR CIKARILMAZ, hepsi modellemeye girer.
#
# NOT: "her alanda bulunan tur" icin ayri bir alt sinir YOKTUR. Gercek yoklugu
# capraz gecerleme kat sayisindan az olan turler icin run_ssdm.R sozde-yokluk
# uretir; boylece o turler de elenmeden modellenebilir.
# =============================================================================

ssdm_filter_species <- function(conv, min_freq_pct = 0) {
  if (is.null(conv) || !isTRUE(conv$ok)) return(NULL)

  n_plots <- conv$n_plots
  ps <- conv$per_species
  ps$frequency_pct <- 100 * ps$n_presence / n_plots

  min_freq_pct <- suppressWarnings(as.numeric(min_freq_pct))
  if (!is.finite(min_freq_pct) || min_freq_pct < 0) min_freq_pct <- 0

  freq_n <- ceiling(min_freq_pct / 100 * n_plots)

  if (min_freq_pct <= 0) {
    ps$keep <- ps$n_presence >= 1
    ps$reason <- ifelse(ps$keep, "", "no presence records")
  } else {
    ps$keep <- ps$frequency_pct >= min_freq_pct & ps$n_presence >= 1
    ps$reason <- ifelse(ps$n_presence < 1, "no presence records",
                 ifelse(ps$keep, "", sprintf("frequency < %g%%", min_freq_pct)))
  }

  # Bilgi amacli: hangi turler sozde-yokluk ile modellenecek (CV.k = 10)
  ps$absence_strategy <- ifelse(ps$n_absence >= 10, "real",
                         ifelse(ps$n_absence > 0, "real+PA", "PA"))

  keep_sp <- ps$species[ps$keep]
  d <- conv$data[conv$data$species %in% keep_sp, , drop = FALSE]
  ps <- ps[order(-ps$frequency_pct, ps$species), , drop = FALSE]

  list(ok = TRUE, data = d, per_species = ps,
       kept = ps[ps$keep, , drop = FALSE],
       dropped = ps[!ps$keep, , drop = FALSE],
       n_total = nrow(ps), n_kept = sum(ps$keep), n_dropped = sum(!ps$keep),
       n_plots = n_plots, n_rows = nrow(d),
       min_freq_pct = min_freq_pct, freq_n = freq_n,
       n_pa = sum(ps$keep & ps$absence_strategy != "real"))
}
