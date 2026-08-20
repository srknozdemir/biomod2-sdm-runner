# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/).

## [1.0.0] - 2026-08-20

First public release.

### Added
- Shiny interface for the biomod2 workflow.
- Three response types: presence/absence, continuous (count / abundance /
  relative) and nominal (multiclass / ordinal).
- Literature-based parameter engine scaled to sample size and predictor count;
  every applied setting exported to `modeling_options_applied.csv`.
- Stacked SDM producing potential species richness for current and future
  conditions, with frequency-based species filtering.
- Threshold-coloured response curves in publication format (PNG + 300 dpi TIFF).
- Automatic absence strategy: real absences when sufficient, pseudo-absences
  otherwise.
- Selectable raster output scale (0-1 or biomod2 native 0-1000); binary maps are
  always derived from the native scale so results are unaffected.
- Bilingual interface (English / Turkish).
- Session and package version capture on every run.
- `install_packages.R` dependency installer and `check_code.R` static checker.
- Runnable synthetic example dataset.
