# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/).

## [1.1.0] - 2026-10-26

### Added
- **Automatic future scenarios.** Point the app at one folder that holds
  `future_<period>_<ssp>_<GCM>` sub-folders (the layout written by common climate
  preparation scripts), press *Scan folder*, and tick the periods, SSPs and GCMs
  to project. The folder-name pattern is editable; folders named `ENSEMBLE` are
  skipped by default. Manual entry of individual folders is still available.
- **GCM ensemble** (`gcm_ensemble.R`). Each GCM is projected separately and the
  predictions are averaged per period and SSP afterwards, keeping the SD among
  GCMs and a GCM-agreement map as explicit uncertainty layers. The ensemble layer
  is selected by algorithm *and* selection metric (`change_em`, default
  `EMwmeanBy<first selection metric>`, e.g. `EMwmeanByRsquared`).
- **Change maps** for continuous responses and habitat suitability: future −
  current, percentage change, GCM agreement on the sign of change, and a
  five-class map (robust / low-agreement decrease, no change, low-agreement /
  robust increase), plus per-GCM change maps.
- **Range change** for presence/absence runs: gain / loss / stable maps per GCM
  and for the GCM consensus (majority vote), coded as in
  `biomod2::BIOMOD_RangeSize`, and `range_change_summary.csv` with Loss, Gain,
  PercLoss, PercGain and SpeciesRangeChange.
- `change_summary.csv`: area-weighted mean change and the area (km²) and share
  of the study area that increases, decreases or stays unchanged, per GCM and
  for the ensemble.
- Optional polygon or raster mask (e.g. a catchment): summaries are repeated
  inside it.
- Optional GCM-exclusion sensitivity ensemble (e.g. without a "hot" model).
- Overview figure of projected change across scenarios (`change_overview.png`).
- `ensemble_members.csv`: the single models kept in the ensemble and their score.
- `ENSEMBLE_BELOW_THRESHOLD.txt` and a log warning when no algorithm reaches the
  ensemble thresholds on average. In that case all algorithms are passed to
  `BIOMOD_EnsembleModeling()`, which keeps only the single cross-validation
  models that pass on their own; the ensemble may therefore rest on very few
  models and should be treated as indicative.

### Changed
- `run_biomod.R`: a future folder that cannot be read is now skipped with a
  warning instead of stopping the whole run.

### Fixed
- Continuous runs (`count`, `abundance`): ensemble projections were rounded to
  whole numbers because `BIOMOD_EnsembleForecasting()` uses `digits = 0` by
  default. The ensemble now uses the same number of decimals as the single-model
  projections (`digits` setting). `relative` responses were not affected.
  Re-run continuous `count`/`abundance` models produced with 1.0.0.
- Threshold-coloured response curves picked the ensemble curve by algorithm name
  only (`EMwmean`), which matched the RMSE-filtered ensemble (`EMwmeanByRMSE`)
  before the R²-filtered one. The curve is now chosen by algorithm and selection
  metric (`rc_em_pick`, default `EMwmeanBy<first selection metric>`).

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
