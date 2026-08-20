# biomod2 SDM Runner

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22036156.svg)](https://doi.org/10.5281/zenodo.22036156)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://cran.r-project.org/)
[![biomod2 >= 4.3](https://img.shields.io/badge/biomod2-%3E%3D4.3-brightgreen.svg)](https://biomodhub.github.io/biomod2/)

**Run biomod2 species distribution models from a graphical interface — no code required.**

Point it at a folder of environmental rasters and a CSV of records, tick the
algorithms you want, and press Run. The application handles data formatting,
parameter selection, cross-validation, ensemble modelling, thresholding and
future projections, then writes every setting it used so the run can be reported.

Interface available in **English** and **Turkish**.

<!-- Ekran goruntusu ekleyin: images/ klasoru olusturup dosyayi yukleyin, sonra:
![Interface](images/screenshot.png)
-->

---

## Getting started

### Step 1 — Download

Click the green **Code** button at the top of this page, then
**Download ZIP**. Unzip it anywhere on your computer, for example
`C:\SDM` or `~/Documents/SDM`.

> **Keep all files in one folder.** The application looks for its companion
> scripts next to `app.R`. Do not move individual files into subfolders.
>
> **Avoid cloud-synced folders** (OneDrive, Google Drive, Dropbox). Files there
> can be stored online only, and R will appear to freeze while it downloads them.

### Step 2 — Point R at that folder

Open R or RStudio and set the working directory to the folder you just unzipped:

```r
setwd("C:/SDM")        # Windows: forward slashes, not backslashes
# setwd("~/Documents/SDM")   # macOS / Linux
```

In RStudio you can also use **Session → Set Working Directory → Choose Directory**.

Check that you are in the right place:

```r
list.files(pattern = "[.]R$")
```

You should see ten `.R` files. If the list is empty, `setwd()` pointed somewhere else.

### Step 3 — Install the required packages

```r
source("install_packages.R")
```

This runs once and may take several minutes the first time. It **skips anything
already installed**, so it is safe to run again later.

When it finishes it prints a summary. Check three things:

| Line | What you want to see |
|---|---|
| `biomod2` | version **4.3** or higher |
| `Rscript` | a file path, not `NOT FOUND` |
| Model availability | `available` next to the algorithms you plan to use |

If biomod2 is older than 4.3, continuous and nominal response types will not
work. Update it with:

```r
install.packages("remotes")
remotes::install_github("biomodhub/biomod2")
```

`java: NOT FOUND` only matters if you want MAXENT. MAXNET does the same job in
pure R and needs no Java.

### Step 4 — Check that the files are consistent (optional)

```r
Rscript check_code.R
```

or, inside R:

```r
source("check_code.R")
```

Expected output:

```
115 fonksiyon denetlendi. Tanimsiz degisken referansi YOK.
```

This takes a couple of seconds and confirms that no script file is missing or
out of date. Useful after updating to a new version. If it reports
`EKSIK DOSYA` (missing file), copy that file into the folder.

### Step 5 — Start the application

```r
shiny::runApp(".")
```

A browser window opens. The panel at the top left repeats the checks from
Step 3 — everything should be green before you continue.

---

## Try it with the example data

The repository ships with a small synthetic dataset so you can see the whole
workflow before preparing your own files.

In the interface, fill in:

| Field | Value |
|---|---|
| Current raster folder | full path ending in `example_data/current` |
| Occurrence CSV | full path ending in `example_data/occurrences.csv` |
| Presence/absence column | `presence` |
| Predictor type for `bedrock` | `categorical` |
| Models | GLM, GAM, RF, MAXNET |

Paths must be **full paths**, for example `C:/SDM/example_data/current`. You can
paste them, or use the browse buttons next to each field.

Press each **Check** button before running: they confirm the folder can be read
and the species name is detected.

Then press **Run biomod2**. Progress appears in the Log tab, and maps and metrics
fill in as they are produced. The example takes a few minutes with four
algorithms.

To try the species richness workflow instead, tick **Enable stacked species
distribution modelling** and supply `example_data/species_matrix.csv` and
`example_data/plot_coords.csv`.

---

## Which script do I run?

**You only ever start `app.R`.** Everything else is either a one-off setup step
or is launched for you in the background.

| Script | When | Do you run it? |
|---|---|---|
| `install_packages.R` | once, before first use | **yes** — Step 3 |
| `check_code.R` | after downloading or updating | optional — Step 4 |
| `app.R` | every session | **yes** — Step 5 |
| `run_biomod.R` | presence/absence runs | no, launched automatically |
| `run_biomod_cont.R` | continuous / nominal runs | no, launched automatically |
| `run_ssdm.R` | stacked SDM runs | no, launched automatically |
| `biomod_opts.R`, `ssdm_convert.R`, `plot_response_threshold.R`, `i18n.R` | loaded as needed | no |

Modelling runs in a **separate background process**, so the interface stays
responsive. Progress appears in the Log tab while it works, and you can close the
browser tab without stopping the run.

---

## If something goes wrong

| Symptom | Cause |
|---|---|
| `could not find function` | Wrong working directory, or a file missing — run `check_code.R` |
| Status stays at one stage | Check the Log tab; large `.asc` files with categorical variables are slow to read |
| Log says "empty or not readable" | Normal in the first seconds; the job is still starting |
| Species name not detected | The CSV needs a species column such as `Tur` or `species` |
| Coordinates rejected | Use a dot as the decimal separator (`30.25`, not `30,25`) |
| Nothing happens after Run | `Rscript` was not found — see the system panel at the top left |

The file `runs/<run_id>/progress.log` records every step with a timestamp and is
the first place to look.

---

## What it can do

The application supports **three ways of describing where a species is**, and a
fourth mode that combines many species at once.

### 1. Presence / absence

The classic case. Your CSV lists coordinates where the species was recorded.
If you also recorded where it was **not** found, select that column and those
real absences are used; otherwise pseudo-absences are generated for you.

Output: habitat suitability maps (0–1) and binary maps thresholded at MaxTSS.

### 2. Abundance and other continuous measurements

Cover percentage, stem counts, biomass, basal area. biomod2 4.3 models these
directly as regression rather than forcing them into presence/absence, and the
application picks the right distribution family for you:

| Your data | Sub-type | Distribution |
|---|---|---|
| Counts (0, 1, 5, 23) | `count` | Poisson |
| Continuous positive values | `abundance` | Gaussian |
| Proportions between 0 and 1 | `relative` | Beta |

Metrics switch to R², RMSE and MAE automatically.

### 3. Nominal classes

Vegetation types, condition classes, ordered categories such as
*sparse → moderate → dense*. Modelled as classification, evaluated with
Accuracy and F1.

### 4. Many species at once → species richness

This is where the pieces come together. Vegetation survey data usually looks
like a **species × plot matrix**:

```
Species,        P1, P2, P3, P4
Cedrus libani,   4,  0, 12,  1
Abies cilicica,  0, 56, 10,  0
Pinus brutia,   75,  0,  0, 25
```

Supply this together with a plot-coordinate file, and the application:

1. **Converts** the matrix into the long format biomod2 needs, turning abundance
   or cover values into 1/0 (Braun-Blanquet codes such as `+`, `r`, `2a` are
   recognised);
2. **Models every species separately** with the algorithms you selected;
3. **Thresholds** each species map at its own MaxTSS cut-off;
4. **Stacks** the binary maps into **potential species richness** — for current
   conditions and for every future scenario you supply, plus change maps.

Zeros in the matrix are treated as real absences. Species with too few records
can be filtered out by frequency (percentage of plots), and species found in
almost every plot are handled with pseudo-absences rather than being discarded.

---

## Preparing your data

### Environmental rasters

One folder per time period, one file per variable, in `.asc`, `.tif` or `.grd`.
All layers must share the same extent and resolution. Future scenario folders
must contain the **same variable names** as the current folder.

After pressing *Check current raster folder* you can mark each variable as
continuous, categorical or nominal.

### Occurrence CSV

```
Tur,x,y,presence
Cedrus libani,312500,4145500,1
Cedrus libani,318500,4132500,0
```

- Coordinate columns do not have to be named `x` and `y`. `X/Y`, `lon/lat`,
  `boylam/enlem`, `POINT_X/POINT_Y` and `utm_x/utm_y` are recognised, and you can
  always pick them manually.
- Coordinates must use a **dot** as the decimal separator.
- Coordinates must be in the same coordinate system as the rasters.
- If you use a species column, every row — including absences — must carry the
  species name, otherwise those rows are filtered out.

### Species matrix (for richness modelling)

Two files: the matrix shown above, and a coordinate file whose first column
matches the matrix header.

```
plot,x,y
P1,312500,4145500
P2,318500,4132500
```

Each plot must appear **once**. The application reports missing and duplicated
plots before you run anything.

---

## What you get

Everything is written to `results/<run_id>/`:

**Maps** (`rasters/`)
- continuous ensemble suitability, current and each future scenario
- binary maps thresholded at MaxTSS
- uncertainty layers (EMcv, EMca)
- species richness and richness change maps, if stacked SDM was used

**Figures** (`plots/`)
- model evaluation, variable importance
- response curves, including a publication-format version where sections above
  the threshold are drawn in red and categorical variables as bars

**Tables**
| File | Contents |
|---|---|
| `auc_tss_summary.csv` | metrics per algorithm, and which entered the ensemble |
| `ensemble_maxTSS_cutoffs.csv` | thresholds on both the 0–1000 and 0–1 scales |
| `modeling_options_applied.csv` | every parameter actually used |
| `species_summary.csv` | per-species results (stacked SDM) |
| `sessionInfo.txt`, `package_versions.csv` | software versions used |
| `progress.log` | timestamped run log |

The last three exist so that a run can be reported under the
[ODMAP protocol](https://doi.org/10.1111/ecog.04960).

---

## Requirements

- R >= 4.1 (tested on 4.6)
- **biomod2 >= 4.3** — required for continuous and nominal response types
- Java and `maxent.jar` only for MAXENT; MAXNET needs neither

All of this is verified by `install_packages.R` and again each time the app starts.

## Good to know

- **Runtime scales with species count.** Each species is a full biomod2 run.
  Test with a few species before launching thirty.
- **Parallel processing is off by default on Windows**, where biomod2 workers can
  stall silently. Confirm a run works with one core before increasing it.
- **Keep your files on a local disk.** Rasters inside OneDrive or Google Drive
  folders can appear to hang while they download.
- **Thresholds stay fixed across scenarios**, so differences between current and
  future maps reflect habitat change rather than a moving cut-off.

---

## What this is — and what it is not

This is a **graphical front end for biomod2**. All modelling, evaluation and
projection is done by biomod2 itself; no new algorithm or statistical method is
introduced here. What the application adds is convenience and reproducibility:
parameters scaled to your data following published SDM guidance, and a complete
record of what was applied.

**If you use this application, please cite biomod2:**

> Gueguen M, Blancheteau H, Lemaire-Patin R, Thuiller W (2026).
> *biomod2: Ensemble Platform for Species Distribution Modeling.*
> R package version 4.3-4-7. https://biomodhub.github.io/biomod2/

Citing this repository as well is welcome — see `CITATION.cff` or the DOI above.

---

## Related software

| Tool | Focus |
|---|---|
| [wallace](https://wallaceecomod.github.io/wallace/) | modular SDM workflow, data acquisition to projection |
| [ShinyBIOMOD](https://gitlab.com/IanOndo/shinybiomod) | biomod2 interface; 2020 GBIF Ebbe Nielsen Challenge winner |
| [SSDM](https://cran.r-project.org/package=SSDM) | stacked SDM with a `gui()` interface |
| [EcoNicheS](https://github.com/armandosunny/EcoNicheS) | niche modelling, overlap and connectivity |
| [eSDM](https://swfsc.github.io/eSDM/) | ensembles of existing SDM predictions |
| [sdm](https://cran.r-project.org/package=sdm) | modelling framework with a `gui()` |

This application differs mainly in exposing biomod2's **continuous and nominal
response types** through a graphical interface, in converting species × plot
matrices for stacked modelling, and in recording applied parameters for reporting.

---

## Files

| File | Purpose |
|---|---|
| `app.R` | Shiny interface |
| `run_biomod.R` | presence/absence pipeline |
| `run_biomod_cont.R` | continuous and nominal pipeline |
| `run_ssdm.R` | stacked SDM (species richness) |
| `biomod_opts.R` | parameter engine and shared helpers |
| `plot_response_threshold.R` | threshold-coloured response curves |
| `ssdm_convert.R` | species matrix conversion |
| `i18n.R` | interface dictionary (EN / TR) |
| `install_packages.R` | dependency installer |
| `check_code.R` | static code check |
| `example_data/` | runnable example dataset |

---

## Parameter sources

Defaults follow published guidance rather than package defaults:
Elith et al. (2008) for GBM; Breiman (2001) and Valavi et al. (2021) for random
forests; Phillips & Dudík (2008) and Radosavljevic & Anderson (2014) for MaxEnt
feature classes and regularisation; Wood (2011) for GAM smoothing;
Austin (2007) for GLM formulation; Leathwick et al. (2006) for MARS;
Barbet-Massin et al. (2012) for pseudo-absence weighting;
Chen & Guestrin (2016) for XGBoost.

Stacked SDM follows Guisan & Rahbek (2011), Calabrese et al. (2014) and
D'Amen et al. (2015). Reporting follows Zurell et al. (2020).

Full references are listed in [`REFERENCES.md`](REFERENCES.md).

## License

GNU General Public License v3.0 or later — see `LICENSE`.

## Contributing

Bug reports welcome through the issue tracker; see `CONTRIBUTING.md`.
