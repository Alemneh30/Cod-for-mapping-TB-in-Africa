Replication code
This repository contains the R code used to reproduce the analyses for:
Liyew AM et al. Mapping tuberculosis prevalence in Africa using a Bayesian geospatial analysis. Communications Medicine. 2025 (5:194).
Folder structure
Run `Runall.R` from the root of the repository; it sources all other scripts
in sequence.
```text
.
├── README.md
├── Runall.R                                 # master script 
├── Mapping\_TB\_Africa\_code.R                 # main analysis
├── Mapping\_TB\_Africa\_robustness\_Ethiopia.R  # exclude-Ethiopia sensitivity
├── Mapping\_TB\_Africa\_robustness.R           # prior- and resolution-sensitivity
├── Analysis\_reported\_vs\_updated.R           # reported vs. re-analysed counts
├── Plots\_TB\_Africa.R                        # main and supplementary figures
├── INPUT/
└── OUTPUT/
```
To access the `INPUT` data access this drive: [storage link]
Unzip the folder and save the unzipped folder`INPUT`in your working directory. The files in `INPUT`include:
```text
TB\_Africa\_v2.xlsx        # TB dataset
vif.R                    # code to run the VIF function
originalcovariates.Rdata # covariate data
supdata3.csv             # reported ADM1-level TB counts

# several files for each country (ISO-3 code given for each country)
gadm36\_ISO-3\_0\_sp.rds     # spatial polygon delimiting ADM-0 border 
gadm36\_ISO-3\_1\_sp.rds     # spatial polygon delimiting ADM-1 border
```
Outputs are organised according to the sensitivity-analysis specification:
```text
OUTPUT/
├── ALL/
│   ├── FILTER/
│   ├── NOFILTER/
│   └── NOFILTER/robustness/
│   └── NOFILTER/reported\_vs\_reanalysis/
├── NOMOZ/
│   ├── FILTER/
│   └── NOFILTER/
└── robustness\_Ethiopia/
└── robustness\_Mozambique/
```
Within these folders, the analysis generates `csv`, `pdf`, `shp`, and `tif` outputs.
R packages
The analysis requires R and R-INLA. The principal packages are:
```r
INLA; raster
viridis; geodata
rnaturalearth; rnaturalearthdata
malariaAtlas; readxl
ggplot2; ggrepel; RColorBrewer
ggmap; tmap; gtools; fmsb
dplyr; sf; scales
```
The original code also used some legacy spatial packages (`rgdal`, `rgeos`, and `maptools`), since superseded by `sf`, `geodata`, and current `terra, raster`functionality.
R-INLA and fmesher versions
Different R-INLA versions may produce small differences in the estimates and WAIC. We recommend the following package versions:
```r
# R-INLA 24.05.10
install.packages("\~/Downloads/INLA\_24.05.10.tgz", repos = NULL, type = "mac.binary")
# from https://inla.r-inla-download.org/R/stable

# fmesher 0.1.7, compatible with INLA 24.05.10
# from https://cran.ms.unimelb.edu.au/src/contrib/Archive/fmesher/
install.packages("\~/Downloads/fmesher\_0.1.7.tar.gz", repos = NULL, type = "source")
```
Deprecated covariate-download functions
Some covariates were originally obtained using `raster::getData()`, since
deprecated in favour of the `geodata` package. Re-downloading covariates
through their current equivalents can introduce small differences from the
original inputs (e.g. updated source data or default resolution).
Running the analysis
The main analysis settings are defined at the top of `Runall.R`:
```r
nn            <- 10000   # posterior samples; use 10000 for final results
mainaggfactor <- 2       # prediction-grid aggregation factor (see below)
mu            <- 0.025   # logistic-regression model parameter
prevunit      <- 1000    # prevalence scaling (per 1,000)
popt          <- 5       # population-density threshold for filtering
allrun        <- TRUE    # run all sensitivity configurations
```
Change these values in `Runall.R` only; the individual scripts read them
from the calling environment and fall back to the same defaults if run
independently.
With `allrun <- TRUE`, the main analysis (`Mapping\_TB\_Africa\_code.R`) runs
the four combinations:
Data                           Population filter     Output folder
---
All observations               No                   `ALL/NOFILTER`
All observations               Yes           `       ALL/FILTER`
Mozambique excluded            No            `       NOMOZ/NOFILTER`
Mozambique excluded            Yes           `       NOMOZ/FILTER`
The population filter excludes areas with population density below 5
persons per square kilometre.
Prediction grid
The covariate raster objects used to gather covariates to fit the model are at 10 arc-minute resolution. For predictions, the grid is aggregated by `mainaggfactor <- 2`. For the main specification, predictions are therefore generated on a **20 × 20 arc-minute grid ** (~37 km north–south at the equator).
Reproducing the analysis
For a complete replication, including all sensitivity analyses:
Clone or download this repository.
Place the required input files in `INPUT/`.
Start a clean R session with the repository root as the working
directory, using the R-INLA/fmesher versions noted above.
Set `allrun <- TRUE` in `Runall.R` (default) to reproduce the main
analysis and all sensitivity checks.
Run `Runall.R` from beginning to end. It sources, in order:
`Mapping\_TB\_Africa\_code.R` (main analysis)
`Mapping\_TB\_Africa\_robustness\_Ethiopia.R` (exclude-Ethiopia check)
`Mapping\_TB\_Africa\_robustness.R` (prior/resolution sensitivity)
`Analysis\_reported\_vs\_updated.R` (reported vs. updated Table 1/Figure 2)
`Plots\_TB\_Africa.R` (figures).
Check the generated files under `OUTPUT/`.
The workflow produces the fitted-model results, model diagnostics,
prevalence predictions, uncertainty summaries, estimated TB cases,
administrative-level summaries, maps, and all sensitivity analyses.
The workflow has been tested (runtime about 12 minutes) on a MacBook Pro, M4 Max (16 processors), 64 GB Memory with R 4.4.3.
Citation
If using this code, please cite:
> Liyew AM et al. Mapping tuberculosis
> prevalence in Africa using a Bayesian geospatial analysis.
> \*Communications Medicine\*. 2025 (5:194).
