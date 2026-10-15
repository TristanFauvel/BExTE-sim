# BExTEsim: simulation studies of Bayesian extrapolation methods

BExTEsim runs simulation studies of the borrowing methods in
[BExTE](https://github.com/TristanFauvel/BExTE). BExTE computes the operating
characteristics of one design. BExTEsim compares methods over a grid of
scenarios:

- It reads YAML configurations of case studies, scenarios and methods.
- It simulates the frequentist and Bayesian operating characteristics of each
  method in each scenario, in parallel if asked.
- It writes the results to files, then summarises them in plots and tables.

## Installation

```r
remotes::install_github("TristanFauvel/BExTE")
remotes::install_github("TristanFauvel/BExTE-sim")
```

`library(BExTEsim)` also attaches BExTE.

## Run a simulation study

The default configurations are in `inst/conf`:

- `full/` holds the scenarios (`scenarios_config.yml`), the methods
  (`methods_config.R`) and the MCMC settings (`mcmc_config.yml`);
- `simulation_config.yml`, `analysis_config.yml` and `outputs_config.yml` set
  the simulation, the analysis steps and the names of the result files;
- `plots_config.R`, `methods_plots_config.R` and `metrics_config.R` set the
  style of the plots and tables.

The case-study configurations are in BExTE (`system.file("conf/case_studies",
package = "BExTE")`).

`run_simulation_env()` runs the simulation study of one environment. It writes
its results to `results/<env>/` and its logs to `logs/<env>/` in the working
directory. `simulation_analysis()` then computes the sweet spots and the
Bayesian operating characteristics.

## Related projects

- [BExTE-app](https://github.com/TristanFauvel/BExTE-app) is a Shiny app to
  configure, run and analyse a simulation study from a browser.
- BExTEpaper holds the scripts that reproduce the paper's simulation study,
  figures and tables.

## Development

Development uses the renv library of the BExTE checkout next to this folder
(`../BExTE/renv`), as BExTEpaper does. Install the BExTE checkout into that
library before you build or test BExTEsim.
