README
================
Stephanie Peacock
06/09/2022

# Testing for broad-scale relationships between freshwater habitat pressure indicators and Pacific salmon population trends

This repository contains code and data associated with an analysis to
quantify the relationships between population status of Pacific salmon
and habitat indicators. The project is a collaboration between DFO and
the Pacific Salmon Foundation, building on the PSF’s recently updated
[Pacific Salmon Explorer](www.salmonexplorer.ca) and associated
datasets.

This work has been published in the journal *Ecological Indicators* and
is available open-access from
<https://www.sciencedirect.com/science/article/pii/S1470160X23000778>.

## Contents of this repo

### Data

Contains *most* datasets in .csv format that are read in by
`code/dataCompilationNuSEDS_CU.R` to compile a single dataset of
population trends, habitat pressure, and other explanatory variables.

Population data were too large to upload to GitHub but can be downloaded
at the following links:

- NuSEDS spawner estimates for each population and year, including data
  quality `Estimate Classification Types` are available
  [here](https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6/resource/f2f34577-600e-46f1-b015-8e9991ced008).

- Conservation_Unit_Data.csv, which organizes NuSEDS data by
  Conservation Units, is available
  [here](https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6/resource/712e87ce-a2b0-4f42-9c17-fdb8d07b8bc8).
  (Note that this dataset does not include some of the fields in the
  above NuSEDS file, hence the need to import both.)

NuSEDS is currently undergoing some restructuring and we cannot
guarantee that the format of the above data will be compatible with the
code in this repo. Email us at the address below if you’d like the exact
files used in our analysis.

The habitat pressure data are also publicly available from the [Salmon
Data Library](https://data.salmonwatersheds.ca/data-library/), as
“spatial datasets” separated by pressure indicator and region.

### Code

#### Files for compiling and checking data

- `dataCompilationNuSEDS_CU.R`: imports [NuSEDS
  data](https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6)
  and several other sources described in the paper (provided in `data/`)
  to yield a filtered dataset of trends in spawner abundance per
  population and habitat pressure values associated with spawning
  watersheds of those populations. The “NuSEDS_CU” part refers to the
  NuSEDS data by CU that were used in this compilation.

- `functions.R`: contains several basic functions used repeatedly in
  data compilation.

#### Files for fitting Bayesian linear model in JAGS and visualizing results

- `fitting.R`: Code for fitting the model above to the
  population/habitat data compiled by `dataCompilationNuSEDS_CU.R` using
  JAGS and rjags. Note that this code involves parallel computation.

- `loadPopData.R`: Called by `fitting.R` and `lookingAtOutput.R` to
  import the population/habitat data and create additional variables for
  factor levels, etc.

- `loadResults.R`: Called by `lookingAtOutput.R` to load MCMC output of
  model fitting and organize parameter output in arrays.

- `lookingAtOutput.R`: Code for plotting the results of the model object
  from `fitting.R` and producing figures in the report.

- `threatened.R`: Code that takes the model object from `fitting.R` and
  samples from the posterior and data to calculate the proportion of
  populations of different salmon species and Freshwater Adaptive Zones
  (FAZs) throughout BC that would be classified as threatened with
  increasing pressure values.

### Simulations

This repo contains code that runs a simple population simulation to
explore how changing harvest and productivity through time may affect
the observed trend in spawner abundance. This simulation was not
described in the paper.

## More information

Visit the [Salmon Watersheds Program project
page](https://salmonwatersheds.ca/projects/assessing-the-vulnerability-of-pacific-salmon-to-freshwater-habitat-pressures/)
and check out the [Pacific Salmon Explorer](www.salmonexplorer.ca).
Questions can be directed at Steph Peacock (speacock at psf dot ca).
