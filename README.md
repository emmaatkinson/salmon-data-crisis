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

This work has been submitted to the journal *Ecological Indicators*.

## Contents of this repo

#### Files for compiling and checking data

-   `dataCompilationNuSEDS_CU.R`: imports [NuSEDS
    data](https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6)
    and several other sources described in the paper to yield a filtered
    dataset of trends in spawner abundance per population and habitat
    pressure values associated with spawning watersheds of those
    populations. The “NuSEDS_CU” part refers to the NuSEDS data by CU
    that were used in this compilation.

-   `functions.R`: contains several basic functions used repeatedly in
    data compilation.

#### Files for fitting Bayesian linear model in JAGS and visualizing results

-   `fitting.R`: Code for fitting the model above to the
    population/habitat data compiled by `dataCompilationNuSEDS_CU.R`
    using JAGS and rjags. Note that this code involves parallel
    computation.

-   `loadPopData.R`: Called by `fitting.R` and `lookingAtOutput.R` to
    import the population/habitat data and create additional variables
    for factor levels, etc.

-   `loadResults.R`: Called by `lookingAtOutput.R` to load MCMC output
    of model fitting and organize parameter output in arrays.

-   `lookingAtOutput.R`: Code for plotting the results of the model
    object from `fitting.R` and producing figures in the report.

-   `threatened.R`: Code that takes the model object from `fitting.R`
    and samples from the posterior and data to calculate the proportion
    of populations of different salmon species and Freshwater Adaptive
    Zones (FAZs) throughout BC that would be classified as threatened
    with increasing pressure values.

## More information

Visit the [Salmon Watersheds Program project
page](https://salmonwatersheds.ca/projects/assessing-the-vulnerability-of-pacific-salmon-to-freshwater-habitat-pressures/)
and check out the [Pacific Salmon Explorer](www.salmonexplorer.ca).
