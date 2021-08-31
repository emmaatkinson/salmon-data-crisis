README
================
Stephanie Peacock
09/08/2021

# Quantifying linkages between habitat and population status to aid recovery of at risk Pacific salmon populations

This repository contains code and data associated with an analysis to
quantify the relationships between population status of Pacific salmon
and habitat indicators. The project is a collaboration between DFO and
the Pacific Salmon Foundation, building on the PSF’s recently updated
[Pacific Salmon Explorer](www.salmonexplorer.ca).

This analysis was undertaken by Stephanie Peacock and all questions
about the code in this repository should be directed to stephanie dot j
dot peacock at gmail dot com.

## Contents

Files for compiling and checking data

-   `dataCompilationNuSEDS_CU.R`: imports NuSEDS data and several other
    sources described in the report to yield a filtered dataset of
    trends in spawner abundance per population and habitat pressure
    values associated with spawning watersheds of those populations. The
    “NuSEDS\_CU” part refers to the NuSEDS data by CU that were used in
    this compilation.

-   `functions.R`: contains several basic functions used repeatedly in
    data compilation.

Files for fitting Bayesian linear model in JAGS and visualizing results

-   `models.R`: Contains multiple JAGS models for different versions of
    the model linking population trends (response) to habitat pressure
    values (explanatory variables). The final model used in the report
    is “model10”.

-   `fitting.R`: Code for fitting the model above to the
    population/habitat data compiled by `dataCompilationNuSEDS_CU.R`
    using JAGS and rjags. Note that this code involves parallel
    computation.

-   `loadPopData.R`: Called by `fitting.R` and `lookingAtOutput.R` to
    import the population/habitat data and create additional variables
    for factor levels, etc.

-   `lookingAtOutput.R`: Code for plotting the results of the model
    object from `fitting.R` and producing figures in the report.

-   `vulnerability.R`: Code that takes the model object from `fitting.R`
    and samples from the posterior and data to calculate the
    sensitivity, exposure, and vulnerability of different salmon species
    and Freshwater Adaptive Zones (FAZs) throughout BC. Produces large
    dot plots for main report showing these metrics.

-   Simulations - contains several files for initial simulation analysis
    of the potential confounding factors of fishing and density
    dependence on trend detection from spawner data alone, as well as
    testing the model fitting to simulated data to ensure known
    parameters could be recovered.

## Background

In British Columbia (BC) 27 Pacific Designatable Units (DU) have been
assessed as either Threatened (7) or Endangered (20) by the Committee on
the Status of Endangered Wildlife in Canada (COSEWIC). These include DUs
from four species: sockeye salmon (Oncorhynchus nerka), Chinook salmon
(O. tshawytscha) and coho salmon (O. kisutch) and steelhead trout (O.
mykiss). The role that degradation of freshwater habitats has played in
contributing to these at risk states is often implied, but is rarely
rigorously quantified, resulting in high uncertainty about which changes
to freshwater systems contribute most to population declines. For
example, recent DFO threats assessments have documented uncertainty in
the degree that threats to freshwater habitat have contributed to
population declines and are preventing recoveries (DFO 2018a, 2018b,
2019a, 2019b, 2020a). There is also substantial uncertainty on the
degree that the impact of changes in habitat quality and quantity have
on life-history parameters (DFO 2020b). These represent substantial
knowledge gaps for Pacific salmonids. This information is critical to
informing area- and threat-based recovery actions and has been
identified as one of the priority strategies within the Wild Salmon
Policy (DFO 2005) and its associated implementation plan (DFO 2018).

Over the past 5 years, DFO and the Pacific Salmon Foundation (PSF) have
collaborated to compile and synthesize information on population and
habitat status for Pacific salmon in BC (Pacific Salmon Foundation
2020). This extensive dataset on population and habitat status presents
a unique opportunity to quantify linkages between the population status
of Pacific salmon and landscape level pressures on their freshwater
habitats. Quantifying these linkages allows for exploration of trends
and differences across species, evaluation of the extent to which
life-history characteristics and/or biogeoclimatic factors are driving
the observed differences across species, and most importantly,
contributions to the identification of recovery actions for mitigating
habitat threats.

## Method overview

Population status was quantified as the trend in spawner abundance using
data from DFO’s new Salmon Escapement Database System (nuSEDS; Fisheries
and Oceans Canada 2020). We removed records prior to 1950, populations
affected by enhancement from hatcheries or artificial channels, and
populations with fewer than 10 years of data. For the remaining
populations, spawner estimates were smoothed using the geometric mean
over the generation length (Pacific Salmon Foundation 2020) and then
log-transformed prior to estimating the linear trend over time
(D’Eon-Eggertson et al. 2015).

Habitat pressure values were calculated by watershed as reported in the
Pacific Salmon Explorer (Pacific Salmon Foundation 2020). We
disaggregated several habitat pressures (e.g., land-cover alteration) to
better reflect pressures that may have different mechanisms of impact
(e.g., agriculture versus urban land cover). The 10 indicators that we
used (Table, below) were not strongly correlated with each other and
represented different types of habitat degradation. The pressure values
for each indicator were standardized by subtracting the mean and
dividing by the standard deviation for that indicator among all
watersheds prior to model fitting.

We used Bayesian hierarchical linear models to assess the relationship
between population trends and habitat pressures within the watersheds in
which those populations spawn. The intercept (i.e., the trend in spawner
abundance in the absence of any habitat pressures) included a random
effect for the population rearing ecotype nested within marine adaptive
zone (MAZ; Holtby and Ciruna 2007) to account for different trends among
populations due to factors operating outside of the spawning watershed.
For example, declines in fishing pressure may result in apparent
increases in spawner abundance (see the files in `Simulation` for an
exploration of these potential confounding factors), and we attempted to
reduce this confounding effect by accounting for shared variation among
populations within the same MAZ. The “rearing ecotypes” we considered
were ocean-type Chinook, stream-type Chinook, chum, coho, pink,
lake-type sockeye, and river-type sockeye. This categorization allowed
for shared variation in survival outside of the spawning watershed,
including marine effects and freshwater residence of juveniles.

The magnitude of effect for each habitat indicator depended on the
“spawning ecotype”, which we defined to group populations that may be
affected similarly by changes to spawning habitat (Chinook, coho,
pink/chum, sockeye). The model also included interactions between each
habitat indicator and the system size. Size was determined based on
Strahler stream order (Hughes et al. 2011). Finally, we included a
random effect of freshwater adaptive zone (FAZ; Holtby and Ciruna 2007)
on the slope of each habitat indicator to account for the potential
mediating effect of general climactic and environmental conditions on
the relationship between population trend and habitat indicators.

| Indicator                        | Pressure value   | Description                                                                                                                                                                                                                                                               |
|----------------------------------|------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Agricultural development         | % watershed area | The percentage of the total watershed area that has been altered by agricultural/rural land use within the total land-cover alteration layer.                                                                                                                             |
| Urban development                | % watershed area | The percentage of the total watershed area that has been altered by urban land use within the total land-cover alteration layer.                                                                                                                                          |
| Riparian disturbance             | % watershed area | The percentage of the total watershed area that has been altered by human activity (forest disturbance, urban land use, agricultural/rural land use, mining development, and other development) within a 30m buffer zone around all streams, rivers, lakes, and wetlands. |
| Linear development               | km/km2           | The density of linear developments within a watershed, excluding roads (considered separately), including railways, utility corridors, pipelines, power lines, telecom cables, right of ways, etc.                                                                        |
| Forestry road density            | km/km2           | The average density of forestry roads within a watershed from the Forest Tenure Road data\*                                                                                                                                                                               |
| Non-forestry road density        | km/km2           | The average density of non-forestry roads (e.g., highways) within a watershed from the Digital Roads Atlas.                                                                                                                                                               |
| Stream crossing density          | \#/km            | The total number of stream crossings per km of the total length of modelled salmon habitat in a watershed. Salmon habitat is defined based on a gradient criterion filtering of the Fish Passage Model (Mount et al. 2011).                                               |
| Forest disturbance\*             | % watershed area | The percentage of total watershed area that has been disturbed by logging and burning in the last 60 years.                                                                                                                                                               |
| Equivalent Clearcut Area (ECA)   | % watershed area | The percentage of total watershed area that is considered functionally and hydrologically comparable to a clearcut forest. Landscapes that have been altered by urban, road, rail, and forestry development as well as crown tenure were considered.                      |
| Mountain pine beetle defoliation | % watershed area | The percentage of pine forests that have been killed by mountain pine beetle in each watershed.                                                                                                                                                                           |

{}\* For watersheds with &gt;50% private ownership, forest disturbance
and associated habitat indicators could not be accurately assessed as
the data were not available, and so we removed those watersheds from the
analysis.
