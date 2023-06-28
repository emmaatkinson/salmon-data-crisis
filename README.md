README
================
Emma Atkinson
27-Jun-2023

# Investigating the state of monitoring data for wild Pacific salmon populations in BC

This repository contains code and data associated with an analysis to
assess the current state of escapement monitoring for wild Pacific salmon
populations in British Columbia

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

