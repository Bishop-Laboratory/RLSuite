# RLSuite Analysis

**Purpose**: The purpose of this repo is to enable replication of the 
figures and analyses from the *RLSuite* manuscript. 

## Quick-start

To re-create the figures from the manuscript, do the following:

1. Open R/RStudio and restore the `renv`

```R
R -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")'
R -e 'renv::restore()'
```

2. Install `RLHub` and `RLSeq` from GitHub

```shell
R -e 'remotes::install_github("Bishop-Laboratory/RLHub", upgrade = "never")'
R -e 'remotes::install_github("Bishop-Laboratory/RLSeq", upgrade = "never")'
```

3. Run the script

```shell
Rscript figures.R
```


