# RLSuite Analysis

**Purpose**: The purpose of this repo is to enable replication of the 
figures and analyses from the *RLSuite* manuscript. 

## Quick-start

To re-create the figures from the manuscript, do the following:

1. Clone the repo

```shell
git clone https://github.com/Bishop-Laboratory/RLSuite.git
cd RLSuite/
```

2. Restore the `renv`

```shell
R -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")'
R -e 'renv::restore()'
```

3. Install `RLHub` and `RLSeq` from GitHub

```shell
R -e 'remotes::install_github("Bishop-Laboratory/RLHub", upgrade = "never")'
R -e 'remotes::install_github("Bishop-Laboratory/RLSeq", upgrade = "never")'
```

4. Run the script

```shell
Rscript figures.R
```

## Other

If any issues arise or anything is unclear, please submit an [issue](https://github.com/Bishop-Laboratory/RLSuite/issues).

