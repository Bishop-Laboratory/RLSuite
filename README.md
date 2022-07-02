# Analysis for Miller et al. 2022

**Purpose**: The purpose of this repo is to enable replication of the 
figures and analyses from the *Miller et al., 2022* manuscript. The published version of the manuscript is found [here](https://doi.org/10.1093/nar/gkac537). 

## Quick-start

To re-create the figures from the manuscript, do the following:

1. Clone the repo

```shell
git clone https://github.com/Bishop-Laboratory/RLoop-QC-Meta-Analysis-Miller-2022.git
cd RLoop-QC-Meta-Analysis-Miller-2022/
```

2. System dependencies

There are some system dependencies which must be installed 
for compiling some packages. For Windows, this requires [RTools](https://cran.r-project.org/bin/windows/Rtools/).

To aid in the discovery of system requirements, use [`getsysreqs`](https://github.com/mdneuzerling/getsysreqs).

Install `getsysreqs`:

```shell
R -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")'
R -e 'remotes::install_github("mdneuzerling/getsysreqs", force=TRUE)'
```

Capture dependencies (replace with your OS/Version):

```shell
OS="ubuntu"
VER="20.04"
REQS=$(Rscript -e 'options(warn = -1); cat(getsysreqs::get_sysreqs("renv.lock", distribution = "'$OS'", release = "'$VER'"))' | sed s/"WARNING: ignoring environment value of R_HOME"//)
```

Install (Ubuntu):

```shell
sudo apt install $REQS -y
```

2. Python dependencies

Install python (only tested on v3.8). Then restore from requirements. 

```shell
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

3. Restore the `renv`

```shell
R -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")'
R -e 'renv::restore()'
```

4. Run the script

```shell
Rscript figures.R
```

## Troubleshooting

If any issues arise or anything is unclear, please submit an [issue](https://github.com/Bishop-Laboratory/RLoop-QC-Meta-Analysis-Miller-2022/issues).

