# Exploratory Data Analysis

This process provides a characteristics summary of the input file [rmap_full_11_25.csv](https://github.com/Bishop-Laboratory/RSeq-supplemental-analysis/blob/main/misc/rmap_full_11_25.csv) including feature statistics and mutual relationships to inform us of any possible confounding factors. This analysis is codified in the R markdown file `eda.Rmd`

**Action**: A [github action](https://docs.github.com/en/actions) will drive this process defined in the workflow file 'eda.yml' located in the `.github/workflows/` directory.

**Trigger**: This process should be executed every time the input file [rmap_full_11_25.csv](https://github.com/Bishop-Laboratory/RSeq-supplemental-analysis/blob/main/misc/rmap_full_11_25.csv) is updated.

The workflow can be [run manually on github]() or [run locally with`act`](https://github.com/nektos/act).

**Returns**: A report in html format named `eda.html`.
