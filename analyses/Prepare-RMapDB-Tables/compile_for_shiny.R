# Compile the data to an rda object for use with the shiny app
DATA_DIR <- "analyses/Prepare-RMapDB-Tables/"
OUT_DATA_DIR <- "~/projects/RMapDB-shiny/data/"

fls <- list.files(path = DATA_DIR, pattern = ".csv.xz", full.names = TRUE)
names(fls) <- gsub(fls, pattern = ".+/([a-zA-Z0-9_]+)\\.csv\\.xz", replacement = "\\1")

dataLst <- fls %>%
  sapply(readr::read_csv)

save(dataLst, file = file.path(OUT_DATA_DIR, "dataLst.rda"), compress = "xz")

