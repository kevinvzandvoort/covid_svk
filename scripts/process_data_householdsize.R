library(data.table)
library(readxl)

setwd("~/workspace/covid_svk2")

#' PROCESS HOUSEHOLD SIZE DATA
#' 
#' * Description: Household size data by age for Slovakia
#' * Source: https://ec.europa.eu/eurostat/databrowser/view/CENS_01RHSIZE__custom_214489/default/table?lang=en
tempfile = "./data/raw/eurostat_household_composition_sk_2011.xlsx"

#' Custom query will likely be expired
#' Changed query in source link by:
#' - including all values under 'Number of persons'
#' - including all values under 'Age class'
if(!file.exists(download.file(
  url="https://ec.europa.eu/eurostat/databrowser-backend/api/query/1.0/LIVE/xlsx/en/download/1e8ac54b-ce97-477e-a3de-c4f5ea489650",
  destfile=temp, mode="ab"
)))

agegrps <- sapply(4:22, function(x) readxl::read_xlsx(tempfile, range = "C7", sheet = x, col_names=F)[[1]])
data <- rbindlist(
  lapply(4:22, function(x, agegrps) as.data.table(
    readxl::read_xlsx(tempfile, range = "B11:L13", sheet = x, col_names=T)
  )[2,][, agegroup := agegrps[x-3]], agegrps)
)
data <- melt(data, id.vars="agegroup")[!is.na(value)]
data[, household_size := as.numeric(gsub("([0-9]+).*$", "\\1", variable))]
agetable <- data.table(
  agegroup = agegrps,
  age_low = seq(0,90,5),
  age_high = c(seq(4,89,5), 100)
)
household_size_age <- merge(data, agetable, by="agegroup")
household_size_age <- household_size_age[, c("household_size", "age_low", "age_high", "value")]
fwrite(household_size_age, "./data/clean/household_size_age.csv")