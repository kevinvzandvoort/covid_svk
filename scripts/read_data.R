#' Age groups to use across in model
#' * Same 16 agegroups as used in contact matrices
#' * 5-year age-bands for u75yo, and one for 75+
age_groups <- data.table(
  age_low = seq(0, 75 ,5),
  age_high = c(seq(4, 74, 5), 100)
)
age_groups[, age_group := seq_len(.N)-1]
age_groups[, name := paste0(age_low,"-",age_high)]
age_groups[age_low == 100, name := "75+"]

#' Prem 2020 synthetic contact matrices for Slovakia
#' source: https://github.com/kieshaprem/synthetic-contact-matrices
#' * rows are age-groups of participants/contactors; columns are age-groups of their contacts
contact_matrices_home <- readRDS("./data/clean/contact_matrices/contact_home.RDS")
contact_matrices_school <- readRDS("./data/clean/contact_matrices/contact_school.RDS")
contact_matrices_work <- readRDS("./data/clean/contact_matrices/contact_work.RDS")
contact_matrices_other <- readRDS("./data/clean/contact_matrices/contact_other.RDS")

dimnames(contact_matrices_home) = list("participants" = age_groups[, name], "contacts" = age_groups[, name])
dimnames(contact_matrices_school) = list("participants" = age_groups[, name], "contacts" = age_groups[, name])
dimnames(contact_matrices_work) = list("participants" = age_groups[, name], "contacts" = age_groups[, name])
dimnames(contact_matrices_other) = list("participants" = age_groups[, name], "contacts" = age_groups[, name])

#' Assume two household members make a single contact per day
contact_prob_household = 1

#' Population data from UNWPP (2019)
#' * Estimates for Slovakia for 2020
#' Source: https://population.un.org/wpp/Download/Standard/Interpolated/
popdata <- fread("./data/clean/unwpp_slovakia_2020.csv")

#' Household size data by age
#' * Needed to construct population size
#' * See ./process_data_householdsize.R
household_size_age <- fread("./data/clean/household_size_age.csv")

#' Allocate model age-groups to population and household size data
for(i in 1:nrow(age_groups)){
  popdata[age >= age_groups[i, age_low] & age <= age_groups[i, age_high], agegroup := factor(age_groups[i, name], age_groups$name, age_groups$name)]
  household_size_age[age_low >= age_groups[i, age_low] & age_high <= age_groups[i, age_high], agegroup := factor(age_groups[i, name], age_groups$name, age_groups$name)]
}
popdata <- popdata[, .(N=sum(popsize)), by=agegroup]
household_size_age <- dcast(household_size_age[, .(N=sum(value)), by=c("agegroup","household_size")], agegroup~household_size, value.var="N")

#' Assume children do not live alone
household_size_age[agegroup %in% age_groups[age_high < 18, name], "1" := 0]
household_size_prob <- rbindlist(lapply(1:nrow(household_size_age), function(x, household_size_age) household_size_age[x, -"agegroup"]/sum(household_size_age[x, -"agegroup"]), household_size_age))

#' Assume contacts in the home are symmetric on population level
contact_matrices_home_total <- contact_matrices_home * popdata[, N]
contact_matrices_home_symmetric <- (contact_matrices_home_total + t(contact_matrices_home_total))/2/popdata[, N]

#' Assumed parameter distributions for Covid-specific parameters
covid_parameters <- list(
  #baseline susceptibility - changed to reach target R0 value
  "base_u" = 0.08,
  #duration of incubation period
  "dE" = list(
    distribution = "gamma",
    values = list("mu" = 2.5, "k" = 4)
  ),
  #duration of pre-clinical period
  "dP" = list(
    distribution = "gamma",
    values = list("mu" = 2.5, "k" = 4)
  ),
  #duration of clinical period
  "dC" = list(
    distribution = "gamma",
    values = list("mu" = 2.5, "k" = 4)
  ),
  #duration of subclinical period
  "dS" = list(
    distribution = "gamma",
    values = list("mu" = 5.0, "k" = 4)
  ),
  #clinical fraction by age; See: getClinicalFraction()
  "clinical_fraction" = list(
    values = list(
      "age_y" = 19, "age_m" = 50, "age_o" = 68,
      "symp_y" = 0.037, "symp_m" = 0.3, "symp_o" = 0.65
    )
  ),
  #relative infectiousness of those that are preclinical
  "fIp" = 1,
  #relative infectiousness of those that are clinical
  "fIc" = 1,
  #relative infectiousness of those that are subclinical
  "fIs" = 0.5,
  #specificity of test (assumed 100%; not currently implemented in model)
  "test_specificity" = 1,
  #sensitivity of test to detect infection
  #preinfectious sensitivity is assumed 0%; not currently implemented in model
  "test_sensitivity_preinfectious" = 0,
  "test_sensitivity_infectious_clinical" = 1,
  "test_sensitivity_infectious_preclinical" = 1,
  "test_sensitivity_infectious_subclinical" = 1
)
