.args <- if (interactive()) c(
  "/nfs/general/covid_svk",
  1
) else commandArgs(trailingOnly = TRUE)

setwd(.args[1])
r = as.numeric(.args[2])

library("Rcpp")
library("RcppGSL")
library("data.table")
library("distcrete")
library("qs")
library("tictoc")

source("./scripts/functions.R")
source("./scripts/read_data.R")

#' where should output files be saved?
out_folder <- "./out"

#' compile c++ code for model using Rcpp
#' compile with openmp flag to increase model speed
Sys.setenv(PKG_CPPFLAGS = "-fopenmp")
sourceCpp(
  "./model/covidIBM.cpp",
  rebuild = T, cacheDir = "./model/build", verbose = T
)

#' target population size to use
n_pop <- 77771
#' observed number who attended test in test-rounds
n_test <- c(48320, 44197, 43983)
#' observed number who tested positive in first test-round
#' used to calculate prevalence when test starts
n_pos <- 1569

#' assume one timestep per day
tstep_day <- 1

#' number of days in which full model is ran
days <- 365

#' target number of infectious individuals at t=1
init_infected <- 10

#' target R0
target_R0 <- 1.5

#' number of days to start lockdown, relative to day when mass-testing is implemented
#' is only used in scenarios where lockdown is implemented
lockdown_start <- -7

#' number of model runs required
runs <- 1

prevalence_trigger <- n_pos/n_test[1]
prop_tested <- (n_test/n_pop) / (sum(popdata[agegroup %in% age_groups[age_high > 10 & age_low < 65, name], N])/sum(popdata[, N]))
if(!dir.exists(out_folder)){
  dir.create(out_folder)
}

#' set a seed for random number generator in R
#' note, does not affect random number generator used in C++ model
set.seed(r)

#' create a new population for this run
bootstrapped_pop <- generatePopulation(n_pop, verbose=T)

#' total household age-pairs for household contact matrix
bootstrap_sample_household_total <- bootstrapped_pop[[2]]

#' population of individuals
bootstrapped_pop = bootstrapped_pop[[1]]
bootstrapped_pop <- setPopDistributions(bootstrapped_pop)
bootstrapped_pop <- setClinicalFraction(bootstrapped_pop, age_groups)

bootstrap_population_popsize <- bootstrapped_pop[, .N, by=age_group][order(age_group), N]

#' this is the matrix for within household mixing (not used in IBM model)
contact_matrix_household <- constructContactMatrix(
  bootstrap_sample_household_total,
  bootstrap_population_popsize,
  bootstrap_population_popsize,
  contact_prob_household
)

#' this is the matrix for outside household mixing
contact_matrix_outside_household <- Reduce("+", list(
  #' assume at-work contacts are reduced by 25% compared to pre-Covid levels
  contact_matrices_work * 0.75,
  #' assume primary schools are open (no change in contacts for u12), and secondary
  #'  schools are closed (no at-school contacts for those 12+)
  contact_matrices_school * matrix(
    c(c(1,1,0.5,rep(0,13)),
      c(1,1,0.5,rep(0,13)),
      c(0.5,0.5,0.5,rep(0,13)),
      rep(0, 16*13)
    ), ncol=16, byrow=T
  ),
  #' assume other contacts are reduced by 75%
  contact_matrices_other * 0.25
))

#' make contact matrix symmetric for this population
#' * ensures total number of contacts between age-groups
#'    i and j == those between j and i
contact_matrix_outside_household <- (
  contact_matrix_outside_household * bootstrap_population_popsize +t(
    contact_matrix_outside_household * bootstrap_population_popsize
  )
)/2/bootstrap_population_popsize

#' calibrate to R0
#' * run model multiple times to validate target R0 level is reached when averaging
#'   over start of multiple simulations
cmm_matrices = list(
  "household" = contact_matrix_household$rate_adjusted * contact_prob_household,
  "outside_household" = t(contact_matrix_outside_household)
)
u <- calculateSusceptibility(cmm_matrices, target_R0)
bootstrapped_pop <- setSusceptibility(bootstrapped_pop, u)

#' baseline parameters for IBM model
ibm_params_base <- list(
  #population of individuals that is used in model
  population = bootstrapped_pop,
  #total number of households
  n_households = length(unique(bootstrapped_pop[, household_id])),
  #number of initial infectious individuals
  n_init_infected = init_infected,
  #total number of agegroups
  agegroups = nrow(age_groups), 
  #days to run the model (if not interrupted prematurely)
  days = 365,
  #time-step used per day
  tstep_day = tstep_day,
  #should quick version of model be used?
  #' if TRUE: faster, does not track who-infects-whom (uses IBM-CMM-hybrid)
  #' if FALSE: slower, does track who-infects-whom (full IBM)
  model_quick = FALSE,
  #allow model to be interrupted after onset of interventions
  iv_stop = FALSE,
  #contact probability by timestep
  contact_household = 1-(1-contact_prob_household)^(1/tstep_day),
  #convert contact matrix to probability of contact between any 2 individuals of ages i and j
  contact_matrix = {
    x = t(t(contact_matrix_outside_household) / bootstrap_population_popsize);
    diag(x) = diag(contact_matrix_outside_household)/(bootstrap_population_popsize-1);
    x
  },
  #should a lockdown be implemented?
  lockdown = list(
    start_value = lockdown_start,
    effect = 1,
    test_after = -lockdown_start
  ),
  #should testing be implemented?
  interventions = list(list(
    #when to start the intervention?
    # type 0: start trigger is prevalence of infectious people
    # type 1: start trigger is days after first intervention
    start_type = 0,
    start_value = prevalence_trigger,
    start_time = -1,
    #duration this intervention is active for
    duration = 10,
    #proportion of eligible people who attend testing
    coverage = 1.00,
    #age groups eligible for testing
    age_eligible = which(!(age_groups[, age_high] < 10 | age_groups[, age_low] >= 65))-1,
    #compliance of self-quarantine with those who test positive
    compliance_positive = 1,
    #compliance of self-quarantine for household members who test positive
    compliance_household = 0,
    #compliance of self-quarantine of eligible people who do not attend testing
    compliance_notest = 0
  )),
  #what parameters should be used for testing?
  #* test specificity and sensitivity
  testpars = list(
    test_specificity = covid_parameters$test_specificity,
    test_sensitivity_preinfectious = covid_parameters$test_sensitivity_preinfectious,
    test_sensitivity_preclinical = covid_parameters$test_sensitivity_infectious_preclinical,
    test_sensitivity_subclinical = covid_parameters$test_sensitivity_infectious_subclinical,
    test_sensitivity_clinical = covid_parameters$test_sensitivity_infectious_clinical 
  )
)

#seed for random number generator in C++ model
model_seed <- r

#scenario 0 - no lockdown or testing
#scen0_baseline <- as.data.table(covidIBM(getIVtype(ibm_params_base, test_coverage=c(), compliance_pos=1, compliance_hh=1), model_seed)$population)

#scenario 1 - no lockdown effect; 
#scenario 1 - full compliance of household members
#scen1_test1_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=1, iv_stop=7*4), model_seed)
#scen1_test2_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=1, iv_stop=7*3), model_seed)
#scen1_test3_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=1, iv_stop=7*2), model_seed)

#scenario 1 - no compliance of household members
#scen1_test1_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0, iv_stop=7*4), model_seed)
#scen1_test2_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0, iv_stop=7*3), model_seed)
#scen1_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0, iv_stop=7*2), model_seed)

#scenario 1 - 50% compliance of household members
scen1_test1_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*4), model_seed)
scen1_test2_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*3), model_seed)
scen1_test3_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*2), model_seed)


#scenario 2 - lockdown brings Rt to 1
#time_test1 <- scen1_test1_compl_full$intervention[[1]]$start_t
time_test1 <- scen1_test1_compl_half$intervention[[1]]$start_t
target_Rt <- 1

#calculate Rt in period before implementation of test
rtHH = getRtHH(
  #as.data.table(scen1_test1_compl_full$population), 
  as.data.table(scen1_test1_compl_half$population), 
  #use longest period of possible infectious period to calculate Rt
  time_test1 + c(-ceiling(qgamma(
    0.999, shape=covid_parameters$dS$values$k,
    scale=covid_parameters$dS$values$mu/covid_parameters$dS$values$k)
  ), 0)
)

if(rtHH["hh"] >= target_Rt){
  stop("HH transmission is too high")
}

#scenario 2 - lockdown without testing
#scen2_test0 <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(0), iv_stop=7*4, compliance_pos=0, compliance_hh=0, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)

#scenario 2 - lockdown with tests and full compliance of household members
#scen2_test1_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=1, iv_stop=7*4, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
#scen2_test2_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=1, iv_stop=7*3, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
#scen2_test3_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=1, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)

#scenario 2 - lockdown with tests and no compliance of household members
#scen2_test1_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0, iv_stop=7*4, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
#scen2_test2_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0, iv_stop=7*3, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
#scen2_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)

#scenario 2 - lockdown with tests and no compliance of household members
scen2_test1_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*4, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
scen2_test2_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*3, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)
scen2_test3_compl_half <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0.5, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)


#scenario 3 - lockdown brings Rt to 0.6, assume compliance with testing
target_Rt <- 0.6
if(rtHH["hh"] >= target_Rt){
  stop("HH transmission is too high")
}
#scen3_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=0, compliance_hh=0, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=(target_Rt - rtHH["hh"])/rtHH["nohh"]), model_seed)

#save results - infected individuals
#qs::qsave(scen0_baseline[part_id > 0], sprintf("%s/run_%s_scen0_baseline.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test0$population)[part_id > 0], sprintf("%s/run_%s_scen2_test0.qs", out_folder, r))

#qs::qsave(as.data.table(scen1_test1_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test1_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen1_test2_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test2_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen1_test3_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test3_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen1_test1_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test1_compl_none.qs", out_folder, r))
#qs::qsave(as.data.table(scen1_test2_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test2_compl_none.qs", out_folder, r))
#qs::qsave(as.data.table(scen1_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test3_compl_none.qs", out_folder, r))

qs::qsave(as.data.table(scen1_test1_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen1_test1_compl_half.qs", out_folder, r))
qs::qsave(as.data.table(scen1_test2_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen1_test2_compl_half.qs", out_folder, r))
qs::qsave(as.data.table(scen1_test3_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen1_test3_compl_half.qs", out_folder, r))

#qs::qsave(as.data.table(scen2_test1_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test1_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test2_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test2_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test3_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test3_compl_full.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test1_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test1_compl_none.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test2_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test2_compl_none.qs", out_folder, r))
#qs::qsave(as.data.table(scen2_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test3_compl_none.qs", out_folder, r))

qs::qsave(as.data.table(scen2_test1_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen2_test1_compl_half.qs", out_folder, r))
qs::qsave(as.data.table(scen2_test2_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen2_test2_compl_half.qs", out_folder, r))
qs::qsave(as.data.table(scen2_test3_compl_half$population)[part_id > 0], sprintf("%s/run_%s_scen2_test3_compl_half.qs", out_folder, r))

#qs::qsave(as.data.table(scen3_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen3_test3_compl_none.qs", out_folder, r))

#' Summarized results - test results in scenario
#tests <- sapply(scen3_test3_compl_none$intervention, "[[", "start_t")
tests <- sapply(scen2_test3_compl_half$intervention, "[[", "start_t")
out <- rbindlist(lapply(c(1:2), function(a){
  rbindlist(lapply(c(1:3), function(b, a){
    #rbindlist(lapply(c("full", "none"), function(d, b, a){
    rbindlist(lapply(c("half"), function(d, b, a){
      iv = get(sprintf("scen%s_test%s_compl_%s", a, b, d))$intervention
      z = as.data.table(iv[b])[, test := b][, compliance := d]
      z[, variable := c("t_start", "test_attend", "test_positive", "in_quarantine", "susceptibles", "preinfectious", "infectious")];
      z[, scenario := a]
      z = dcast(z, scenario+compliance+test~variable, value.var="V1")
      return(z)
    }, b, a))
  }, a))
}))

#out_baseline <- rbindlist(lapply(1:3, function(x, tests){
#  rbindlist(lapply(1:2, function(d, dats, x, tests){
#    dat <- get(dats[d])
#    if(!is.data.table(dat)){
#      dat <- as.data.table(dat$population)
#    }
#    return(data.table(
#      "scenario" = ifelse(d==1, 0, d),
#      "compliance" = NA,
#      "test" = 0,
#      "in_quarantine" = NA,
#      "infectious" = dat[part_id > 0 & infectious_at <= tests[x] & recovered_at > tests[x], .N],
#      "preinfectious" = dat[part_id > 0 & infected_at <= tests[x] & infectious_at > tests[x], .N],
#      "susceptibles" = sum(bootstrap_population_popsize) - scen0_baseline[part_id > 0 & infected_at <= tests[x], .N],
#      "t_start" = tests[x],
#      "test_attend" = NA,
#      "test_positive" = NA
#    ))
#  }, c("scen0_baseline", "scen2_test0"), x, tests))
#}, tests))

#iv = scen3_test3_compl_none$intervention
#out_scen3 = rbindlist(lapply(1:3, function(x, iv){
#  z <- as.data.table(iv[[x]])
#  colnames(z) <- c("t_start", "test_attend", "test_positive", "in_quarantine", "susceptibles", "preinfectious", "infectious")
#  z[, test := x][, compliance := "none"][, scenario := 3]
#  return(z)
#}, iv))

#' combine summarized results
#out <- rbindlist(list(
#  out, out_baseline, out_scen3
#), use.names = TRUE)
out[, popsize := sum(bootstrap_population_popsize)]
out[, rt := sum(rtHH)]
setorder(out, scenario, compliance, test)

#' save data
fwrite(out, sprintf("%s/run_%s_half.csv", out_folder, r))
  