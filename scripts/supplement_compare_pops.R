setwd("~/workspace/covid_svk/")

library("Rcpp")
library("data.table")
library("distcrete")
library("qs")
library("tictoc")

source("./scripts/functions.R")
source("./scripts/read_data.R")

out_folder <- "/media/lsh1604011/writable/home/kevin/workspace/slovak_runs/"

Sys.setenv(PKG_CPPFLAGS = "-fopenmp") #compile with openmp to reduce model speed
sourceCpp(
  "~/workspace/covid_svk/model/covidIBMslovakia_quick.cpp",
  rebuild = T,
  cacheDir = "./model/build",
  verbose = T
)

n_pop <- 77771
n_test <- c(48320, 44197, 43983)
n_pos <- 1569

set.seed(101)
tstep_day <- 1
days <- 365
init_infected <- 10
target_R0 <- 1.5

lockdown_start <- -7

prevalence_trigger <- n_pos/n_test[1]
prop_tested <- (n_test/n_pop) / (sum(popdata[agegroup %in% age_groups[age_high > 10 & age_low < 65, name], N])/sum(popdata[, N]))

getRtHH <- function(data, period){
  #primary cases
  infectious <- data[part_id > 0 & recovered_at > -1 & recovered_at < period[2] & infectious_at >= period[1], part_id]
  if(length(infectious) == 0){ warning("no infectious individuals in period"); return(NA) }
  #secondary cases
  
  rt_hh = mean(sapply(infectious, function(x){data[infected_by == x & household_id == data[part_id==x, household_id], .N]}))
  rt_ohh = mean(sapply(infectious, function(x){data[infected_by == x & household_id != data[part_id==x, household_id], .N]}))
  
  return(c("rt_household"=rt_hh, "rt_outside"=rt_ohh))
}

outRt <- list()

out <- list()
outbak <- out
runs <- 100
for(r in 1:runs){
  message(sprintf("Run %s/%s", r, runs))
  tictoc::tic()
  
  bootstrapped_pop <- generatePopulation(n_pop, verbose=T)
  bootstrap_sample_household_total <- bootstrapped_pop[[2]]
  
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
    contact_matrices_work * 0.75,
    contact_matrices_school * matrix(c(c(1,1,0.5,rep(0,13)),c(1,1,0.5,rep(0,13)),c(0.5,0.5,0.5,rep(0,13)),rep(0, 16*13)), ncol=16, byrow=T), #secondary schools closed
    contact_matrices_other * 0.25
  ))
  
  #' make symmetric for this population
  contact_matrix_outside_household = (contact_matrix_outside_household * bootstrap_population_popsize +
                                        t(contact_matrix_outside_household * bootstrap_population_popsize)
  )/2/bootstrap_population_popsize
  
  cmm_matrices = list(
    "household" = contact_matrix_household$rate_adjusted * contact_prob_household,
    "outside_household" = t(contact_matrix_outside_household)
  )
  
  out[[length(out)+1]] <- list(
    pop = bootstrapped_pop,
    cmm_matrices = cmm_matrices
  )
}

#average household size
quantile(sapply(out, function(x) mean(x$pop[, .N, by="household_id"][, N])), c(0.5, 0.025, 0.975))

#population distribution  
x = rbindlist(lapply(out, function(x){
  return(x$pop[, .N, by="age_group"][, N := N/nrow(x$pop)][order(age_group)])
}))[, .(median=median(N), low95=quantile(N, 0.025), high95=quantile(N, 0.975)), by="age_group"][order(age_group)]

eig = sapply(out, function(x){
  return(eigen(x$cmm_matrices$household)$values[1])
})

eig_i <- which(abs(eig - median(eig)) == min(abs(eig - median(eig))))[1]
med_cm <- out[[eig_i]]$cmm_matrices$household
med_cm = as.data.table(med_cm)
colnames(med_cm)=age_groups$name
med_cm[, "contactor"] = age_groups$name
med_cm = melt(med_cm, id.vars="contactor", variable.name="contactee")

med_cm_o <- out[[eig_i]]$cmm_matrices$outside_household
med_cm_o = as.data.table(med_cm_o)
colnames(med_cm_o)=age_groups$name
med_cm_o[, "contactor"] = age_groups$name
med_cm_o = melt(med_cm_o, id.vars="contactor", variable.name="contactee")

syn_cm_home = contact_matrices_home_symmetric#contact_matrices_home
syn_cm_home = as.data.table(syn_cm_home)
colnames(syn_cm_home)=age_groups$name
syn_cm_home[, "contactor"] = age_groups$name
syn_cm_home = melt(syn_cm_home, id.vars="contactor", variable.name="contactee")

syn_cm_o = contact_matrix_outside_household
syn_cm_o = as.data.table(syn_cm_o)
colnames(syn_cm_o)=age_groups$name
syn_cm_o[, "contactor"] = age_groups$name
syn_cm_o = melt(syn_cm_o, id.vars="contactor", variable.name="contactee")


cm_plot = rbindlist(list(
  med_cm[, setting := "home"][, type := "model"],
  med_cm_o[, setting := "outside home"][, type := "model"],
  syn_cm_home[, setting := "home"][, type := "synthetic"],
  syn_cm_o[, setting := "outside home"][, type := "synthetic"]
))

lshcols <- c("#000000", "#0D5257", "#00BF6F", "#00AEC7", "#A7A8AA")

library(patchwork)

plot_home = ggplot(
  cm_plot[setting == "home"],
  aes(x=factor(contactor, age_groups$name), y=factor(contactee, age_groups$name), fill=value)
)+geom_tile(
)+facet_grid(setting~type)+scale_fill_gradientn(
  colours=c(lshcols[2], high=lshcols[3])
)+theme_classic()+labs(
  x="contactor (age group)", y="contactee (age group)", fill="avg contacts\nper day"
)+theme(axis.text.x = element_text(angle=45, hjust=1))


plot_nohome = ggplot(cm_plot[setting == "outside home"],
  aes(x=factor(contactor, age_groups$name), y=factor(contactee, age_groups$name), fill=value)
)+geom_tile(
)+facet_grid(setting~type)+facet_grid(setting~type)+scale_fill_gradientn(
  colours=c(lshcols[2], lshcols[3], "#FFFFFF", "#00AEC7", "#FFB81C", "#FE5000")
)+theme_classic()+labs(
  x="contactor (age group)", y="contactee (age group)", fill="avg contacts\nper day"
)+theme(axis.text.x = element_text(angle=45, hjust=1))

plot_home/plot_nohome + plot_annotation(tag_levels = 'A')
ggsave("./suppl_compare_matrices.png")

plot_pop = ggplot(
  x[, group := "model"],
  aes(
    x=factor(
      as.character(age_group),
      as.character(0:15),
      age_groups$name
    ),
    y=median, ymin=low95, ymax=high95, colour=group, fill=group
  )
)+geom_ribbon(
  aes(group=group)
)+geom_line(
  colour=lshcols[2],
  aes(group=group)
)+geom_point(
  data=data.table(
    age_group=0:15,
    median=popdata$N/sum(popdata$N),
    low95=rep(0,16),
    high95=rep(0,16),
    group=rep("UNWPP", 16)
  ),
  shape=15,
  size=3
)+coord_cartesian(
  ylim=c(0.04,0.085)
)+theme_classic()+labs(
  x="age group",
  y="relative population size (%)",
  colour = "population",
  fill = "population"
)+scale_y_continuous(labels=scales::percent)+scale_colour_manual(
  values=c("UNWPP" = lshcols[3], "model" = lshcols[2])
)+scale_fill_manual(
  values=c("UNWPP" = lshcols[3], "model" = lshcols[2])
)

plot_pop/plot_home/plot_nohome + plot_annotation(tag_levels = 'A')
ggsave("./compare_model_pop.png", units="mm", width=210, height=297*0.75)
  #' calibrate to R0
  u <- calculateSusceptibility(cmm_matrices, target_R0)
  bootstrapped_pop <- setSusceptibility(bootstrapped_pop, u)
  
  ibm_params_base <- list(
    population = bootstrapped_pop,
    n_households = length(unique(bootstrapped_pop[, household_id])),
    n_init_infected = init_infected,
    agegroups = nrow(age_groups), 
    days = 365,
    tstep_day = tstep_day,
    model_quick = TRUE,
    iv_stop = FALSE,
    contact_household = 1-(1-contact_prob_household)^(1/tstep_day),
    contact_matrix = {
      x = t(t(contact_matrix_outside_household) / bootstrap_population_popsize);
      diag(x) = diag(contact_matrix_outside_household)/(bootstrap_population_popsize-1);
      x
    },
    lockdown = list(
      start_value = lockdown_start,
      effect = 1,
      test_after = -lockdown_start
    ),
    interventions = list(list(
      start_type = 0, #"prevalence",
      start_value = prevalence_trigger,
      start_time = -1,
      duration = 10, #days
      coverage = 0.86,
      age_eligible = which(!(age_groups[, age_high] < 10 | age_groups[, age_low] >= 65))-1, #age groups tested
      compliance_positive = 1, #compliance of those testing positive
      compliance_household = 0, #compliance of household members of those who test positive
      compliance_notest = 0 #compliance of those eligible but not tested
    )),
    testpars = list(
      test_specificity = covid_parameters$test_specificity,
      test_sensitivity_preinfectious = covid_parameters$test_sensitivity_preinfectious,
      test_sensitivity_preclinical = covid_parameters$test_sensitivity_infectious_preclinical,
      test_sensitivity_subclinical = covid_parameters$test_sensitivity_infectious_subclinical,
      test_sensitivity_clinical = covid_parameters$test_sensitivity_infectious_clinical 
    )
  )
  
  model_seed <- r
  
  #scenario 0 - no lockdown or testing
  scen0_baseline <- as.data.table(covidIBM(getIVtype(ibm_params_base, test_coverage=c(), compliance_pos=1, compliance_hh=1), model_seed)$population)
  
  #scenario 1 - no lockdown effect; 
  scen1_test1_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=1, iv_stop=7*4), model_seed)
  scen1_test2_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=1, iv_stop=7*3), model_seed)
  scen1_test3_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=1, iv_stop=7*2), model_seed)
  
  scen1_test1_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0, iv_stop=7*4), model_seed)
  scen1_test2_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0, iv_stop=7*3), model_seed)
  scen1_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0, iv_stop=7*2), model_seed)
  
  #scenario 2 - lockdown brings Rt to 1
  time_test1 <- scen1_test1_compl_full$intervention[[1]]$start_t
  target_Rt <- 1
  
  rt = getRt(
    as.data.table(scen1_test1_compl_full$population), 
    time_test1 + c(-ceiling(qgamma(
      0.999, shape=covid_parameters$dS$values$k,
      scale=covid_parameters$dS$values$mu/covid_parameters$dS$values$k)
    ), 0)
  )
  
  scen2_test0 <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(0), iv_stop=7*4, compliance_pos=0, compliance_hh=0, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  
  scen2_test1_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=1, iv_stop=7*4, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  scen2_test2_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=1, iv_stop=7*3, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  scen2_test3_compl_full <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=1, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  
  scen2_test1_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1)]), compliance_pos=1, compliance_hh=0, iv_stop=7*4, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  scen2_test2_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2)]), compliance_pos=1, compliance_hh=0, iv_stop=7*3, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  scen2_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=1, compliance_hh=0, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  
  #scenario 3 - lockdown brings Rt to 0.6, no compliance with testing
  target_Rt <- 0.6
  scen3_test3_compl_none <- covidIBM(getIVtype(ibm_params_base, test_coverage=c(prop_tested[c(1,2,3)]), compliance_pos=0, compliance_hh=0, iv_stop=7*2, lockdown_start = lockdown_start + time_test1, lockdown_effect=target_Rt/rt), model_seed)
  
  qs::qsave(scen0_baseline[part_id > 0], sprintf("%s/run_%s_scen0_baseline.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test0$population)[part_id > 0], sprintf("%s/run_%s_scen2_test0.qs", out_folder, r))
  
  qs::qsave(as.data.table(scen1_test1_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test1_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen1_test2_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test2_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen1_test3_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen1_test3_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen1_test1_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test1_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen1_test2_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test2_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen1_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen1_test3_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test1_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test1_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test2_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test2_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test3_compl_full$population)[part_id > 0], sprintf("%s/run_%s_scen2_test3_compl_full.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test1_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test1_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test2_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test2_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen2_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen2_test3_compl_none.qs", out_folder, r))
  qs::qsave(as.data.table(scen3_test3_compl_none$population)[part_id > 0], sprintf("%s/run_%s_scen3_test3_compl_none.qs", out_folder, r))
  
  tests <- sapply(scen3_test3_compl_none$intervention, "[[", "start_t")
  
  out <- rbindlist(lapply(c(1:2), function(a){
    rbindlist(lapply(c(1:3), function(b, a){
      rbindlist(lapply(c("full", "none"), function(d, b, a){
        iv = get(sprintf("scen%s_test%s_compl_%s", a, b, d))$intervention
        z = as.data.table(iv[b])[, test := b][, compliance := d]
        z[, variable := c("t_start", "test_attend", "test_positive", "in_quarantine", "susceptibles", "preinfectious", "infectious")];
        z[, scenario := a]
        z = dcast(z, scenario+compliance+test~variable, value.var="V1")
        return(z)
      }, b, a))
    }, a))
  }))
  
  out_baseline <- rbindlist(lapply(1:3, function(x, tests){
    rbindlist(lapply(1:2, function(d, dats, x, tests){
      dat <- get(dats[d])
      if(!is.data.table(dat)){
        dat <- as.data.table(dat$population)
      }
      return(data.table(
        "scenario" = ifelse(d==1, 0, d),
        "compliance" = NA,
        "test" = 0,
        "in_quarantine" = NA,
        "infectious" = dat[part_id > 0 & infectious_at <= tests[x] & recovered_at > tests[x], .N],
        "preinfectious" = dat[part_id > 0 & infected_at <= tests[x] & infectious_at > tests[x], .N],
        "susceptibles" = sum(bootstrap_population_popsize) - scen0_baseline[part_id > 0 & infected_at <= tests[x], .N],
        "t_start" = tests[x],
        "test_attend" = NA,
        "test_positive" = NA
      ))
    }, c("scen0_baseline", "scen2_test0"), x, tests))
  }, tests))
  
  iv = scen3_test3_compl_none$intervention
  out_scen3 = rbindlist(lapply(1:3, function(x, iv){
    z <- as.data.table(iv[[x]])
    colnames(z) <- c("t_start", "test_attend", "test_positive", "in_quarantine", "susceptibles", "preinfectious", "infectious")
    z[, test := x][, compliance := "none"][, scenario := 3]
    return(z)
  }, iv))
  
  out <- rbindlist(list(
    out, out_baseline, out_scen3
  ), use.names = TRUE)
  out[, popsize := sum(bootstrap_population_popsize)]
  out[, rt := rt]
  setorder(out, scenario, compliance, test)
  
  fwrite(out, sprintf("%s/run_%s.csv", out_folder, r))
  
  tictoc::toc()
}
