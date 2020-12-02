#' Calculate probability of clinical disease by age group
#' * following posterior estimates by Davies et al
#' * source: https://www.nature.com/articles/s41591-020-0962-9
getClinicalFraction <- function(age_groups){
  #' smoothly interpolate between points (x0, y0) and (x1, y1) using cosine interpolation.
  #' for x < x0, returns y0; for x > x1, returns y1; for x0 < x < x1, returns the cosine interpolation between y0 and y1
  interpolate_cos = function(x, x0, y0, x1, y1)
  {
    ifelse(x < x0, y0, ifelse(x > x1, y1, y0 + (y1 - y0) * (0.5 - 0.5 * cos(pi * (x - x0) / (x1 - x0)))))
  }
  
  age_groups[, mid := mean(c(age_low, age_high)), by=seq_len(nrow(age_groups))]
  
  age_y = covid_parameters[["clinical_fraction"]][["values"]][["age_y"]]
  age_m = covid_parameters[["clinical_fraction"]][["values"]][["age_m"]]
  age_o = covid_parameters[["clinical_fraction"]][["values"]][["age_o"]]
  
  # definition of "young", "middle", and "old"
  young  = interpolate_cos(age_groups[, mid], age_y, 1, age_m, 0);
  old    = interpolate_cos(age_groups[, mid], age_m, 0, age_o, 1);
  middle = 1 - young - old;
  
  symp_y = covid_parameters[["clinical_fraction"]][["values"]][["symp_y"]]
  symp_m = covid_parameters[["clinical_fraction"]][["values"]][["symp_m"]]
  symp_o = covid_parameters[["clinical_fraction"]][["values"]][["symp_o"]]
  
  return(young * symp_y + middle * symp_m + old * symp_o)
}

#' Set age-specific clinical fraction in population
setClinicalFraction <- function(population_sample, age_groups){
  init_colnames <- colnames(population_sample)
  age_groups[, "y"] <- getClinicalFraction(age_groups)
  population_sample <- merge(population_sample, age_groups[, c("age_group", "y")], by="age_group")
  setorder(population_sample, participant_id)
  
  return(population_sample[, c(init_colnames, "y"), with=FALSE])
}

#' Set susceptibility to infection
setSusceptibility <- function(population_sample, u=covid_parameters$base_u){
  population_sample[, "u"] <- u
  
  return(population_sample)
}

#' Get random values from discrete gamma functions
#' * paramerised using mu and k
rgammaAltDiscrete <- function(n, mu, k, tstep=1){
  distcrete::distcrete("gamma", tstep, shape = k, scale = mu/k)$r(n)
}

#' Get bootstrapped estimate of population
#' * requires global dataset: population_data (source UNWPP)
bootstrapPopulation <- function(n_households){
  household_sample <- sample(unique(population_data[, household_id]), n_households, TRUE)
  population_sample <- rbindlist(
    lapply(
      1:n_households,
      function(x) population_data[household_id == household_sample[x]][, household_id := x]
    )
  )
  
  setorder(population_sample, household_id)
  population_sample[, "participant_id"] <- c(1:nrow(population_sample))
  
  return(population_sample[, c("participant_id", "household_id", "age_group")])
}

#' Get distribution to use
getRandDist <- function(parameter, n=1, tstep=1){
  p <- covid_parameters[[parameter]]
  if("distribution" %in% names(p)){
    
    return(switch(
      p[["distribution"]],
      "gamma" = rgammaAltDiscrete(n, p[["values"]][["mu"]], p[["values"]][["k"]], tstep),
      NULL
    ))
  } else {
    stop(sprintf("No distribution set for parameter %s", parameter))
  }
}

#' Sample random values for distributions
#' * This will be the duration an individual will be in a compartment, if they would get there
#' * dE: incubation period
#' * dP: duration of pre-clinical period
#' * dC: duration of clinical period
#' * dS: duration of subclinical period
setPopDistributions <- function(population_sample, tstep=1){
  population_sample[, "dE"] <- tstep * getRandDist("dE", nrow(population_sample), 1/tstep)
  population_sample[, "dP"] <- tstep * getRandDist("dP", nrow(population_sample), 1/tstep)
  population_sample[, "dC"] <- tstep * getRandDist("dC", nrow(population_sample), 1/tstep)
  population_sample[, "dS"] <- tstep * getRandDist("dS", nrow(population_sample), 1/tstep)
  
  return(population_sample)
}

#' Calculate total number of age-pairs in each household included
#' * used to create within-household contact matrix corresponding with household structures
#' * account for not being able to contact yourself
totalAgePairs <- function(households){
  households <- households[, c("household_id", "age_group")]
  out <- rbindlist(lapply(
    unique(households[, household_id]),
    function(x){
      y = households[household_id == x];
      return(rbindlist(lapply(
        1:nrow(y),
        function(i){y[-i][, age_group2 := y[i, age_group]]}
      )))
    }
  ))
  out <- as.matrix(dcast(out, age_group2 ~ age_group, fun.aggregate = length)[, -"age_group2"])
  
  return(out)
}

#' Construct contact matrix from total number of contacts by age
constructContactMatrix <- function(total_contacts, sample_popsize, population_popsize, adjust_ratio=1){
  rate_raw <- total_contacts
  for(j in 1:ncol(rate_raw)){
    rate_raw[,j] <- total_contacts[,j]/sample_popsize[j]
  }
  rate_adjusted <- rate_raw
  for(j in 1:ncol(rate_adjusted)){
    for(i in 1:nrow(rate_adjusted)){
      rate_adjusted[i,j] <- (rate_raw[i,j] + rate_raw[j,i] * population_popsize[i]/population_popsize[j])/2
    }
  }
  rate_adjusted * adjust_ratio
  prob <- rate_adjusted/population_popsize
  #' cannot contact oneself, so need to subtract N with 1 for i == j
  diag(prob) <- diag(rate_adjusted)/(population_popsize-1)
  
  return(list(
    total = total_contacts,
    rate_raw = rate_raw,
    rate_adjusted = rate_adjusted,
    prob = prob
  ))
}

#' Calculate R0 if this would be a compartmental model
#' * Create Next Generation Matrix
#' * Differs with IBM as it does not constrain contacts within the household, but uses
#'   the same average number of contacts as used in the IBM
calculateR0cmm <- function(matrices, covid_parameters) {
  dIp = covid_parameters$dP$values$mu
  dIs = covid_parameters$dC$values$mu
  dIa = covid_parameters$dS$values$mu
  
  y = getClinicalFraction(age_groups)
  
  cm = Reduce("+", matrices)
  
  ngm = covid_parameters$base_u * t(t(cm) * (
    y * (covid_parameters$fIp * dIp + covid_parameters$fIc * dIs) + 
      (1 - y) * covid_parameters$fIs * dIa)
  )
  
  return(abs(eigen(ngm)$values[1]))
}

#' Calculate new susceptibility to reach target R0
#' * using compartmental model structure
calculateSusceptibility <- function(matrices, target_R0, covid_parameters.=covid_parameters){
  return( covid_parameters$base_u * (target_R0/calculateR0cmm(matrices, covid_parameters)) )
}

#' Create a new population for simulation
#' * Uses global datasets: popdata; age_groups; contacts_home; household_size_prob
generatePopulation <- function(target_size, contacts_home=contact_matrices_home_symmetric, popdata.=popdata, household_size_age.=household_size_age, verbose=FALSE){
  contact_home_prob <- contacts_home/rowSums(contacts_home)
  agegroups <- nrow(popdata)
  adult_agegroups <- age_groups[, age_high] >= 18
  households <- matrix(0, nrow=agegroups, ncol=target_size)
  
  h <- 0
  target_popdist <- popdata[, N]/sum(popdata[, N]) * target_size
  households_rsums <- rowSums(households)
  
  total_agepairs <- matrix(0, nrow=agegroups, ncol=agegroups)
  
  if(verbose) tictoc::tic()
  while(sum(households_rsums) < target_size){
    h <- h+1
    
    #sample agegroup
    agegrp <- sample(1:nrow(popdata), 1, prob = pmax((target_popdist - households_rsums)/sum(target_popdist - households_rsums), 0))
    
    #sample household size
    hhsize <- sample(1:6, 1, prob=household_size_prob[agegrp,])
    
    #sample other household member agegroups
    members <- matrix(0, nrow=agegroups, ncol=1)
    
    #ensure at least 1 adult
    while(sum(members[adult_agegroups, 1]) == 0){
      #probabilities of household member age based on at home contact matrix
      members <- rmultinom(
        1:agegroups,
        hhsize-1,
        contact_home_prob[agegrp,] * pmax((target_popdist - households_rsums)/sum((target_popdist - households_rsums)),0)
      )
      members[agegrp,1] <- members[agegrp,1]+1  
    }
    
    for(i in 1:agegroups){
      total_agepairs[i, ] <- total_agepairs[i, ] + members[i]*(members - sapply(1:agegroups, function(a,i) a==i, i))  
    }
    
    households[,h] <- members
    households_rsums <- households_rsums + members
    if(verbose & h %% 500 == 0){
      print(sprintf("households: %s; popsize: %s/%s",h, sum(households_rsums), target_size))
    }
  }
  if(verbose) tictoc::toc()
  
  households <- households[,1:h]
  
  if(verbose) message("Convert to population structure...")
  population <- rbindlist(lapply(1:agegroups, function(a) data.table(
    household_id=unlist(
      sapply(
        seq(1, ncol(households)),
        function(i, hhid, hhmem) rep(hhid[i], hhmem[i]),
        seq(1, ncol(households)),
        unlist(households)[seq(a, by=agegroups, length.out=ncol(households))]
      )
    ),
    age_group=unlist(
      sapply(
        seq(1, ncol(households)),
        function(i, age, hhmem) rep(age-1, hhmem[i]),
        a,
        unlist(households)[seq(a, by=agegroups, length.out=ncol(households))]
      )
    )
  )))
  
  setorder(population, household_id, age_group)
  population[, participant_id := 1:.N]
  
  return(list(
    population[, c("participant_id", "household_id", "age_group")],
    total_agepairs
  ))
}

#' Function to overwrite basic model-options to apply interventions
getIVtype = function(params, test_coverage, compliance_pos, compliance_hh, model_quick=FALSE, iv_stop=999, lockdown_start=-1, lockdown_effect=1){
  params[["model_quick"]] = model_quick
  params[["iv_stop"]] = iv_stop
  params[["lockdown"]][["start_value"]] = lockdown_start
  params[["lockdown"]][["effect"]] = lockdown_effect
  params[["interventions"]][[1]][["compliance_positive"]] = compliance_pos
  params[["interventions"]][[1]][["compliance_household"]] = compliance_hh
  params[["interventions"]][[1]][["compliance_notest"]] = compliance_hh
  if(length(test_coverage) > 0){
    params[["interventions"]][[1]][["coverage"]] = test_coverage[1]
  }
  if(length(test_coverage) > 1){
    params[["interventions"]] = c(
      params[["interventions"]][1],
      lapply(2:length(test_coverage), function(w, iv, test_coverage){
        iv[["start_type"]] = 1
        iv[["start_value"]] = 7
        iv[["coverage"]] = test_coverage[w]
        return(iv)
      }, params[["interventions"]][[1]], test_coverage)
    )  
  } else if(length(test_coverage) == 0){
    params[["interventions"]] = list()
  }
  return(params)
}

#' Function to calculate Rt in a given period
#' * Returns average number of secondary infections for a cohort of infectious individuals, who 
#'   were infectious in the period provided and who completed their infectious period
#' * Requires large enough period so those with long infectious period are not censored
getRt <- function(data, period){
  #primary cases
  infectious <- data[part_id > 0 & recovered_at > -1 & recovered_at < period[2] & infectious_at >= period[1], part_id]
  if(length(infectious) == 0){ warning("no infectious individuals in period"); return(NA) }
  #secondary cases
  infected <- data[part_id > 0 & infected_by %in% infectious]
  return(nrow(infected)/length(infectious))
}