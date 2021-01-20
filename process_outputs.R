library("qs")
library("data.table")
library("ggplot2")
library("patchwork")

setwd("~/workspace/covid_svk")

#' TODO: still need to clean this script

lshcols <- c("#000000", "#0D5257", "#00BF6F", "#00AEC7", "#A7A8AA", "#32006E", "#1E22AA", "#FE5000", "#FFB81C")
out_folder <- "./out/"
runs <- 400
#how many weeks to plot before the first test and after the last test?
plotdata_weeks <- c(3,2)

#' summarize output data
#' returns incidence and prevalence over time
processOutput <- function(data, plot_times, test_times, r, scen, test, compliance){
  prevalence <- sapply(plot_times, function(x) data[status > 1 & infectious_at <= x & (recovered_at > x | recovered_at == -1), .N])
  incidence <- sapply(plot_times, function(x) data[infectious_at == x, .N])
  out <- rbindlist(list(
    data.table(
      time=plot_times, relative_time=plot_times - test_times[1],
      value=prevalence, type="prevalence"
    ),
    data.table(
      time=plot_times, relative_time=plot_times - test_times[1],
      value=incidence, type="incidence"
    )
  ))[, run := r][, scen := scen][, test := test][, compliance := compliance]
  return(out)
}

#' get observed summary model results
#'  observed prevalence during test
#'  return median and 95% interval
summary_data_obs <- rbindlist(lapply(1:runs, function(x){
  if(!file.exists(sprintf("%s/run_%s.csv", out_folder, x))) return(NULL)
  smrydata <- fread(sprintf("%s/run_%s.csv", out_folder, x))
  smrydata[, prevalence := test_positive/test_attend][, actual_prevalence := infectious/popsize]
  return(smrydata[, run := x])
}))[!is.na(prevalence), .(median=median(prevalence), low95=quantile(prevalence, 0.025), high95=quantile(prevalence, 0.975)), by=c("scenario","test","compliance")]

#' add actual model results
#'  actual prevalence at time test would be implemented
#'  for scenarios without testing (baseline, lockdown)
#'  returns prevalence ratio between test 2 and 1 (t2t1) and test 3 and 1 (t3t1)
summary_data <- rbindlist(lapply(1:runs, function(x){
  if(!file.exists(sprintf("%s/run_%s.csv", out_folder, x))) return(NULL)
  smrydata <- fread(sprintf("%s/run_%s.csv", out_folder, x))
  smrydata[, prevalence := test_positive/test_attend][, actual_prevalence := infectious/popsize][is.na(prevalence), prevalence := actual_prevalence]
  rundata = dcast(smrydata, scenario+compliance~t_start, value.var = "prevalence")
  colnames(rundata) <- c("scenario", "compliance", "test1", "test2", "test3")
  return(melt(rundata[, t2t1 := test2/test1][, t3t1 := test3/test1][, -c("test1","test2","test3")],
       measure.vars=c("t2t1","t3t1"))[, run := x])
}))[, .(median=median(value), low95=quantile(value, 0.025), high95=quantile(value, 0.975)), by=c("scenario","compliance","variable")]

#' add observed reduction between tests
summary_data <- rbindlist(
  list(summary_data, data.table(
    scenario=c("o","o"),
    compliance=c(NA,NA),
    variable=c("t2t1","t3t1"),
    median=c(1-0.56, 1-0.82),
    low95=c(1-0.58, 1-0.83),
    high95=c(1-0.54, 1-0.81)
  ))
)

#' what scenarios are available
scenarios <- list(
  list(name = "%s/run_%s_scen0_baseline.qs",
       scen = 0, test = NA, compliance = NA),
  list(name = "%s/run_%s_scen1_test1_compl_full.qs",
       scen = 1, test = 1, compliance = "full"),
  list(name = "%s/run_%s_scen1_test1_compl_full.qs",
       scen = 1, test = 1, compliance = "full"),
  list(name = "%s/run_%s_scen1_test2_compl_full.qs",
       scen = 1, test = 2, compliance = "full"),
  list(name = "%s/run_%s_scen1_test3_compl_full.qs",
       scen = 1, test = 3, compliance = "full"),
  list(name = "%s/run_%s_scen1_test1_compl_none.qs",
       scen = 1, test = 1, compliance = "none"),
  list(name = "%s/run_%s_scen1_test2_compl_none.qs",
       scen = 1, test = 2, compliance = "none"),
  list(name = "%s/run_%s_scen1_test3_compl_none.qs",
       scen = 1, test = 3, compliance = "none"),
  list(name = "%s/run_%s_scen2_test0.qs",
       scen = 2, test = 0, compliance = NA),
  list(name = "%s/run_%s_scen2_test1_compl_full.qs",
       scen = 2, test = 1, compliance = "full"),
  list(name = "%s/run_%s_scen2_test2_compl_full.qs",
       scen = 2, test = 2, compliance = "full"),
  list(name = "%s/run_%s_scen2_test3_compl_full.qs",
       scen = 2, test = 3, compliance = "full"),
  list(name = "%s/run_%s_scen2_test1_compl_none.qs",
       scen = 2, test = 1, compliance = "none"),
  list(name = "%s/run_%s_scen2_test2_compl_none.qs",
       scen = 2, test = 2, compliance = "none"),
  list(name = "%s/run_%s_scen2_test3_compl_none.qs",
       scen = 2, test = 3, compliance = "none"),
  list(name = "%s/run_%s_scen3_test3_compl_none.qs",
       scen = 3, test = 3, compliance = "none")
)

#read output data files
rundata <- list()
for(r in 1:runs){
  message(sprintf("run: %s/%s", r, runs))
  if(!file.exists(sprintf("%s/run_%s.csv", out_folder, r))) next;
  smrydata <- fread(sprintf("%s/run_%s.csv", out_folder, r))
  test_times <- smrydata[scenario == 0, t_start]
  plot_times <- (test_times[1]-plotdata_weeks[1]*7):(test_times[3]+plotdata_weeks[2]*7)
  
  plotdata <- rbindlist(lapply(scenarios, function(x){
    data <- qs::qread(sprintf(x[["name"]], out_folder, r))
    out <- processOutput(data, plot_times, test_times, r, x[["scen"]], x[["test"]], x[["compliance"]])  
  }))
  plotdata[, "popsize"] <- unique(smrydata$popsize)
  rundata[[length(rundata)+1]] <- plotdata
}

plotdata <- rbindlist(rundata)
setorder(plotdata, test)

#colours to use in plot
testcols <- c(
  "Baseline - No test" = lshcols[1],
  "Lockdown (Re: 1) - No test" = lshcols[5],
  "Lockdown (Re: 0.6) - No test" = lshcols[6],
  "Only pilot round" = lshcols[2],
  "Pilot and 1st round" = lshcols[3],
  "Pilot, 1st, and 2nd round" = lshcols[4]
)

#observed prevalence data
observed_prevalence <- summary_data_obs[compliance == "full"]
observed_prevalence[test==1, relative_time := 0]
observed_prevalence[test==2, relative_time := 7]
observed_prevalence[test==3, relative_time := 14]
observed_prevalence[, scen := scenario]

#summarize data across runs
plotdata2 <- plotdata[compliance == "full" | is.na(compliance) | scen==3, .(
  median=median(value/popsize*1000),
  low95=quantile(value/popsize*1000, 0.025),
  high95=quantile(value/popsize*1000, 0.975),
  low50=quantile(value/popsize*1000, 0.25),
  high50=quantile(value/popsize*1000, 0.75)
), by=c("relative_time", "type", "scen", "test", "compliance")]

plotdata2[scen==0, collab := "Baseline - No test"]
plotdata2[scen==2 & test==0, collab := "Lockdown (Re: 1) - No test"]
plotdata2[scen==3, collab := "Lockdown (Re: 0.6) - No test"]
plotdata2[scen!=3 & test==1, collab := "Only pilot round"]
plotdata2[scen!=3 & test==2, collab := "Pilot and 1st round"]
plotdata2[scen!=3 & test==3, collab := "Pilot, 1st, and 2nd round"]
plotdata2[, "collab"] <- factor(plotdata2[, collab], levels=names(testcols), labels=names(testcols))

#plot prevalence over time
#TODO ggplot is not plotting in correct order - currently specifying lines layer by layer
plot_prevtime <- ggplot(
  data=NULL,
  aes(x = relative_time, y = median, group=test)
)+facet_grid(
  factor(
    as.character(scen),
    as.character(1:2),
    c("No lockdown", "Lockdown")
  )~paste0("Compliance household: ", compliance),
  scales="free_y"
)+geom_vline(
  data=data.table(
    relative_time=c(0,7,14)
  ), aes(xintercept=relative_time)
)+geom_vline(
  data=data.table(
    relative_time=c(-7),
    scen=2
  ), aes(xintercept=relative_time), linetype=2
)+geom_line(
  data=plotdata2[type=="prevalence" & scen==3][, -c("compliance")][, test := 3][, scen := 2],
  aes(colour=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen==3][, -c("compliance")][, test := 3][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 3],
  aes(colour=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 3],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 1],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen > 0 & scen < 3 & test == 1],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen==2 & test == 0][, -c("compliance", "scen")][, test := 0][, scen := 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen==2 & test == 0][, -c("compliance", "scen")][, test := 0][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen==0][, -c("compliance", "scen")][, test := 0][, scen := 1],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen==0][, -c("compliance", "scen")][, test := 0][, scen := 1],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="prevalence" & scen==0 & relative_time <= -7][, -c("compliance", "scen")][, test := 0][, scen := 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="prevalence" & scen==0 & relative_time <= -7][, -c("compliance", "scen")][, test := 0][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_errorbar(
  data = observed_prevalence[scen==1],
  aes(x=relative_time, ymin=low95*1000, ymax=high95*1000),
  colour=lshcols[8], width=1
)+geom_point(
  data = observed_prevalence[scen==1],
  aes(x=relative_time, y=median*1000),
  fill=lshcols[8],
  colour="#000000",
  shape=23,
  size=3
)+geom_errorbar(
  data = observed_prevalence[scen==2],
  aes(x=relative_time, ymin=low95*1000, ymax=high95*1000),
  colour=lshcols[9], width=1
)+geom_point(
  data = observed_prevalence[scen==2],
  aes(x=relative_time, y=median*1000),
  fill=lshcols[9],
  colour="#000000",
  shape=23,
  size=3
)+theme_classic(
)+scale_colour_manual(
  values=testcols
)+scale_fill_manual(
  values=testcols
)+labs(
  colour = "Mass tests",
  fill = "Mass tests",
  x = "time relative to pilot testing (days)",
  y = "prevalence (per 1000)"
)+theme(legend.position = "bottom")

plot_inctime <- ggplot(
  data=NULL,
  aes(x = relative_time, y = median, group=test)
)+facet_grid(
  factor(
    as.character(scen),
    as.character(1:2),
    c("No lockdown", "Lockdown")
  )~paste0("Compliance household: ", compliance),
  scales="free_y"
)+geom_vline(
  data=data.table(
    relative_time=c(0,7,14)
  ), aes(xintercept=relative_time)
)+geom_vline(
  data=data.table(
    relative_time=c(-7),
    scen=2
  ), aes(xintercept=relative_time), linetype=2
)+geom_line(
  data=plotdata2[type=="incidence" & scen==3][, -c("compliance")][, test := 3][, scen := 2],
  aes(colour=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen==3][, -c("compliance")][, test := 3][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 3],
  aes(colour=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 3],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 1],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen > 0 & scen < 3 & test == 1],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen==2 & test == 0][, -c("compliance", "scen")][, test := 0][, scen := 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen==2 & test == 0][, -c("compliance", "scen")][, test := 0][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen==0][, -c("compliance", "scen")][, test := 0][, scen := 1],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen==0][, -c("compliance", "scen")][, test := 0][, scen := 1],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+geom_line(
  data=plotdata2[type=="incidence" & scen==0 & relative_time <= -7][, -c("compliance", "scen")][, test := 0][, scen := 2],
  aes(colour=collab,fill=collab)
)+geom_ribbon(
  data=plotdata2[type=="incidence" & scen==0 & relative_time <= -7][, -c("compliance", "scen")][, test := 0][, scen := 2],
  colour=NA, alpha=0.5,
  aes(ymin=low95, ymax=high95, fill=collab)
)+theme_classic(
)+scale_colour_manual(
  values=testcols
)+scale_fill_manual(
  values=testcols
)+labs(
  colour = "Mass tests",
  fill = "Mass tests",
  x = "time relative to pilot testing (days)",
  y = "incidence (per 1000)"
)+theme(legend.position = "bottom")

#colours for prevalence plot
prevcols = c(
  "Observed" = lshcols[7],
  "Baseline - No test" = lshcols[1],
  "Lockdown (Re: 0.6) - No test" = lshcols[6],
  "No lockdown - Mass testing" = lshcols[8],
  "Lockdown (Re: 1)\n+ Mass-testing" = lshcols[9],
  "Lockdown (Re: 1) - No test" = lshcols[5]
)

#plot effectiveness
summary_data2 <- summary_data[compliance != "none" | scenario %in% c(0,3,"o")]
summary_data2[scenario == 2 & compliance == "", scenario := "2o"]

scenario_labels <- rev(c(
  "o" = "Observed",
  "0" = "Baseline - No test",
  "3" = "Lockdown (Re: 0.6) - No test",
  "2o" = "Lockdown (Re: 1) - No test",
  "1" = "No lockdown - Mass testing",
  "2" = "Lockdown (Re: 1)\n+ Mass-testing"
))

plot_effectiveness <- ggplot(
  summary_data2,
  aes(
    x=factor(
      scenario,
      names(scenario_labels),
      scenario_labels
    ),
    colour=factor(
      scenario,
      names(scenario_labels),
      scenario_labels
    ),
    y=median
  )
)+geom_segment(
  aes(
    x=factor(
      scenario,
      names(scenario_labels),
      scenario_labels
    ),
    xend=factor(
      scenario,
      names(scenario_labels),
      scenario_labels
    ),
    y=low95, yend=high95
  ),
  size=1
)+geom_point()+facet_grid(
  .~factor(
    variable,
    c("t2t1","t3t1"),
    c("Test round 1 vs pilot round", "Test round 2 vs test round 1")
  )
)+theme_classic(
)+geom_hline(
  yintercept=1,
  size=0.25,
  colour="#777777"
)+geom_hline(
  data=summary_data2[scenario=="o"],
  linetype=2,
  aes(yintercept=median, colour=factor(
    scenario,
    names(scenario_labels),
    scenario_labels
  ))
)+labs(
  x="",
  y="prevalence ratio"
)+scale_y_continuous(
  trans="log10", breaks=c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2)
)+theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")+coord_flip(
)+scale_colour_manual(values=prevcols)

(plot_effectiveness)/(plot_prevtime)+ plot_annotation(tag_levels = 'A')+plot_layout(
  heights = c(2,3)
)
ggsave("./results20201227.png", units="mm", width=200, height=200)

plot_inctime
ggsave("./results_incidence20201227.png", units="mm", width=200, height=200)
