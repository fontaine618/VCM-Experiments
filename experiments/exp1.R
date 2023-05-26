library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

theme_set(theme_minimal())


# ==============================================================================
# Setup batchtools registry
DIR = "./experiments/exp1/"
# removeRegistry()
registry = makeExperimentRegistry(
  file.dir=DIR, 
  seed=1,
  packages=c("dplyr", "magrittr"),
  source=c("./batchtools/utils/metrics.R")
)
registry$cluster.functions = makeClusterFunctionsSocket(ncpus = 10)
export = list()
# ------------------------------------------------------------------------------



# ==============================================================================
# Setup problem
source("./batchtools/problems/synthetic.R")
addProblem(
  name="synthetic",
  fun=synthetic,
  data=NULL
)
export[["synthetic"]] = synthetic

# for debugging
instance = synthetic(NULL, NULL)
# ------------------------------------------------------------------------------


# ==============================================================================
# Setup algorithms
# source("./batchtools/utils/metrics.R")
# export[["evaluate"]] = evaluate
source("./batchtools/algorithms/ttest.R")
addAlgorithm(
  name="T_test",
  fun=ttest_wrapper
)
source("./batchtools/algorithms/spfda.R")
addAlgorithm(
  name="SPFDA",
  fun=spfda_wrapper
)
source("./batchtools/algorithms/lsvcmm.R")
addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_wrapper
)
source("./batchtools/algorithms/gam.R")
addAlgorithm(
  name="GAMM",
  fun=gam_wrapper
)
export[["ttest_wrapper"]] = ttest_wrapper
export[["spfda_wrapper"]] = spfda_wrapper
export[["lsvcmm_wrapper"]] = lsvcmm_wrapper
export[["gam_wrapper"]] = gam_wrapper
# ------------------------------------------------------------------------------


# ==============================================================================
# Experimental design
problems = list(
  synthetic=CJ(
    n_subjects=c(100),
    n_timepoints=c(11),
    prop_observed=c(0.5, 1.), # with and without missing values
    observation_variance=c(1.),
    random_effect_variance=c(1., 0.), # with and without random effect
    effect_size=c(1.),
    seed=seq(3)
  )
)

algorithms = list(
  LSVCMM=CJ(
    selection=c("aic"),
    boot=c(T),
    adaptive=c(1.),
    estimate_variance_components=c(F, T),
    random_effect=c(F, T)
  ),
  T_test=CJ(method=c("hommel")),
  SPFDA=CJ(),
  GAMM=CJ()
)

addExperiments(
  prob.designs=problems,
  algo.designs=algorithms,
  repls=1
)
# ------------------------------------------------------------------------------



# ==============================================================================
# Get current status
summarizeExperiments()
batchExport(export)

# send one job to see if we coded things correctly
testJob(1, external=T)

submitJobs(resources=list(walltime=1000))
getStatus()
# ------------------------------------------------------------------------------




# ==============================================================================
# Get results by timepoint
loadRegistry(DIR)

result = loadResult(1)
b1true = result$vc$true_value
t0 = result$vc$time
b1_excludes_zero = 1*(abs(b1true)>1e-10)

decision = function(result) result$vc$excludes_zero
decisions = reduceResultsDataTable(fun = decision) %>% unwrap()

estimate = function(result) result$vc$estimate
estimates = reduceResultsDataTable(fun = estimate) %>% unwrap()

mse = function(result) (result$vc$estimate - b1true)^2
mses = reduceResultsDataTable(fun = mse) %>% unwrap()

parameters = getJobPars() %>% unwrap()
parameters %<>% mutate(
  algo_display = paste(
    algorithm, 
    ifelse(
      algorithm=="LSVCMM", 
      paste0(" (", selection, ", a=", adaptive, ", ", ifelse(boot, "Boot", "Est"), ")"),
      ""
    ), 
    sep=""
  )
)

# merge together
estimates = parameters %>% select(job.id, algo_display, prop_observed) %>% left_join(estimates, by="job.id")
decisions = parameters %>% select(job.id, algo_display, prop_observed) %>% left_join(decisions, by="job.id")
mses = parameters %>% select(job.id, algo_display, prop_observed) %>% left_join(mses, by="job.id")

# plot proportion included

# [1] "LSVCMM (aic, a=0, Est)"  "LSVCMM (aic, a=1, Est)"  "LSVCMM (aic, a=0, Boot)" "LSVCMM (aic, a=1, Boot)"
# [5] "LSVCMM (bic, a=0, Est)"  "LSVCMM (bic, a=1, Est)"  "LSVCMM (bic, a=0, Boot)" "LSVCMM (bic, a=1, Boot)"
# [9] "T_test"                  "SPFDA"                   "GAMM" 
algos = c(
  "LSVCMM (aic, a=0, Est)",  "LSVCMM (aic, a=1, Est)",  "LSVCMM (aic, a=0, Boot)", "LSVCMM (aic, a=1, Boot)",
  "LSVCMM (bic, a=0, Est)",  "LSVCMM (bic, a=1, Est)",  "LSVCMM (bic, a=0, Boot)", "LSVCMM (bic, a=1, Boot)"
)
algos = c(
  "LSVCMM (aic, a=0, Boot)", "LSVCMM (aic, a=1, Boot)",
  "T_test",                  "SPFDA",                   "GAMM" 
)
prop = 0.5

colnames(decisions) = c("job.id", "algo", "prop_observed", t0)
decisions_long = decisions %>% 
  filter(algo %in% algos) %>% 
  filter(prop_observed==prop) %>%
  select(-job.id, -prop_observed) %>% 
  pivot_longer(-algo, names_to="time", values_to="prop_included") %>% 
  group_by(algo, time) %>% summarize(prop_included=mean(prop_included, na.rm=T))
decisions_long %<>% mutate(time=as.numeric(time))

gprop = ggplot(
)  + 
  xlab("Time") + ylab("Proportion non-zero") + 
  geom_bar(
    data=data.frame(time=t0, prop_included=b1_excludes_zero),
    mapping=aes(x=time, y=prop_included),
    stat="identity",
    position="dodge",
    alpha=0.2,
    width=0.05
  ) + geom_bar(
    data=decisions_long,
    mapping=aes(x=time, y=prop_included, fill=algo),
    stat="identity",
    position="dodge",
    width=0.05
  ) + scale_x_continuous(breaks=0:5/5) + 
  ggtitle(paste0("Proportion observed: ", prop*100, "%"))

# plot error

colnames(mses) =c("job.id", "algo", "prop_observed", t0)
mses_long = mses %>% 
  filter(algo %in% algos) %>% 
  filter(prop_observed==prop) %>%
  select(-job.id, -prop_observed) %>% 
  pivot_longer(-algo, names_to="time", values_to="mse") %>% 
  group_by(algo, time) %>% summarize(mse=mean(mse, na.rm=T))
mses_long %<>% mutate(time=as.numeric(time))

gmse = ggplot(
)  + 
  xlab("Time") + ylab("RMSE") + geom_bar(
    data=mses_long,
    mapping=aes(x=time, y=sqrt(mse), fill=algo),
    stat="identity",
    position="dodge",
    width=0.05
  ) + scale_x_continuous(breaks=0:5/5) + 
  geom_line(
    data=data.frame(time=t0, mse=b1true*max(sqrt(mses_long$mse))/max(b1true)),
    mapping=aes(x=time, y=mse),
    linewidth=3,
    alpha=0.2
  )



g = cowplot::plot_grid(
  gprop, gmse, nrow=2
)

ggsave(paste0(DIR, "comparison_metrics_", prop, ".pdf"), g, width=20, height=10)

# ------------------------------------------------------------------------------















# ==============================================================================
# Process results
loadRegistry(DIR)

result = loadResult(33)

reduce = function(result) {
  c(
    rmse=result$metrics$rmse,
    overall=result$metrics$cmat$overall,
    byClass=result$metrics$cmat$byClass,
    global=result$metrics$global,
    time=result$time[1]
  )
}
results = reduceResultsDataTable(fun = reduce) %>% unwrap()

# get the job parameters (algo parameters and problem parameters)
parameters = getJobPars() %>% unwrap()

# merge together
out = parameters %>% left_join(results, by="job.id")
head(out)

res_table = out %>% group_by(algorithm, effect_size) %>% summarize(
  recall_sd=sd(byClass.Recall, na.rm=T),
  fpr_sd=sd(1-byClass.Specificity, na.rm=T),
  fdr_sd=sd(1-`byClass.Pos Pred Value`, na.rm=T),
  fwer_sd=sd(byClass.Specificity<1, na.rm=T),
  global_rejection_rate_sd=sd(global, na.rm=T),
  rmse.global_sd=sd(rmse.global, na.rm=T),
  rmse.nonnull_sd=sd(rmse.nonnull, na.rm=T),
  rmse.null_sd=sd(rmse.null, na.rm=T),
  time_sd=sd(time.user.self),
  
  recall=mean(byClass.Recall, na.rm=T),
  fpr=mean(1-byClass.Specificity, na.rm=T),
  fdr=mean(1-`byClass.Pos Pred Value`, na.rm=T),
  fwer=mean(byClass.Specificity<1, na.rm=T),
  global_rejection_rate=mean(global, na.rm=T),
  rmse.global=mean(rmse.global, na.rm=T),
  rmse.nonnull=mean(rmse.nonnull, na.rm=T),
  rmse.null=mean(rmse.null, na.rm=T),
  time=mean(time.user.self)
)

write_csv(res_table, paste0(DIR, "metrics.csv"))
# ------------------------------------------------------------------------------



# ==============================================================================
# Plot
theme_set(theme_minimal())

g_recall = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=recall, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=recall - recall_sd/10, 
    ymax=recall + recall_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("Recall") 

g_fdr = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=fdr, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=fdr - fdr_sd/10, 
    ymax=fdr + fdr_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("FDR") 


g_spec = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=fpr, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=fpr - fpr_sd/10, 
    ymax=fpr + fpr_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("1-Specificity") + 
  geom_hline(yintercept=0.05, linetype="dotted", color="black")

g_rmse.global = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=rmse.global, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=rmse.global - rmse.global_sd/10, 
    ymax=rmse.global + rmse.global_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("RMSE (global)") 

g_rmse.null = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=rmse.null, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=rmse.null - rmse.null_sd/10, 
    ymax=rmse.null + rmse.null_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("RMSE (null)") 


g_rmse.nonnull = ggplot(
  data=res_table,
  mapping=aes(
    x=effect_size, 
    y=rmse.nonnull, 
    color=algorithm, 
    group=algorithm,
    fill=algorithm,
    ymin=rmse.nonnull - rmse.nonnull_sd/10, 
    ymax=rmse.nonnull + rmse.nonnull_sd/10
  )
) + geom_line() + geom_ribbon(alpha=0.2, linewidth=0) +
  xlab("Effect Size") + ylab("RMSE (non-null)") 

g = cowplot::plot_grid(
  g_spec, g_recall, g_rmse.global, g_rmse.null, g_rmse.nonnull,
  nrow=5
)

ggsave(paste0(DIR, "metrics.pdf"), g, width=8, height=12)

# ------------------------------------------------------------------------------