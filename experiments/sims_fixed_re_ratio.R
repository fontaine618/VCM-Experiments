library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)


# ==============================================================================
# Setup batchtools registry
DIR = "./experiments/re_ratio/"

# registry = loadRegistry(DIR, writeable=T)

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
source("./batchtools/problems/synthetic_ar.R")
addProblem(
  name="synthetic_ar",
  fun=synthetic_ar,
  data=NULL
)
export[["synthetic_ar"]] = synthetic_ar

# for debugging
instance = synthetic_ar(NULL, NULL)
# ------------------------------------------------------------------------------


# ==============================================================================
# Setup algorithms
source("./batchtools/algorithms/lsvcmm.R")
addAlgorithm(
  name="LSVCMM",
  fun=lsvcmm_wrapper
)
export[["lsvcmm_wrapper"]] = lsvcmm_wrapper
# ------------------------------------------------------------------------------


# ==============================================================================
# Experimental design
n_reps=500
problems = list(
  synthetic_ar=CJ(
    n_subjects=c(60),
    n_timepoints=c(11),
    prop_observed=c(0.8), 
    observation_variance=c(0.2),
    random_effect_variance=c(0., 1., 4.), # with and without random effect
    random_effect_ar1_correlation=c(0.95), # with and without random effect
    effect_size=c(1.),
    seed=seq(n_reps)
  )
)

algorithms = list(
  LSVCMM=data.table(
    selection=c("aic"),
    boot=c(F),
    adaptive=c(1.),
    estimate_variance_components=c(F, F, F, F, F, F, F, T),
    random_effect=c(T),
    kernel_scale=c(0.3),
    lambda=c(-1),
    re_ratio=c(0, 0.01, 0.1, 1, -1, 10, 100, -1),
    name=c("no RE", "v=0.01", "v=0.1", "v=1", "proxy", "v=10", "v=100", "estimated")
  )
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
submitJobs(resources=list(walltime=1000))
getStatus()

# # in case walltime was not enough, run this to resume
not_done = findNotDone()
submitJobs(not_done$job.id, resources=list(walltime=1000))
# job = makeJob(2000)
# job$prob.pars
# instance = job$instance
# df = instance$data_long
# out = VCMM::lsvcmm(
#   data=df,
#   response="response",
#   subject="subject_id",
#   time="time",
#   vc_covariates="group",
#   cv=0,
#   estimate_variance_components=T,
#   random_effect=T,
#   kernel_scale=0.3,
#   lambda=0.02,
#   adaptive=1,
#   re_ratio=-1
# )

# ------------------------------------------------------------------------------




# ==============================================================================
# Get results by timepoint
library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

DIR = "./experiments/re_ratio/"
loadRegistry(DIR)

source("./batchtools/problems/synthetic_ar.R")
instance = synthetic_ar(NULL, NULL)
t0 = instance$true_values$time

decision = function(result) result$vc$excludes_zero
decisions = reduceResultsDataTable(fun = decision) %>% unwrap()
colnames(decisions) = c("job.id", paste0("T_", t0))

estimate = function(result) result$vc$estimate
estimates = reduceResultsDataTable(fun = estimate) %>% unwrap()
colnames(estimates) = c("job.id", paste0("T_", t0))

write.csv(decisions, paste0(DIR, "decisions.csv"), row.names=F)
write.csv(estimates, paste0(DIR, "estimates.csv"), row.names=F)

parameters = getJobPars() %>% unwrap()
write.csv(parameters, paste0(DIR, "parameters.csv"), row.names=F)


# estpars = reduceResultsDataTable(fun = function(result) result$log$re_ratio) %>% unwrap()
# estpars %<>% mutate(job.id=as.character(job.id)) %>% left_join(parameters%>% 
#                          rownames_to_column(var="job.id") %>% 
#                          select(job.id, algo_display, prob_display), by="job.id")
# 
# with(estpars, table(result.1, algo_display, prob_display))
# ------------------------------------------------------------------------------




# ==============================================================================
# Get results by timepoint
library(batchtools)
library(data.table)
library(tidyverse)
library(magrittr)

DIR = "./experiments/re_ratio/"
loadRegistry(DIR)

theme_set(theme_minimal())

decisions = read.csv(paste0(DIR, "decisions.csv")) %>% column_to_rownames(var="job.id")
estimates = read.csv(paste0(DIR, "estimates.csv")) %>% column_to_rownames(var="job.id")
parameters = read.csv(paste0(DIR, "parameters.csv")) %>% column_to_rownames(var="job.id")

parameters %<>% mutate(
  algo_display = paste0(algorithm, " (", name, ")")
)

parameters %<>% mutate(
  prob_display = paste0("RE variance: ", random_effect_variance)
)

b1true = instance$true_values$b1
b1nz = 1*(abs(b1true) > 1e-10)

tnames = paste0("T_", t0)
tnamesz = tnames[b1nz==0]
tnamesnz = tnames[b1nz==1]
probs = parameters$prob_display %>% unique

accmat = decisions %>% apply(1, function(row) 1*(row==b1nz)) %>% t() %>% data.frame()
accmat %<>% rownames_to_column(var="job.id") %>% 
  left_join(parameters %>% 
              rownames_to_column(var="job.id") %>% 
              select(job.id, algo_display, prob_display), 
            by="job.id")
accmat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="acc")
accmat %<>% mutate(time2=as.numeric(stringr::str_split_i(time, "_", 2)))

pred0mat = decisions %>% apply(1, function(row) 1*(row==0)) %>% t() %>% data.frame()
pred0mat %<>% rownames_to_column(var="job.id") %>% 
  left_join(parameters %>% 
              rownames_to_column(var="job.id") %>% 
              select(job.id, algo_display, prob_display), 
            by="job.id")
pred0mat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="pred0")
pred0mat %<>% mutate(time2=as.numeric(stringr::str_split_i(time, "_", 2)))

pred1mat = decisions %>% apply(1, function(row) 1*(row==1)) %>% t() %>% data.frame()
pred1mat %<>% rownames_to_column(var="job.id") %>% 
  left_join(parameters %>% 
              rownames_to_column(var="job.id") %>% 
              select(job.id, algo_display, prob_display), 
            by="job.id")
pred1mat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="pred1")
pred1mat %<>% mutate(time2=as.numeric(stringr::str_split_i(time, "_", 2)))

sqerrmat = estimates %>% apply(1, function(row) (row-b1true)^2) %>% t() %>% data.frame()
sqerrmat %<>% rownames_to_column(var="job.id") %>% 
  left_join(parameters %>% 
              rownames_to_column(var="job.id") %>% 
              select(job.id, algo_display, prob_display), 
            by="job.id")
sqerrmat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="sqerr")
sqerrmat %<>% mutate(time2=as.numeric(stringr::str_split_i(time, "_", 2)))

estmat = estimates %>% apply(1, function(row) row) %>% t() %>% data.frame()
estmat %<>% rownames_to_column(var="job.id") %>% 
  left_join(parameters %>% 
              rownames_to_column(var="job.id") %>% 
              select(job.id, algo_display, prob_display), 
            by="job.id")
estmat %<>% tidyr::pivot_longer(cols=all_of(tnames), names_to="time", values_to="estimate")
estmat %<>% mutate(time2=as.numeric(stringr::str_split_i(time, "_", 2)))


gs = list()
for(i in seq_along(probs)){
  prob = probs[i]
  
  gest = ggplot() + 
    geom_line(
      data=data.frame(y=b1true, x=t0),
      mapping=aes(x=x, y=y),
      linewidth=5, alpha=0.2
    ) + 
    geom_boxplot(
      data=estmat %>% filter(prob_display==prob),
      mapping=aes(x=time2, y=estimate, 
                  color=algo_display, 
                  fill=algo_display,
                  group=paste(algo_display, time2)),
      outlier.alpha=0.05, alpha=0.2
    ) + theme(
      legend.position="none",
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    ) + ylab("Estimated difference") + ggtitle(prob) + ylim(-1, 2)
  
  g0 = ggplot() + 
    geom_bar(
      data=pred0mat %>% filter(prob_display==prob),
      mapping=aes(x=time2, y=pred0, fill=algo_display),
      stat="summary",
      position="dodge",
      fun="mean"
    ) + theme(
      legend.position="none",
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    ) + ylab("Prop. est. zero") + ylim(0, 1)
  
  g1 = ggplot() + 
    geom_bar(
      data=pred1mat %>% filter(prob_display==prob),
      mapping=aes(x=time2, y=pred1, fill=algo_display),
      stat="summary",
      position="dodge",
      fun="mean"
    ) + theme(
      legend.position="none",
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    ) + ylab("Prop. est. non-zero") + ylim(0, 1)
  
  gerr = ggplot() + 
    geom_boxplot(
      data=sqerrmat %>% filter(prob_display==prob),
      mapping=aes(x=time2, y=sqerr, 
                  color=algo_display, 
                  fill=algo_display,
                  group=paste(algo_display, time2)),
      outlier.alpha=0.05, alpha=0.2
    ) + scale_y_sqrt(limits=c(0, 2)) + theme(
      legend.position="none"
    ) + ylab("Squared error") + xlab("Time")
  
  if(i>1){
    gest = gest + theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
    g0 = g0 + theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
    g1 = g1 + theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
    gerr = gerr + theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
  }
  
  gs  = c(gs, list(gest, g0, g1, gerr, NULL))
  
}


glegend = ggplot(
  data=sqerrmat %>% filter(prob_display==probs[1], algo_display!="GAMM", algo_display!="T_test"), 
  aes(fill=algo_display, x=time2, y=sqerr)
) + 
  geom_bar(alpha=0, stat="identity") + 
  guides(fill = guide_legend(override.aes = list(alpha=1)))+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill='NA'),
        panel.grid = element_blank()
  )

gs[[15]] = glegend

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(probs),
  nrow=5, 
  align="v", axis="tblr",
  byrow=F
)

ggsave(
  paste0(DIR, "results.pdf"),
  g,
  width=12, height=12
)





gs = list()
for(i in seq_along(probs)){
  prob = probs[i]
  
  gest = ggplot() + 
    geom_line(
      data=data.frame(y=b1true, x=t0),
      mapping=aes(x=x, y=y),
      linewidth=1, alpha=1
    ) + 
    geom_boxplot(
      data=estmat %>% filter(prob_display==prob),
      mapping=aes(x=time2, y=estimate, 
                  color=algo_display, 
                  fill=algo_display,
                  group=paste(algo_display, time2)),
      outlier.alpha=0.05, alpha=0.2
    ) + theme(
      legend.position="none"
    ) + ylab("Estimated group difference") + ggtitle(prob) + ylim(-1, 2) + xlab("Time")
  
  if(i>1){
    gest = gest + theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank()
    )
  }
  
  gs  = c(gs, list(gest, NULL))
  
}


glegend = ggplot(
  data=sqerrmat %>% filter(prob_display==probs[1]),
  aes(fill=algo_display, x=time2, y=sqerr)
) + 
  geom_bar(alpha=0, stat="identity") + 
  guides(fill = guide_legend(override.aes = list(alpha=1), nrow=2))+ 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key = element_rect(fill='NA'),
        panel.grid = element_blank()
  )

gs[[4]] = glegend

g = cowplot::plot_grid(
  plotlist=gs,
  ncol=length(probs),
  nrow=2, 
  align="v", axis="tblr",
  byrow=F, 
  rel_heights=c(5, 1)
)

ggsave(
  paste0(DIR, "results2.pdf"),
  g,
  width=12, height=4
)


# ------------------------------------------------------------------------------