# Function to generate data with AR correlation
synthetic_ar = function(
    n_subjects=100,
    n_timepoints=11,  # default to equidistant times on [0,1]
    prop_observed=0.7,  # SRS proportion
    observation_variance=1.,
    random_effect_variance=1.,
    random_effect_ar1_correlation=0.9,
    effect_size=1, # multiplies the group difference function 
    seed=1
){
  require(dplyr)
  require(magrittr)
  set.seed(seed)
  
  corr = random_effect_ar1_correlation^seq(0, n_timepoints-1)
  corrmat = toeplitz(corr)
  thetamat = mvtnorm::rmvnorm(n_subjects, sigma=corrmat) * sqrt(random_effect_variance)
  t0 = seq(0, n_timepoints-1) / (n_timepoints-1)
  timemat = matrix(t0, n_subjects, n_timepoints, byrow=T)
  # intercept function
  f0 = function(t) -t
  # difference function
  f1raw = function(t) 1/(1+exp((0.6-t)*20))
  f1 = function(t) ifelse(abs(f1raw(t)) < 0.1, 0, f1raw(t))  # small values a shrunk to 0
  group = sample(0:1, n_subjects, T)
  groupmat = matrix(group, n_subjects, n_timepoints)
  term0mat = f0(timemat)
  term1mat = f1(timemat * groupmat)
  errormat = matrix(rnorm(n_timepoints*n_subjects), n_subjects, n_timepoints) * sqrt(observation_variance)
  ymat = term0mat + term1mat * effect_size + thetamat + errormat
  smat = matrix(seq(n_subjects), n_subjects, n_timepoints)
  omat = matrix(runif(n_timepoints*n_subjects) < prop_observed, n_subjects, n_timepoints)
  
  # all observations
  data_full = data.frame(
    response=as.vector(ymat),
    time=as.vector(timemat),
    subject_id=as.vector(smat),
    group=as.vector(groupmat),
    observed=as.vector(omat)
  )
  
  # subset to partially observed
  data_long = data_full %>% filter(observed) %>% select(-observed)
  
  # wide format for other methods (will include NAs)
  data_wide = data_long %>% tidyr::pivot_wider(
    id_cols=c("subject_id", "group"),
    names_from="time",
    values_from="response",
    names_prefix="t",
    names_sort=TRUE
  ) %>% arrange(subject_id)
  
  # wide format, with imputation using FPCA
  Y = data_wide %>% select(starts_with("t")) %>% as.matrix
  Yhat = refund::fpca.sc(Y=Y, argvals=t0)$Yhat
  Yhat[!is.na(Y)] = Y[!is.na(Y)]
  
  data_wide_imputed = bind_cols(
    data_wide %>% select(subject_id, group),
    data.frame(Yhat)
  )
  
  # store generating values
  true_values = data.frame(
    time=t0, 
    b0=f0(t0),
    b1=f1(t0) * effect_size
  )
  
  # output
  instance = list(
    data_long=data_long,
    data_wide=data_wide, 
    data_wide_imputed=data_wide_imputed,
    true_values=true_values,
    times=t0,
    colnames=list(
      subject_id="subject_id",
      group="group",
      long_index="time",
      long_response="response",
      wide_times=paste0("t", t0)
    )
  )
  
  return(instance)
}

# simulate data
instance = synthetic_ar(
  n_subjects=100,
  n_timepoints=11,
  prop_observed=0.7,
  observation_variance=1.,
  random_effect_variance=1.,
  random_effect_ar1_correlation=0.9,
  effect_size=1,
  seed=1
)

# how to select the best tuning parameters
selection = "aic"

out = VCMM::lsvcmm(
  data=instance$data_long, 
  response=instance$colnames$long_response,
  subject=instance$colnames$subject_id,
  time=instance$colnames$long_index,
  vc_covariates=instance$colnames$group,
  cv=ifelse(selection=="cv", 5, 0), # number of CV folds (0 means no CV)
  estimate_variance_components=F, # whether to estimate variance parameters (F=proxy)
  random_effect=T, # whether to model random effects (F=Independent)
  sgl=1., # Sparse group Lasso weight (1=Lasso, 0=Group Lasso)
  lambda=NULL, # regularization parameter, NULL=start from max value and decrease
  lambda_factor=0.005, # lambda_min = lambda_factor * lambda_max
  n_lambda=100, # how many lambda to estimate (log sequence)
  adaptive=1., # adaptive power parameter, 0=no adaptation
  kernel="squared_exponential",  # only one implemented yet
  kernel_scale=NULL, # NULL=estimate at multiple values
  kernel_scale_factor=10, # range will be [h/10, h*10], h=n^-0.2
  n_kernel_scale=21 # how many kernel scale to estimate (log sequence)
)

best_idx = switch(selection,
                  "aic"=which.min(out$models_path$aic),
                  "bic"=which.min(out$models_path$bic),
                  "cv"=which.min(out$models_path$cv_score)
)
lambda = out$models_path$lambda[best_idx]
kernel_scale = out$models_path$kernel_scale[best_idx]
re_ratio = out$models_path$re_ratio[best_idx]
b1est = out$vc_path[2,,best_idx]

# bootstrap
boot_obj = VCMM::lsvcmm.boot(
  data=instance$data_long, 
  response=instance$colnames$long_response,
  subject=instance$colnames$subject_id,
  time=instance$colnames$long_index,
  vc_covariates=instance$colnames$group,
  estimate_variance_components=F,
  random_effect=T, 
  sgl=1.,
  lambda=lambda, 
  adaptive=1., 
  kernel="squared_exponential",
  kernel_scale=kernel_scale, 
  n_samples=1000
)

band_df = VCMM::confidence_band(boot_obj, 0.95, "quantile", 2)