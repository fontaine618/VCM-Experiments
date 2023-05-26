synthetic = function(
    data, 
    job, 
    n_subjects=100,
    n_timepoints=11,
    prop_observed=0.7,
    observation_variance=1.,
    random_effect_variance=1.,
    effect_size=1.,
    seed=1
){
  library(dplyr)
  library(magrittr)
  set.seed(seed)
  complete_design = data.frame(
    subject_id=rep(seq(n_subjects), n_timepoints),
    time=ceiling(seq(n_subjects * n_timepoints) / n_subjects)
  )
  data = complete_design %>% sample_n(ceiling(prop_observed * n_subjects * n_timepoints))
  data %<>% arrange(subject_id, time)
  group_design = data.frame(
    subject_id=seq(n_subjects),
    group=sample(0:1, n_subjects, T),
    random_effect=rnorm(n_subjects, 0, sqrt(random_effect_variance))
  )
  data %<>% left_join(group_design, by="subject_id")
  data %<>% mutate(time=(time-1)/(n_timepoints-1))
  f0 = function(t) -t
  f1raw = function(t) 1/(1+exp((0.6-t)*20))
  f1 = function(t) ifelse(abs(f1raw(t)) < 0.1, 0, f1raw(t))
  data %<>% mutate(
    b0=f0(time),
    b1=f1(time)
  )
  data %<>% mutate(mean=b0 + b1*group*effect_size + random_effect)
  data %<>% mutate(response=mean+rnorm(nrow(.), 0, sqrt(observation_variance)))
  times = (seq(n_timepoints)-1)/(n_timepoints-1)
  true_values = data.frame(
    time=times, 
    b0=f0(times),
    b1=f1(times) * effect_size
  )
  
  data_wide = data %>% tidyr::pivot_wider(
    id_cols=c("subject_id", "group"),
    names_from="time",
    values_from="response",
    names_prefix="t",
    names_sort=TRUE
  )
  
  Y = data_wide %>% select(starts_with("t")) %>% as.matrix
  Yhat = refund::fpca.sc(Y=Y, argvals=times)$Yhat
  Yhat[!is.na(Y)] = Y[!is.na(Y)]
  
  data_wide_imputed = bind_cols(
    data_wide %>% select(subject_id, group),
    data.frame(Yhat)
  )
  
  instance = list(
    data_long=data %>% select(subject_id, time, group, response),
    data_wide=data_wide, 
    data_wide_imputed=data_wide_imputed,
    true_values=true_values,
    times=times,
    colnames=list(
      subject_id="subject_id",
      group="group",
      long_index="time",
      long_response="response",
      wide_times=paste0("t", times)
    )
  )
  
  return(instance)
}

