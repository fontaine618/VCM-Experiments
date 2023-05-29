lsvcmm_wrapper = function(data, job, instance, selection="aic", boot=F, 
                          estimate_variance_components=F,
                          random_effect=T, ...){
  t0 = proc.time()
  df = instance$data_long
  out = VCMM::lsvcmm(
    data=df, 
    response="response",
    subject="subject_id",
    time="time",
    vc_covariates="group",
    cv=ifelse(selection=="cv", 5, 0),
    estimate_variance_components=estimate_variance_components,
    random_effect=random_effect,
    ...
  )
  best_idx = switch(selection,
                "aic"=which.min(out$models_path$aic),
                "bic"=which.min(out$models_path$bic),
                "cv"=which.min(out$models_path$cv_Score)
  )
  lambda = out$models_path$lambda[best_idx]
  kernel_scale = out$models_path$kernel_scale[best_idx]
  
  b1true = instance$true_values$b1
  
  if(!boot){
    b = out$vc_path[2,,best_idx]
    vc = data.frame(
      time=out$estimated_time,
      estimate=b,
      lower=NA,
      upper=NA,
      pvalue=NA,
      adj_pvalue=NA,
      excludes_zero=1*(abs(b) > 1e-10),
      true_value=b1true
    )
  }else{
    out_boot = VCMM::lsvcmm.boot(
      data=df, 
      response="response",
      subject="subject_id",
      time="time",
      vc_covariates="group",
      lambda=lambda, 
      kernel_scale=kernel_scale,
      estimate_variance_components=estimate_variance_components,
      random_effect=random_effect,
      ...
    )
    band_df = VCMM::confidence_band(out_boot, 0.95, "quantile", 2)
    b = out_boot$vc[2, ]
    # b = apply(out_boot$vc_boot[2,,], 1, median)
    
    vc = data.frame(
      time=out$estimated_time,
      estimate=b,
      lower=band_df$L,
      upper=band_df$U,
      pvalue=NA,
      adj_pvalue=NA,
      excludes_zero=band_df$excludes_zero,
      true_value=b1true
    )
  }
  
  log = list(
    lambda=lambda,
    kernel_scale=kernel_scale
  )
  metrics = evaluate(vc, b1true)
  
  return(list(
    metrics=metrics,
    vc=vc, 
    log=log,
    time=proc.time()-t0
  ))
}