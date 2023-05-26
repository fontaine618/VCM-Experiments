gam_wrapper = function(data, job, instance, df=6, ...){
  require(mgcv)
  require(gratia)
  
  t0 = proc.time()
  dat = instance$data_long
  dat$group = as.factor(dat$group)
  
  out = mgcv::gamm(
    response ~ s(time, k=df, bs="cr", by=group),
    family=gaussian(),
    random=list(subject_id=~1),
    data=dat
  )
  
  pred_dat = data.frame(
    time=instance$times,
    group=factor(1, levels=0:1)
  )
  
  ci = confint(out, parm="s(time):group1", data=pred_dat, type="simultaneous", partial_match=T)
  
  b1true = instance$true_values$b1
  
  vc = data.frame(
    time=instance$times,
    estimate=ci$est,
    lower=ci$lower,
    upper=ci$upper,
    pvalue=NA,
    adj_pvalue=NA,
    excludes_zero=(ci$lower>0) | (ci$upper<0),
    true_value=b1true
  )
  
  log = list()
  metrics = evaluate(vc, b1true)
  
  return(list(
    metrics=metrics,
    vc=vc, 
    log=log,
    time=proc.time()-t0
  ))
}