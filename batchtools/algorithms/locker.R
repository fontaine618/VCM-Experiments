locker_wrapper = function(data, job, instance, ...){
  require(dplyr)
  
  t0 = proc.time()
  dat = instance$data_long  
  
  subjects = unique(dat$subject_id)
  
  X = lapply(subjects, function(i) dat %>% filter(subject_id==i) %>% pull(group))
  X_obser = lapply(subjects, function(i) dat %>% filter(subject_id==i) %>% pull(time))
  X_obser_num = sapply(X, length)
  Y = lapply(subjects, function(i) dat %>% filter(subject_id==i) %>% pull(response))
  Y_obser = lapply(subjects, function(i) dat %>% filter(subject_id==i) %>% pull(time))
  Y_obser_num = sapply(Y, length)
  
  out = LocKer::LocKer(
    X=X, X_obser=X_obser, X_obser_num=X_obser_num,
    Y=Y, Y_obser=Y_obser, Y_obser_num=Y_obser_num,
    family="Gaussian", 
    timeint=c(-0.1, 1.1),
    L_list=c(11), 
    absTol_list=rep(1e-4),
    roupen_para_list=10 ^ seq(-5, -2, length.out=21),
    lambda_list=10 ^ seq(-5, -2, length.out=21)
  )
  
  # b0 = fda::eval.fd(instance$times, out$beta0fd_est)
  b1 = fda::eval.fd(instance$times, out$betafd_est)
  
  b1true = instance$true_values$b1
  
  vc = data.frame(
    time=instance$times,
    estimate=b1,
    lower=NA,
    upper=NA,
    pvalue=NA,
    adj_pvalue=NA,
    excludes_zero=1*(abs(b1) > 1e-10),
    true_value=b1true
  )
  
  log = list(K=out$L, alpha=out$roupen_select, lambda=out$lambda_select)
  metrics = evaluate(vc, b1true)
  
  return(list(
    metrics=metrics,
    vc=vc, 
    log=log,
    time=proc.time()-t0
  ))
}