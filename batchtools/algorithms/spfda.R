spfda_wrapper = function(data, job, instance, random_effect=T, ...){
  require(dplyr)
  
  t0 = proc.time()
  dat = instance$data_wide_imputed
  X = dat %>% select(group) %>% as.matrix()
  X = cbind(1, X)
  Y = dat %>% select(starts_with("t")) %>% as.matrix()
  t = instance$times
  
  if(random_effect){
    W = spfda::spfda_weight(
      X=X,
      Y=Y,
      bandwidth=2,
      part=list(seq(length(t)))
    )
  }else{
    W = NULL
  }
  
  alphas = seq(0.1, 0.9, by = 0.05)
  lambdas = 10^seq(log10(0.1), log10(10), length.out = 20)
  Ks = seq(5, 10)
  all_params = expand.grid(lambdas, alphas, Ks)
  names(all_params) = c("lambda", "alpha", "K")
  
  BICs = lapply(seq_len(nrow(all_params)), function(ii){
    param = all_params[ii,]
    lambda = param$lambda
    alpha = param$alpha
    K = param$K
    res = spfda::spfda(
      Y=Y,
      X=X,
      time=instance$times,
      lambda=lambda,
      alpha=alpha,
      nsp=K,
      ord=3,
      CI=F,
      W=W
    )
    BIC(res)
  })
  
  BICfd = cbind(all_params, BIC=unlist(BICs))
  parms = BICfd[which.min(BICfd$BIC),  ]
  
  out = spfda::spfda(
    Y=Y,
    X=X,
    time=instance$times,
    lambda=parms$lambda,
    alpha=parms$alpha,
    nsp=parms$K,
    ord=3,
    CI=T,
    W=W
  )
  
  est = out$get_coef()[2,]
  se = out$get_se()[2,]
  lower = est - 2*se
  upper = est + 2*se
  
  b1true = instance$true_values$b1
  
  vc = data.frame(
    time=instance$times,
    estimate=est,
    lower=lower,
    upper=upper,
    pvalue=NA,
    adj_pvalue=NA,
    excludes_zero=(lower>0) | (upper<0),
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