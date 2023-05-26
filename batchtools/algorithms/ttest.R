ttest_wrapper = function(data, job, instance, method="hommel", ...){
  require(dplyr)
  
  t0 = proc.time()
  dat = instance$data_long  
  
  by_time = sapply(instance$times, function(t){
    datt = dat %>% filter(time == t)
    y0 = datt %>% filter(group==0) %>% pull(response)
    y1 = datt %>% filter(group==1) %>% pull(response)
    out = t.test(y0, y1)
    c(
      estimate=diff(out$estimate),
      stderr=out$stderr,
      pvalue=out$p.value
    )
  })
  out = data.frame(t(by_time))
  colnames(out) = c("est", "se", "pvalue")
  out$time=instance$times
  mult = qnorm(0.025/length(instance$times), lower.tail=F)
  adj_pvalue = p.adjust(out$pvalue, method=method)
  
  b1true = instance$true_values$b1
  
  vc = data.frame(
    time=instance$times,
    estimate=out$est,
    lower=out$est - mult*out$se,
    upper=out$est + mult*out$se,
    pvalue=out$pvalue,
    adj_pvalue=adj_pvalue,
    excludes_zero=adj_pvalue < 0.05,
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