evaluate = function(output, b){
  p = abs(b) > 1e-10
  pos = factor(as.numeric(p), levels=c(0, 1))
  ppos = factor(as.numeric(output$excludes_zero), levels=c(0, 1))
  estpost = factor(as.numeric(abs(output$est)>1e-10), levels=c(0, 1))
  cmat = caret::confusionMatrix(reference=pos, data=ppos, mode="everything", positive="1")
  cmatest = caret::confusionMatrix(reference=pos, data=estpost, mode="everything", positive="1")
  best = output$estimate
  rmse = sqrt(mean((best-b)^2))
  rmse0 = sqrt(mean((best-b)[!p]^2))
  rmse1 = sqrt(mean((best-b)[p]^2))
  rmse = c(global=rmse, null=rmse0, nonnull=rmse1)
  b0 = any(p)
  best0 = any(output$excludes_zero == 1)
  metrics = list(
    cmat=cmat, 
    cmatest=cmatest, 
    rmse=rmse,
    global=1*(b0==best0)
  )
  return(metrics)
}