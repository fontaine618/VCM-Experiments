library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(VCMM)
theme_set(theme_minimal())

setwd("/home/simon/Documents/VCM-Experiments/")

dir_figures = "./figures/adaptive1_nt6/"
load("./data/BWmissing_tensor.Rdata")
demo = read.csv("./data/demo.csv")


subjects = dimnames(clrX)[[1]]
taxas = dimnames(clrX)[[2]]
timepoints = dimnames(clrX)[[3]]
subject_df = data.frame(id=as.integer(subjects))
tobs = c(0, 4, 8, 12, 16, 22)
# t0 = seq(0, 22, 0.5)
t0 = tobs
nt = length(t0)
breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1., 3.)
adaptive = 1
kernel_scale = pracma::logseq(1/20, 5, 20)
cv = 0

for(taxa in taxas){
  
  taxa_df = data.frame(clrX[,taxa,])
  taxa_df %<>% bind_cols(subject_df)
  taxa_df %<>% left_join(demo, by=c("id"="ID"))
  taxa_df %<>% select(-Diagnosis)
  taxa_df %<>% pivot_longer(-c(id, Group, Gender), names_to="time", values_to="clr", values_drop_na=T)
  taxa_df %<>% mutate(time=stringr::str_sub(time, 2))
  taxa_df %<>% mutate(time=as.integer(time))
  
  taxa_df %<>% mutate(
    group=ifelse(Group=="K", 1, 0),
    sex=ifelse(Gender=="M", 1, 0)
  )
  
  out = lsvcmm(
    data=taxa_df,
    response="clr",
    subject="id",
    vc_covariates="group",
    nvc_covariates="sex",
    time="time",
    estimated_time=t0,
    kernel_scale=kernel_scale,
    n_lambda=100,
    lambda_factor=0.0005,
    cv=cv,
    cv_seed=0,
    adaptive=adaptive
  )
  
  
  which_bic = out$models_path$bic %>% which.min
  which_aic = out$models_path$aic %>% which.min
  # which_cv = out$models_path$predparss %>% which.min
  
  h_aic = out$models_path$kernel_scale[which_aic]
  h_bic = out$models_path$kernel_scale[which_bic]
  # h_cv = out$models_path$kernel_scale[which_cv]
  
  l_aic = out$models_path$lambda[which_aic]
  l_bic = out$models_path$lambda[which_bic]
  # l_cv = out$models_path$lambda[which_cv]
  
  b0_aic = out$vc_path[1,,which_aic]
  b0_bic = out$vc_path[1,,which_bic]
  # b0_cv = out$vc_path[1,,which_cv]
  
  b1_aic = out$vc_path[2,,which_aic]
  b1_bic = out$vc_path[2,,which_bic]
  # b1_cv = out$vc_path[2,,which_cv]
  
  curve_df = bind_rows(
    data.frame(b1=b1_aic, b0=b0_aic, criterion=paste0("aic(h=", round(h_aic, 3), ",l=", round(l_aic, 5), ")"), t=t0),
    data.frame(b1=b1_bic, b0=b0_bic, criterion=paste0("bic(h=", round(h_bic, 3), ",l=", round(l_bic, 5), ")"), t=t0)
    # data.frame(b1=b1_cv, b0=b0_cv, criterion=paste0("cv(h=", round(h_cv, 3), ",l=", round(l_cv, 5), ")"), t=t0)
  )
  yrange = quantile(taxa_df$clr, c(0.025, 0.975))
  drange = yrange - mean(yrange)
  
  g0 = ggplot(
    data=curve_df,
    mapping=aes(x=t, y=b0, group=criterion, color=criterion)
  ) + geom_line() + 
    ggtitle(paste0(taxa, " (", nt, "time points)")) + 
    ylab("Group W estimate") + coord_cartesian(ylim=yrange) + 
    geom_vline(xintercept=tobs, linetype="dashed", alpha=0.5)
  
  
  gdiff = ggplot(
    data=curve_df,
    mapping=aes(x=t, y=b1, group=criterion, color=criterion)
  ) + geom_line() + 
    ylab("K-W estimate") + coord_cartesian(ylim=drange) + 
    geom_vline(xintercept=tobs, linetype="dashed", alpha=0.5)
  
  
  g1 = ggplot(
    data=curve_df,
    mapping=aes(x=t, y=b0+b1, group=criterion, color=criterion)
  ) + geom_line() + 
    ylab("Group K estimate") + coord_cartesian(ylim=yrange) + 
    geom_vline(xintercept=tobs, linetype="dashed", alpha=0.5)
  
  gx = ggplot(
    data=taxa_df,
    mapping=aes(x=time, y=clr, color=Group, group=Group)
  ) + geom_smooth() + coord_cartesian(ylim=yrange) + 
    geom_vline(xintercept=tobs, linetype="dashed", alpha=0.5)
  
  # gcv = ggplot(
  #   data=out$models_path,
  #   mapping=aes(x=lambda, y=predparss, group=kernel_scale, color=kernel_scale)
  # ) + geom_line() +
  #   scale_x_log10() + 
  #   scale_color_gradientn(breaks=breaks, trans="log", colors=rainbow(5)) +
  #   ylab("CV Error")
  
  gaic = ggplot(
    data=out$models_path,
    mapping=aes(x=lambda, y=aic, group=kernel_scale, color=kernel_scale)
  ) + geom_line() +
    scale_x_log10() + 
    scale_color_gradientn(breaks=breaks, trans="log", colors=rainbow(5)) +
    ylab("AIC")
  
  gbic = ggplot(
    data=out$models_path,
    mapping=aes(x=lambda, y=bic, group=kernel_scale, color=kernel_scale)
  ) + geom_line() +
    scale_x_log10() + 
    scale_color_gradientn(breaks=breaks, trans="log", colors=rainbow(5)) +
    ylab("BIC")
  
  g = ggpubr::ggarrange(g0, g1, gdiff, gx, gaic, gbic, ncol=1, align="hv")
  ggsave(paste0(dir_figures, stringr::str_replace_all(taxa, "/", "_"), ".pdf"), g, width=8, height=20)
  
}