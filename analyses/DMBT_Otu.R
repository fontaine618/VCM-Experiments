library(tidyverse)
library(magrittr)
library(VCMM)
theme_set(theme_minimal())
extrafont::font_import(paths = "~/.local/share/fonts/", prompt=F)



# ==============================================================================
# PREPARE DATA
source("/home/simon/Documents/VCM-Experiments/process_data.R")

# filter Otus
prevalent_otus = microbiome::core(pseq, detection=0, prevalence=0.05) %>% phyloseq::taxa_names()
pseq %<>% microbiome::transform(transform="clr")
otus = pseq %>% phyloseq::taxa_names()
pseq %<>% phyloseq::subset_taxa(otus %in% prevalent_otus)

# Settings
FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_Otu/"
RESULT_PATH = "/home/simon/Documents/VCM-Experiments/results/DMBT_Otu/"
t0 = c(0, 4, 8, 12, 16, 22)
nt = length(t0)
kernel_scale = pracma::logseq(1/20, 1, 20)
adaptive = 0
sgl = 1.
n_samples = 1000

# Display values
vc_vars = c("scc", "group", "g_scc")
nvc_vars = c("sex")

display_varnames = c(
  "(Intercept)"="(Intercept)",
  "scc"="SCC",
  "group"="KO",
  "g_scc"="KO:SCC"
  )

design_diffs = matrix(c(
  0, 1, 0, 0, 
  0, 0, 1, 0, 
  0, 1, 0, 1, 
  0, 0, 1, 1
), nrow=4, ncol=4, byrow=T)


display_diffs = c(
  "scc"="WT: HP/CIS->SCC",
  "group"="HP/CIS: WT->KO",
  "scc_g_scc"="KO: HP/CIS->SCC",
  "group_g_scc"="SCC: WT->KO"
)

design_groups = matrix(c(
  1, 0, 0, 0, 
  1, 1, 0, 0, 
  1, 0, 1, 0, 
  1, 1, 1, 1
), nrow=4, ncol=4, byrow=T)

display_groupnames = c(
    "(Intercept)"="WT, HP/CIS",
    "scc"="WT, SCC",
    "group"="KO, HP/CIS",
    "g_scc"="KO, SCC"
  )

display_colors = c(
    "(Intercept)"="HP/CIS",
    "scc"="SCC",
    "group"="HP/CIS",
    "g_scc"="SCC"
  )

display_linetype = c(
    "(Intercept)"="WT",
    "scc"="WT",
    "group"="KO",
    "g_scc"="KO"
  )

# ------------------------------------------------------------------------------




# ==============================================================================
# PREPARE DATA
tax = phyloseq::tax_table(pseq)
clr = phyloseq::otu_table(pseq) %>% data.frame()
taxas = colnames(clr)
meta = phyloseq::sample_data(pseq)
results = list()
# ------------------------------------------------------------------------------



for(otu in taxas){
  
  # ==============================================================================
  # PREPARE DATA
  df = bind_cols(data.frame(meta), clr[, otu])
  colnames(df) = c("subject_id", "time", "Type", "Gender", "Diagnosis", "Diagnosis2", "CLR")
  df %<>% mutate(
    group=ifelse(Type=="KO", 1, 0),
    sex=ifelse(Gender=="F", 1, 0),
    hp=ifelse(Diagnosis2=="hyperplasia", 1, 0),
    cis=ifelse(Diagnosis2=="CIS", 1, 0),
    scc=ifelse(Diagnosis2=="SCC", 1, 0)
  ) %>% mutate(
    cis_or_scc=pmax(cis, scc)
  ) %>% mutate(
    g_hp=group*hp,
    g_cis=group*cis, 
    g_scc=group*cis,
    g_cis_or_scc=group*cis_or_scc,
    g_sex=group*sex
  )
  varnames = c("(Intercept)", vc_vars)
  
  taxonomy = as.vector(tax[otu, ])
  while(is.na(tail(taxonomy, 1))) taxonomy = head(taxonomy, -1)
  best_taxonomy = tail(taxonomy, 1)
  # ------------------------------------------------------------------------------
  
  
  
  # ==============================================================================
  # RUN CV
  fit = lsvcmm(
    data=df,
    response="CLR",
    subject="subject_id",
    vc_covariates=vc_vars,
    nvc_covariates=nvc_vars,
    time="time",
    estimated_time=t0,
    kernel_scale=kernel_scale,
    n_lambda=100,
    lambda_factor=0.0005,
    adaptive=adaptive,
    sgl=sgl
  )
  
  which_aic = fit$models_path$aic %>% which.min
  h_aic = fit$models_path$kernel_scale[which_aic]
  l_aic = fit$models_path$lambda[which_aic]
  b_aic = fit$vc_path[,,which_aic]
  a_aic = fit$nvc_path[,which_aic]
  sig2_aic = fit$models_path$sig2[which_aic]
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # RUN Bootstrap
  obj = lsvcmm.boot(
    data=df,
    response="CLR",
    subject="subject_id",
    vc_covariates=vc_vars,
    nvc_covariates=nvc_vars,
    time="time",
    estimated_time=t0,
    kernel_scale=h_aic,
    lambda=l_aic,
    adaptive=adaptive,
    sgl=sgl,
    n_samples=n_samples
  )
  grp_obj = transform_obj(obj, design_groups)
  diff_obj = transform_obj(obj, design_diffs)
  coef_ci = confidence_band(obj, 0.95, "quantile", 1:4)
  grp_ci = confidence_band(grp_obj, 0.95, "quantile", 1:4)
  diff_ci = confidence_band(diff_obj, 0.95, "quantile", 1:4)
  sex_ci = obj$nvc_boot %>% quantile(c(0.5, 0.025, 0.975))
  # ------------------------------------------------------------------------------
  
  
  
  
  
  
  # ==============================================================================
  # Plot data
  df %<>% mutate(Diagnosis3 = ifelse(scc, "SCC", "HP/CIS"))
  df %<>% mutate(Group = paste(Type, Diagnosis3, sep=", "))
  g_data = ggplot(
    data=df %>% mutate(CLR=robustbase::huberize(CLR, k=3)),
    mapping=aes(x=time, y=CLR, color=Type, fill=Type, group=paste(Group, time),
                linetype=Diagnosis3)
  ) + 
    geom_boxplot(width=1, alpha=0.5) +
    ylab("CLR (trimmed)") + 
    xlab("Week") +
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none", 
          text=element_text(family="Roboto"),
          plot.title=element_text(size=20)) +
    ggtitle(otu)
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Plot Groups
  names(grp_ci) = display_groupnames
  grp_ci_df = grp_ci %>% bind_rows(.id="group")
  grp_ci_df$Type = sapply(grp_ci_df$group, function(g) stringr::str_split(g, ", ")[[1]][1])
  grp_ci_df$Diagnosis = sapply(grp_ci_df$group, function(g) stringr::str_split(g, ", ")[[1]][2])
  g_grps = ggplot(
    data=grp_ci_df,
    mapping=aes(x=estimated_time, y=median, ymin=L, ymax=U, 
                color=Type, group=group, fill=Type, linetype=Diagnosis)
  ) + 
    geom_line() + 
    geom_ribbon(alpha=0.1, color=NA) +
    # geom_errorbar(width=0.5) + 
    ylab("Fitted mean (CLR)") + 
    xlab("Week") +
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) +
    labs(color="Group", fill="Group") + 
    theme(text=element_text(family="Roboto"), legend.position="bottom")
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Plot Differences
  names(diff_ci) = names(display_diffs)
  yrange = range(diff_obj$vc_boot)
  g_diffs = lapply(seq_along(display_diffs), function(i){
    diff_name = names(display_diffs)[i]
    g = ggplot() +
      geom_hline(yintercept=0, linewidth=0.5, linetype="dashed") +
      geom_line(
        data=diff_ci[[diff_name]],
        mapping=aes(x=estimated_time, y=median),
        linewidth=1
      ) + 
      geom_errorbar(
        data=diff_ci[[diff_name]],
        mapping=aes(x=estimated_time, ymin=L, ymax=U, color=excludes_zero),
        linewidth=1,
        width=0.5
      ) + 
      ggtitle(display_diffs[diff_name]) + 
      ylab("") + xlab("Week") + 
      scale_color_manual(breaks=c(F, T), values=c("black", "red")) + 
      scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
      theme(text=element_text(family="Roboto"), legend.position="none") +
      coord_cartesian(ylim=yrange)
    if(i %% 2 == 0) g = g +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    if(i <3) g = g +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    return(g)
  })
  
  g_diff = cowplot::plot_grid(plotlist=g_diffs, ncol=2, align="v", axis="lr", byrow=T)
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Additional info
  sex = paste0("Sex (M->F): ", round(sex_ci[1], 2), " [", round(sex_ci[2], 2), ",",round(sex_ci[3], 2), "]")
  selection = paste0(
    "AIC selection: ", 
    paste("h", round(h_aic, 4), sep="="), ", ",
    paste("lambda", round(l_aic, 4), sep="=")
  )
  parms = paste0(
    "Parameters: ",
    "sgl=", sgl, ", adaptive=", adaptive
  )
  info = paste(parms, selection, sex, sep="\n")
  g_text = ggplot() + 
    ggtitle(info) + 
    theme_void()  + 
    theme(text=element_text(family="Roboto", size=8))
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # All together
  g = cowplot::plot_grid(g_data, g_grps, g_diff, g_text, nrow=4, align="v", axis="lr", rel_heights=c(1, 1, 2, 0.3))
  
  
  # ggsave(paste0(FIG_PATH, taxa_name, ".png"), g, width=8, height=10, bg="white")
  ggsave(paste0(FIG_PATH, otu, ".png"), g, width=10, height=16, bg="white")
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Save results
  names(diff_ci) = display_diffs
  diff_ci_df = diff_ci %>% bind_rows(.id="diff")
  diff_ci_df$Within = sapply(diff_ci_df$diff, function(g) stringr::str_split(g, ": ")[[1]][1])
  diff_ci_df$Change = sapply(diff_ci_df$diff, function(g) stringr::str_split(g, ": ")[[1]][2])
  names(coef_ci) = display_varnames
  coef_ci_df = coef_ci %>% bind_rows(.id="varname")
  
  coef_ci_df$grouping = "coefs"
  grp_ci_df$grouping = "groups"
  diff_ci_df$grouping = "diffs"
  sex_ci_df = data.frame(
    grouping="sex",
    median=sex_ci[1],
    mean=obj$nvc_boot %>% mean,
    estimate=a_aic,
    L=sex_ci[2],
    U=sex_ci[3],
    excludes_zero=(sex_ci[3]<0) || (sex_ci[2]>0)
  )
  
  outdf = bind_rows(sex_ci_df, grp_ci_df, diff_ci_df, coef_ci_df)
  outdf$otu = otu
  outdf$taxonomy = best_taxonomy
  results[[otu]] = outdf
  # ------------------------------------------------------------------------------
  
}

write_csv(bind_rows(results), paste0(RESULT_PATH, "coefs.csv"))
