library(tidyverse)
library(magrittr)
library(VCMM)
theme_set(theme_minimal())
extrafont::font_import(paths = "~/.local/share/fonts/", prompt=F)



# ==============================================================================
# PREPARE DATA
source("/home/simon/Documents/VCM-Experiments/process_data.R")
pseq = pseq
# pseq = pseq_genus
FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_2way_genotype_sex/"

otus = c(1, 2, 6, 7, 11, 13, 14, 15, 17, 21, 45)
selected_taxas = c(
  "g__Rothia",
  "g__Cupriavidus",
  "g__Bradyrhizobium",
  "g__Atopostipes",
  "g__Eisenbergiella",
  "g__Anaerofustis",
  "g__Clostridium_XlVb",
  "g__Caulobacter",
  "g__Akkermansia",
  "Unknown"
)

t0 = c(0, 4, 8, 12, 16, 22)
nt = length(t0)
kernel_scale = pracma::logseq(1/20, 1, 20)
adaptive = 0
n_samples = 1000

vc_vars = c("sex", "group", "g_sex")

display_varnames = c(
    "(Intercept)"="WT, M",
    "sex"="WT, F-M",
    "group"="KO-WT, M",
    "g_sex"="KO-WT, F-M"
  )

design = matrix(c(
  1, 0, 0, 0, 
  1, 1, 0, 0, 
  1, 0, 1, 0, 
  1, 1, 1, 1
), nrow=4, ncol=4, byrow=T)


display_groupnames = c(
    "(Intercept)"="WT, M",
    "scc"="WT, F",
    "group"="KO, M",
    "g_sex"="KO, F"
  )

display_colors = c(
    "(Intercept)"="M",
    "sex"="F",
    "group"="M",
    "g_sex"="F"
  )

display_linetype = c(
    "(Intercept)"="WT",
    "sex"="WT",
    "group"="KO",
    "g_sex"="KO"
  )


# ------------------------------------------------------------------------------




# ==============================================================================
# PREPARE DATA
tax = phyloseq::tax_table(pseq)
otu = phyloseq::otu_table(pseq)
clr = otu %>% data.frame()
clr[clr==0] = 1
clr = log(clr)
clr = scale(clr, scale=F)
clr = t(clr)
taxas = colnames(clr)
meta = phyloseq::sample_data(pseq)
# ------------------------------------------------------------------------------



for(i in otus){
# for(taxa_name in selected_taxas){
  
  # taxa = taxa_name
  otu_name = paste0("Otu", stringr::str_pad(i, 4, pad="0"))
  taxa = which.max(colnames(otu) == otu_name)
  taxa_name = tax[taxa, ]
  taxa_name = taxa_name[1, max(1, which.max(is.na(taxa_name))-1)]
  
  # ==============================================================================
  # PREPARE DATA
  df = bind_cols(data.frame(meta), clr[taxa, ])
  # df = bind_cols(data.frame(meta), clr[, taxa])
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
  # ------------------------------------------------------------------------------
  
  
  
  # ==============================================================================
  # RUN CV
  fit = lsvcmm(
    data=df,
    response="CLR",
    subject="subject_id",
    vc_covariates=vc_vars,
    nvc_covariates=NULL,
    time="time",
    estimated_time=t0,
    kernel_scale=kernel_scale,
    n_lambda=100,
    lambda_factor=0.0005,
    adaptive=adaptive
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
    nvc_covariates=NULL,
    time="time",
    estimated_time=t0,
    kernel_scale=h_aic,
    lambda=l_aic,
    adaptive=adaptive,
    n_samples=n_samples
  )
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Plot curves
  px = length(varnames)
  gs = list()
  yrange = range(obj$vc_boot)
  
  for(j in seq(px)){
    is_intercept = j==1
    band_df = confidence_band(obj, 1-0.05/px, "quantile", j)
    boot_samples = data.frame(obj$vc_boot[j, , ])
    boot_samples$t = t0
    boot_samples %<>% tidyr::pivot_longer(-t)
    
    g = ggplot() + 
      geom_hline(yintercept=0, linewidth=0.5, linetype="dashed") + 
      geom_line(
        data=boot_samples,
        mapping=aes(x=t, y=value, group=name),
        alpha = 0.01
      ) + 
      geom_line(
        data=band_df,
        mapping=aes(x=estimated_time, y=estimate),
        linewidth=2
      ) + 
      geom_errorbar(
        data=band_df, 
        mapping=aes(x=estimated_time, ymin=L, ymax=U, color=excludes_zero),
        width=0.5,
        linewidth=1
      ) + xlab("") + ggtitle(display_varnames[varnames[j]]) + ylab("") + 
      scale_color_manual(breaks=c(F, T), values=c("black", "red")) + 
      scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
      theme(text=element_text(family="Roboto"), legend.position="none") +
      coord_cartesian(ylim=yrange)
    gs[[j]] = g
  }
  
  gb = cowplot::plot_grid(plotlist=gs, ncol=2, align="v", axis="lr", byrow=F)
  
  df_groups = design %*% obj$vc
  rownames(df_groups) = display_groupnames
  colnames(df_groups) = t0
  
  display_df = data.frame(
    var=varnames,
    name=display_groupnames,
    color=display_colors,
    linetype=display_linetype
  )
  
  df_groups_long = data.frame(
    b = as.vector(df_groups %>% t),
    group = as.vector(sapply(display_groupnames, function(x) rep(x, 6))),
    week = t0
  )
  
  df_groups_long %<>% left_join(display_df, by=c("group"="name"))
  
  gg = ggplot(
    data=df_groups_long,
    mapping=aes(x=week, y=b, group=group, color=color, linetype=linetype)
  ) + geom_line(linewidth=1) +
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
    ylab("") + 
    ggtitle(paste0(otu_name, " (", taxa_name, ")")) +
    # ggtitle(paste("Genus:", taxa_name)) + 
    theme(legend.position="bottom") + 
    guides(
      color=guide_legend(title="Sex"),
      linetype=guide_legend(title="Type")
    )
  
  g = cowplot::plot_grid(gg, gb, nrow=2, align="v", axis="lr", rel_heights=c(1, 2))
  
  
  # ggsave(paste0(FIG_PATH, taxa_name, ".png"), g, width=8, height=10, bg="white")
  ggsave(paste0(FIG_PATH, otu_name, ".png"), g, width=8, height=10, bg="white")
  # ------------------------------------------------------------------------------
  
  
  
}
