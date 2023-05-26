library(tidyverse)
library(magrittr)
library(VCMM)
theme_set(theme_minimal())
extrafont::font_import(paths = "~/.local/share/fonts/", prompt=F)



# ==============================================================================
# PREPARE DATA
source("/home/simon/Documents/VCM-Experiments/process_data.R")
# pseq = pseq
pseq = pseq_genus
FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_2way_genotype_histo/"

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

method = "monotone"


# HP vs CIS/SCC: use cis_or_scc
# HP/CIS vs SCC: use SCC
# HP -> CIS -> SCC: use cis_or_scc, scc
vc_vars = switch(
  method,
  "monotone"=c("cis_or_scc", "scc", "group", "g_cis_or_scc", "g_scc"),
  "hp_vs_cis_or_scc"=c("cis_or_scc", "group", "g_cis_or_scc"),
  "hp_or_cis_vs_scc"=c("scc", "group", "g_scc")
)

display_varnames = switch(
  method,
  "monotone"=c(
    "(Intercept)"="WT, HP",
    "cis_or_scc"="WT, CIS-HP",
    "scc"="WT, SCC-CIS",
    "group"="KO-WT, HP",
    "g_cis_or_scc"="KO-WT, CIS-HP",
    "g_scc"="KO-WT, SCC-CIS"
  ),
  "hp_vs_cis_or_scc"=c(
    "(Intercept)"="WT, HP",
    "cis_or_scc"="WT, CIS/SCC-HP",
    "group"="KO-WT, HP",
    "g_cis_or_scc"="KO-WT, CIS/SCC-HP"
  ),
  "hp_or_cis_vs_scc"=c(
    "(Intercept)"="WT, HP/CIS",
    "scc"="WT, SCC-HP/CIS",
    "group"="KO-WT, HP/CIS",
    "g_scc"="KO-WT, SCC-HP/CIS"
  )
)

design4 = matrix(c(
  1, 0, 0, 0, 
  1, 1, 0, 0, 
  1, 0, 1, 0, 
  1, 1, 1, 1
), nrow=4, ncol=4, byrow=T)


design6 = matrix(c(
  1, 0, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 
  1, 1, 1, 0, 0, 0, 
  1, 0, 0, 1, 0, 0, 
  1, 1, 0, 1, 1, 0, 
  1, 1, 1, 1, 1, 1
), nrow=6, ncol=6, byrow=T)

design = switch(
  method,
  "monotone"=design6, 
  "hp_vs_cis_or_scc"=design4,
  "hp_or_cis_vs_scc"=design4
)

display_groupnames = switch(
  method,
  "monotone"=c(
    "(Intercept)"="WT, HP",
    "cis_or_scc"="WT, CIS",
    "scc"="WT, SCC",
    "group"="KO, HP",
    "g_cis_or_scc"="KO, CIS",
    "g_scc"="KO, SCC"
  ),
  "hp_vs_cis_or_scc"=c(
    "(Intercept)"="WT, HP",
    "cis_or_scc"="WT, CIS/SCC",
    "group"="KO, HP",
    "g_cis_or_scc"="KO, CIS/SCC"
  ),
  "hp_or_cis_vs_scc"=c(
    "(Intercept)"="WT, HP/CIS",
    "scc"="WT, SCC",
    "group"="KO, HP/CIS",
    "g_scc"="KO, SCC"
  )
)

display_colors = switch(
  method,
  "monotone"=c(
    "(Intercept)"="HP",
    "cis_or_scc"="CIS",
    "scc"="SCC",
    "group"="HP",
    "g_cis_or_scc"="CIS",
    "g_scc"="SCC"
  ),
  "hp_vs_cis_or_scc"=c(
    "(Intercept)"="HP",
    "cis_or_scc"="CIS/SCC",
    "group"="HP",
    "g_cis_or_scc"="CIS/SCC"
  ),
  "hp_or_cis_vs_scc"=c(
    "(Intercept)"="HP/CIS",
    "scc"="SCC",
    "group"="HP/CIS",
    "g_scc"="SCC"
  )
)

display_linetype = switch(
  method,
  "monotone"=c(
    "(Intercept)"="WT",
    "cis_or_scc"="WT",
    "scc"="WT",
    "group"="KO",
    "g_cis_or_scc"="KO",
    "g_scc"="KO"
  ),
  "hp_vs_cis_or_scc"=c(
    "(Intercept)"="WT",
    "cis_or_scc"="WT",
    "group"="KO",
    "g_cis_or_scc"="KO"
  ),
  "hp_or_cis_vs_scc"=c(
    "(Intercept)"="WT",
    "scc"="WT",
    "group"="KO",
    "g_scc"="KO"
  )
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



# for(i in otus){
for(taxa_name in selected_taxas){
  
  taxa = taxa_name
  # otu_name = paste0("Otu", stringr::str_pad(i, 4, pad="0"))
  # taxa = which.max(colnames(otu) == otu_name)
  # taxa_name = tax[taxa, ]
  # taxa_name = taxa_name[1, max(1, which.max(is.na(taxa_name))-1)]
  
  # ==============================================================================
  # PREPARE DATA
  # df = bind_cols(data.frame(meta), clr[taxa, ])
  df = bind_cols(data.frame(meta), clr[, taxa])
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
    g_cis_or_scc=group*cis_or_scc
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
    nvc_covariates="sex",
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
    nvc_covariates="sex",
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
    # ggtitle(paste0(otu_name, " (", taxa_name, ")")) + 
    ggtitle(paste("Genus:", taxa_name)) + 
    theme(legend.position="bottom") + 
    guides(
      color=guide_legend(title="Diagnosis"),
      linetype=guide_legend(title="Type")
    )
  
  g = cowplot::plot_grid(gg, gb, nrow=2, align="v", axis="lr", rel_heights=c(1, 2))
  
  
  ggsave(paste0(FIG_PATH, method, "/", taxa_name, ".png"), g, width=8, height=10, bg="white")
  # ggsave(paste0(FIG_PATH, method, "/", otu_name, ".png"), g, width=8, height=10, bg="white")
  # ------------------------------------------------------------------------------
  
  
  
}
