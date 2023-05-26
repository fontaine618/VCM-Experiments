library(tidyverse)
library(magrittr)
library(VCMM)
theme_set(theme_minimal())
extrafont::font_import(paths = "~/.local/share/fonts/", prompt=F)

DATA_PATH = "/home/simon/Documents/VCM-Experiments/data/DMBT_Processed/"
FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_Genus/"
aggregation = "none"
filename = c(
  "genus"="pseq_genus.Rdata",
  "family"="pseq_family.Rdata",
  "order"="pseq_order.Rdata",
  "unique"="pseq_unique.Rdata",
  "none"="pseq.Rdata"
)[aggregation]
name = c(
  "genus"="Genus",
  "family"="Family",
  "order"="Order",
  "unique"="Highest Known",
  "none"="OTU"
)[aggregation]
t0 = c(0, 4, 8, 12, 16, 22)
# t0 = seq(0, 22, 1)
nt = length(t0)
kernel_scale = pracma::logseq(1/20, 1, 20)
adaptive = 0
n_samples = 1000

# ==============================================================================
# LOAD DATA
# TODO: save pseq objects better
pseq = pseq_genus
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





# ==============================================================================
# Cross-sectional DA
pval = list()
lfc = list()
for( tt in t0){
  pseq_t = pseq %>% phyloseq::subset_samples(visit_id == tt)
  x = data.frame(Type=phyloseq::get_variable(pseq_t, "Type"))
  rownames(x) = phyloseq::sample_names(pseq_t)
  o = phyloseq::otu_table(pseq_t) %>% as.data.frame()
  deseq = DESeq2::DESeqDataSetFromMatrix(
    countData=as.matrix(o),
    colData=x,
    design=~Type
  )
  test_res = DESeq2::DESeq(deseq, test="Wald")
  test_res = DESeq2::results(test_res)
  pval[[as.character(tt)]] = test_res$pvalue
  lfc[[as.character(tt)]] = test_res$log2FoldChange
}
pvalmat = bind_cols(pval)
lfcmat = bind_cols(lfc)
rownames(pvalmat) = taxas
rownames(lfcmat) = taxas
colnames(pvalmat) = t0
colnames(lfcmat) = t0
# ------------------------------------------------------------------------------


out = list()
for(taxa in taxas){
  # ==============================================================================
  # PREPARE DATA
  df = bind_cols(data.frame(meta), clr[, taxa])
  colnames(df) = c("subject_id", "time", "Type", "Gender", "Diagnosis", "CLR")
  df %<>% mutate(
    group=ifelse(Type=="KO", 1, 0),
    sex=ifelse(Gender=="F", 1, 0)
  )
  # ------------------------------------------------------------------------------
  
  
  
  # ==============================================================================
  # RUN CV
  fit = lsvcmm(
    data=df,
    response="CLR",
    subject="subject_id",
    vc_covariates="group",
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
  b0_aic = fit$vc_path[1,,which_aic]
  b1_aic = fit$vc_path[2,,which_aic]
  sig2_aic = fit$models_path$sig2[which_aic]
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # RUN Bootstrap
  obj = lsvcmm.boot(
    data=df,
    response="CLR",
    subject="subject_id",
    vc_covariates="group",
    nvc_covariates="sex",
    time="time",
    estimated_time=t0,
    kernel_scale=h_aic,
    lambda=l_aic,
    adaptive=adaptive,
    n_samples=n_samples
  )
  band_df = confidence_band(obj, 0.95, "quantile", 2)
  
  boot_samples = data.frame(obj$vc_boot[2, , ])
  boot_samples$t = t0
  boot_samples %<>% tidyr::pivot_longer(-t)
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Plot
  
  gdata = 
    ggplot() + 
    geom_boxplot(
      data=df %>% mutate(CLR=robustbase::huberize(CLR, k=3)), 
      mapping=aes(x=time, y=CLR, fill=Type, group=paste(Type, time)),
      width=1
    ) + ylab("CLR (trimmed)") + labs(color="Type") + xlab("Week") + 
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom", 
          text=element_text(family="Roboto")) +
    ggtitle(taxa) 
  
  
  
  gb = ggplot() + 
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
    ) + xlab("") + ylab("KO-WT CLR") + guides(color=guide_legend(title="Locally DA")) + 
    scale_color_manual(breaks=c(F, T), values=c("black", "red")) + 
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
    ggtitle("LSVCMM") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom", 
          text=element_text(family="Roboto")) 
  
  g0 = ggplot(
    data=band_df,
    mapping=aes(x=estimated_time, y=prop0)
  ) + geom_bar(alpha=0.5, stat="identity", width=1) + xlab("Week") + ylab("Prop. non-zero") + 
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
    ylim(0, 1)  + 
    theme(text=element_text(family="Roboto"))
  
  csda_df = data.frame(
    time=t0,
    pval=pvalmat[taxa, ] %>% t,
    lfc=-lfcmat[taxa, ] %>% t,
    padj=p.adjust(pvalmat[taxa, ], method="bonferroni")
  ) %>% mutate(
    lda=padj < 0.05
  )
  
  gcsda = ggplot(
    data=csda_df,
    mapping=aes(x=t0, y=lfc, fill=lda)
  ) + geom_bar(stat="identity", width=1)+ xlab("Week") + ylab("Log2 Fold Change") + 
    scale_x_continuous(breaks=t0, limits=c(-1, 23), minor_breaks=NULL) + 
    theme(text=element_text(family="Roboto"), legend.position="none") + 
    scale_fill_manual(breaks=c(F, T), values=c("black", "red")) + 
    geom_hline(yintercept=0, linewidth=0.5, linetype="dashed")+
    ggtitle("DESeq2")
    
  
  
  g = cowplot::plot_grid(gdata, gb, gcsda, nrow=3, ncol=1, align="v", axis="lr", rel_heights=c(2, 2, 1))
  
  ggsave(paste0(FIG_PATH, "taxawise/", taxa, ".png"), g, width=8, height=8, bg="white")
  # ------------------------------------------------------------------------------
  
  
  
  
  
  # ==============================================================================
  # Store output
  out[[taxa]] = list(
    kernel_scale=h_aic, 
    lambda=l_aic, 
    b1=b1_aic,
    sig2=sig2_aic,
    band=band_df,
    csda=csda_df
  )
  # ------------------------------------------------------------------------------

}

save(out, file=paste0(FIG_PATH, "results.RData"))



FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_Genus/"
load(file=paste0(FIG_PATH, "results.RData"))






# ==============================================================================
# Process results
results = data.frame(
  taxa=names(out),
  kernel_scale=sapply(out, function(x) x$kernel_scale),
  lambda=sapply(out, function(x) x$lambda),
  sig2=sapply(out, function(x) x$sig2),
  normb=sapply(out, function(x) sqrt(sum(x$b1^2))),
  n_non_zero=sapply(out, function(x) sum(x$band$excludes_zero)),
  max_prop0=sapply(out, function(x) max(x$band$prop0)),
  n_non_zero_deseq=sapply(out, function(x) sum(x$csda$lda))
)
results %<>% mutate(score=normb/sqrt(sig2))
results %<>% arrange(desc(n_non_zero), desc(score))

non_zero = results %>% filter(n_non_zero > 0 | n_non_zero_deseq > 0)
non_zero_b = sapply(non_zero$taxa, function(x) out[[x]]$b1) %>% t()
colnames(non_zero_b) = t0
non_zero_0 = sapply(non_zero$taxa, function(x) out[[x]]$band$excludes_zero) %>% t() %>% data.frame()
colnames(non_zero_0) = t0
non_zero_0_long = non_zero_0 %>% mutate(taxa=non_zero_0 %>% rownames) %>%
  tidyr::pivot_longer(-taxa, names_to="week", values_to="value")

non_zero_clipped = data.frame(non_zero_b)
non_zero_cliped_long = non_zero_clipped %>% mutate(taxa=non_zero_clipped %>% rownames) %>%
  tidyr::pivot_longer(-taxa, names_to="week", values_to="value")

non_zero_deseq = sapply(non_zero$taxa, function(x) out[[x]]$csda$lda) %>% t() %>% data.frame()
colnames(non_zero_deseq) = t0
non_zero_deseq_long = non_zero_deseq %>% mutate(taxa=non_zero_deseq %>% rownames) %>%
  tidyr::pivot_longer(-taxa, names_to="week", values_to="value")

non_zero_cliped_long$cs = non_zero_deseq_long$value
non_zero_cliped_long$CSDA = ifelse(non_zero_cliped_long$cs, "stripe", "none")
non_zero_cliped_long %<>% mutate(CSDA=ifelse(is.na(CSDA), "none", CSDA))

non_zero_cliped_long$nz = non_zero_0_long$value
non_zero_cliped_long %<>% mutate(nztxt = ifelse(nz, "*", ""))

g = ggplot(
  data=non_zero_cliped_long,
  mapping=aes(x=week, y=taxa)
) + ggpattern::geom_tile_pattern(
  mapping=aes(fill=value, pattern=CSDA),
  pattern_fill="grey80",
  pattern_size=0,
  pattern_density=0.5
) + 
  geom_text(mapping=aes(label=nztxt), size=7, vjust="center", color="black") + 
  scale_fill_gradient2(midpoint=0) + 
  scale_x_discrete(breaks=paste0("X", t0), limits=paste0("X", t0), labels=t0) + 
  xlab("Week") + ylab("") + ggtitle("Locally DA taxa") + labs(fill="KO-WT") + 
  theme(
    text=element_text(family="Roboto"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  scale_y_discrete(limits=rev)

ggsave(paste0(FIG_PATH, "non_zero.png"), g, width=6, height=4, bg="white")


# ------------------------------------------------------------------------------