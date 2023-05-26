library(tidyverse)
library(magrittr)
library(VCMM)
theme_set(theme_bw())
extrafont::font_import(paths = "~/.local/share/fonts/", prompt=F)



# ==============================================================================
# Settings
FIG_PATH = "/home/simon/Documents/VCM-Experiments/figures/DMBT_Otu/"
RESULT_PATH = "/home/simon/Documents/VCM-Experiments/results/DMBT_Otu/"
results = read_csv(paste0(RESULT_PATH, "coefs.csv"))
results %<>% mutate(
  masked_estimate=factor(
    ifelse(excludes_zero, sign(estimate), "NA"),
    levels=c("-1", "1", "NA"), labels=c("Decrease", "Increase", "Non sig.")
  ),
  display=paste0(otu, "\n(", taxonomy, ")")
)

coefs_nz = results %>% filter(grouping == "diffs") %>% 
  group_by(otu) %>% summarize(nz=any(excludes_zero)) %>% filter(nz) %>% pull(otu)
results_coef_nz = results %>% filter(otu %in% coefs_nz)

ordering = results_coef_nz %>% select(otu, taxonomy, display) %>% group_by(otu) %>% 
  summarize(otu=head(otu, 1), taxonomy=head(taxonomy, 1), display=head(display, 1)) %>%
  arrange(taxonomy, otu) %>% pull(display) %>% rev()

results_coef_nz %<>% mutate(display=factor(display, levels=ordering, labels=ordering, ordered=T))
# ------------------------------------------------------------------------------




# ==============================================================================
# First test

gko = ggplot(
  data=results_coef_nz %>% filter(grouping == "coefs", varname=="KO"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("KO") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    legend.position="none"
  ) + coord_fixed()

gscc = ggplot(
  data=results_coef_nz %>% filter(grouping == "coefs", varname=="SCC"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("SCC") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    legend.position="none",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )  + coord_fixed()

gkoscc = ggplot(
  data=results_coef_nz %>% filter(grouping == "coefs", varname=="KO:SCC"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("KO:SCC") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  ) + guides(alpha="none") + coord_fixed()

g = cowplot::plot_grid(gko, gscc, gkoscc, ncol=3, align="v", axis="tb", rel_widths=c(1.68, 0.952, 1.4))

ggsave(paste0(FIG_PATH, "coefs.png"), g, width=10, height=12, bg="white")
# ------------------------------------------------------------------------------





# ==============================================================================
# First test
g1 = ggplot(
  data=results_coef_nz %>% filter(grouping == "diffs", diff=="WT: HP/CIS->SCC"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("WT: HP/CIS->SCC") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    legend.position="none"
  ) + coord_fixed()

g2 = ggplot(
  data=results_coef_nz %>% filter(grouping == "diffs", diff=="KO: HP/CIS->SCC"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("KO: HP/CIS->SCC") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    legend.position="none",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )  + coord_fixed()

g3 = ggplot(
  data=results_coef_nz %>% filter(grouping == "diffs", diff=="HP/CIS: WT->KO"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("HP/CIS: WT->KO") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    legend.position="none",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )  + coord_fixed()

g4 = ggplot(
  data=results_coef_nz %>% filter(grouping == "diffs", diff=="SCC: WT->KO"),
  mapping=aes(x=as.factor(estimated_time), y=display, fill=masked_estimate, alpha=masked_estimate)
) + geom_tile() + 
  xlab("Week") + ggtitle("SCC: WT->KO") + labs(fill="Direction") + ylab("") + 
  scale_fill_manual(values=c("red", "blue", "white")) + 
  scale_alpha_manual(values=c(1, 1, 0)) + 
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank()
  )  + coord_fixed()  + guides(alpha="none")

g = cowplot::plot_grid(g1, g2, g3, g4, ncol=4, align="v", axis="tb", rel_widths=c(2, 1, 1, 1.5))

ggsave(paste0(FIG_PATH, "diffs.png"), g, width=12, height=12, bg="white")
# ------------------------------------------------------------------------------


