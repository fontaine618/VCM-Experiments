library(tidyverse)
library(magrittr)

DATA_PATH = "/home/simon/Documents/VCM-Experiments/data/DMBT/"
GroupGeneTime_PATH = paste0(DATA_PATH, "GroupGeneTime.csv")
SALIVA_OTU_PATH = paste0(DATA_PATH, "SALIVA_OTU.csv")
SALIVA_TAXO_PATH = paste0(DATA_PATH, "SALIVA_TAXO.csv")
DEMO_PATH = paste0(DATA_PATH, "demo.csv")
META_PATH = paste0(DATA_PATH, "meta.csv")
OUT_PATH = "/home/simon/Documents/VCM-Experiments/data/DMBT_Processed/"



# ==============================================================================
# LOAD OTU TABLE & SAMPLE METADATA
otu = read.table(file=SALIVA_OTU_PATH, header = TRUE, sep=',')
# Format:
# - label (always 0.03)
# - Group (sample descriptor)
# - numOtus (always 1726)
# - OtuXXXX, between 0001 and 1866 (only 1726)
group = otu$Group
subject_id = as.numeric(gsub(".*?([0-9]+).*","\\1", group))
visit_id = as.numeric(gsub("^[^_]*[_]?([0-9]*).*",'\\1', group))
visit_id[is.na(visit_id)] = 0
otu = otu[, -c(1, 2, 3)]  # drop non-otu columns
meta_sample = data.frame(subject_id=subject_id, visit_id=visit_id)
rm(subject_id, visit_id, group)
# ------------------------------------------------------------------------------



# ==============================================================================
# LOAD SUBJECT METADATA
meta_subject = read.table(file=GroupGeneTime_PATH, header = TRUE, sep=',')
# Format:
# - ID
# - Type (WT, KO)
demo = read.table(file=DEMO_PATH, header=T, sep=",")
# Format:
# - ID
# - Gender
# - Group (W, K)
# - Diagnosis (1, 2, 3)
meta_subject %<>% left_join(demo %>% select(ID, Gender, Diagnosis), by="ID")
# ------------------------------------------------------------------------------



# ==============================================================================
# MERGE
meta = meta_sample %>% left_join(meta_subject, by=c("subject_id"="ID"))
# Format:
# - subject_id
# - visit_id
# - Type
# - Gender (M/F)
# - Diagnosis (1/2/3)
# ------------------------------------------------------------------------------



# ==============================================================================
# UPDATED DIAGNOSIS
meta2 = read.table(file=META_PATH, header = TRUE, sep=',')
meta2 %<>% mutate(
  subject_id = as.numeric(gsub(".*?([0-9]+).*","\\1", Group))
)
meta2 %>% group_by(subject_id) %>% summarise(
  CIS=sum(Diagnosis2=="CIS"),
  SCC=sum(Diagnosis2=="SCC"),
  HP=sum(Diagnosis2=="hyperplasia"),
) %>% select(-subject_id) %>% apply(1, function(x) sum(x>0)) # looks fine!
meta2 %<>% group_by(subject_id) %>% summarize(Diagnosis2=head(Diagnosis2, 1))
meta = meta %>% left_join(meta2, by="subject_id")

with(meta, table(Diagnosis, Diagnosis2)) # turns out, we had the same
# ------------------------------------------------------------------------------



# ==============================================================================
# LOAD TAXONOMY
taxonomy = read.table(SALIVA_TAXO_PATH, header=TRUE, sep=",")
# Format:
# - OTU
# - Size
# - Domain, Phylum, Class, Order, Family, Genus
otus = colnames(otu)
rownames(taxonomy) = taxonomy$OTU
taxonomy %<>% filter(OTU %in% otus)
taxonomy = taxonomy[match(otus, taxonomy$OTU), ]
taxonomy %<>% mutate(
  Domain = paste0("d__", Domain),
  Phylum = paste0("p__", Phylum),
  Class = paste0("c__", Class),
  Order = paste0("o__", Order),
  Family = paste0("f__", Family),
  Genus = paste0("g__", Genus)
)
taxonomy %<>% select(-c(Size, OTU))
for(j in 1:ncol(taxonomy)){
  taxonomy[grepl('unclassified', taxonomy[, j], fixed = TRUE), j] = NA
}
# ------------------------------------------------------------------------------



# ==============================================================================
# TO PHYLOSEQ
otu = phyloseq::otu_table(otu, taxa_are_rows=F)
tax = phyloseq::tax_table(as.matrix(taxonomy))
sam = phyloseq::sample_data(meta)
pseq = phyloseq::phyloseq(otu, tax, sam)
save(pseq, file=paste0(OUT_PATH, "pseq.RData"))
# ------------------------------------------------------------------------------



# ==============================================================================
# AGGREGATE TO FINEST KNOWN TAXONOMY
tax = phyloseq::tax_table(pseq) %>% as.data.frame()
tax$unique_name = apply(format(tax), 1, paste, collapse=":")
tax = phyloseq::tax_table(as.matrix(tax))
pseq_unique = phyloseq::phyloseq(otu, tax, sam)
pseq_unique = microbiome::aggregate_taxa(pseq_unique, "unique_name")
save(pseq_unique, file=paste0(OUT_PATH, "pseq_unique.RData"))
# ------------------------------------------------------------------------------



# ==============================================================================
# AGGREGATE TO HIGHER LEVELS
pseq_genus = microbiome::aggregate_taxa(pseq, "Genus")
save(pseq_genus, file=paste0(OUT_PATH, "pseq_genus.RData"))
pseq_family = microbiome::aggregate_taxa(pseq, "Family")
save(pseq_family, file=paste0(OUT_PATH, "pseq_family.RData"))
pseq_order = microbiome::aggregate_taxa(pseq, "Order")
save(pseq_order, file=paste0(OUT_PATH, "pseq_order.RData"))
# ------------------------------------------------------------------------------