# this file processes DMBT data (longitudinal microbiome) (GL, 6/8/2022)
setwd('C:\\Users\\ligen\\Dropbox (Personal)\\Research File\\Collaborative Projects\\WithNisha\\DMBT longitudinal Study\\Data')



library(phyloseq)
library(microbiome)
library(ggplot2)


# aggregate_top_taxa <- function (x, top, level) {
#   
#   .Deprecated("aggregate_rare",
#               "The microbiome::aggregate_top_taxa function is deprecated.")
#   
#   x <- aggregate_taxa(x, level)
#   
#   tops <- top_taxa(x, top)
#   tax <- tax_table(x)
#   
#   inds <- which(!rownames(tax) %in% tops)
#   
#   tax[inds, level] <- "Other"
#   
#   tax_table(x) <- tax
#   
#   tt <- tax_table(x)[, level]
#   tax_table(x) <- tax_table(tt)
#   
#   aggregate_taxa(x, level)
#   
# }







#### LOAD DATA ######################

## load OTU table
microbiome=read.table(file='SALIVA_OTU.csv',header = TRUE,sep=',')
dim(microbiome) # 294 samples, 1726 OTU
microbiome[1:5,1:5] # read counts
Micro.mat=microbiome[,-c(1,2,3)]
rownames(Micro.mat)=microbiome$Group
samp.lab=microbiome$Group
# get sample id and visit id
sample.id=gsub(".*?([0-9]+).*","\\1",samp.lab)
length(unique(sample.id)) # 65 unique subjects
visit.id=as.numeric(gsub("^[^_]*[_]?([0-9]*).*",'\\1',samp.lab))
visit.id[is.na(visit.id)]=0
table(visit.id) # all have baseline, all but one have endpoint
#
# load gene type cvrt
gene=read.table(file='GroupGeneTime.csv',header = TRUE,sep=',')
mch=match(as.numeric(sample.id),gene$ID)
gene.type=gene$Type[mch]
# metadata
ID=data.frame(sample=sample.id,visit=visit.id,gene=gene.type)



## load taxonomy
taxonomy=read.table(file='SALIVA_TAXO.csv',header = TRUE,sep=',')
dim(taxonomy) # 1866*8
colnames(taxonomy)
rownames(taxonomy)=taxonomy$OTU
mch=match(colnames(Micro.mat),rownames(taxonomy))
sum(is.na(mch)) # 0
taxo=taxonomy[mch,c(-1,-2)]
# add taxonomic level identifier
taxo$Domain=paste('d__',taxo$Domain,sep='')
taxo$Phylum=paste('p__',taxo$Phylum,sep='')
taxo$Class=paste('c__',taxo$Class,sep='')
taxo$Order=paste('o__',taxo$Order,sep='')
taxo$Family=paste('f__',taxo$Family,sep='')
taxo$Genus=paste('g__',taxo$Genus,sep='')
# replace "unclassified" by NA
for(j in 1:ncol(taxo)){
  taxo[grepl('unclassified',taxo[,j],fixed = TRUE),j]=NA
}
TAX=taxo # 1726*6, a taxonomic tree with various depths






## create phyloseq data structure
rownames(ID)=rownames(Micro.mat)
OTU = otu_table(Micro.mat, taxa_are_rows = FALSE) 
TAX = tax_table(as.matrix(TAX))
samples = sample_data(ID)

mydata <- phyloseq(OTU, TAX, samples)
mydata 




## save data
save(Micro.mat, TAX, ID, mydata, file="DMBT_processed_data.Rdata")













#### aggregate data to different levels
##################################################################################################
load("DMBT_processed_data.Rdata")
 

mydata 


## aggregate to finest known taxonomy for each taxon 
# (taxa may be at different level, but the tree is kept at the deepest equal level)
TAX=tax_table(mydata)
TAX=as.data.frame(TAX)
UniqueTax=apply(format(TAX), 1, paste, collapse="") # unique identifier (by concatenating names)
length(unique(UniqueTax)) # p=262 unique leaf nodes in the tax tree
ASV=rownames(TAX) 
full_TAX=as.matrix(cbind(TAX,UniqueTax,ASV))
sum(is.na(full_TAX)) #2603, indeed with NA in the matrix
 
mydata.temp=phyloseq(otu_table(mydata),tax_table(full_TAX),sample_data(mydata))
mydata.temp
mydata.full= aggregate_taxa(mydata.temp, "UniqueTax", verbose = FALSE)
mydata.full # 262 taxa, 294 samples
sum(is.na(tax_table(mydata.full))) #0, aggregate_taxa always replace NA by "unknown", and thus create a equal-depth full tree (which is good)
# check
sum(otu_table(mydata))
sum(otu_table(mydata.full)) # do not lose any OTU data




#### aggregate to genus level
taxlevel="Genus" # desired level of data aggregation 
aggreg.data= aggregate_taxa(mydata, taxlevel, verbose = FALSE)
aggreg.data   # 294 samples, 169 taxa at genus level (incl last row: unknown)

# remove last row in tax_table and otu_table (i.e. data lost for those OTUs with unknown genus names)
new.TAX=tax_table(aggreg.data)
new.TAX=new.TAX[-nrow(new.TAX),]
new.OTU=otu_table(aggreg.data)
dim(new.OTU) # 169*294
new.OTU=new.OTU[-nrow(new.OTU),]

# get usable phyloseq data
mydata.aggreg=phyloseq(otu_table(new.OTU),tax_table(new.TAX),sample_data(mydata))
mydata.aggreg # 168 taxa, 294 samples


save(mydata, mydata.full, mydata.aggreg, file="DMBT_processed_data1.Rdata")



# #### aggregate to Family level
# taxlevel="Family" # desired level of data aggregation 
# aggreg.data= aggregate_taxa(mydata, taxlevel, verbose = FALSE)
# aggreg.data   # 294 samples, 169 taxa at genus level (incl last row: unknown)
# 
# # remove last row in tax_table and otu_table (i.e. data lost for those OTUs with unknown genus names)
# new.TAX=tax_table(aggreg.data)
# new.TAX=new.TAX[-nrow(new.TAX),]
# new.OTU=otu_table(aggreg.data)
# dim(new.OTU) # 169*294
# new.OTU=new.OTU[-nrow(new.OTU),]
# 
# # get usable phyloseq data
# mydata.aggreg=phyloseq(otu_table(new.OTU),tax_table(new.TAX),sample_data(mydata))
# mydata.aggreg # 168 taxa, 294 samples
# 
# 
# save(mydata, mydata.full, mydata.aggreg, file="DMBT_processed_data2.Rdata")




##### 7/25/2022 GL: process data to get tensor data for Eric Lock
load("DMBT_processed_data1.Rdata")
mydata.aggreg # genus-level data
temp=otu_table(mydata.aggreg)
currX=t(temp)
currID=sample_data(mydata.aggreg)
samplename=unique(currID$sample)
visitname=sort(unique(currID$visit))
# create tensor (sample*taxon*visit)
tensorX=array(NA,dim=c(length(samplename),ncol(currX),length(visitname) ), 
              dimnames=list(samplename,colnames(currX),visitname) )
for(i in 1:nrow(currX)){
  mch1=match(currID$sample[i],samplename)
  mch2=match(currID$visit[i],visitname)
  tensorX[mch1,,mch2]=currX[i,]
}
# check missing proportion
sum(is.na(tensorX))/prod(dim(tensorX))  # 24.6% block-wise missing
# CLR transformation per microbiome vector
clrX=tensorX
for(i in 1:length(samplename)){
  for(j in 1:length(visitname)){
    if(all(!is.na(tensorX[i,,j]))){
      temp=tensorX[i,,j]
      temp[temp==0]=1
      clrX[i,,j]=scale(log(temp),center=TRUE,scale=FALSE)
    }
  }
}
#check
max(clrX,na.rm=TRUE)
min(clrX,na.rm=TRUE)
save(tensorX,clrX,file='BWmissing_tensor.Rdata') 




#### Preliminary analysis of Microbiome data
##################################################################################################
load("DMBT_processed_data1.Rdata")
curr.data=mydata.aggreg   # genus level data
otu=t(otu_table(curr.data)) # n*p
# check meta data
ID=sample_data(curr.data)
length(unique(ID$sample)) # 65 unique mice


pdf(paste0('Meta_data_info.pdf'), width = 10, height = 5)
par(mfrow=c(2,2))
out=table(ID$sample)
hist(out,xlab='# of visits',main=paste0('Histogram of repeated measurements (among ', length(unique(ID$sample)),' mice)'))
out=table(ID$visit)
barplot(out~as.numeric(names(out)),xlab='week',main='# of mice with measurements in each week')

# check library size
hist(rowSums(otu),xlab='library size',main=paste0('Histogram of library sizes (among ',nrow(otu),' samples)')) # library sizes for 294 samples

# check zeros
sort(colSums(otu==0)) # each OTU has as least one non-zero
hist(colSums(otu==0)/nrow(otu),20,xlab='zero percentage per taxon',main=paste0('Histogram of zero percentages (among ',ncol(otu),' taxa)'))
dev.off()




### convert reads to relative abundances
relabn=otu/rowSums(otu)
curr.data1=phyloseq(otu_table(relabn,taxa_are_rows = FALSE),tax_table(curr.data),sample_data(curr.data))





############### VISUALIZATION ##########################
## longitudinal composition plot

# # get top 10 taxa at genus level
# curr.data2=aggregate_top_taxa(curr.data1, 10, 'Genus')
# # plot
# plot_composition(curr.data2) + theme(legend.position = "bottom") +
#   scale_fill_brewer("Genus", palette = "Paired") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle("Relative abundance") +  theme(legend.title = element_text(size = 18))


## KO group
KO.data <- subset_samples(curr.data1, gene == "KO")
# merge all taxa that are detected rare
KO.data1 <- aggregate_rare(KO.data, level="Genus", detection = 0.1/100, prevalence = .31) # detection=P/A threshold, prevalence=non-zero percentage cutoff
KO.data1 # aggregate to 11 taxa (including one other)
# Remove the "g__" before the taxa names
taxa_names(KO.data1) <- gsub("g__", "", taxa_names(KO.data1) )

# group by subj
plot_composition(KO.data1, 
                          sample.sort = 'visit',
                          otu.sort = NULL,
                          plot.type = "barplot",
                          group_by='sample',
                          verbose = FALSE) +
  labs(x = "Sample",
       y = "Relative abundance",
       title = "Longitudinal Compositional Plot",
       subtitle = 'KO Mice') +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer("Genus", palette = "Paired")
ggsave(paste0('BarPlot_KO_samp_genus.pdf'),width=24, height=8) 

# group by week
plot_composition(KO.data1, 
                 sample.sort = 'sample',
                 otu.sort = NULL,
                 plot.type = "barplot",
                 group_by='visit',
                 verbose = FALSE) +
  labs(x = "Week",
       y = "Relative abundance",
       title = "Longitudinal Compositional Plot",
       subtitle = 'KO Mice') +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer("Genus", palette = "Paired")
ggsave(paste0('BarPlot_KO_visit_genus.pdf'),width=24, height=8) 



## WT group
WT.data <- subset_samples(curr.data1, gene == "WT")
# merge all taxa that are detected rare
WT.data1 <- aggregate_rare(WT.data, level="Genus", detection = 0.1/100, prevalence = .3) # detection=P/A threshold, prevalence=non-zero percentage cutoff
WT.data1 # aggregate to 10 taxa (including one other)
# Remove the "g__" before the taxa names
taxa_names(WT.data1) <- gsub("g__", "", taxa_names(WT.data1) )

# group by subj
plot_composition(WT.data1, 
                 sample.sort = 'visit',
                 otu.sort = NULL,
                 plot.type = "barplot",
                 group_by='sample',
                 verbose = FALSE) +
  labs(x = "Sample",
       y = "Relative abundance",
       title = "Longitudinal Compositional Plot",
       subtitle = 'WT Mice') +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer("Genus", palette = "Paired")
ggsave(paste0('BarPlot_WT_samp_genus.pdf'),width=24, height=8) 

# group by week
plot_composition(WT.data1, 
                 sample.sort = 'sample',
                 otu.sort = NULL,
                 plot.type = "barplot",
                 group_by='visit',
                 verbose = FALSE) +
  labs(x = "Week",
       y = "Relative abundance",
       title = "Longitudinal Compositional Plot",
       subtitle = 'WT Mice') +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer("Genus", palette = "Paired")
ggsave(paste0('BarPlot_WT_visit_genus.pdf'),width=24, height=8) 



## together, average across samples within each genotype group
#
curr.data2 <- aggregate_rare(curr.data1, level="Genus", detection = 0.1/100, prevalence = .3) # detection=P/A threshold, prevalence=non-zero percentage cutoff
curr.data2 # aggregate to 10 taxa (including one other)
# Remove the "g__" before the taxa names
taxa_names(curr.data2) <- gsub("g__", "", taxa_names(curr.data2) )
# create a index column for sample average
sample_data(curr.data2)$avg.ind = paste0(sample_data(curr.data2)$visit,sample_data(curr.data2)$gene)
# average otu
avg.otu=aggregate(t(otu_table(curr.data2)), list(sample_data(curr.data2)$avg.ind), mean)
rownames(avg.otu)=avg.otu[,1]
temp=avg.otu[,1]
avg.meta=data.frame(visit=as.numeric(gsub("([0-9]+).*$", "\\1",temp)),gene=gsub('[[:digit:]]+', '', temp))
rownames(avg.meta)=temp
curr.data3=phyloseq(otu_table(avg.otu[,-1],taxa_are_rows = FALSE), tax_table(curr.data2),sample_data(avg.meta))


plot_composition(curr.data3, 
                 sample.sort = 'gene',
                 x.label = 'gene',
                 group_by='visit',
                 otu.sort = NULL,
                 plot.type = "barplot",
                 verbose = FALSE) +
  labs(y = "Relative abundance", x='') + 
      #x = "Week",
      #title = "Longitudinal Compositional Plot") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.ticks.x=element_blank(),
        text = element_text(size = 36)) +
  scale_fill_brewer("Genus", palette = "Paired")
#ggsave(paste0('BarPlot_All_avg_genus.pdf'),width=24, height=8) 
ggsave(paste0('OSCC_TDA.pdf'),width=24, height=10) 




##################### Ordination Analysis ########################
### PCoA plot (of all samples) 
# color by visit; shape by gene type (pool all mice)
CF.ord <- ordinate(curr.data1, "PCoA", "bray")
shapes=sample_data(curr.data1)$gene 
shapevalue=c(17,19)
plot_ordination(curr.data1, CF.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(color=sample_data(curr.data1)$visit, shape=shapes),size=2) + 
  scale_shape_manual(values = shapevalue) + 
  ggtitle(paste("Genus-Level Ordination Plot (Bray-Curtis) for All Samples (Colored by Visit)"))
ggsave(paste0('PCoA_All_genus.pdf'),width=8, height=8) 


## PCoA (of KO/WT at T0)
subset1=subset_samples(curr.data1,visit==0) 
CF.ord <- ordinate(subset1, "PCoA", "bray")
shapes=sample_data(subset1)$gene 
shapevalue=c(17,19)
plot_ordination(subset1, CF.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(shape=shapes),size=3) + 
  scale_shape_manual(values = shapevalue) + 
  ggtitle(paste("Genus-Level Ordination Plot (Bray-Curtis) for T0 Samples"))
ggsave(paste0('PCoA_T0_genus.pdf'),width=8, height=8) 


## PCoA (of KO/WT at T22)
subset2=subset_samples(curr.data1,visit==22) 
CF.ord <- ordinate(subset2, "PCoA", "bray")
shapes=sample_data(subset2)$gene 
shapevalue=c(17,19)
plot_ordination(subset2, CF.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(shape=shapes),size=3) + 
  scale_shape_manual(values = shapevalue) + 
  ggtitle(paste("Genus-Level Ordination Plot (Bray-Curtis) for T22 Samples"))
ggsave(paste0('PCoA_T22_genus.pdf'),width=8, height=8) 



## PCoA (of KO)
subset3=subset_samples(curr.data1,gene=='KO') 
CF.ord <- ordinate(subset3, "PCoA", "bray")
plot_ordination(subset3, CF.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(color=sample_data(subset3)$visit),shape=17,size=2) + 
  ggtitle(paste("Genus-Level Ordination Plot (Bray-Curtis) for KO Samples (Colored by Visit)"))
ggsave(paste0('PCoA_KO_genus.pdf'),width=8, height=8) 


## PCoA (of WT)
subset4=subset_samples(curr.data1,gene=='WT') 
CF.ord <- ordinate(subset4, "PCoA", "bray")
plot_ordination(subset4, CF.ord, type="samples") +
  geom_point(size=2,color='white',shape=16)+ # mask up default points from ordination plot
  geom_point(aes(color=sample_data(subset4)$visit),shape=19,size=2) + 
  ggtitle(paste("Genus-Level Ordination Plot (Bray-Curtis) for WT Samples (Colored by Visit)"))
ggsave(paste0('PCoA_WT_genus.pdf'),width=8, height=8) 



################### FORMAL TESTING ##################################
## PERMANOVA (entire profile)
library(vegan)
# T0: KO vs WT
adonis(distance(subset1, method="bray") ~ sample_data(subset1)$gene)  # p=0.712, no significance between KO and WT
# T22: KO vs WT
adonis(distance(subset2, method="bray") ~ sample_data(subset2)$gene)  # p=0.632, no significance between KO and WT



### Differential abundance analysis (individual taxon)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")
library(ANCOMBC)
# T0: KO vs WT
out = ancombc(phyloseq = subset1,formula='gene', p_adj_method = "holm", 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)
res = out$res
res$p_val # raw pvalue
sum(res$diff_abn)  # 0 DA taxa (out of 32)



# T22 KO vs WT
out = ancombc(phyloseq = subset2,formula='gene', p_adj_method = "holm", 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)
res = out$res
res$p_val # raw pvalue
sum(res$diff_abn)  # 0 DA taxa (out of 30)
