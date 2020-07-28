---
title: "Ulbrich et al. 2020 SG cultivar microbiome OTU filtering & normalization"
author: "Tayler Ulbrich"
date: "7/1/2020"
---
*OTU filtering notes:*
Three microbiome files are processed in this script, each of which were filtered for non-bacterial or non-fungal reads, singeltons (<2 reads), and rarefied.
- bacterial root subset of 4 cultivars (endo) rarefied to 1984 reads (dropped from downstream analysis; used combined root & soil dataset instead) 
- bacterial soil & root subset of 4 cultivars (rhizoendo) rarefied to 2026 reads (used for downstream analysis)
- bacterial soil (rhizo) for 12 cultivars rarefied to 4694 reads (used for downstream analysis)
- fungal soil (rhizo) for 12 cultivars rarefied to 4153 reads (used for downstream analysis)


# Set working directory and load library
```{r setup, include=FALSE}

# Load Phyloseq: 
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE)
# Install Microbiome 
#library(BiocInstaller)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("microbiome")

library(magrittr)
library(phyloseq)
library(microbiome)
library(scales)
library(grid)
library(ggplot2)
library(vegan)
library(multcompView)
library(multcomp)
library(lme4)
library(lmerTest)
library(ggpubr)
library(car)
library(lsmeans)
library(dplyr)
library(ape)
library(MASS)
library(stats)
library(tidyr)
library(tibble)

theme_set(theme_classic())
```


# BACTERIA OTU FILTERING & NORMALIZATION 
## 1) Make phyloseq object with otu table, taxa file, and tree file 
```{r}

# load metadata file with X.SampleID as the column with the same sample ids in otu table 
metadata <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
SAMP <- sample_data(metadata) # this is the metadata file that you merge with phyloseq
rownames(SAMP) <- SAMP$X.SampleID # row names must match OTU table


# load otu table and otu table 
# Note: Usearch output has a #OTU ID in cell A1 - you need to remove this for it to work.
otufile <- read.table("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/Final_otu_tax_tree_files/TMC_SGVar16S_RhizoEndo_combined_merged_ALLruns_OTU_table.txt", header = TRUE)
otu <- otu_table(otufile, taxa_are_rows = TRUE)

# load taxa file 
    # tax file - from export, separated by ; delimitation, removed taxonomy headings "k:", cleaned to have only taxonomys with 80perc confidence; changed all blanks to "unclassified"

taxafile_otu <- as.matrix(read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Final_otu_tax_tree_files/TMC_SGVar16S_RhizoEndo_combined_merged_ALLruns_otus_taxonomy_SilvaV123_sintax_80perc_editted.csv", header = TRUE, row.names =1))
tax_otu <- tax_table(taxafile_otu)


# load tree file 
# note: I used figtree program to convert the .tre file to a newick file format - required for read_tree
tree_otu_file <- "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Final_otu_tax_tree_files/2018.06.22_SGvarietyRhizoEndo16s_pasta.nwk"
tree_otu <- read_tree(tree_otu_file)
tree_otu # 48759 tips and 48757 internal nodes 
#plot(tree_otu)

# combine data into phyloseq-object (include map, otu, tax, tree)
myphy_otus<- merge_phyloseq(otu, tax_otu, SAMP, tree_otu)
sampledata <- as(sample_data(myphy_otus), "data.frame")

# check the dimensions of OTU 
dim(otu_table(myphy_otus)) #48759, 192 (taxa, samples)

```

## 2) Rename column names (taxa) and clean up OTUS (remove 0 abundance OTUS and remove any non-bacterial reads)
```{r}

#what are the column names of our taxonomy file? 
colnames(tax_table(myphy_otus)) # good to go

# remove otus not in our samples 
otus.p <- prune_taxa(taxa_sums(myphy_otus) > 0, myphy_otus)
ntaxa(myphy_otus) - ntaxa(otus.p) # removed 27787 otus that weren't present in my samples (present in other samples on the same sequencing run)
ntaxa(myphy_otus) #48759 taxa in my samples 


# How many total sequences for endo vs. rhizo vs. rhizoendo? 
otus_endo <- subset_samples(otus.p, sampledata$SampleSite == "Endosphere")
  otus_endo <- prune_taxa(taxa_sums(otus_endo) > 0, otus_endo)
  sum(sample_sums(otus_endo))#1028968
  ntaxa(otus_endo) #3297
otus_rhizo <- subset_samples(otus.p, sampledata$SampleSite == "Rhizosphere")
  otus_rhizo <- prune_taxa(taxa_sums(otus_rhizo) > 0, otus_rhizo)
  sum(sample_sums(otus_rhizo))#2294871
  ntaxa(otus_rhizo) # 20278
otus_rhizoendo <- otus.p
  sum(sample_sums(otus_rhizoendo))#3323839
  ntaxa(otus_rhizoendo) # 20972


# subset for only bacteria and remove any mitochondria and chloroplast
bact.p <- otus.p %>%
  subset_taxa(
    Kingdom == "Bacteria" & # removed 114 non-bact kingdom 
    Class != "Chloroplast" & # removed 167 chloro
    Family != "Mitochondria" # removed 65 mitochondria
    )


# I wrote the taxaonomy file into csv to open in excel to confirm how many non-bacterial samples should be removed 
#tax <- as(tax_table(otus.p), "matrix")
# write.csv(tax, "C:/Users/Owner/Documents/01_Grad School/06_Sequencing/SGvarietyRhizo/16S/final_files/RhizoEndo/Final_filtered_phyloseq_objects/2018.10.30_80perc_taxtable_ONLYmysamples.csv")
# found 114 non-bact, 167 chloroplast,65 mitochondria - 

# how many non-bact otus were removed?
ntaxa(otus.p)- ntaxa(bact.p) # removed 346 non-bacterial otus in whole dataset -- this is what I found in my excel file 

# What % of sequences were bacterial?
sum(otu_table(bact.p))/sum(otu_table(otus.p))*100 # 80.93274


# subset by rhizosphere and endosphere to see where the non-bacterial samples are coming from 
# remove otus not present in the rhizosphere, endopshere 
otus_rhizo <- subset_samples(myphy_otus, sampledata$SampleSite == "Rhizosphere")
dim(otu_table(otus_rhizo)) # 48759 144
rhizo.p<-prune_taxa(taxa_sums(otus_rhizo) > 0, otus_rhizo)
dim(otu_table(rhizo.p)) #20278 otus not in rhizo removed 
rhizo.p.bact <- rhizo.p %>%
  subset_taxa(
    Kingdom == "Bacteria" &# removed 103 non-bact kingdom 
    Class != "Chloroplast" & #removed 107 chloro
    Family != "Mitochondria" # removed 56 mito
    )

dim(otu_table(rhizo.p.bact)) #20012 144 
ntaxa(rhizo.p) - ntaxa(rhizo.p.bact) #226 non-bact/chloro/mito samples removed from rhizosphere 

#endosphere
otus_endo <- subset_samples(myphy_otus, sampledata$SampleSite == "Endosphere")
dim(otu_table(otus_endo)) #48759 48

endo.p<-prune_taxa(taxa_sums(otus_endo) > 0, otus_endo)
dim(otu_table(endo.p)) #3297 otus not in endosphere removed
endo.p.bact <- endo.p %>%
  subset_taxa(
   Kingdom == "Bacteria" & # removed 14 non-bact kingdom 
   Class != "Chloroplast" & #removed 90 chloro
   Family != "Mitochondria" # removed 24 mito
    )

dim(otu_table(endo.p.bact)) #3169 48
ntaxa(endo.p) - ntaxa(endo.p.bact) #128 non-bact/chloro/mito samples removed from endosphere 


# where did the 346 non-bacterial samples we removed come from? 
ntaxa(rhizo.p) - ntaxa(rhizo.p.bact) # removed 266 non-bacterial otus in rhizosphere
ntaxa(endo.p) - ntaxa(endo.p.bact) # removed 128 otus non-bacterial otus in endosphere


# % read number remaining after removing non-bacterial/chloro/mito samples for rhizo?
# what % of total rhizo reads were bacterial 
sum(otu_table(rhizo.p.bact))/sum(otu_table(rhizo.p))*100 # 99% of sample reads maintained - only lost <1% of reads 

# More specifically, how many reads were lost to non-bacterial samples in rhizo? 
sum(otu_table(rhizo.p)) - sum(otu_table(rhizo.p.bact)) # 18098 reads lost in rhizosphere samples

# # % read number remaining after removing removing non-bacterial/chloro/mito samples for endo?
sum(otu_table(endo.p.bact))/sum(otu_table(endo.p))*100 # = 40.17% of the endosphere reads were bacterial :-/ 
sum(otu_table(endo.p)) - sum(otu_table(endo.p.bact)) # 615667 reads lost in endosphere samples 

# NOTE: files to continue with: bact.p; endo.p.bact; rhizo.p.bact (these were filtered for only otus in our samples and all non-bacterial reads/otus were removed)

```


## 3) Remove Singeltons
```{r}

# Remove singletons 
endo.pp.bact <- prune_taxa(taxa_sums(endo.p.bact) > 1, endo.p.bact)
ntaxa(endo.p.bact)- ntaxa(endo.pp.bact) # removed 499 endo singletons
ntaxa(endo.pp.bact) #2670

rhizo.pp.bact <- prune_taxa(taxa_sums(rhizo.p.bact) > 1, rhizo.p.bact)
ntaxa(rhizo.p.bact)- ntaxa(rhizo.pp.bact) # removed 2134 rhizo singletons
ntaxa(rhizo.pp.bact) #17878

bact.pp <- prune_taxa(taxa_sums(bact.p) > 1, bact.p)
ntaxa(bact.p)- ntaxa(bact.pp) # removed 2091 total singletons
ntaxa(bact.pp) #18535


``` 

## 4) Remove low sequence coverage samples
```{r}
# all samples 
sort(sample_sums(bact.pp))  
bact.pp_pruned <- prune_samples(sample_sums(bact.pp) > 2000, bact.pp) # set at 2000
sort(sample_sums(bact.pp_pruned)) 
# OTUS: removed samples G2, H2, SGendoG2, G11, SGendoJ2, SGendoK11, SGendoJ7, SGendoJ3, H5, SGendoG10
nsamples(bact.pp_pruned) # 182
mean(sample_sums(bact.pp_pruned)) # average read # = 14726.79
median(sample_sums(bact.pp_pruned)) # median read # = 14698
max(sample_sums(bact.pp_pruned)) #55045 
min(sample_sums(bact.pp_pruned)) #2026
sum(sample_sums(bact.pp_pruned)) #2680275
ntaxa(bact.pp_pruned) #18535

## rhizosphere 
# remove low sequence coverage samples
sort(sample_sums(rhizo.pp.bact))  
rhizo.pp_pruned <- prune_samples(sample_sums(rhizo.pp.bact) > 4600, rhizo.pp.bact) # set at 4600 
sort(sample_sums(rhizo.pp_pruned)) 
# OTUS: removed samples G2, H2, H5, G11, E9, E10
nsamples(rhizo.pp_pruned) # 138
# confirm average/median/min read # after pruning
mean(sample_sums(rhizo.pp_pruned)) # average read # = 16430.12
median(sample_sums(rhizo.pp_pruned)) # median read # = 16348.5
min(sample_sums(rhizo.pp_pruned)) # min read # = 4694
max(sample_sums(rhizo.pp_pruned))
sum(sample_sums(rhizo.pp_pruned)) # 2267356
ntaxa(rhizo.pp_pruned) #17878

## endosphere 
# remove low sequence coverage samples
sort(sample_sums(endo.pp.bact))  
endo.pp_pruned <- prune_samples(sample_sums(endo.pp.bact) > 1900, endo.pp.bact) # set at 1900 
sort(sample_sums(endo.pp_pruned)) 
# OTUS: removed samples TMCSGendoG2, TMCSGendoJ2, TMCSGendoK11, TMCSGendoJ7, TMCSGendoJ3
nsamples(endo.pp_pruned) # 43
mean(sample_sums(endo.pp_pruned)) # average read # =  9499.93
median(sample_sums(endo.pp_pruned)) # median read # = 8651
min(sample_sums(endo.pp_pruned)) # min read # = 1984
max(sample_sums(endo.pp_pruned)) #24356
sum(sample_sums(endo.pp_pruned)) # 408497
ntaxa(endo.pp_pruned) #2670

```
### Save raw data 
```{r}
# files with no singletons and low-sequence depth samples removed; no normalization.

#save(bact.pp_pruned, file = "16S/PhyloseqObjects/20200501_SGVar16Sbact_pp_nosingles_raw.Rdata")

#save(rhizo.pp_pruned, file = "16S/PhyloseqObjects/20200501_SGVar16Sbact_RHIZO_pp_nosingles_raw.Rdata")

#save(endo.pp_pruned, file = "16S/PhyloseqObjects/20200501_SGVar16Sbact_ENDO_pp_nosingles_raw.Rdata")


```
## 5) Investigate Sequencing Depth
```{r}

# What does the total reads per sample distribution look like? 
sample_sum_df <- data.frame(sum = sample_sums(bact.pp))

ggplot(sample_sum_df, aes(x = sum)) +   # Histogram of sample read counts 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) + 
  ggtitle("Distribution of sample sequencing depth of otus (bact.pp)") + 
   xlab("Read counts") + 
  theme(axis.title.y = element_blank())
# extreme outlier at 50000 read counts! Average around 1500 reads

# mean, max, min, median of sample read counts to inform normalization techniques 
min(sample_sums(bact.pp)) # 6 reads
mean(sample_sums(bact.pp)) # 13999.91 reads
max(sample_sums(bact.pp)) # 55045 reads
median(sample_sums(bact.pp)) # 13724 reads

################################
# PLOT read # by sample - variation within variety 

sampledata_bact.pp_pruned <- as.data.frame(as(sample_data(bact.pp_pruned),"matrix"))
sampledata_bact.pp_pruned$samplesums <- sample_sums(bact.pp_pruned)


ggplot(sampledata_bact.pp_pruned, aes(x = X.SampleID, y = samplesums, fill = Variety)) +
  geom_bar(stat = "identity", width = 0.85) + 
  #scale_fill_manual(values = speciesPalette) +
  facet_grid(.~SampleSite, drop = TRUE, space = "free", scales = "free") +
  # use to plot multiple panels if you want to look at sample variation within your treatment.
  ylab("Total Read # ") + 
  xlab("Soil Sample Replicates") +
  theme_classic() + 
theme(plot.title = element_text(size = 16), 
  axis.title.x = element_text(size = 40), 
  axis.title.y = element_text(size = 40), 
  axis.text.x = element_blank(), #element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 30, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 30),
  legend.title = element_text(size = 30),
  legend.position="right", # remove legend "none"
  strip.text = element_text(size = 30))  #text in boxes for facet grid


```
####  Rarefaction Curves 
```{r, echo=FALSE}
############
# ALL SAMPLES 
#############

# Rarefaction curve for every sample 
otus.tab=t(as(otu_table(bact.pp_pruned), "matrix"))
otus.map=sample_data(bact.pp_pruned)
row.names(otus.tab)=otus.map$X.SampleID 
raremax <- min(rowSums(otus.tab)) #  2026 - this is the min read number 

rarecurve <- rarecurve(otus.tab, step = 1000, sample = raremax,label = FALSE)

#####################

sampledata_bact.pp_pruned <- as.data.frame(as(sample_data(bact.pp_pruned),"matrix"))
sampledata_bact.pp_pruned$samplesums <- sample_sums(bact.pp_pruned)


Observed <- estimate_richness(bact.pp_pruned, measures = "Observed")
library(tibble)
Observed <- rownames_to_column(Observed, "X.SampleID")
sampledata_bact.pp_pruned <- merge(sampledata_bact.pp_pruned, Observed, by = "X.SampleID" )

# find min read # that you could rarefy to 
raremax = min(sampledata_bact.pp_pruned$samplesums) 

# Plot rarefcation curve, richness * read # 
require(ggplot2)
ggplot(sampledata_bact.pp_pruned, aes(x= samplesums, y = Observed, color = SampleSite)) + 
    geom_point(size = 4) + 
  ggtitle("Rarefaction curve for all bacterial reads") + 
  labs(x = "Total Reads", y = "Observed Richness")+
 # geom_text(label = sampledata_bact.pp_pruned$X.SampleID) +
  geom_vline(xintercept = raremax, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #


# Compare rarefied richness to observed richness 

bact.pp_pruned.otu =t(as(otu_table(bact.pp_pruned), "matrix"))
row.names(bact.pp_pruned.otu)=sampledata_bact.pp_pruned$X.SampleID 

richness_rarefied <- rarefy(bact.pp_pruned.otu, raremax)
bact.spp = specnumber(bact.pp_pruned.otu) #specnumber finds frequences of species

rare_plot <- plot(bact.spp, richness_rarefied, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "bact.pp_pruned, rarefied to 2026")
abline(0,1) # this line should nicely match the points

# Or plot with ggplot
bact.richness.comp <- merge(richness_rarefied, bact.spp, by = 0)
colnames(bact.richness.comp) <- c("X.SampleID","richness_rarefied", "richness_raw")

require(ggplot2)
ggplot(bact.richness.comp, aes(x= richness_rarefied, y = richness_raw)) + 
    geom_point(size = 4) + 
  ggtitle("") + 
  labs(x = "richness_rarefied", y = "richness_raw")+
 # geom_text(label = sampledata_bact.pp_pruned$X.SampleID) +
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #

##################
# Root and Soil 4 cultivar subset 
##################


bact.pp_pruned_4sub <- subset_samples(bact.pp_pruned, sample_data(bact.pp_pruned)$Variety == "Alamo" |  sample_data(bact.pp_pruned)$Variety == "Kanlow" | sample_data(bact.pp_pruned)$Variety == "Cave-in-Rock" | sample_data(bact.pp_pruned)$Variety == "Southlow")
nsamples(bact.pp_pruned_4sub) # 88

sampledata_bact.pp_pruned_4sub <- as.data.frame(as(sample_data(bact.pp_pruned_4sub),"matrix"))
sampledata_bact.pp_pruned_4sub$samplesums <- sample_sums(bact.pp_pruned_4sub)


# Rarefaction curve for every sample 
otus.tab=t(as(otu_table(bact.pp_pruned_4sub), "matrix"))
otus.map=sample_data(bact.pp_pruned_4sub)
row.names(otus.tab)=otus.map$X.SampleID 
raremax <- min(rowSums(otus.tab)) #  2110 - this is the min read number 

rarecurve <- rarecurve(otus.tab, step = 1000, sample = raremax,label = FALSE)


#####################

Observed <- estimate_richness(bact.pp_pruned_4sub, measures = "Observed")
library(tibble)
Observed <- rownames_to_column(Observed, "X.SampleID")
sampledata_bact.pp_pruned_4sub <- merge(sampledata_bact.pp_pruned_4sub, Observed, by = "X.SampleID" )

# find min read # that you could rarefy to 
raremax = min(sampledata_bact.pp_pruned_4sub$samplesums) 
raremax

# Plot rarefcation curve, richness * read # 
require(ggplot2)
ggplot(sampledata_bact.pp_pruned_4sub, aes(x= samplesums, y = Observed, color = SampleSite)) + 
    geom_point(size = 4) + 
  ggtitle("Rarefaction curve for root & soil of 4 cultivars") + 
  labs(x = "Total Reads", y = "Observed Richness")+
 # geom_text(label = sampledata_bact.pp_pruned_4sub$X.SampleID) +
  geom_vline(xintercept = raremax, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #


# Compare rarefied richness to observed richness 

bact.pp_pruned_4sub.otu =t(as(otu_table(bact.pp_pruned_4sub), "matrix"))
row.names(bact.pp_pruned_4sub.otu)=sampledata_bact.pp_pruned_4sub$X.SampleID 

richness_rarefied <- rarefy(bact.pp_pruned_4sub.otu, raremax)
bact.spp = specnumber(bact.pp_pruned_4sub.otu) #specnumber finds frequences of species

rare_plot <- plot(bact.spp, richness_rarefied, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "bact.pp_pruned_4sub, rarefied to 2026")
abline(0,1) # this line should nicely match the points

# Or plot with ggplot
bact.richness.comp <- merge(richness_rarefied, bact.spp, by = 0)
colnames(bact.richness.comp) <- c("X.SampleID","richness_rarefied", "richness_raw")

require(ggplot2)
ggplot(bact.richness.comp, aes(x= richness_rarefied, y = richness_raw)) + 
    geom_point(size = 4) + 
  ggtitle("") + 
  labs(x = "richness_rarefied", y = "richness_raw")+
 # geom_text(label = sampledata_bact.pp_pruned_4sub$X.SampleID) +
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #



#################
## Rhizosphere 
################

# Rarefaction curve for every sample 
otus.tab=t(as(otu_table(rhizo.pp_pruned), "matrix"))
otus.map=sample_data(rhizo.pp_pruned)
row.names(otus.tab)=otus.map$X.SampleID 
raremax <- min(rowSums(otus.tab)) #  
raremax # 4964
rarecurve <- rarecurve(otus.tab, step = 1000, sample = raremax,label = FALSE)


###################

sampledata_rhizo.pp_pruned <- as.data.frame(as(sample_data(rhizo.pp_pruned),"matrix"))
sampledata_rhizo.pp_pruned$samplesums <- sample_sums(rhizo.pp_pruned)


Observed <- estimate_richness(rhizo.pp_pruned, measures = "Observed")
library(tibble)
Observed <- rownames_to_column(Observed, "X.SampleID")
sampledata_rhizo.pp_pruned <- merge(sampledata_rhizo.pp_pruned, Observed, by = "X.SampleID" )

# find min read # that you could rarefy to 
raremax = min(sampledata_rhizo.pp_pruned$samplesums) 
raremax # 4694 

# Plot rarefcation curve, richness * read # 
require(ggplot2)
ggplot(sampledata_rhizo.pp_pruned, aes(x= samplesums, y = Observed, color = Variety)) + 
    geom_point(size = 4) + 
  ggtitle("Rarefaction curve for all Rhizosphere reads") + 
  labs(x = "Total Reads", y = "Observed Richness")+
 # geom_text(label = sampledata_rhizo.pp_pruned$X.SampleID) +
  geom_vline(xintercept = raremax, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #


# Compare rarefied richness to observed richness 

rhizo.pp_pruned.otu =t(as(otu_table(rhizo.pp_pruned), "matrix"))
row.names(rhizo.pp_pruned.otu)=sampledata_rhizo.pp_pruned$X.SampleID 

richness_rarefied <- rarefy(rhizo.pp_pruned.otu, raremax)
rhizo.spp = specnumber(rhizo.pp_pruned.otu) #specnumber finds frequences of species

rare_plot <- plot(rhizo.spp, richness_rarefied, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "rhizo.pp_pruned, rarefied to 2026")
abline(0,1) # this line should nicely match the points

# Or plot with ggplot
rhizo.richness.comp <- merge(richness_rarefied, rhizo.spp, by = 0)
colnames(rhizo.richness.comp) <- c("X.SampleID","richness_rarefied", "richness_raw")

require(ggplot2)
ggplot(rhizo.richness.comp, aes(x= richness_rarefied, y = richness_raw)) + 
    geom_point(size = 4) + 
  ggtitle("") + 
  labs(x = "richness_rarefied", y = "richness_raw")+
 # geom_text(label = sampledata_rhizo.pp_pruned$X.SampleID) +
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #

###############
## Endosphere
################

# Rarefaction curve for every sample 
otus.tab=t(as(otu_table(endo.pp_pruned), "matrix"))
otus.map=sample_data(endo.pp_pruned)
row.names(otus.tab)=otus.map$X.SampleID 
raremax <- min(rowSums(otus.tab)) #  
raremax # 1984
rarecurve <- rarecurve(otus.tab, step = 1000, sample = raremax,label = FALSE)


###################

sampledata_endo.pp_pruned <- as.data.frame(as(sample_data(endo.pp_pruned),"matrix"))
sampledata_endo.pp_pruned$samplesums <- sample_sums(endo.pp_pruned)


Observed <- estimate_richness(endo.pp_pruned, measures = "Observed")
library(tibble)
Observed <- rownames_to_column(Observed, "X.SampleID")
sampledata_endo.pp_pruned <- merge(sampledata_endo.pp_pruned, Observed, by = "X.SampleID" )

# find min read # that you could rarefy to 
raremax = min(sampledata_endo.pp_pruned$samplesums) 
raremax # 1984

# Plot rarefcation curve, richness * read # 
require(ggplot2)
ggplot(sampledata_endo.pp_pruned, aes(x= samplesums, y = Observed, color = Variety)) + 
    geom_point(size = 4) + 
  ggtitle("Rarefaction curve for all Endosphere reads") + 
  labs(x = "Total Reads", y = "Observed Richness")+
 # geom_text(label = sampledata_endo.pp_pruned$X.SampleID) +
  geom_vline(xintercept = raremax, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #


# Compare rarefied richness to observed richness 

endo.pp_pruned.otu =t(as(otu_table(endo.pp_pruned), "matrix"))
row.names(endo.pp_pruned.otu)=sampledata_endo.pp_pruned$X.SampleID 

richness_rarefied <- rarefy(endo.pp_pruned.otu, raremax)
endo.spp = specnumber(endo.pp_pruned.otu) #specnumber finds frequences of species

rare_plot <- plot(endo.spp, richness_rarefied, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "endo.pp_pruned, rarefied to 2026")
abline(0,1) # this line should nicely match the points

# Or plot with ggplot
endo.richness.comp <- merge(richness_rarefied, endo.spp, by = 0)
colnames(endo.richness.comp) <- c("X.SampleID","richness_rarefied", "richness_raw")

require(ggplot2)
ggplot(endo.richness.comp, aes(x= richness_rarefied, y = richness_raw)) + 
    geom_point(size = 4) + 
  ggtitle("") + 
  labs(x = "richness_rarefied", y = "richness_raw")+
 # geom_text(label = sampledata_endo.pp_pruned$X.SampleID) +
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.5) + 
  theme(
  plot.title = element_text(size = 30), 
  axis.title.x = element_text(size = 25), 
  axis.title.y = element_text(size = 25), 
  axis.text.x = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20))# add visual line of where you could rarefy to min. read #

```
##6a) Rarefy Normalization
```{r warning=TRUE}
## ALL data 
set.seed(28132) # set this so your data is reproducable
bact.rfy2026 <- rarefy_even_depth(bact.pp_pruned, sample.size = 2026, replace = FALSE, rngseed = TRUE)
  #6338OTUs were removed because they are no longer present in any sample after random subsampling 
ntaxa(bact.rfy2026) # 12197
sum(sample_sums(bact.rfy2026)) # 368732 total reads

## Rhizosphere
set.seed(28132) # set this so your data is reproducable
rhizo.rfy4694 <- rarefy_even_depth(rhizo.pp_pruned, sample.size = 4694, replace = FALSE, rngseed = TRUE)
  #3288OTUs were removed because they are no longer present in any sample after random subsampling
 ntaxa(rhizo.rfy4694) # 14590
sum(sample_sums(rhizo.rfy4694)) # 647772 total reads

## Endosphere
set.seed(28132) # set this so your data is reproducable
endo.rfy1984 <- rarefy_even_depth(endo.pp_pruned, sample.size = 1984, replace = FALSE, rngseed = TRUE)
  #578OTUs were removed because they are no longer present in any sample after random subsampling 
ntaxa(endo.rfy1984) # 1944
sum(sample_sums(endo.rfy1984)) # 85312 total reads

###########################

# 8) Save files as phyloseq objects to open easier later 

# Rarified Data 
bact.rfy <- bact.rfy2026
#save(bact.rfy, file = "SGvar_16S_RhizoEndo_rfy.Rdata")

rhizo.rfy <- rhizo.rfy4694
#save(rhizo.rfy, file = "SGvar_16S_Rhizo_rfy.Rdata")

endo.rfy <- endo.rfy1984
#save(endo.rfy, file = "SGvar_16S_EndoOnly_rfy.Rdata")


```


## 6b) Deseq VST normalization 

```{r}


#https://joey711.github.io/phyloseq-extensions/DESeq2.html
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")
#biocLite("vsn")
#biocLite("phyloseq")
library(DESeq2)#bioconductor package
library(ggplot2)
library(vsn)#bioconductor package

#############
# All samples (rhizo & endo) 
#############


#nomalizing using DESEQ2
bact.pp.pruned.deseq=phyloseq_to_deseq2(bact.pp_pruned, ~X.SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
# Blind = TRUE means you make the transformation blind to the experimental design. 
vsd.bact<-varianceStabilizingTransformation(bact.pp.pruned.deseq)
# If no error - use this for your final transformation!
#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

plotSparsity(assay(vsd.bact))
meanSdPlot(assay(vsd.bact))

bact.matrix_obj_vst=assay(vsd.bact)#variance stabilized transformation OTU table in matrix form
 
#now we need to turn it back into a phyloseq object
bact.OTU.table.vst = otu_table(bact.matrix_obj_vst, taxa_are_rows = TRUE)

## add the minimum value back to the dataframe to make negative values == 0 
min(otu_table(bact.OTU.table.vst)) # -0.1154438
bact.vst_otu <- as.data.frame(as(otu_table(bact.OTU.table.vst),"matrix"))
bact.vst_otu_int <- bact.vst_otu + abs(min(otu_table(bact.OTU.table.vst)))
min(bact.vst_otu_int) # 0

bact.vst_otu_int = otu_table(bact.vst_otu_int, taxa_are_rows = TRUE)


bact.vst=phyloseq(bact.vst_otu_int,tax_table(bact.pp_pruned), sample_data(bact.pp_pruned), phy_tree(bact.pp_pruned))

#save(bact.vst, file = "SGvar_16S_RhizoEndo_vst.Rdata")


############
# Rhizosphere (Soil) samples 
#############

#nomalizing using DESEQ2
rhizo.pp.pruned.deseq=phyloseq_to_deseq2(rhizo.pp_pruned, ~X.SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
# Blind = TRUE means you make the transformation blind to the experimental design. 
vsd.rhizo<-varianceStabilizingTransformation(rhizo.pp.pruned.deseq)
# If no error - use this for your final transformation!
#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

plotSparsity(assay(vsd.rhizo))
meanSdPlot(assay(vsd.rhizo))

rhizo.matrix_obj_vst=assay(vsd.rhizo)#variance stabilized transformation OTU table in matrix form

#now we need to turn it back into a phyloseq object
rhizo.OTU.table.vst = otu_table(rhizo.matrix_obj_vst, taxa_are_rows = TRUE)

## add the minimum value back to the dataframe to make negative values == 0 
min(otu_table(rhizo.OTU.table.vst)) # 3.37 - no negative values so it's good as is.


rhizo.vst=phyloseq(rhizo.OTU.table.vst,tax_table(rhizo.pp_pruned), sample_data(rhizo.pp_pruned),phy_tree(rhizo.pp_pruned))

#save(rhizo.vst, file = "SGvar_16S_Rhizo_vst.Rdata")


############
# Endosphere (Root) samples
#############


#nomalizing using DESEQ2
endo.pp.pruned.deseq=phyloseq_to_deseq2(endo.pp_pruned, ~X.SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
# Blind = TRUE means you make the transformation blind to the experimental design. 
vsd.endo<-varianceStabilizingTransformation(endo.pp.pruned.deseq)
# If no error - use this for your final transformation!
#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

plotSparsity(assay(vsd.endo))
meanSdPlot(assay(vsd.endo))

endo.matrix_obj_vst=assay(vsd.endo)#variance stabilized transformation OTU table in matrix form
 
#now we need to turn it back into a phyloseq object
endo.OTU.table.vst = otu_table(endo.matrix_obj_vst, taxa_are_rows = TRUE)

## add the minimum value back to the dataframe to make negative values == 0 
min(otu_table(endo.OTU.table.vst)) # -5.786465
endo.vst_otu <- as.data.frame(as(otu_table(endo.OTU.table.vst),"matrix"))
endo.vst_otu_int <- endo.vst_otu + abs(min(endo.vst_otu))
min(endo.vst_otu_int) # 0

endo.OTU.table.vst_int = otu_table(endo.vst_otu_int, taxa_are_rows = TRUE)


endo.vst=phyloseq(endo.OTU.table.vst_int,tax_table(endo.pp_pruned), sample_data(endo.pp_pruned),phy_tree(endo.pp_pruned))

#save(endo.vst, file = "SGvar_16S_EndOnly_vst.Rdata")

```
---
# FUNGI OTU FILTERING & NORMALIZATION 

## 1) Make phyloseq object with otu table, taxa file, and tree file (don't have the tree file yet)
```{r}
# load metadata file with X.SampleID as the column with the same sample ids in otu table 
metadata <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
SAMP <- sample_data(metadata) # this is the metadata file that you merge with phyloseq
rownames(SAMP) <- SAMP$X.SampleID # row names must match OTU table


# load otu table and otu table 
# Note: Usearch output has a #OTU ID in cell A1 - you need to remove this for it to work.
otufile <- read.table("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Final_otu_tax_tree_files/TMC_SGVariety_ITS_OTU97_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.txt", header = TRUE)
otu <- otu_table(otufile, taxa_are_rows = TRUE)

# load taxa file 
    # tax file - from export, separated by ; delimitation, removed taxonomy headings "k:", cleaned to have only taxonomys with 80perc confidence; changed all blanks to "unclassified"
taxafile_otu <- as.matrix(read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Final_otu_tax_tree_files/TMC_SGVariety_ITS_80percTaxonomy_unite7.2_editted.csv", header = TRUE, row.names =1))
tax_otu <- tax_table(taxafile_otu)

# combine data into phyloseq-object (include map, otu, tax, tree)
myphy_otus<- merge_phyloseq(otu, tax_otu, SAMP)
sampledata <- as(sample_data(myphy_otus), "data.frame")

# check the dimensions of OTU 
dim(otu_table(myphy_otus)) #5790 144 (taxa, samples)

```

## 2) Rename column names (taxa) 
```{r}
#what are the column names of our taxonomy file? 
colnames(tax_table(myphy_otus)) # good to go

# remove otus not in our samples 
otus.p <- prune_taxa(taxa_sums(myphy_otus) > 0, myphy_otus)
ntaxa(myphy_otus) - ntaxa(otus.p) # removed 1054 otus that weren't present in my samples (only in Robert's and MMPRNTS) 
ntaxa(otus.p) #4736 taxa in my samples 
sum(sample_sums(otus.p)) # 2202804 sequences before filtering


```

## 3) Remove Singletons
```{r}
# Now remove singletons 
otus.pp <- prune_taxa(taxa_sums(otus.p) > 1, otus.p)
ntaxa(otus.p)- ntaxa(otus.pp) # removed 97 singletons
ntaxa(otus.pp) #4639

sum(sample_sums(otus.pp)) #2202707 total reads
```

## 4) Remove low sequence samples 
```{r}

# all samples 
sort(sample_sums(otus.pp))  
  # Filtered samples below 4100 reads; Lost samples G5, E10, C4, G3, G1, J4, H5, B2, A9
otus.pp_prune <- prune_samples(sample_sums(otus.pp) > 4100, otus.pp) 
sort(sample_sums(otus.pp_prune)) 
nsamples(otus.pp_prune) # 135
mean(sample_sums(otus.pp_prune)) # average read # = 16268.73
median(sample_sums(otus.pp_prune)) # median read # = 16415
min(sample_sums(otus.pp_prune)) #4163
max(sample_sums(otus.pp_prune)) # 38470
sum(sample_sums(otus.pp_prune)) #2196278
ntaxa(otus.pp_prune)# 4639
 
max(sample_sums(otus.pp_prune))/min(sample_sums(otus.pp_prune)) # 9 fold difference in reads
``` 

## 5) Investigate Sequencing Depth
```{r}

# What does the total reads per sample distribution look like? 
sample_sum_df <- data.frame(sum = sample_sums(otus.pp))

ggplot(sample_sum_df, aes(x = sum)) +   # Histogram of sample read counts 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) + 
  ggtitle("Distribution of sample sequencing depth of otus") + 
   xlab("Read counts") + 
  theme(axis.title.y = element_blank())

# mean, max, min, median of sample read counts to inform normalization techniques 
min(sample_sums(otus.pp)) # 12 reads
mean(sample_sums(otus.pp)) # 15296.58 reads
max(sample_sums(otus.pp)) # 38370 reads
median(sample_sums(otus.pp)) # 15715 reads

# plot sequencing depth 
require(data.table)
# make a histogram of total counts 
tdt <- data.table(tax_table(otus.pp),
                  TotalCounts = taxa_sums(otus.pp),
                  OTU = taxa_names(otus.pp))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts for otus")


################
# all samples rarefaction curve 
richness <- estimate_richness(otus.pp, measures = "Observed")
samplesums <- sample_sums(otus.pp)
richness_sums <- merge(richness,samplesums, by = "row.names")
# y is samplesums
require(ggplot2)
ggplot(richness_sums, aes(x=y, y = Observed)) + 
    geom_point() + 
  ggtitle("otu Observed Richness by Total Reads") + 
  labs(x = "Total Reads", y = "Observed Richness")+
  geom_text(label = richness_sums$Row.names)

``` 
### Rarefaction Curves 
```{r, echo=FALSE}

richness <- estimate_richness(otus.pp, measures = "Observed")
samplesums <- sample_sums(otus.pp)
richness_sums <- merge(richness,samplesums, by = "row.names")
# y is samplesums
require(ggplot2)
ggplot(richness_sums, aes(x=y, y = Observed)) + 
    geom_point() + 
  ggtitle("otu Observed Richness by Total Reads") + 
  labs(x = "Total Reads", y = "Observed Richness")+
  geom_text(label = richness_sums$Row.names)


# Plot a rarefaction curve to determine cut-off point for read depth
otus.tab=t(as(otu_table(otus.pp_prune), "matrix"))
otus.map=sample_data(otus.pp_prune)
row.names(otus.tab)=otus.map$X.SampleID 
raremax <- min(rowSums(otus.tab)) #  4153 - this is the min read number 

rarecurve <- rarecurve(otus.tab, step = 1000, sample = raremax,label = FALSE)

otus.rare = rarefy(otus.tab, raremax)
otus.spp = specnumber(otus.tab) #specnumber finds frequences of species

rare_plot <- plot(otus.spp, otus.rare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "otus.pp_prune, rarefied to 4153")
abline(0,1) # this line should nicely match the points

```
## 6a) Rarefy Normalization
```{r}
## ALL data 
set.seed(28132) # set this so your data is reproducable
fungi.rfy4153 <- rarefy_even_depth(otus.pp_prune, sample.size = 4153, replace = FALSE, rngseed = TRUE)
  #575OTUS were removed because they are no longer present in any sample after random subsampling 
ntaxa(fungi.rfy4153) # 4064
sum(sample_sums(fungi.rfy4153)) # 560655 total reads


# Save files as phyloseq objects to open easier later 

# Rarified Data 
fungi.rfy <- fungi.rfy4153
#save(fungi.rfy, file = "SGvar_ITS_Rhizo_rfy.Rdata")


```
## 6b) Deseq VST normalization 
```{r}
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")
#biocLite("vsn")
#biocLite("phyloseq")
library(DESeq2)#bioconductor package
library(ggplot2)
library(vsn)#bioconductor package

#nomalizing using DESEQ2
fungi.pp.pruned.deseq=phyloseq_to_deseq2(otus.pp_prune, ~X.SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
# Blind = TRUE means you make the transformation blind to the experimental design. 
vsd.fungi<-varianceStabilizingTransformation(fungi.pp.pruned.deseq)
# If no error - use this for your final transformation!
#may result in the below error
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

plotSparsity(assay(vsd.fungi))
meanSdPlot(assay(vsd.fungi))

fungi.matrix_obj_vst=assay(vsd.fungi)#variance stabilized transformation OTU table in matrix form
 
#now we need to turn it back into a phyloseq object
fungi.OTU.table.vst = otu_table(fungi.matrix_obj_vst, taxa_are_rows = TRUE)

# Now make the min. negative values == 0 
min(otu_table(fungi.OTU.table.vst)) # -4.42
fungi.vst_otu <- as.data.frame(as(otu_table(fungi.OTU.table.vst),"matrix"))
fungi.vst_otu_int <- fungi.vst_otu + abs(min(fungi.vst_otu))
min(fungi.vst_otu_int) # 0

fungi.vst_otu_int = otu_table(fungi.vst_otu_int, taxa_are_rows = TRUE)

fungi.vst=phyloseq(fungi.vst_otu_int,tax_table(otus.pp_prune), sample_data(otus.pp_prune))


#save(fungi.vst, file = "SGvar_ITS_Rhizo_vst.Rdata")

```