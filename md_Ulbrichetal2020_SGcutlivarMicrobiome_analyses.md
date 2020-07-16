---
title: "Ulbrich et al. 2020 SG cultivar analyses"
author: "Tayler Ulbrich"
date: "7/1/2020"
output: 
  html_document: 
    keep_md: yes
    theme: cosmo
    toc: yes
editor_options: 
  chunk_output_type: console
  
---

```{r}
library(knitr)
knitr::opts_chunk$set(
  warning = FALSE, # show warnings
  message = TRUE, # show messages
  error = TRUE, # do not interrupt generation in case of errors,
  echo = TRUE  # show R code
)
```
*Resources Used:*
EDAMAME Tutorial: https://github.com/raleva/Edamame_phyloseq/blob/master/Edamame_phyloseq.Rmd
Phyloseq Tutorial: http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
https://rpubs.com/maddieSC/R_SOP_UCR_Jan_2018
Microbiome course tutorial: https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html#subset-by-taxonomy

*metadata file:* C:/Users/Owner/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/Final_metadata/
2018.10.25_SGvariety_RhizoEndo_metadata_SroleyFinal_RootsFinal_allpseudoreps_coords_TOCdata.csv

*Bacterial data files:*
- rhizo.rfy - rarified soil-associated samples (12 cultivars)
- bact.rfy - combined and rarified root and soil samples (4 cultivars); subset into 4 cultivars = rhizoendo; or rhizoendo_root for the 4 cultivars' root samples; or rhizoendo_soil for the 4 cultivars' soil samples 

*Fungal data files:* 
fungi.rfy - rarified soil-associated samples (12 cultivars)


*OTU filtering notes:*
Three microbiome files are processed in this script, each of which were filtered for non-bacterial or non-fungal reads, singeltons (<2 reads), and rarefied.
- bacterial root subset of 4 cultivars (endo) rarefied to 1984 reads
- bacterial soil & root subset of 4 cultivars (rhizoendo) rarefied to 2026 reads
- bacterial soil (rhizo) for 12 cultivars rarefied to 4694 reads
- fungal soil (rhizo) for 12 cultivars rarefied to 4153 reads

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
metadata <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
SAMP <- sample_data(metadata) # this is the metadata file that you merge with phyloseq
rownames(SAMP) <- SAMP$X.SampleID # row names must match OTU table


# load otu table and otu table 
# Note: Usearch output has a #OTU ID in cell A1 - you need to remove this for it to work.
otufile <- read.table("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/Final_otu_tax_tree_files/TMC_SGVar16S_RhizoEndo_combined_merged_ALLruns_OTU_table.txt", header = TRUE)
otu <- otu_table(otufile, taxa_are_rows = TRUE)

# load taxa file 
    # tax file - from export, separated by ; delimitation, removed taxonomy headings "k:", cleaned to have only taxonomys with 80perc confidence; changed all blanks to "unclassified"

taxafile_otu <- as.matrix(read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Final_otu_tax_tree_files/TMC_SGVar16S_RhizoEndo_combined_merged_ALLruns_otus_taxonomy_SilvaV123_sintax_80perc_editted.csv", header = TRUE, row.names =1))
tax_otu <- tax_table(taxafile_otu)


# load tree file 
# note: I used figtree program to convert the .tre file to a newick file format - required for read_tree
tree_otu_file <- "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Final_otu_tax_tree_files/2018.06.22_SGvarietyRhizoEndo16s_pasta.nwk"
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
ntaxa(myphy_otus) - ntaxa(otus.p) # removed 27787 otus that weren't present in my samples (only in Robert's and MMPRNTS) 
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
#save(bact.rfy, file = "C:/Users/Owner/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20190212_SGVar16S_pp_nosingles_rfy_RhizoEndo.Rdata")

rhizo.rfy <- rhizo.rfy4694
#save(rhizo.rfy, file = "C:/Users/Owner/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20190212_SGVar16S_pp_nosingles_rfy_RhizoONLY.Rdata")

endo.rfy <- endo.rfy1984
#save(endo.rfy, file = "C:/Users/Owner/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20190212_SGVar16S_pp_nosingles_rfy_EndoONLY.Rdata")


```


## 6b) Deseq VST normalization 
### Bact.pp VST
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

#save(bact.vst, file = "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20200501_SGVar16Sbact_pp_nosingles_vst.Rdata")


```
### Rhizo.pp VST
```{r}

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

#save(rhizo.vst, file = "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20200501_SGVar16Srhizo_pp_nosingles_vst_rhizoOnly.Rdata")




```
### Endo.pp VST
```{r}


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

#save(endo.vst, file = "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/PhyloseqObjects/20200501_SGVar16Sendo_pp_nosingles_vst_endoOnly.Rdata")





```
---
# FUNGI OTU FILTERING & NORMALIZATION 
## 1) Make phyloseq object with otu table, taxa file, and tree file (don't have the tree file yet)
```{r}
# load metadata file with X.SampleID as the column with the same sample ids in otu table 
metadata <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
SAMP <- sample_data(metadata) # this is the metadata file that you merge with phyloseq
rownames(SAMP) <- SAMP$X.SampleID # row names must match OTU table


# load otu table and otu table 
# Note: Usearch output has a #OTU ID in cell A1 - you need to remove this for it to work.
otufile <- read.table("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Final_otu_tax_tree_files/TMC_SGVariety_ITS_OTU97_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.txt", header = TRUE)
otu <- otu_table(otufile, taxa_are_rows = TRUE)

# load taxa file 
    # tax file - from export, separated by ; delimitation, removed taxonomy headings "k:", cleaned to have only taxonomys with 80perc confidence; changed all blanks to "unclassified"
taxafile_otu <- as.matrix(read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Final_otu_tax_tree_files/TMC_SGVariety_ITS_80percTaxonomy_unite7.2_editted.csv", header = TRUE, row.names =1))
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
#save(fungi.rfy, file = "C:/Users/Owner/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/ITS/PhyloseqObjects/2019.02.14_TMC_SGVariety_ITS_rfy_nosingles.Rdata")


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



#save(fungi.vst, file = "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/ITS/PhyloseqObjects/20200501_SGVarITS_fungi_pp_nosingles_vst.Rdata")

```
---
# ANALYSIS 
# Load metadata 
```{r}
metadata <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
dim(metadata)

# subset for only rhizosphere samples so there aren't duplicate environmental data with the root samples (should be n = 144)
metadata <- filter(metadata, SampleSite == "Rhizosphere")
dim(metadata)


```

# Load phyloseq objects
```{r}
# Do not use the separately rarefied root samples (endo.rfy or endo.vst) for analysis 
# Only use soil (n = 12 cultivars) and root or soil (rhizoendo, n = 4 cultivars)
### Soil-associated (rhizosphere) Only 
#  rarified data (rhizo.rfy)
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_16S_Rhizo_rfy.Rdata")
# extract metadata for future analysis 
sampledata_rhizo.rfy <- as(sample_data(rhizo.rfy), "data.frame")


### Combined rhizosphere and endosphere  
#  rarified data (bact.rfy)
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_16S_RhizoEndo_rfy.Rdata")
# extract metadata for future analysis 
sampledata_bact.rfy <- as(sample_data(bact.rfy), "data.frame")


### Rhizosphere Fungi 
#  rarified data (fungi.rfy)
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_ITS_Rhizo_rfy.Rdata")

# extract metadata for future analysis 
sampledata_fungi.rfy <- as(sample_data(fungi.rfy), "data.frame")
# reorder cultivars so that they always appear in the same order 
#sampledata_fungi.rfy$Variety <- factor(sampledata_fungi.rfy$Variety, c("Alamo", "Kanlow", "EG1101", "EG1102", "Blackwell", "Cave-in-Rock", "Dacotah", "NE28", "EG2101", "Shelter", "Southlow","Trailblazer"))

###########################################
# DESEQ Variance Stablizing Transformation data 

### Add the lowest minimum value to the matrix so you don't have neg values. 

#Rhizo.vst
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_16S_Rhizo_vst.Rdata")
sampledata_rhizo.vst <- as(sample_data(rhizo.vst), "data.frame")

# Bact (Rhizoendo, bact.vst)
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_16S_RhizoEndo_vst.Rdata")
sampledata_bact.vst <- as(sample_data(bact.vst), "data.frame")


# Fung.vst (soil fungi)
load("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PhyloseqObjects_forAnalysis/SGvar_ITS_Rhizo_vst.Rdata")
sampledata_fungi.vst <- as(sample_data(fungi.vst),"data.frame")


```

## Split bact.rfy into rhizoendo for 4 cultivars analyzed for rhizo and endo
```{r}
# 4 cultivars with endo and rhizo samples 
rhizoendo <- subset_samples(bact.rfy, sampledata_bact.rfy$Variety == "Alamo" | sampledata_bact.rfy$Variety =="Kanlow" | sampledata_bact.rfy$Variety =="Cave-in-Rock" |sampledata_bact.rfy$Variety == "Southlow")
  # remove otus not in rhizo or endo subset
rhizoendo<- prune_taxa(taxa_sums(rhizoendo) >0, rhizoendo)
#ntaxa(rhizoendo1) - ntaxa(rhizoendo) #3557 otus removed that weren't present in the 4 cultivars

# extract metadata 
sampledata_rhizoendo <- as(sample_data(rhizoendo), "data.frame")

# Subset for soil and root (4 cultivars each)

rhizoendo.soil <- subset_samples(rhizoendo, sampledata_rhizoendo$SampleSite == "Rhizosphere")
sampledata_rhizoendo.soil <- as(sample_data(rhizoendo.soil), "data.frame")
nsamples(rhizoendo.soil) # 46 

rhizoendo.root <- subset_samples(rhizoendo, sampledata_rhizoendo$SampleSite == "Endosphere")
sampledata_rhizoendo.root <- as(sample_data(rhizoendo.root), "data.frame")
nsamples(rhizoendo.root)#42 

```
---
# Compare Normalization techniques 
*Protest to compare data normalized with DESEQ2 VST or rarefaction*
```{r}
# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 

###################
### bacterial soil community
###################
# rfy
set.seed(2)
rhizo.rfy_w.unifrac <- phyloseq::distance(rhizo.rfy, "wunifrac" )
rhizo.rfy_w.unifrac.m <- as.matrix(rhizo.rfy_w.unifrac)
# stress = 0.18  

# VST 
set.seed(2)
rhizo.vst_w.unifrac <- phyloseq::distance(rhizo.vst, "wunifrac")
rhizo.vst_w.unifrac.m <- as.matrix(rhizo.vst_w.unifrac)

# protest compare ordinations
## RFY vs. VST 
set.seed(2)
protest( X = rhizo.rfy_w.unifrac, Y = rhizo.vst_w.unifrac, scores = "sites", permutations = 999)

#Procrustes Sum of Squares (m12 squared):        0.1661 
#Correlation in a symmetric Procrustes rotation: 0.9132 
#Significance:  0.001 

#######################
### Bact - Soil &  Root (4 cultivars)
#######################

# 4 cultivars with endo and rhizo samples 
rhizoendo.rfy <- subset_samples(bact.rfy, sampledata_bact.rfy$Variety == "Alamo" | sampledata_bact.rfy$Variety =="Kanlow" | sampledata_bact.rfy$Variety =="Cave-in-Rock" |sampledata_bact.rfy$Variety == "Southlow")
  # remove otus not in rhizo or endo subset
rhizoendo.rfy <- prune_taxa(taxa_sums(rhizoendo.rfy) >0, rhizoendo.rfy)
nsamples(rhizoendo.rfy) #88 

rhizoendo.vst <- subset_samples(bact.vst, sampledata_bact.vst$Variety == "Alamo" | sampledata_bact.vst$Variety =="Kanlow" | sampledata_bact.vst$Variety =="Cave-in-Rock" |sampledata_bact.vst$Variety == "Southlow")
  # remove otus not in rhizo or endo subset
rhizoendo.vst <- prune_taxa(taxa_sums(rhizoendo.vst) >0, rhizoendo.vst)
nsamples(rhizoendo.vst) #88 

# rfy
set.seed(2)
rhizoendo.rfy_w.unifrac <- phyloseq::distance(rhizoendo.rfy, "wunifrac" )
rhizoendo.rfy_w.unifrac.m <- as.matrix(rhizoendo.rfy_w.unifrac)

# VST 
set.seed(2)
rhizoendo.vst_w.unifrac <- phyloseq::distance(rhizoendo.vst, "wunifrac")
rhizoendo.vst_w.unifrac.m <- as.matrix(rhizoendo.vst_w.unifrac)

# protest compare ordinations
# RFY vs. VST 
set.seed(2)
protest( X = rhizoendo.rfy_w.unifrac.m, Y = rhizoendo.vst_w.unifrac.m, scores = "sites", permutations = 999)

# after adding zero to get rid of neg. values in VST transformation 
#Procrustes Sum of Squares (m12 squared):        0.8324 
#Correlation in a symmetric Procrustes rotation: 0.4094 
#Significance:  0.001 

 

##################
### Fungal soil community
##################
# rfy
set.seed(2)
fungi.rfy_bray <- phyloseq::distance(fungi.rfy, "bray" )
fungi.rfy_bray.m <- as.matrix(fungi.rfy_bray)

# VST 

set.seed(2)
fungi.vst_bray <- phyloseq::distance(fungi.vst, "bray")
fungi.vst_bray.m <- as.matrix(fungi.vst_bray)
# Warning data may be meaningless because of negative entries in bray 

# protest compare ordinations
# RFY vs. vst 
# after addin constant to get rid of vst negatives, bray can't run with negatives 
set.seed(2)
protest( X = fungi.rfy_bray.m, Y = fungi.vst_bray.m, scores = "sites", permutations = 999)

#Procrustes Sum of Squares (m12 squared):        0.3223 
#Correlation in a symmetric Procrustes rotation: 0.8232 
#Significance:  0.001 

```
### Permanova compare rfy and vst 
```{r}

#################
### Bacteria SOIL  
################
# rfy
set.seed(2)
rhizo.rfy.wunifrac<- phyloseq::distance(rhizo.rfy,"wunifrac")

# vst
set.seed(2) 
rhizo.vst.wunifrac<- phyloseq::distance(rhizo.vst,"wunifrac")

# prop
rhizo.raw_ra  = transform_sample_counts(rhizo.pp_pruned, function(x) x / sum(x))
set.seed(2) 
rhizo.raw_ra.wunifrac <- phyloseq::distance(rhizo.raw_ra, "wunifrac")


## RFY 
set.seed(2)
adonis(rhizo.rfy.wunifrac~sampledata_rhizo.rfy$Date_Julian + sampledata_rhizo.rfy$Variety, permutations = 999)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#sampledata_rhizo.rfy$Date_Julian   1   0.09919 0.099195 15.4432 0.09515  0.001 ***
#sampledata_rhizo.rfy$Variety      10   0.13396 0.013396  2.0856 0.12850  0.001 ***
#Residuals                        126   0.80933 0.006423         0.77634           
#Total                            137   1.04248                  1.00000                                  42   2.02617                  1.00000    

#VST 
set.seed(2)
adonis(rhizo.vst.wunifrac~sampledata_rhizo.rfy$Date_Julian + sampledata_rhizo.rfy$Variety, permutations = 999)
#Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)    
#sampledata_rhizo.rfy$Date_Julian   1 0.00008258 8.2584e-05 15.2577 0.09358  0.001 ***
#sampledata_rhizo.rfy$Variety      10 0.00011797 1.1797e-05  2.1795 0.13367  0.001 ***
#Residuals                        126 0.00068199 5.4130e-06         0.77276           
#Total                            137 0.00088254                    1.00000                                42  0.056671                   #1.00000     


#################
### Bacteria Rhizo vs. Endo 
################

# 4 cultivars with endo and rhizo samples 
rhizoendo.rfy <- subset_samples(bact.rfy, sampledata_bact.rfy$Variety == "Alamo" | sampledata_bact.rfy$Variety =="Kanlow" | sampledata_bact.rfy$Variety =="Cave-in-Rock" |sampledata_bact.rfy$Variety == "Southlow")
  # remove otus not in rhizo or endo subset
rhizoendo.rfy <- prune_taxa(taxa_sums(rhizoendo.rfy) >0, rhizoendo.rfy)
nsamples(rhizoendo.rfy) #88 

rhizoendo.vst <- subset_samples(bact.vst, sampledata_bact.vst$Variety == "Alamo" | sampledata_bact.vst$Variety =="Kanlow" | sampledata_bact.vst$Variety =="Cave-in-Rock" |sampledata_bact.vst$Variety == "Southlow")
  # remove otus not in rhizo or endo subset
rhizoendo.vst <- prune_taxa(taxa_sums(rhizoendo.vst) >0, rhizoendo.vst)
nsamples(rhizoendo.vst) #88 


# rfy
set.seed(2)
rhizoendo.rfy.wunifrac<- phyloseq::distance(rhizoendo.rfy,"wunifrac")

# vst
set.seed(2) 
rhizoendo.vst.wunifrac<- phyloseq::distance(rhizoendo.vst,"wunifrac")

## RFY 
set.seed(2)
adonis(rhizoendo.rfy.wunifrac~sample_data(rhizoendo.rfy)$Variety* sample_data(rhizoendo.rfy)$SampleSite, permutations = 999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sample_data(rhizoendo.rfy)$Variety                                        3    0.1415  0.0472   1.762 0.02454  0.100 .  
#sample_data(rhizoendo.rfy)$SampleSite                                     1    3.3858  3.3858 126.458 0.58718  0.001 ***
#sample_data(rhizoendo.rfy)$Variety:sample_data(rhizoendo.rfy)$SampleSite  3    0.0969  0.0323   1.207 0.01681  0.289    
#Residuals                                                                80    2.1419  0.0268         0.37146           
#Total                                                                    87    5.7662                 1.00000           

#VST 
set.seed(2)
adonis(rhizoendo.vst.wunifrac~sample_data(rhizoendo.vst)$Variety*sample_data(rhizoendo.vst)$SampleSite, permutations = 999)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#sample_data(rhizoendo.vst)$Variety                                        3     13.56  4.5209 1.13150 0.03852  0.443
#sample_data(rhizoendo.vst)$SampleSite                                     1      1.66  1.6580 0.41497 0.00471  0.571
#sample_data(rhizoendo.vst)$Variety:sample_data(rhizoendo.vst)$SampleSite  3     17.28  5.7599 1.44160 0.04907  0.267
#Residuals                                                                80    319.64  3.9955         0.90771       
#Total                                                                    87    352.14                 1.00000    

###################
# Fungi soil 
##################

# rfy
set.seed(2)
fungi.rfy.bray<- phyloseq::distance(fungi.rfy,"bray")

# vst
set.seed(2) 
fungi.vst.bray<- phyloseq::distance(fungi.vst,"bray")


## RFY 
set.seed(2)
adonis(fungi.rfy.bray~sampledata_fungi.rfy$Date_Julian + sampledata_fungi.rfy$Variety, permutations = 999)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sampledata_fungi.rfy$Date_Julian   1    0.4994 0.49942  2.2491 0.01618  0.001 ***
#sampledata_fungi.rfy$Variety      10    3.0607 0.30607  1.3783 0.09914  0.001 ***
#Residuals                        123   27.3132 0.22206         0.88468           
#Total                            134   30.8734                 1.00000       

#VST 
set.seed(2)
adonis(fungi.vst.bray~sampledata_fungi.rfy$Date_Julian + sampledata_fungi.rfy$Variety, permutations = 999)
#Df SumsOfSqs MeanSqs F.Model R2 Pr(>F)
#sampledata_fungi.rfy$Date_Julian   1         0       0                  
#sampledata_fungi.rfy$Variety      10         0       0                  
#Residuals                        123         0       0                  
#Total                            134         0          

```
---
# Summarize dominant taxonomic groups 
*what are the dominant phyla in each microbiome?*
```{r}
# Create summary data on the dominate phylum, etc.https://github.com/joey711/phyloseq/issues/418


library("phyloseq")
library("data.table")
library("ggplot2")

fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
Nsamples = nsamples(physeq)
  summarydt = mdt[, list(meanRA = sum(RelativeAbundance)/Nsamples,
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]

  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}

# Relative abundance of soil bacteria 
plot_taxa_summary(rhizo.rfy, "Phylum") # mean relative abundance of allphyla
rhizo.rfy_phyla_summary <- summarize_taxa(rhizo.rfy, "Phylum")
rhzo.rfy_phyla_var <- summarize_taxa(rhizo.rfy, "Phylum", "Variety")

# Relative abundance of soil fungi  
plot_taxa_summary(fungi.rfy, "Phylum") # mean relative abundance of allphyla
fungi.rfy_phyla_summary <- summarize_taxa(fungi.rfy, "Phylum")
fungi.rfy_phyla_var <- summarize_taxa(fungi.rfy, "Phylum", "Variety")


# Relative abundance of dominant phyla for root and soil (n = 4 cultivars) 

# roots 
rhizoendo.root_phyla_summary <- summarize_taxa(rhizoendo.root, "Phylum")

# soils
rhizoendo.soil_phyla_summary <- summarize_taxa(rhizoendo.soil, "Phylum")

```

# In-text Figures 
*Note: figures saved as metafiles and text size and colors editted in powerpoint, then exported as .tiff for final manuscript*
## Figure 1: 6 panel plot with Root traits, MBC, MBN, bact div, Nfix
```{r}

Variety.colors.bw <- c("grey","grey","grey","grey","white","white","white","white","white","white","white","white","grey","white")

####################
# A) SRL (volume) boxplot 
metadata_var <- metadata[,c("Variety", "SRL_giaroots")]
metadata_eco <- metadata[,c("Ecotype", "SRL_giaroots")]
colnames(metadata_eco) <- colnames(metadata_var)
SRLgia_data <-rbind(metadata_var, metadata_eco)

# reorder cultivars 
SRLgia_data$Variety <- factor(SRLgia_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

SRLgia_plot <- ggplot(data = SRLgia_data, aes(x = Variety, y = SRL_giaroots, fill = Variety)) + 
       geom_boxplot()+
    scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  stat_summary(geom = 'text', label =c("ab","a","ab","b","ab","a","ab","ab","b","ab","ab","ab","0","0") , fun.y = max, vjust = -1) +
  theme_classic() + 
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "Specific Root Length (volume-weighted) (cm/cm3)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_blank(), 
  axis.text.x = element_text(size = 15, color = "black"), 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "right" , 
  #plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 
SRLgia_plot

###############
# B) Avg_root_width_diam boxplot 
metadata_var <- metadata[,c("Variety", "Avg_root_width_diam")]
metadata_eco <- metadata[,c("Ecotype", "Avg_root_width_diam")]
colnames(metadata_eco) <- colnames(metadata_var)
RootDiam_data <-rbind(metadata_var, metadata_eco)

# reorder cultivars 
RootDiam_data$Variety <- factor(RootDiam_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

RootDiamplot <- ggplot(data = RootDiam_data, aes(x = Variety, y = Avg_root_width_diam, fill = Variety)) + 
         geom_boxplot()+
    scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  stat_summary(geom = 'text', label =c("cd","d","bcd","a","cd","d","bcd","d","ab","bcd","abc","abc","0","0") , fun.y = max, vjust = -1) +
  theme_classic() + 
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "Average root diameter (cm)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_blank(), 
  axis.text.x = element_text(size = 15, color = "black") , 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "none" , 
 # plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 
RootDiamplot

########
# C) MBC boxplot 

# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 

metadata_noSMCout <- filter(metadata, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
dim(metadata_noSMCout)


metadata_var <- metadata_noSMCout[,c("Variety", "ugC_MicBiomass_g_dry_soil")]
metadata_eco <- metadata_noSMCout[,c("Ecotype", "ugC_MicBiomass_g_dry_soil")]
colnames(metadata_eco) <- colnames(metadata_var)
MBC_data <-rbind(metadata_var, metadata_eco)
dim(MBC_data)

# reorder cultivars 
MBC_data$Variety <- factor(MBC_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

MBCplot <- ggplot(data = MBC_data, aes(x = Variety, y = ugC_MicBiomass_g_dry_soil, fill = Variety)) + 
    geom_boxplot()+
      scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  theme_classic() + 
  stat_summary(geom = 'text', label =c("cde","e","ab","bcde","abcd","a","abc","de","cde","cde","a","a","0","0") , fun.y = max, vjust = -1) +
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "MBC (C(ug)/ soil(g))") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_blank(),#text(size = 15, color = "black"), 
  axis.text.x = element_text(size = 15, color = "black",angle=0, vjust=0.6) , 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "none" , 
  plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 


MBCplot

####################
# D) MBN boxplot 

 #soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 

metadata_noSMCout <- filter(metadata, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
dim(metadata_noSMCout)



metadata_var <- metadata_noSMCout[,c("Variety", "ugN_MicBiomass_g_dry_soil")]
metadata_eco <- metadata_noSMCout[,c("Ecotype", "ugN_MicBiomass_g_dry_soil")]
colnames(metadata_eco) <- colnames(metadata_var)
MBN_data <-rbind(metadata_var, metadata_eco)

# reorder cultivars 
MBN_data$Variety <- factor(MBN_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

MBNplot <- ggplot(data = MBN_data, aes(x = Variety, y = ugN_MicBiomass_g_dry_soil, fill = Variety)) + 
    geom_boxplot()+
    scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  stat_summary(geom = 'text', label =c("cd","c","c","d","b","ab","ab","c","ab","c","ab","a","0","0") , fun.y = max, vjust = -1) +
  theme_classic() + 
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "MBN (N(ug)/ soil(g))") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_blank(),#(size = 15, color = "black"), 
  axis.text.x = element_text(size = 15, color = "black",angle=0, vjust=0.6) , 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "none" , 
  plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 
MBNplot


###############################
# E) Bacterial Shannon Diversity (data table with shannon diversity data created with 'estimate_richness' function, see "shannon diversity" section below) 

# read in shannon diversity data 
rhizo.rfy_diversity <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Figures/rhizo.rfy_diversity.csv", header = TRUE)

# MERGE variety and ecotype 
metadata_var <- rhizo.rfy_diversity[,c("Variety", "Shannon")]
metadata_eco <- rhizo.rfy_diversity[,c("Ecotype", "Shannon")]
colnames(metadata_eco) <- colnames(metadata_var)
Shannon_data <-rbind(metadata_var, metadata_eco)

 
Shannon_data$Variety <- factor(Shannon_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

# Shanno Div
shannon_rhizo_plot <- ggplot(data = Shannon_data, aes(x = Variety, y = Shannon, fill = Variety)) +geom_boxplot()+
      scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  theme_classic() + 
  stat_summary(geom = 'text', label =c("ab","ab","a","a","bc","ab","c","ab","a","ab","ab","abc","0","0") , fun.y = max, vjust = -1) +
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "Shannon diversity") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_text(size = 15, color = "black"), 
  axis.text.x = element_text(size = 15, color = "black",angle=0, vjust=0.6) , 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "none" , 
  plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 
  
shannon_rhizo_plot
#################################################
# F) Prop. N fixers

#read in rhizosphere and endosphere proportion of OTUS with Nfix genes/normalized 16S OTU abundance datatables 
Rhizo_Nfix_Prop <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/Picrust/Filtered_greengenes_PicrustRESULTS/16SRhizo/2019.09.23_Rhizo_NfixOTUcounts_NormSeqSum_meta.csv", header = TRUE)


# boxplot 
metadata_var <- Rhizo_Nfix_Prop[,c("Variety", "Sum_NfixCount_Normalized")]
metadata_eco <- Rhizo_Nfix_Prop[,c("Ecotype", "Sum_NfixCount_Normalized")]
colnames(metadata_eco) <- colnames(metadata_var)
NfixRhizo_data <-rbind(metadata_var, metadata_eco)

# reorder cultivars 
NfixRhizo_data$Variety <- factor(NfixRhizo_data$Variety, levels = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"))

NfixRhizo_data_plot <- ggplot(data = NfixRhizo_data, aes(x = Variety, y = Sum_NfixCount_Normalized, fill = Variety)) + 
    geom_boxplot()+
      scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  theme_classic() + 
  stat_summary(geom = 'text', label =c("bcd","cde","a","bc","de","cde","e","bcde","ab","cde","bcde","e","0","0") , fun.y = max, vjust = -1) +
  scale_fill_manual(values = Variety.colors.bw)+
  labs(x = "Cultivar", y = "Predicted proportion of N-fixers") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
  axis.title.y = element_text(size = 15, color = "black"), 
  axis.text.y = element_text(size = 15, color = "black"),
  axis.title.x =  element_text(size = 15, color = "black"), 
  axis.text.x = element_text(size = 15, color = "black",angle=0, vjust=0.6) , 
  axis.ticks.x = element_blank(),
  legend.title = element_text(size = 25, color = "black") , 
  legend.text = element_text(size = 20, color = "black"), legend.position = "none" , 
  plot.margin=unit(c(.5,1,.5,.5),"cm") , # add margins
  panel.border = element_rect(colour = "black", fill=NA, size=1)) + # add border
      guides(fill = guide_legend(title = "Cultivar")) 

NfixRhizo_data_plot

############################################################################
#### Arrange all these plots together 
library('ggpubr')

 ggarrange(SRLgia_plot, RootDiamplot,MBCplot, MBNplot, shannon_rhizo_plot, NfixRhizo_data_plot,
                    labels = c("A", "B", "C","D","E","F"),
                    ncol = 2, nrow = 3, 
                    widths = c(1,1,1,1), heights = c(1,1,1,1),
                      common.legend = TRUE, legend = "none",
            align = "hv")
 
 

```
## Figure 2A: Ordination - RhizoEndo 
```{r}

set.seed(2)
w.unifrac <- ordinate(rhizoendo, method = "NMDS", distance = "unifrac",weighted = TRUE)
# w.unifrac - stress = 0.08 - no convergence 

# extract out x and y dimensions from ordination 
w.unifrac.points <- as.data.frame(w.unifrac$points)
w.unifrac.points <- rownames_to_column(w.unifrac.points, "X.SampleID")

# merge with metadata 
w.unifrac.points.m <- merge(w.unifrac.points, sampledata_rhizoendo[,c(2,3,7,10,14)], by = "X.SampleID")


# Plot every point 
# reorder by ecotype 
w.unifrac.points.m$Variety <- factor(w.unifrac.points.m$Variety, levels = c("Alamo","Kanlow", "Cave-in-Rock", "Southlow"))

GreyWhite_colors <- c("gray33","gray65","white","gray91")

ggplot(data=w.unifrac.points.m,aes(x=MDS1, y= MDS2 , shape = SampleSite, fill = Variety))+
geom_point(aes(shape = SampleSite, color = Variety, fill = Variety),size = 5, stroke = 2)+
  theme(plot.title = element_text(size = 16)) + 
  labs(x = "NMDS1", y = "NMDS2") +
  scale_fill_manual(values = GreyWhite_colors, name = "Cultivar") +
  scale_colour_manual(values = c("black","black","black","black"))+
  scale_shape_manual(values = c(21,24))+
  theme(axis.title.x = element_text(size = 30), 
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 25, color = "black"),
  axis.text.y = element_text(size = 25, color = "black"), 
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 25), 
  legend.title = element_text(size = 40), 
  legend.position="right",  #remove legend  
  plot.margin=unit(c(.5,1,.5,.5),"cm"))  # add margins

# Opened in powerpoint and corrected legend 


```
## Figure 2B: Soil vs. Root Bacteria Phyla Stacked Bar plot 
```{r}
# Relative Abundance of Soil and Root PHyla (and Proteobacteria classes)

# Populate Proteobacteria classes in phylum level 
#https://github.com/joey711/phyloseq/issues/659
rhizoendo2 = rhizoendo #create a new phyloseq object, so the old isn't disrupted
plot.tax = as.data.frame(tax_table(rhizoendo2))
plot.tax = data.frame(lapply(plot.tax, as.character), stringsAsFactors = F)
plot.tax$Phylum[plot.tax$Phylum=="Proteobacteria"] = plot.tax$Class[plot.tax$Phylum=="Proteobacteria"] 
plot.tax[] = lapply(plot.tax, factor)
plot.tax.table = tax_table(plot.tax)
rownames(plot.tax.table) = rownames(tax_table(rhizoendo2))
colnames(plot.tax.table) = colnames(tax_table(rhizoendo2))

tax_table(rhizoendo2) = plot.tax.table
rhizoendo2

# glom by Phylum 
rhizoendo_phylum <- tax_glom(rhizoendo2, taxrank = "Phylum")
ntaxa(rhizoendo_phylum)


# relative abundance
rhizoendo_phylum.ra = transform_sample_counts(rhizoendo_phylum, function(x) 100 * x/sum(x))

# Put into a dataframe
rhizoendo_phylum.ra.df<- psmelt(rhizoendo_phylum.ra) #create dataframe from phyloseq


################# PLOT ###############################

# Summarize the mean relative abundance of Class by Treatment (to plot just by the treatment, not individual samples) 
### #NOTE: If you don't take the mean, it will just plot the total Rel.Abund across all samples - I think mean is more informative.
rhizoendo_phylum.ra.df_mean <- rhizoendo_phylum.ra.df %>% group_by(SampleSite,Phylum) %>%  summarize_at("Abundance",mean)


# For meaned data
rhizoendo_phylum.ra.df_mean$Phylum <- as.character(rhizoendo_phylum.ra.df_mean$Phylum) # convert to character
rhizoendo_phylum.ra.df_mean$Phylum[rhizoendo_phylum.ra.df_mean$Abundance < 1 ] <- "< 1% abund."

## first set colors so that you have consistant colors in your plots 
library("RColorBrewer")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getPalette = colorRampPalette(brewer.pal(12, "Paired")) # colorblind safe colors

# edit this to standardize colors by grouping factor 
  # prune is the filtered phyloseq object
speciesList = unique(rhizoendo_phylum.ra.df_mean$Phylum)
speciesPalette = getPalette(length(speciesList))
names(speciesPalette) = speciesList


# change names of samplesites to root- and soil-associated
rhizoendo_phylum.ra.df_mean$SampleSite <- factor(rhizoendo_phylum.ra.df_mean$SampleSite, levels = c("Rhizosphere","Endosphere"), 
                  labels = c("Soil", "Root"))

p = ggplot(rhizoendo_phylum.ra.df_mean, aes(x = SampleSite, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.85) + 
  scale_fill_manual(values = speciesPalette, name = "Bacterial Phyla") +
  facet_grid(.~SampleSite, drop = TRUE, space = "free", scales = "free") + #renames facet titles
  # use to plot multiple panels if you want to look at sample variation within your treatment.
  ylab("Mean Relative Abundance of Bacterial Phyla in Roots and Soils") + 
  xlab("Cultivar") +
  theme_classic() + 
theme(plot.title = element_text(size = 16), 
  axis.title.x = element_text(size = 20), 
  axis.title.y = element_text(size = 15), 
  axis.text.x = element_text(size = 20, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  legend.position="right", # remove legend "none"
  strip.text = element_text(size = 20)) #text in boxes for facet grid

p

```
## Figure 3: Soil bacteria & fungi centroid NMDS 
```{r}
# SOIL BACTERIA CENTROID PLOT 

# weighted unifrac distance matrix for soil samples 
set.seed(2)
w.unifrac <- ordinate(rhizo.rfy, method = "NMDS", distance = "unifrac", weighted = TRUE )
# stress = 0.18 

# extract out x and y dimensions from ordination 
w.unifrac.points <- as.data.frame(w.unifrac$points)
w.unifrac.points <- rownames_to_column(w.unifrac.points, "X.SampleID")

# merge with metadata 
w.unifrac.points.m <- merge(w.unifrac.points, metadata[,c(2,3,7,10)], by = "X.SampleID")


# calculate the mean and se of the points from the ordination by treatment 
w.unifrac.points.mean <- w.unifrac.points.m %>%
      group_by(Variety) %>% 
  summarise(Mean_MDS1 = mean(MDS1), 
            SD_MDS1 = sd(MDS1), 
            N = n(),
            SE_MDS1 = SD_MDS1/sqrt(N), 
            Mean_MDS2 = mean(MDS2), 
            SD_MDS2 = sd(MDS2), 
            N = n(),
            SE_MDS2 = SD_MDS2/sqrt(N), 
                        ) 


# add a column with number identifiers for each variety 
w.unifrac.points.mean$VarietyNumber <- w.unifrac.points.mean$Variety
w.unifrac.points.mean$VarietyNumber<- c("1","5","6","7","2","3","8","4","9","10","11","12")

# add ecotype back into df
w.unifrac.points.mean <- merge(w.unifrac.points.mean, metadata[,c(3,10,15)], by = "Variety")


 # reorder 
w.unifrac.points.mean$Variety <- factor(w.unifrac.points.mean$Variety, levels = c("Alamo","EG1101","EG1102","Kanlow", "Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer"))

# set colors 
BW_colors <- c("grey","grey","grey","grey","white","white","white","white","white","white","white","white")

# By Variety and Ecotype 
library(ggplot2)

bact_centroid <- ggplot(data=w.unifrac.points.mean,aes(x=Mean_MDS1, y= Mean_MDS2))+
  geom_errorbar(aes(ymax = Mean_MDS2 + SE_MDS2, ymin = Mean_MDS2 - SE_MDS2, width = .001))+ 
 geom_errorbarh(aes(xmax = Mean_MDS1 + SE_MDS1, xmin = Mean_MDS1 - SE_MDS1, height = .001)) +
  geom_point(aes(color = Ecotype, fill = Variety),size = 5, stroke = 2) + 
   #geom_text(aes(label = VarietyNumber),vjust = 0, nudge_y = 0.002, colour = "black", size = 6)+# puts numbers in middle without circle
 geom_label(aes(label=VarietyNumber, fontface = "bold", fill =Variety),   size = 10)+
      #scale_color_manual(values = c("goldenrod1","steelblue4"), name = "Ecotype") +
      scale_color_manual(values = c("black","black"), name = "Ecotype") +
       scale_fill_manual(values = BW_colors, name = "Cultivar") +
  scale_shape_manual(values = c(22,21)) + 
  labs(x = "NMDS1", y = "NMDS2") +
  theme(axis.title.x = element_text(size = 30), 
  axis.title.y = element_text(size = 30), 
  axis.text.x = element_text(size = 25, color = "black"), 
  axis.text.y = element_text(size = 25, color = "black"), 
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 25), 
  legend.title = element_text(size = 40), 
  legend.position="none",  #remove legend  
  plot.margin=unit(c(.5,1,.5,.5),"cm"))  # add margins



################################
### SOIL FUNGI CENTROID PLOT ###
set.seed(2)
bray.its <- ordinate(fungi.rfy, method = "NMDS", distance = "bray" )
# stress = 0.26 

# extract out x and y dimensions from ordination 
bray.its.points <- as.data.frame(bray.its$points)
library(tibble)
bray.its.points <- rownames_to_column(bray.its.points, "X.SampleID")

# merge with metadata 
bray.its.points.m <- merge(bray.its.points, sampledata_fungi.rfy[,c(2,3,7,10,15)], by = "X.SampleID")


# calculate the mean and se of the points from the ordination by treatment 
bray.its.points.mean <- bray.its.points.m %>%
      group_by(Variety) %>% 
  summarise(Mean_MDS1 = mean(MDS1), 
            SD_MDS1 = sd(MDS1), 
            N = n(),
            SE_MDS1 = SD_MDS1/sqrt(N), 
            Mean_MDS2 = mean(MDS2), 
            SD_MDS2 = sd(MDS2), 
            N = n(),
            SE_MDS2 = SD_MDS2/sqrt(N), 
                        ) 


# add a column with number identifiers for each variety 
bray.its.points.mean$VarietyNumber <- bray.its.points.mean$Variety
bray.its.points.mean$VarietyNumber<- c("1","5","6","7","2","3","8","4","9","10","11","12")

# add ecotype back into df
bray.its.points.mean <- merge(bray.its.points.mean,  sampledata_fungi.rfy[,c(3,7,10,15)], by = "Variety")


 # reorder 
bray.its.points.mean$Variety <- factor(bray.its.points.mean$Variety, levels = c("Alamo","EG1101","EG1102","Kanlow", "Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer"))

# set colors 
Eco_colors <- c("darkgoldenrod","gold","lightgoldenrod1","yellow","blue4","cornflowerblue","cyan3","deepskyblue4","lightseagreen","turquoise4","steelblue4","royalblue3")

BW_colors <- c("grey","grey","grey","grey","white","white","white","white","white","white","white","white")
BW_colors_outline <- c("black","black","black","black","black","black","black","black","black","black","black","black")


WarmCool_colors <- c("red","orangered3","firebrick3","orange1","mediumblue","darkcyan","springgreen4","greenyellow","aquamarine","cornflowerblue","deepskyblue","plum4")


# By Variety and Ecotype 
library(ggplot2)
fungi_centroid <-ggplot(data=bray.its.points.mean,aes(x=Mean_MDS1, y= Mean_MDS2))+
  geom_errorbar(aes(ymax = Mean_MDS2 + SE_MDS2, ymin = Mean_MDS2 - SE_MDS2, width = .001))+ 
 geom_errorbarh(aes(xmax = Mean_MDS1 + SE_MDS1, xmin = Mean_MDS1 - SE_MDS1, height = .001)) +
  geom_point(aes(color = Ecotype, fill = Variety),size = 5, stroke = 2) + 
   #geom_text(aes(label = VarietyNumber),vjust = 0, nudge_y = 0.002, colour = "black", size = 6)+# puts numbers in middle without circle
 geom_label(aes(label=VarietyNumber, fontface = "bold", fill =Variety),   size = 10)+
      #scale_color_manual(values = c("goldenrod1","steelblue4"), name = "Ecotype") +
      scale_color_manual(values = c("black","black"), name = "Ecotype") +
       scale_fill_manual(values = BW_colors, name = "Cultivar") +
  scale_shape_manual(values = c(22,21)) + 
  labs(x = "NMDS1", y = "NMDS2") +
  theme(axis.title.x = element_text(size = 30)) + 
  theme(axis.title.x = element_text(size = 30), 
  axis.title.y = element_text(size = 30), 
  axis.text.x = element_text(size = 25, color = "black"), 
  axis.text.y = element_text(size = 25, color = "black"), 
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 25), 
  legend.title = element_text(size = 40), 
  legend.position="none",  #remove legend  
  plot.margin=unit(c(.5,1,.5,.5),"cm"))  # add margins



############## PLOT TOGETHER ####################
# first run fungi.rfy centroid ordination code in fungi.rfy_OTUanalysis Rcode

#require(devEMF)
#emf(file = "bactfung_centroid.emf", width = 7, height = 7,bg = "white", fg = "black", pointsize = 12,family = "Arial", coordDPI = 300)

library('ggpubr')

 ggarrange(bact_centroid, fungi_centroid,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, 
                    widths = c(1,1,1,1), heights = c(1,1,1,1),
                      common.legend = TRUE, legend = "none",
            align = "hv")
#dev.off()

# opened in powerpoint and cleaned (e.g. made text box and arrows for hidden centroids )


```
## Figure 4: Variation in Abundant Bacterial Phyla plot 
```{r}

# First, I found phylum that drive variation in the communities with MVAbund 


# Populate Proteobacteria classes in phylum level 
#https://github.com/joey711/phyloseq/issues/659

rhizo.rfy2 = rhizo.rfy #create a new phyloseq object, so the old isn't disrupted
plot.tax = as.data.frame(tax_table(rhizo.rfy2))
plot.tax = data.frame(lapply(plot.tax, as.character), stringsAsFactors = F)
plot.tax$Phylum[plot.tax$Phylum=="Proteobacteria"] = plot.tax$Class[plot.tax$Phylum=="Proteobacteria"] 
plot.tax[] = lapply(plot.tax, factor)
plot.tax.table = tax_table(plot.tax)
rownames(plot.tax.table) = rownames(tax_table(rhizo.rfy2))
colnames(plot.tax.table) = colnames(tax_table(rhizo.rfy2))
identical(plot.tax.table[1:14590,3:7], tax_table(rhizo.rfy2)[1:14590,3:7]) #final check; should be true. Ensure the number of rows is adjusted to the dimensions of your tax table. 
tax_table(rhizo.rfy2) = plot.tax.table


# glom at phyla level (with proteo classes)
rhizo.rfy_phyla <- tax_glom(rhizo.rfy2, "Phylum")
ntaxa(rhizo.rfy_phyla) #43 

#########################################
# list of significantly different phyla (found with MVAbund function)

sig_phyla <- as.vector(c("Acidobacteria","Actinobacteria","Verrucomicrobia","Planctomycetes","Firmicutes","Bacteroidetes","Gemmatimonadetes", "Deltaproteobacteria", "Betaproteobacteria")) 


# what % of sequences do these dominant phyla make up ?
rhizo.rfy2_phy <- tax_glom(rhizo.rfy2, "Phylum")
taxnames <- as.data.frame(as(tax_table(rhizo.rfy2_phy),"matrix"))
taxa_names(rhizo.rfy2_phy) <- taxnames$Phylum
rhizo.rfy2_phy_sig <- prune_taxa(sig_phyla, rhizo.rfy2_phy)
ntaxa(rhizo.rfy2_phy_sig) # 9 

sum(sample_sums(rhizo.rfy2_phy_sig))/sum(sample_sums(rhizo.rfy2_phy))*100 # 74.3



####################################################

#### #PLOT 

# list of significantly different phyla 
sig_phyla <- as.vector(c("Acidobacteria","Actinobacteria","Verrucomicrobia","Planctomycetes","Firmicutes","Bacteroidetes","Gemmatimonadetes", "Deltaproteobacteria", "Betaproteobacteria")) 

# change otu names to phylumn names to subset by sig phylums
taxnames <- as.data.frame(as(tax_table(rhizo.rfy_phyla),"matrix"))
taxa_names(rhizo.rfy_phyla) <- taxnames$Phylum

# calculate relative abundance 
library(microbiome)
# convert all abundance data to rel.abund.
# use the phyloseq object that was manipulated to include proteo class in phyla
rhizo.rfy_phyla_ra <- transform_sample_counts(rhizo.rfy_phyla, function(x) 100 * x/sum(x))
sample_sums(rhizo.rfy_phyla_ra) # check - all should = 1 


# change otu names to phylumn names to subset by sig phylums
taxnames <- as.data.frame(as(tax_table(rhizo.rfy_phyla_ra),"matrix"))
taxa_names(rhizo.rfy_phyla_ra) <- taxnames$Phylum

rhizo.rfy_phy_sig <- prune_taxa(sig_phyla, rhizo.rfy_phyla_ra)
ntaxa(rhizo.rfy_phy_sig) # 9

# convert to dataframe 
rhizo.rfy_phy_sig.df <- psmelt(rhizo.rfy_phy_sig)

# Sum each phyla by sample 
rhizo.rfy_phy_sig.df_sums <- rhizo.rfy_phy_sig.df %>% group_by(X.SampleID, Variety, Ecotype, Phylum) %>%  summarize_at("Abundance",sum)


# Summarize the mean relative abundance of phylum by Treatment (to plot just by the treatment, not individual samples) 
### #NOTE: If you don't take the mean, it will just plot the total Rel.Abund across all samples

library(dplyr)
rhizo.rfy_phy_sig.df_mean <- rhizo.rfy_phy_sig.df_sums %>% group_by(Variety, Phylum) %>%    summarise(Mean_Abundance = mean(Abundance), 
            SD_Abundance = sd(Abundance), 
            N = n(),
            SE_Abundance = SD_Abundance/sqrt(N))

rhizo.rfy_phy_sig.df_eco_mean <- rhizo.rfy_phy_sig.df_sums %>% group_by(Ecotype, Phylum) %>%    summarise(Mean_Abundance = mean(Abundance), 
            SD_Abundance = sd(Abundance), 
            N = n(),
            SE_Abundance = SD_Abundance/sqrt(N))

# merge the upland and lowland means to the variety table
colnames(rhizo.rfy_phy_sig.df_eco_mean) <- colnames(rhizo.rfy_phy_sig.df_mean)
rhizo.rfy_phy_sig.df_mean_var_eco <- rbind(rhizo.rfy_phy_sig.df_mean,rhizo.rfy_phy_sig.df_eco_mean)


#### for horizontal plot ####
rhizo.rfy_phy_sig.df_mean_var_eco$Variety <- factor(rhizo.rfy_phy_sig.df_mean_var_eco$Variety, levels = c("Upland","Lowland","Trailblazer","Southlow","Shelter","NE28","EG2101","Dacotah","Cave-in-Rock","Blackwell","Kanlow","EG1102","EG1101","Alamo"))

WarmCool_colors <- c("dodgerblue4","red2","darkcyan","greenyellow","cornflowerblue","aquamarine","deepskyblue","plum","springgreen4","mediumpurple2","orange1","orangered3","goldenrod","red")

BW_colors_outline <- c("black","black","black","black","black","black","black","black","black","black","black","black","black","black")

 # change order of grouping factor by relative abundance 
rhizo.rfy_phy_sig.df_mean_var_eco$Phylum <- factor(rhizo.rfy_phy_sig.df_mean_var_eco$Phylum, levels = c("Acidobacteria","Verrucomicrobia", "Planctomycetes","Actinobacteria","Betaproteobacteria","Deltaproteobacteria","Bacteroidetes","Gemmatimonadetes","Firmicutes"))

# Stacked bar plot faceted by phyla  
p = ggplot(rhizo.rfy_phy_sig.df_mean_var_eco, aes(x = Variety, y = Mean_Abundance, fill = Variety, colour = Variety)) +
  geom_bar(stat = "identity", width = 1) + # width = 1 removes space in between bars 
  geom_errorbar(aes(ymin = Mean_Abundance-SE_Abundance, ymax = Mean_Abundance+SE_Abundance), width = .2, position = position_dodge(.9)) + 
  scale_fill_manual(values = WarmCool_colors, name = "Cultivar") +
  scale_colour_manual(values = BW_colors_outline, name = "Cultivar") +  facet_grid(.~Phylum, drop = TRUE, space = "fixed", scales = "free") +
  scale_x_discrete(breaks = c("Alamo","EG1101", "EG1102", "Kanlow","Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer", "Lowland","Upland"),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","L","U")) +
  ylab("Mean Relative Abundance (%)") + 
  xlab("Cultivar") +
  coord_flip()+
  theme_classic() + 
theme(plot.title = element_text(size = 16), 
  axis.title.x = element_text(size =10, color = "black"), 
  axis.title.y = element_text(size = 15), 
  axis.text.x = element_text(size =10, color = "black"),
  axis.text.y = element_text(size = 20, color = "black"), 
  axis.ticks.x = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  panel.spacing=unit(0, "lines"),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  legend.position="none", # remove legend "none"
  strip.text = element_text(size = 10)) #text in boxes for facet grid

p

# File saved as a metafile and modified in powerpoint (added significance letters from Anova tests, adjusted colors and font sizes)

```
## Figure S1 - Ordination with all samples
```{r}

### SOIL BACTERIA 

# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 
rhizo.rfy_noSMCout <- subset_samples(rhizo.rfy, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
nsamples(rhizo.rfy_noSMCout) #136 

set.seed(2)
w.unifrac <- ordinate(rhizo.rfy_noSMCout, method = "NMDS", distance = "unifrac", weighted = TRUE )
# stress = 0.18

# extract out x and y dimensions from ordination 
w.unifrac.points <- as.data.frame(w.unifrac$points)
w.unifrac.points <- rownames_to_column(w.unifrac.points, "X.SampleID")

# merge with metadata 
w.unifrac.points.m <- merge(w.unifrac.points, metadata[,c(2,3,7,10)], by = "X.SampleID")


# Plot every point 
# reorder by ecotype 
w.unifrac.points.m$Variety <- factor(w.unifrac.points.m$Variety, levels = c("Alamo","EG1101","EG1102","Kanlow", "Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer"))

WarmCool_colors <- c("red","goldenrod","orangered3","orange1","mediumblue","darkcyan","springgreen4","greenyellow","aquamarine","cornflowerblue","deepskyblue","plum4")


library(ggplot2)
ggplot(data=w.unifrac.points.m,aes(x=MDS1, y= MDS2 ,colour = Variety, shape = Ecotype))+
      geom_point(size = 5)+ 
  theme(plot.title = element_text(size = 16)) + 
  labs(x = "NMDS1", y = "NMDS2") +
 # geom_errorbar(aes(ymax = Mean_MDS2 + SE_MDS2, ymin = Mean_MDS2 - SE_MDS2, width = .001))+ 
 #geom_errorbarh(aes(xmax = Mean_MDS1 + SE_MDS1, xmin = Mean_MDS1 - SE_MDS1, height = .001)) +
  scale_color_manual(values = WarmCool_colors, name = "Cultivar")+
  scale_shape_manual(values = c(17,16))+
  theme(axis.title.x = element_text(size = 30), 
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 25, color = "black"),
  axis.text.y = element_text(size = 25, color = "black"), 
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 25), 
  legend.title = element_text(size = 40), 
  legend.position="none",  #remove legend  
  plot.margin=unit(c(.5,1,.5,.5),"cm"))  # add margins

## SOIL FUNGI 
fungi.rfy_noSMCout <- subset_samples(fungi.rfy, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
nsamples(fungi.rfy_noSMCout) #133
set.seed(2)
bray <- ordinate(fungi.rfy_noSMCout, method = "NMDS", distance = "bray" )
# stress = 0.26

# extract out x and y dimensions from ordination 
bray.points <- as.data.frame(bray$points)
bray.points <- rownames_to_column(bray.points, "X.SampleID")

# merge with metadata 
bray.points.m <- merge(bray.points, sampledata_fungi.rfy[,c(1,2,3,7,10,15)], by = "X.SampleID")


# Plot every point 
# plot
bray.points.m$Variety <- factor(bray.points.m$Variety, levels = c("Alamo","EG1101","EG1102","Kanlow", "Blackwell","Cave-in-Rock","Dacotah","EG2101","NE28","Shelter", "Southlow", "Trailblazer"))

WarmCool_colors <- c("red","goldenrod","orangered3","orange1","mediumblue","darkcyan","springgreen4","greenyellow","aquamarine","cornflowerblue","deepskyblue","plum4")


library(ggplot2)
ggplot(data=bray.points.m,aes(x=MDS1, y= MDS2 ,colour = Variety, shape = Ecotype))+
      geom_point(size = 5)+ 
  theme(plot.title = element_text(size = 16)) + 
  labs(x = "NMDS1", y = "NMDS2") +
 # geom_errorbar(aes(ymax = Mean_MDS2 + SE_MDS2, ymin = Mean_MDS2 - SE_MDS2, width = .001))+ 
 #geom_errorbarh(aes(xmax = Mean_MDS1 + SE_MDS1, xmin = Mean_MDS1 - SE_MDS1, height = .001)) +
  scale_color_manual(values = WarmCool_colors, name = "Cultivar")+
  scale_shape_manual(values = c(17,16))+
  theme(axis.title.x = element_text(size = 30), 
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 25, color = "black"),
  axis.text.y = element_text(size = 25, color = "black"), 
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.text = element_text(size = 25), 
  legend.title = element_text(size = 40), 
  legend.position="none",  #remove legend  
  plot.margin=unit(c(.5,1,.5,.5),"cm"))  # add margins




```
---
# Alpha Diversity
*effect of cultivar and ecotype on microbiome alpha diversity*

## Soil Fungi - 12 cultivars
```{r}

## create datatable with diversity data 
alpha_meas <- c("Observed", "Chao1", "Shannon", "InvSimpson")
# calculate diversity and merge with metadata, then alculate the averages for diversity 
require(tidyr)

its.diversity <- estimate_richness(fungi.rfy, measures = alpha_meas)
fungi.rfy_diversity<- merge(its.diversity, sampledata_fungi.rfy, by = "row.names")

# calculate Pielou's eveness (shannon/log(richness)
fungi.rfy_diversity$Evenness_Pielou <- fungi.rfy_diversity$Shannon/log(fungi.rfy_diversity$Observed)

```
### Shannon - Soil Fungi 
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(fungi.rfy_diversity$Shannon)
hist(fungi.rfy_diversity$Shannon) # Very left skewed

fungi_shannon = lmer(Shannon ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(fungi_shannon) 
shapiro.test(res) 
# shannon p <0.001 
# log(shannon) p <0.001 - can't normalize
# try to remove the low outlier at 1.96 to see if normality improves 
fungi.rfy_diversity$Shannon_OutlierRmv <-  fungi.rfy_diversity$Shannon
fungi.rfy_diversity$Shannon_OutlierRmv[fungi.rfy_diversity$Shannon_OutlierRmv < 2] <- NA
# data still not normal 


# 2) homoegeneity of variances
bartlett.test(log(Shannon) ~ Variety , data = fungi.rfy_diversity) 
# shannon p <0.001
# log(shannon), p <0.001

# check if anova is sig. with mixed effects despite normality
fungi_shannon = lmer(Shannon ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)

fungi_shannon_SMC = lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = fungi.rfy_diversity)

fungi_shannon_date <- lmer(Shannon ~ Date_Julian + Variety + (1|Block/Variety), data =   fungi.rfy_diversity)
# rank deficient modle 

AIC(fungi_shannon,fungi_shannon_date)
  #      df      AIC
#fungi_shannon      15 216.8383
#fungi_shannon_date 15 220.7301

# ANOVA without SMC as covariate
Anova(fungi_shannon, type = "III")
#Response: Shannon
#              Chisq Df Pr(>Chisq)    
#(Intercept) 909.310  1     <2e-16 ***
#Variety       8.081 11      0.706   

### Non-parametric stats for Shannon diversity (don't know how to add Block)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
fungi.shannon_Kruskal <- kruskal.test(Shannon ~ Variety , data = fungi.rfy_diversity)
fungi.shannon_Kruskal
#Kruskal-Wallis chi-squared = 7.2202, df = 11, p-value = 0.781
#use Dunn test for Kruskal post-hoc 
library(FSA)
fungi.shannon_dunn <- dunnTest(Shannon ~ Variety, data = fungi.rfy_diversity, method= "bh") # bh is the same as fdr 
fungi.shannon_dunn
# no sig.


### Non-parametric t-test (Wilcoxon Signed-ranked test)

wilcox.test(Shannon ~ Ecotype, data = fungi.rfy_diversity) 
#Wilcoxon rank sum test with continuity correction
#data:  Shannon by Ecotype
#W = 2177, p-value = 0.6163



#> model = lmer(Shannon ~ Ecotype + (1|Block/Plot), data = fungi.rfy_diversity)
#> emmeans(model,pairwise~Ecotype, adjust = "fdr")
#$emmeans
 #Ecotype emmean     SE   df lower.CL upper.CL
 #Lowland   4.42 0.0799 9.58     4.24     4.60
 #Upland    4.47 0.0657 4.57     4.30     4.65

```
### Observed richness - Fungal Soil 
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(fungi.rfy_diversity$Observed)
hist(fungi.rfy_diversity$Observed)

fungi_rich = lmer(Observed ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(fungi_rich) 
shapiro.test(res) 
# observed p = 0.27

# 2) homoegeneity of variances
bartlett.test(Observed ~ Variety , data = fungi.rfy_diversity)
# observed p =0.46 - equal variance

# Check if the model is stronger with SMC as covariate 
Obs_Variety_model <- lmer(Observed ~ Variety +(1|Block/Variety), data = fungi.rfy_diversity)
#boundary (singular) fit: see ?isSingular
Obs_Variety_model

# model with SMC as covariate
Obs_Variety_model_SMC = lmer(Observed ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = fungi.rfy_diversity)

# check AIC values
AIC(Obs_Variety_model,Obs_Variety_model_SMC)
#df      AIC
#Obs_Variety_model     15 1380.103
#Obs_Variety_model_SMC 16 1355.605 # YES 

anova(Obs_Variety_model_SMC, type = "III")

#Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Variety                         15198.2  1381.7    11 34.009  0.6289 0.7914
#SMC_perc_GravWaterCont_OutlierRmv  3082.3  3082.3     1 45.079  1.4030 0.2424



########################################
# By Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Observed_EcoBl_model <- lmer(Observed ~ Ecotype + (1|Block/Plot), data = fungi.rfy_diversity)

plot(Observed_EcoBl_model)


# model with SMC as covariate
Observed_EcoBl_model_SMC <- lmer(Observed ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Plot), data = fungi.rfy_diversity) 

# check AIC values
AIC(Observed_EcoBl_model,Observed_EcoBl_model_SMC)
#    df      AIC
#Observed_EcoBl_model     5 1447.150
#Observed_EcoBl_model_SMC  6 1424.182 #YES

anova(Observed_EcoBl_model_SMC, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
#Ecotype                          85.605  85.605     1  45.068  0.0383 0.8457
#SMC_                           228.086 228.086     1 111.641  0.1020 0.7500

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")

#$emmeans
# Ecotype emmean   SE    df lower.CL upper.CL
# Lowland    355 9.26 19.93      336      374
# Upland     351 6.70  6.82      336      367

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate   SE   df t.ratio p.value
# Lowland - Upland      3.5 11.4 41.6 0.306   0.7610 


```

### Pielou's Evenness - Soil Fungi 
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(fungi.rfy_diversity$Evenness_Pielou)
hist((fungi.rfy_diversity$Evenness_Pielou)^2) # Very left skewed

model = lmer(sqrt(Evenness_Pielou) ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model) 
shapiro.test(res) 
# log p <0.001
# sqrt p<0.001 


# mixed effects model despite non-normal data 
Fungi_Even = lmer(log(Evenness_Pielou) ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)

Anova(Fungi_Even, type = "III")
#Chisq Df Pr(>Chisq)    
#(Intercept) 59.8396  1  1.029e-14 ***
#Variety      9.6255 11     0.5644    


### Non-parametric stats for Evenness diversity (don't know how to add Block)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
fungi.even_Kruskal <- kruskal.test(Evenness_Pielou ~ Variety , data = fungi.rfy_diversity)
fungi.even_Kruskal
#data:  Evenness_Pielou by Variety
#Kruskal-Wallis chi-squared = 8.9806, df = 11, p-value = 0.6237 
library(FSA)
fungi.even_dunn <- dunnTest(Evenness_Pielou ~ Variety, data = fungi.rfy_diversity, method= "bh") # bh is the same as fdr 
fungi.even_dunn
# no sig.


### Non-parametric t-test (Wilcoxon Signed-ranked test)

wilcox.test(Evenness_Pielou ~ Ecotype, data = fungi.rfy_diversity) 
#Wilcoxon rank sum test with continuity correction
#data:  Evenness_Pielou by Ecotype
#W = 2109, p-value = 0.8516
#alternative hypothesis: true location shift is not equal to 0


 model = lmer(Evenness_Pielou ~ Ecotype + (1|Block/Plot), data = fungi.rfy_diversity)
 emmeans(model,pairwise~Ecotype, adjust = "fdr")
#Ecotype emmean     SE   df lower.CL upper.CL
# Lowland  0.753 0.0124 9.58    0.726    0.781
# Upland   0.764 0.0102 4.58    0.737    0.791

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate     SE   df t.ratio p.value
# Lowland - Upland  -0.0107 0.0128 40.8 -0.835  0.4088 
```
## Soil Bacteria - 12 cultivars
```{r}
###Resource: https://rpubs.com/dillmcfarlan/R_microbiota_BRC2017

## create datatable with diversity data 
alpha_meas <- c("Observed", "Chao1", "Shannon", "InvSimpson")
# calculate diversity and merge with metadata, then alculate the averages for diversity 
require(tidyr)

r.diversity <- estimate_richness(rhizo.rfy, measures = alpha_meas)
rhizo.rfy_diversity<- merge(r.diversity, sampledata_rhizo.rfy, by = "row.names")

# calculate Pielou's eveness (shannon/log(richness)
rhizo.rfy_diversity$Evenness_Pielou <- rhizo.rfy_diversity$Shannon/log(rhizo.rfy_diversity$Observed)


################################
# Remove Dacotah Cultivar and see if trends hold# 
################################
rhizo.rfy_NoDac <- subset_samples(rhizo.rfy, Variety != "Dacotah")
nsamples(rhizo.rfy_NoDac) #126 - good removed Dacotah samples 

rhizo.rfy_NoDac_div <- estimate_richness(rhizo.rfy_NoDac, measures = alpha_meas)
rhizo.rfy_NoDac_diversity<- merge(rhizo.rfy_NoDac_div,sample_data(rhizo.rfy_NoDac), by = "row.names")

# calculate Pielou's eveness (shannon/log(richness)
rhizo.rfy_NoDac_diversity$Evenness_Pielou <- rhizo.rfy_NoDac_diversity$Shannon/log(rhizo.rfy_NoDac_diversity$Observed)

```

###  Shannon Diversity - Soil Bacteria
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# normal residuals?
ggqqplot(rhizo.rfy_diversity$Shannon)
hist(rhizo.rfy_diversity$Shannon)

# shapiro test - null hypothesis: residuals are normally distributed 
model = lmer(Shannon ~ Variety +  (1|Block/Variety) , data =   rhizo.rfy_diversity)
res = resid(model)
shapiro.test(res) 
# shannon ~ Variety, p = 0.6

# homoegeneity of variances
bartlett.test(Shannon ~ Variety, data = rhizo.rfy_diversity) 
# Shannon ~ Variety, p = 0.9


# Does Shannon differ by Variety? (Shannon is normally distributed)
require(lme4)
require(lmerTest)

# first check if SMC as covariate improves model fit 
Shannon_VarietyBl_model <- lmer(Shannon ~ Variety + (1|Block/Variety), data = rhizo.rfy_diversity)

Shannon_VarietyBl_model_SMC <- lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data =   rhizo.rfy_diversity)

Shannon_VarietyBl_model_date <- lmer(Shannon ~ Date_Julian + Variety + (1|Block/Variety), data =   rhizo.rfy_diversity)
# rank deficient model 

Shannon_VarietyBl_model_NO3 <- lmer(Shannon ~ Variety +NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+ (1|Block/Variety), data =   rhizo.rfy_diversity)
# rank deficient model 

AIC(Shannon_VarietyBl_model,  Shannon_VarietyBl_model_NO3)
#  df       AIC
#Shannon_VarietyBl_model      15 -117.5242
#Shannon_VarietyBl_model_date 15 -113.6324
#Shannon_VarietyBl_model_SMC  16 -106.3809 
#Shannon_VarietyBl_model_NO3 16 -110.3907


anova(Shannon_VarietyBl_model_date, type = "III")
#Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
#Variety 0.57019 0.051836    11 32.309  4.4052 0.0004719 ***


# pairwise test (gives you p.values for every combination)
emmeans(Shannon_VarietyBl_model, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Shannon_VarietyBl_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD
 #Variety      lsmean     SE   df lower.CL upper.CL .group
#  Kanlow         6.36 0.0590 8.63     6.13     6.59  a    
# EG1102         6.36 0.0590 8.63     6.13     6.59  a    
# NE28           6.38 0.0590 8.63     6.15     6.61  a    
# EG2101         6.45 0.0590 8.63     6.22     6.67  ab   
# Alamo          6.45 0.0590 8.63     6.22     6.67  ab   
# Southlow       6.45 0.0590 8.63     6.23     6.68  ab   
# Shelter        6.46 0.0608 9.71     6.23     6.68  ab   
# Cave-in-Rock   6.47 0.0608 9.71     6.25     6.70  ab   
# EG1101         6.48 0.0590 8.63     6.25     6.71  ab   
# Trailblazer    6.52 0.0608 9.71     6.29     6.74  abc  
# Blackwell      6.54 0.0590 8.63     6.31     6.77   bc  
# Dacotah        6.67 0.0590 8.63     6.44     6.89    c  


################################
# Remove Dacotah and see if trend holds# 
################################

# confirm normality 
ggqqplot(rhizo.rfy_NoDac_diversity$Shannon)
hist(rhizo.rfy_NoDac_diversity$Shannon)

# shapiro test - null hypothesis: residuals are normally distributed 
model = lmer(Shannon ~ Variety +  (1|Block/Variety), rhizo.rfy_NoDac_diversity)
plot(model)
res = resid(model)
shapiro.test(res) 
# shannon ~ Variety, p = 0.58

# homoegeneity of variances
bartlett.test(Shannon ~ Variety, data = rhizo.rfy_NoDac_diversity) 
# Shannon ~ Variety, p = 0.8


# Does Shannon differ by Variety? (Shannon is normally distributed)
require(lme4)
require(lmerTest)
Shannon_VarietyBl_model_noDac <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_NoDac_diversity)

Shannon_VarietyBl_model_noDac_SMC <- lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = rhizo.rfy_NoDac_diversity)

AIC(Shannon_VarietyBl_model_noDac, Shannon_VarietyBl_model_noDac_SMC)
# df        AIC
#Shannon_VarietyBl_model_noDac     14 -102.15589
#Shannon_VarietyBl_model_noDac_SMC 15  -92.79755#

anova(Shannon_VarietyBl_model_noDac, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Variety 0.26685 0.026685    10 29.325   2.144 0.05295 .

# pairwise test (gives you p.values for every combination)
emmeans(Shannon_VarietyBl_model_noDac, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Shannon_VarietyBl_model_noDac,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

#Variety      lsmean     SE    df lower.CL upper.CL .group
# Kanlow         6.36 0.0560  9.67     6.15     6.57  a    
# EG1102         6.36 0.0560  9.67     6.16     6.57  a    
# NE28           6.38 0.0560  9.67     6.17     6.58  a    
# EG2101         6.45 0.0560  9.67     6.24     6.65  a    
## Alamo          6.45 0.0560  9.67     6.24     6.65  a    
# Southlow       6.45 0.0560  9.67     6.25     6.66  a    
## Shelter        6.46 0.0581 11.07     6.25     6.66  a    
# Cave-in-Rock   6.47 0.0581 11.07     6.27     6.68  a    
# EG1101         6.48 0.0560  9.67     6.27     6.68  a    
## Trailblazer    6.52 0.0581 11.07     6.31     6.72  a    
# Blackwell      6.54 0.0560  9.67     6.33     6.75  a  

########################################
# Does Shannon Diversity differ by Ecotype? Upland > Lowland (but driven by Dacotah)
#########################################
require(lme4)
require(lmerTest)
Shannon_EcoBl_model <- lmer(Shannon ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_diversity)

Shannon_EcoBl_model_SMC <- lmer(Shannon ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Shannon_EcoBl_model, Shannon_EcoBl_model_SMC)
# df      AIC
#Shannon_EcoBl_model      5 -148.263
#Shannon_EcoBl_model_SMC  6 -143.107

plot(Shannon_EcoBl_model)

anova(Shannon_EcoBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
  #        Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Ecotype 0.071919 0.071919     1 42.443  6.1455 0.01723 *

# pairise test
library(emmeans)
emmeans(Shannon_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")
#contrast           estimate         SE    df t.ratio p.value
 #Lowland - Upland -0.0796252 0.03212222 42.16  -2.479  0.0173
# Upland > diversity than lowland

################ Is this driven by Dacotah? 
require(lme4)
require(lmerTest)
Shannon_EcoBl_noDac_model <- lmer(Shannon ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

anova(Shannon_EcoBl_noDac_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
   #      Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#Ecotype 0.04606 0.04606     1 38.146  3.7205 0.06121

emmeans(Shannon_EcoBl_noDac_model, pairwise~Ecotype, p.adjust = "fdr")
#Ecotype emmean     SE   df lower.CL upper.CL
# Lowland   6.41 0.0449 4.16     6.29     6.53
# Upland    6.47 0.0426 3.40     6.34     6.59

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate    SE df t.ratio p.value
# Lowland - Upland   -0.054 0.028 38 -1.928  0.0613 


##################################################
# Does Shannon Div. differ among cultivars sampled on the same date? (Julian dates: 180, 195, 202, 209)
##################################################

# subset diversity data by date 
rhizo.rfy_diversity_180 <- filter(rhizo.rfy_diversity, Date_Julian == "180")
rhizo.rfy_diversity_195 <- filter(rhizo.rfy_diversity, Date_Julian == "195")
rhizo.rfy_diversity_202 <- filter(rhizo.rfy_diversity, Date_Julian == "202")
rhizo.rfy_diversity_209 <- filter(rhizo.rfy_diversity, Date_Julian == "209")


Shannon_VarietyBl_model_180 <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_diversity_180)
anova(Shannon_VarietyBl_model_180, type = "III") # p = 0.002

Shannon_VarietyBl_model_195 <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_diversity_195)
anova(Shannon_VarietyBl_model_195, type = "III") # p = 0.808 

Shannon_VarietyBl_model_202 <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_diversity_202)
# model failed to converge 
anova(Shannon_VarietyBl_model_202, type = "III") # p = 0.04 

Shannon_VarietyBl_model_209 <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_diversity_209)
anova(Shannon_VarietyBl_model_209, type = "III") # p = 0.206 


```

### Observed richness - Soil Bacteria
```{r}
# Does Richness differ by Variety or Ecotype? 

# test for normal residuals and equal variances (ANOVA assumptions)
# normal residuals
ggqqplot(rhizo.rfy_diversity$Observed)
hist(rhizo.rfy_diversity$Observed)

model = lmer(Observed ~ Variety + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizo.rfy_diversity)
res = resid(model)
shapiro.test(res)


# 2) homoegeneity of variances
bartlett.test(Observed ~ Variety, data = rhizo.rfy_diversity) 
# Richness by Var, p  0.718

## Is model fit improved with SMC as covariate?
Observed_VarietyBl_model <- lmer(Observed ~ Variety + (1|Block/Variety), data = rhizo.rfy_diversity)

Observed_VarietyBl_SMC_model <- lmer(Observed ~ SMC_perc_GravWaterCont_OutlierRmv + Variety  +  (1|Block/Variety), data = rhizo.rfy_diversity)


Observed_VarietyBl_date <- lmer(Observed ~ Date_Julian + Variety  +  (1|Block/Variety), data = rhizo.rfy_diversity)



AIC(Observed_VarietyBl_date, Observed_VarietyBl_SMC_model)
 # df      AIC
#Observed_VarietyBl_model     15 1563.189
#Observed_VarietyBl_SMC_model 16 1535.468 # YES 


anova(Observed_VarietyBl_SMC_model, type = "III")
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
#Variety                       184306 16755.1    11  33.352  2.1775 0.04129
#SMC_perc_OutlierRmv           5003  5002.7     1 121.594  0.6501 0.42163



emmeans(Observed_VarietyBl_SMC_model, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Observed_VarietyBl_SMC_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD
#Variety      lsmean   SE    df lower.CL upper.CL .group
# EG1102         1339 44.9  8.51     1165     1513  a    
# NE28           1406 44.4  8.19     1231     1581  ab   
# EG2101         1410 45.3  8.82     1237     1584  ab   
# Kanlow         1421 45.7  9.01     1247     1595  ab   
# EG1101         1441 44.4  8.20     1266     1615  ab   
# Shelter        1442 46.2  9.48     1269     1615  ab   
# Southlow       1443 44.5  8.24     1269     1617  ab   
# Trailblazer    1468 46.1  9.45     1295     1641  ab   
# Cave-in-Rock   1477 45.9  9.27     1304     1650  ab   
# Alamo          1478 44.5  8.25     1304     1653  ab   
# Blackwell      1479 44.8  8.41     1304     1653  ab   
# Dacotah        1527 47.4 10.34     1353     1700   b   


####################################
# Is it driven by dakota? 

## By Variety? 
Obs_Variety_noDak_SMC <- lmer(Observed ~ Variety +SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizo.rfy_NoDac_diversity)


anova(Obs_Variety_noDak_SMC, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
#Variety                    161604 16160.4    10  30.301  2.0842 0.05839 .
#SMC_perc                       2055  2054.5     1 110.144  0.2650 0.60775 
                              
                              
emmeans(Obs_Variety_noDak_SMC, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Obs_Variety_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

#Variety      lsmean   SE    df lower.CL upper.CL .group
# EG1102         1333 39.1  9.21     1187     1479  a    
#NE28           1403 39.1  9.21     1257     1549  ab   
# Kanlow         1412 39.1  9.21     1266     1558  ab   
# EG2101         1425 39.1  9.21     1279     1571  ab   
# Southlow       1440 39.1  9.21     1294     1585  ab   
# EG1101         1444 39.1  9.21     1298     1589  ab   
# Shelter        1445 40.9 10.89     1299     1591  ab   
# Trailblazer    1465 40.9 10.89     1319     1611   b   
# Alamo          1475 39.1  9.21     1329     1621   b   
# Cave-in-Rock   1479 40.9 10.89     1333     1624   b   
# Blackwell      1484 39.1  9.21     1338     1630   b   

########################################
# By Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Observed_EcoBl_model <- lmer(Observed ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_diversity)

Observed_EcoBl_model_SMC <- lmer(Observed ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Observed_EcoBl_model_SMC, Observed_EcoBl_model)
#df      AIC
#Observed_EcoBl_model_SMC  6 1626.683 # YES
#Observed_EcoBl_model      5 1657.701


plot(Observed_EcoBl_model_SMC)

anova(Observed_EcoBl_model_SMC, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
      #                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)  
#Ecotype                        13725   13725     1  42.985  1.7703 0.1904  
#SMC_perc_GravWaterCont        32438   32438     1 119.689  4.1841 0.0430 *

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")
#Ecotype emmean   SE   df lower.CL upper.CL
# Lowland   1416 35.5 4.23     1320     1512
# Upland    1460 33.3 3.31     1359     1561

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate   SE   df t.ratio p.value
# Lowland - Upland    -43.9 21.5 41.8 -2.043  0.0474 

############### 
# Without dacotah 
###############
Observed_EcoBl_model <- lmer(Observed ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

plot(Observed_EcoBl_model)

anova(Observed_EcoBl_model, type = "III")
 #   Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Ecotype  19544   19544     1 39.159  2.4665 0.1243

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")

#Ecotype emmean   SE   df lower.CL upper.CL
# Lowland   1416 31.7 4.21     1330     1502
# Upland    1448 30.1 3.42     1358     1537

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate   SE   df t.ratio p.value
# Lowland - Upland    -31.7 20.2 37.8 -1.570  0.1247 

```

### Pielou's Evenness - Soil Bacteria
```{r}
# Does Evenness_Pielou differ by Variety or Ecotype? 

# test for normal residuals and equal variances (ANOVA assumptions)
# normal residuals
ggqqplot(rhizo.rfy_diversity$Evenness_Pielou)
hist(rhizo.rfy_diversity$Evenness_Pielou)

model = lmer(Evenness_Pielou ~ Variety + (1|Block/Variety), data = rhizo.rfy_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# Evenness_Pielou by Var, p=0.1

# 2) homoegeneity of variances
bartlett.test(Evenness_Pielou ~ Variety, data = rhizo.rfy_diversity) 
# Evenness_Pielou by Var, p = 0.7156

## Analysis


## Will SMC as covariate improve model fit? 
Even_Variety_model <- lmer(Evenness_Pielou ~ Variety +(1|Block/Variety), data = rhizo.rfy_diversity)

Even_Variety_model_SMC <- lmer(Evenness_Pielou ~ Variety + SMC_perc_GravWaterCont_OutlierRmv+(1|Block/Variety), data = rhizo.rfy_diversity)

Even_Variety_model_date <- lmer(Evenness_Pielou ~ Variety + Date_Julian+(1|Block/Variety), data = rhizo.rfy_diversity)

AIC(Even_Variety_model, Even_Variety_model_date)
# df       AIC
#Even_Variety_model    15 -723.4241
#Even_Variety_mode_SMC 16 -696.9520 # NO



anova(Even_Variety_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)    
#Variety 0.0049035 0.00044577    11 32.036  4.7091 0.0002802 ***

emmeans(Even_Variety_model, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Even_Variety_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

#Variety      lsmean      SE   df lower.CL upper.CL .group
# Kanlow        0.877 0.00517 10.1    0.858    0.896  a    
# NE28          0.880 0.00517 10.1    0.861    0.899  ab   
# Alamo         0.884 0.00517 10.1    0.865    0.903  abc  
# EG1102        0.884 0.00517 10.1    0.865    0.903  abc  
# Cave-in-Rock  0.887 0.00534 11.5    0.868    0.906  abc  
# Shelter       0.888 0.00534 11.5    0.869    0.907  abc  
# Southlow      0.888 0.00517 10.1    0.869    0.907  abc  
# EG2101        0.888 0.00517 10.1    0.869    0.907  abc  
# EG1101        0.891 0.00517 10.1    0.872    0.910  abc  
# Trailblazer   0.894 0.00534 11.5    0.875    0.913   bcd 
# Blackwell     0.896 0.00517 10.1    0.877    0.915    cd 
# Dacotah       0.908 0.00517 10.1    0.889    0.927     d 



########### Is it driven by dakota? 

## By Variety? 
Even_Variety_model <- lmer(Evenness_Pielou ~ Variety +(1|Block/Variety), data = rhizo.rfy_NoDac_diversity)
Even_Variety_model

anova(Even_Variety_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
 #          Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#Variety 0.0020957 0.00020956    10  29.2  2.0648 0.06217 .

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Even_Variety_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

# Variety      lsmean      SE   df lower.CL upper.CL .group
# Kanlow        0.877 0.00507 11.5    0.859    0.895  a    
# NE28          0.880 0.00507 11.5    0.862    0.898  a    
# Alamo         0.884 0.00507 11.5    0.866    0.902  a    
# EG1102        0.884 0.00507 11.5    0.866    0.902  a    
## Cave-in-Rock  0.887 0.00526 13.2    0.869    0.905  a    
# Shelter       0.888 0.00526 13.2    0.870    0.906  a    
# Southlow      0.888 0.00507 11.5    0.870    0.906  a    
# EG2101        0.888 0.00507 11.5    0.870    0.906  a    
# EG1101        0.891 0.00507 11.5    0.873    0.909  a    
# Trailblazer   0.894 0.00526 13.2    0.877    0.912  a    
# Blackwell     0.896 0.00507 11.5    0.878    0.914  a   



########################################
# By Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Even_Eco_model <- lmer(Evenness_Pielou ~ Ecotype +(1|Block/Plot), data = rhizo.rfy_diversity)

Even_Eco_model_SMC <- lmer(Evenness_Pielou ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Even_Eco_model, Even_Eco_model_SMC)
#df       AIC
#Even_Eco_model      5 -799.3514
#Even_Eco_model_SMC  6 -779.0489

plot(Even_Eco_model)

anova(Even_Eco_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq    Mean Sq NumDF  DenDF F value  Pr(>F)  
#Ecotype 0.00050872 0.00050872     1 42.288  5.4062 0.02494 *

# pairise test
library(emmeans)
emmeans(Even_Eco_model, pairwise ~ Ecotype, adjust = "fdr")

#Ecotype emmean      SE   df lower.CL upper.CL
# Lowland  0.884 0.00427 4.97    0.873    0.895
# Upland   0.891 0.00390 3.48    0.880    0.903

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate      SE   df t.ratio p.value
# Lowland - Upland -0.00713 0.00306 42.2 -2.325  0.0250 



## without dacotah 
Even_Eco_model <- lmer(Evenness_Pielou ~ Ecotype +(1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

plot(Even_Eco_model)

anova(Even_Eco_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq    Mean Sq NumDF  DenDF F value  Pr(>F)  
#Ecotype 0.00029972 0.00029972     1 37.874  2.9646 0.09326 .

#emmeans(Even_Eco_model, pairwise ~ Ecotype, adjust = "fdr")
# Ecotype emmean      SE   df lower.CL upper.CL
# Lowland  0.884 0.00392 4.47    0.874    0.895
# Upland   0.889 0.00369 3.49    0.878    0.900#

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate      SE   df t.ratio p.value
# Lowland - Upland -0.00461 0.00268 38.1 -1.722  0.0933 

```

## Soil and Root (Rhizoendo) - 4 cultivars
```{r}
## create datatable with diversity data 
alpha_meas <- c("Observed", "Chao1", "Shannon", "InvSimpson")
# calculate diversity and merge with metadata, then alculate the averages for diversity 
require(tidyr)

#Bact_EndoRhizo only 
diversity <- estimate_richness(rhizoendo, measures = alpha_meas)
rhizoendo_diversity <- merge(diversity, sampledata_rhizoendo, by = "row.names")

# calculate Pielou's eveness (shannon/log(richness)
rhizoendo_diversity$Evenness_Pielou <- rhizoendo_diversity$Shannon/log(rhizoendo_diversity$Observed)
dim(rhizoendo_diversity)
```
###Shannon - RhizoEndo
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(rhizoendo_diversity$Shannon)
hist(rhizoendo_diversity$Shannon)

model = lmer((Shannon)^2 ~ Variety*SampleSite + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizoendo_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# shannon p <0.001
# log_shannon, p<0.001
# sq shannon, p = 0.12

# 2) homoegeneity of variances
bartlett.test(Shannon^2 ~ Variety , data = rhizoendo_diversity)
# shannon p =0.6256 - equal variance
# sq shannon p = 0.71407

# check if SMC as covariate improves model fit 
rhizoendo_shannon = lmer((Shannon)^2 ~ Variety + (1|Block/Variety), data = rhizoendo_diversity)
# singular fit 

rhizoendo_shannon_SMC = lmer((Shannon)^2 ~ Variety +  SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

AIC(rhizoendo_shannon,rhizoendo_shannon_SMC)
  #    df      AIC
#rhizoendo_shannon      7 690.9829 # NO 
#rhizoendo_shannon_SMC  8 692.0643


# Anova - Does Shannon differ by Variety? 
require(lme4)
require(lmerTest)
Shannon_Variety_model <- lmer(Shannon^2 ~ Variety*SampleSite + (1|Block/Variety), data = rhizoendo_diversity)
Shannon_Variety_model
anova(Shannon_Variety_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF  DenDF   F value Pr(>F)    
#Variety               54.8    18.3     3  8.478    1.7311 0.2338    
#SampleSite         12421.5 12421.5     1 69.162 1178.0998 <2e-16 ***
#Variety:SampleSite    55.1    18.4     3 69.038    1.7414 0.1666 

emmeans(Shannon_Variety_model, pairwise ~ SampleSite, adjust = "fdr")
#SampleSite  emmean    SE   df lower.CL upper.CL
# Endosphere    14.2 0.602 7.46     12.8     15.6
# Rhizosphere   38.1 0.579 6.55     36.7     39.5

## Plot 
# consistant colors 
library("RColorBrewer")
colourCount = length(unique(rhizoendo_diversity$Variety))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# reorder cultivars 
rhizoendo_diversity$Variety <- factor(rhizoendo_diversity$Variety, levels = c("Alamo", "Kanlow", "Cave-in-Rock", "Southlow"))

ggplot(data = rhizoendo_diversity, aes(x = SampleSite, y = Shannon, fill = getPalette(colourCount))) + 
  geom_boxplot(aes(factor(SampleSite), fill = factor(Variety))) + theme_classic() + 
  scale_fill_manual(values = getPalette(colourCount))+
  labs(x = "SampleSite", y = "Shannon Diversity") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  theme(axis.title = element_text(size = 25)) + 
  theme(axis.text.y = element_text(size = 25)) +
  theme(legend.title = element_text(size = 25)) + 
  theme(legend.text = element_text(size = 20)) + 
  theme(axis.title.x = element_text(size = 25)) + 
  theme(axis.text.x = element_text(size = 20)) + 
  guides(fill = guide_legend(title = "Cultivar")) 




```

### Observed richness - RhizoEndo
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(rhizoendo_diversity$Observed)
hist(rhizoendo_diversity$Observed)

model = lmer(Observed ~ Variety*SampleSite + SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Variety), data = rhizoendo_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# observed p =0.66

# 2) homoegeneity of variances
bartlett.test(Observed ~ Variety , data = rhizoendo_diversity) 
# observed p =0.9914 - equal variance


# Anova - Does Richness differ by Variety?  AND does SMC include model fit?
rhizoendo_obs = lmer(Observed ~ Variety*SampleSite + (1|Block/Variety), data = rhizoendo_diversity)


rhizoendo_obs_SMC = lmer(Observed ~ Variety*SampleSite +  SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

AIC(rhizoendo_obs,rhizoendo_obs_SMC)
#    df      AIC
#rhizoendo_obs     11 925.4898
#rhizoendo_obs_SMC 12 921.4109 #YES

require(lme4)
require(lmerTest)

anova(rhizoendo_obs_SMC, type = "III")
#  Sum Sq  Mean Sq NumDF  DenDF   F value Pr(>F)    
#Variety                        15796     5265     3  9.011    1.5736 0.2625 #SampleSite              11155060 11155060     1 67.541 3333.8986 <2e-16 ***
#SMC_perc_GravWaterCont_    7192     7192     1 42.333    2.1495 0.1500    
#Variety:SampleSite          6298     2099     3 67.226    0.6274 0.5998 

emmeans(rhizoendo_obs_SMC, pairwise ~ SampleSite, adjust = "fdr")
#SampleSite  emmean   SE   df lower.CL upper.CL
# Endosphere     171 12.7 4.94      138      203
 #Rhizosphere    889 12.4 4.47      856      922



## Analysis
# does Observed difer by ecotype? 
Obs_EcoSampleSite_model <- lmer(Observed ~ Ecotype*SampleSite + (1|Block/Plot), data = rhizoendo_diversity)
anova(Obs_EcoSampleSite_model)
#Type III Analysis of Variance Table with Satterthwaite's method
 #                    Sum Sq  Mean Sq NumDF  DenDF   F value Pr(>F)    
#Ecotype                 513      513     1 11.399    0.1604 0.6962    
#SampleSite         11222706 11222706     1 71.359 3509.8013 <2e-16 ***
#Ecotype:SampleSite     3386     3386     1 71.359    1.0588 0.3070   


```

### Pielou's Evenness
```{r}

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot((rhizoendo_diversity$Evenness_Pielou))
hist(rhizoendo_diversity$Evenness_Pielou)

model = lmer((Evenness_Pielou) ~ SampleSite + (1|Block/Variety), data = rhizoendo_diversity)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# Evenness_Pielou p <0.001 
# log_Evenness_Pielou, p <0.001
# sq Evenness_Pielou, p = 0.09

# 2) homoegeneity of variances
bartlett.test(Evenness_Pielou^2 ~ Variety , data = rhizoendo_diversity)
# Evenness_Pielou p =0.6256 - equal variance
# sq Evenness_Pielou p = 0.71407

emmeans(model, pairwise~SampleSite, adjust = "fdr")
#SampleSite  emmean     SE   df lower.CL upper.CL
 #Endosphere   0.729 0.0103 12.3    0.707    0.752
 #Rhizosphere  0.909 0.0098 10.6    0.887    0.931

# Test with proper mixed-effects model despite non-normal data 
rhizoendo_even = lmer(Evenness_Pielou^2 ~ Variety*SampleSite+(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

Anova(rhizoendo_even, Type = "III")
#Response: Evenness_Pielou^2
#                      Chisq Df Pr(>Chisq)    
#Variety              4.9359  3     0.1765    
#SampleSite         261.4840  1     <2e-16 ***
#Variety:SampleSite   6.1541  3     0.1044 

### Non-parametric stats for Shannon diversity (don't know how to add Block)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
rhizoendo.even <- kruskal.test(Evenness_Pielou~ SampleSite, rhizoendo_diversity)
rhizoendo.even
#Kruskal-Wallis rank sum test
#data:  Evenness_Pielou by Variety
#Kruskal-Wallis chi-squared = 64.989, df = 1, p-value = 7.532e-16

wilcox.test(Evenness_Pielou ~ SampleSite, data = rhizoendo_diversity) 
#Wilcoxon rank sum test with continuity correction
#data:  Shannon by Ecotype
#W = 2177, p-value = 0.6163


```
---
# BetaDiversity 
*effect of cultivar and ecotype on microbiome structure*
## Export Distance Matrices for Primer analysis
```{r}
# Export distance matrices to analyse in Primer 
  #Anderson, M.L., R.N. Gorley, and K.R. Clarke. 2008. “PERMANOVA+ for PRIMER.”

## BACTERIA 
#Weighted and unweigthed unifrac distance

# Soil
# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 

rhizo.rfy_noSMCout <- subset_samples(rhizo.rfy, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
nsamples(rhizo.rfy_noSMCout)#136


set.seed(2)
 rhizo.rfy_uni<- phyloseq::distance(rhizo.rfy_noSMCout, "unifrac")
set.seed(2)
rhizo.rfy_w.uni <- phyloseq::distance(rhizo.rfy_noSMCout, "wunifrac" )

# Root & Soil (combined)

# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 

## Remove samples B10 and C8 because soil moisture values were outliers for these. 
rhizoendo_noSMCout <- subset_samples(rhizoendo, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
nsamples(rhizoendo_noSMCout)# 88

# w.unifrac distance matrix
set.seed(2)
rhizoendo_noSMCout_w.uni <- phyloseq::distance(rhizoendo_noSMCout, "wunifrac" )

#write.csv(as.matrix(rhizoendo_noSMCout_w.uni), "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_NoSMC_Wuni_primer.csv")

# extract metadata in same order as distance matrix for primer analysis
names <- as.data.frame(row.names(as.matrix(rhizoendo_noSMCout_w.uni)))
colnames(names)[1] <- "X.SampleID"
names$id  <- 1:nrow(names) # create an id for row order so we can keep the same order as the distance matrix
rhizoendo_noSMCout_w.uni_primer_meta <- merge(names, sampledata_rhizoendo, by = "X.SampleID")
rhizoendo_noSMCout_w.uni_primer_meta <- rhizoendo_noSMCout_w.uni_primer_meta[order(rhizoendo_noSMCout_w.uni_primer_meta$id), ]
dim(rhizoendo_noSMCout_w.uni_primer_meta)


#write.csv(rhizoendo_noSMCout_w.uni_primer_meta, "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_NoSMC_Wuni_primer_metadata.csv")

#############
# Root and soil subset of 4 cultivars

rhizoendo_noSMCout_root <- subset_samples(rhizoendo_noSMCout, SampleSite == "Endosphere")
rhizoendo_noSMCout_soil <- subset_samples(rhizoendo_noSMCout, SampleSite == "Rhizosphere")

# w.unifrac distance matrix
set.seed(2)
rhizoendo_noSMCout_root_wuni <- phyloseq::distance(rhizoendo_noSMCout_root, "wunifrac" )

set.seed(2)
rhizoendo_noSMCout_soil_wuni <- phyloseq::distance(rhizoendo_noSMCout_soil, "wunifrac" )

#write.csv(as.matrix(rhizoendo_noSMCout_root_wuni), "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_ROOTonly_NoSMC_Wuni_primer.csv")

#write.csv(as.matrix(rhizoendo_noSMCout_soil_wuni), "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_SOILonly_NoSMC_Wuni_primer.csv")


# extract metadata in same order as distance matrix for primer analysis
names <- as.data.frame(row.names(as.matrix(rhizoendo_noSMCout_root_wuni)))
colnames(names)[1] <- "X.SampleID"
names$id  <- 1:nrow(names) # create an id for row order so we can keep the same order as the distance matrix
rhizoendo_noSMCout_root_wuni_meta <- merge(names, sampledata_rhizoendo, by = "X.SampleID")
rhizoendo_noSMCout_root_wuni_meta <- rhizoendo_noSMCout_root_wuni_meta[order(rhizoendo_noSMCout_root_wuni_meta$id), ]
dim(rhizoendo_noSMCout_root_wuni_meta)

#write.csv(rhizoendo_noSMCout_root_wuni_meta, "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_ROOTonly_NoSMC_Wuni_meta.csv")

names <- as.data.frame(row.names(as.matrix(rhizoendo_noSMCout_soil_wuni)))
colnames(names)[1] <- "X.SampleID"
names$id  <- 1:nrow(names) # create an id for row order so we can keep the same order as the distance matrix
rhizoendo_noSMCout_soil_wuni_meta <- merge(names, sampledata_rhizoendo, by = "X.SampleID")
rhizoendo_noSMCout_soil_wuni_meta <- rhizoendo_noSMCout_soil_wuni_meta[order(rhizoendo_noSMCout_soil_wuni_meta$id), ]
dim(rhizoendo_noSMCout_soil_wuni_meta)

#write.csv(rhizoendo_noSMCout_soil_wuni_meta, "16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2020.05.15_rhizoendo.rfy_4Var_soilonly_NoSMC_Wuni_meta.csv")



#######################
## FUNGI 
# presence absence for jaccard

# soil moisture content will be included as a covariate in the models 
## Remove samples B10 and C8 because soil moisture values were outliers for these. 
fungi.rfy_noSMCout <- subset_samples(fungi.rfy, X.SampleID != "TMCB10" & X.SampleID != "TMCC8")
nsamples(fungi.rfy_noSMCout) #133


pa <- function(x)(ifelse(x>0,1,0))
fungi.rfy.pa <- transform_sample_counts(otu_table(fungi.rfy_noSMCout),pa)
fungi.rfy.pa <- merge_phyloseq(fungi.rfy.pa,sample_data(fungi.rfy_noSMCout))
set.seed(2)
fungi.rfy.pa_jac<- phyloseq::distance(fungi.rfy.pa,"jaccard", binary = TRUE)

# bray curtis 
set.seed(2)
fungi.rfy_bray <- phyloseq::distance(fungi.rfy_noSMCout,"bray")

### write distance matrix as csv for primer analysis

#write.csv(as.matrix(matrix), "filedestination.csv")


```
## Final paiwise padjust
```{r}
### BACTERIA 

# read in list of pvalues from primer pairwise analysis 
Variety_pairwise_Bact_wuni <-read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PrimerAnalysis/Bacteria/2019.12.17_SMCoutliersRmvd_primer_pairwise_BACT.csv", header = TRUE)

# extract pvalues & confirm that read list is a vector
Variety_pairwise_Bact_wuni_pvalues <- Variety_pairwise_Bact_wuni$P_perm
is.vector(Variety_pairwise_Bact_wuni_pvalues)
Variety_pairwise_Bact_wuni_pvalues <- as.numeric(Variety_pairwise_Bact_wuni_pvalues)

library(stats)
Variety_pairwise_Bact_wuni_pvalues_FDR <- as.matrix(p.adjust(Variety_pairwise_Bact_wuni_pvalues, method = "fdr", length(Variety_pairwise_Bact_wuni_pvalues)))


Variety_pairwise_Bact_wuni_pvalues2 <- cbind(Variety_pairwise_Bact_wuni,Variety_pairwise_Bact_wuni_pvalues_FDR)
View(Variety_pairwise_Bact_wuni_pvalues2)
colnames(Variety_pairwise_Bact_wuni_pvalues2)[4] <- "pvalue_fdr"

# Let's reshape this into a table to more easily compare which pairwise diffs. are sig. 
library(tidyr)
Variety_pairwise_W.uni_pvalues3 <- separate(Variety_pairwise_Bact_wuni_pvalues2, Variety_Comparison, c("Comparison1","Comparison2"), sep ="_" )

#write.csv(Variety_pairwise_W.uni_pvalues3, "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2019.12.17_16Srhizo.rfy_w.uni_primer_pairwiseP_FDRadjust.csv")


# Reshape into a table 
## NOT perfect, couldn't get it into the lower triangle of the matrix.
library(reshape2)
pvalue_table <- dcast(Variety_pairwise_W.uni_pvalues3, Comparison2~Comparison1, value.var = "pvalue_fdr")

#write.csv(pvalue_table, "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2019.12.17_16Srhizo.rfy_w.uni_primer_pairwiseP_FDRadjust_TABLE.csv")

### FUNGI 

Variety_pairwise_Fungi_Bray<-read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/PrimerAnalysis/Fungi/2019.12.17_SMCoutliersRmvd_primer_pairwise_FUNGI.csv", header = TRUE)


# extract pvalues & confirm that read list is a vector
Variety_pairwise_Fungi_Bray_pvalues <- Variety_pairwise_Fungi_Bray$P_perm
is.vector(Variety_pairwise_Fungi_Bray_pvalues)
Variety_pairwise_Fungi_Bray_pvalues <- as.numeric(Variety_pairwise_Fungi_Bray_pvalues)

library(stats)
Variety_pairwise_Fungi_Bray_pvalues_FDR <- as.matrix(p.adjust(Variety_pairwise_Fungi_Bray_pvalues, method = "fdr", length(Variety_pairwise_Fungi_Bray_pvalues)))


Variety_pairwise_Fungi_Bray_pvalues2 <- cbind(Variety_pairwise_Fungi_Bray,Variety_pairwise_Fungi_Bray_pvalues_FDR)
View(Variety_pairwise_Fungi_Bray_pvalues2)
colnames(Variety_pairwise_Fungi_Bray_pvalues2)[4] <- "pvalue_fdr"

# Let's reshape this into a table to more easily compare which pairwise diffs. are sig. 
library(tidyr)
Variety_pairwise_Fungi_pvalues3 <- separate(Variety_pairwise_Fungi_Bray_pvalues2, Variety_Comparison, c("Comparison1","Comparison2"), sep ="_" )

#write.csv(Variety_pairwise_Fungi_pvalues3, "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/R_output/PrimerAnalysis/Final/noSMCoutliers/2019.12.17_fungi.bray_primer_pairwiseP_FDRadjust.csv")

```
---
# dbRDA analysis
*correlate edaphic conditions & root trait predictor variables (Table 3)*

## First check for collinearity among envirnomental variables 
```{r}

library(Hmisc)

# subset metadata for the env. variables of interest 
metadata_env <- select(metadata,SRL_giaroots_OutlierRmv, SRL, SMC_perc_GravWaterCont_OutlierRmv, Avg_root_width_diam, Network_length, NO3_ugN_g_drysoil_K2SO4_SMCOutRmv,  NH4_ugN_g_drysoil_K2SO4_SMCOutRmv,DryRootWt_total_g )
dim(metadata_env) # 144

# remove rows with NAs 
metadata_env <- metadata_env[complete.cases(metadata_env), ]
dim(metadata_env) # 139


env_corr <- rcorr(as.matrix(metadata_env),type = "pearson")

# flattenCorrMatrix to nicely visualize results 
  #http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
##################

env_corr_results <- flattenCorrMatrix(env_corr$r, env_corr$P)
#View(env_corr_results)

##### Removed variables that are correlated with -0.5<r>0.05
# SRL & avg. root diam
# SRL & dry root weight
# rooth length & dry root weight 
# Remaining predictor variables for dbRDA analyses: soil moisture content, soil ammonium, soil nitrate, average root diameter, average root length 

```

## 12 cultivars Soil Bacteria 
```{r}
#http://finzi.psych.upenn.edu/R/library/vegan/html/capscale.html

# remove any samples with missing metadata from phyloseq object
rhizo.rfy_noNA <- rhizo.rfy %>%
  subset_samples(
      !is.na(SMC_perc_GravWaterCont_OutlierRmv) & 
      !is.na(Avg_root_width_diam) &
      !is.na(Network_length) & 
      !is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv  ) & 
      !is.na(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv ) 
  )

set.seed(2)
rhizo.rfy_wuni <- phyloseq::distance(rhizo.rfy_noNA, "wunifrac" )


# Distance matrix - should be as matrix for dbrda
rhizo.rfy_w.uni.m <- as.matrix(rhizo.rfy_wuni)

# sampledata as dataframe 
sampledata_rhizo.rfy <- as.data.frame(sample_data(rhizo.rfy_noNA))
sampledata_rhizo.rfy_sub <- select(sampledata_rhizo.rfy, SMC_perc_GravWaterCont_OutlierRmv,Avg_root_width_diam,Network_length,NO3_ugN_g_drysoil_K2SO4_SMCOutRmv,NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)

# scale all variables by z score (#/mean of variable/sd variable)
#https://www.r-bloggers.com/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/

sampledata_rhizo.rfy_sub_scale <- scale(sampledata_rhizo.rfy_sub, center = TRUE, scale = TRUE)
sampledata_rhizo.rfy_sub_scale <- merge(sampledata_rhizo.rfy_sub_scale, sampledata_rhizo.rfy[,c(7,16)], by = 0) # add block back to df

# run RDA (distance matrix ~ dataframe$factor1 + dataframe$factor2... + Condition(dataframe$Block)) 
  # Condition is a random factor that is controlle for first before evaluating effect of factors 

set.seed(2)
dbRDA <- dbrda(rhizo.rfy_w.uni.m ~ sampledata_rhizo.rfy_sub_scale$NO3_ugN_g_drysoil_K2SO4 +sampledata_rhizo.rfy_sub_scale$NH4_ugN_g_drysoil_K2SO4 +  sampledata_rhizo.rfy_sub_scale$Network_length + sampledata_rhizo.rfy_sub_scale$Avg_root_width_diam + Condition(sampledata_rhizo.rfy_sub_scale$Block))

dbRDA
# Inertia Proportion Rank RealDims
#Total         1.02143    1.00000              
#Conditional   0.07830    0.07665    3         
#Constrained   0.10333    0.10117    5        5
#Unconstrained 0.83980    0.82218  127      108

summary(dbRDA)

plot(dbRDA)

anova(dbRDA) # overal significance 
#   Df SumOfSqs      F Pr(>F)    
#Model      5  0.10333 3.1253  0.001 ***
#Residual 127  0.83980                  
           
anova(dbRDA, by = "axis", perm.max = 500)

set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 
x
#                                               Df SumOfSqs      F Pr(>F)    
#sampledata_rhizo.rfy_sub_scale$NO3_ugN_g_drysoil_K2SO4   1  0.05995 9.0661  0.001 ***
#sampledata_rhizo.rfy_sub_scale$NH4_ugN_g_drysoil_K2SO4   1  0.01094 1.6551  0.066 .  
#sampledata_rhizo.rfy_sub_scale$SMC_perc_GravWaterCont    1  0.01705 2.5782  0.009 ** 
#sampledata_rhizo.rfy_sub_scale$Avg_root_width_diam       1  0.00774 1.1708  0.247    
#sampledata_rhizo.rfy_sub_scale$Network_length            1  0.00765 1.1565  0.262    
#Residual                                               127  0.83980    

# Calculate R2 
x$R2 <- (x$SumOfSqs/sum(x$SumOfSqs))*100

```

## 12 cultivars Soil Fungi
```{r}
#http://finzi.psych.upenn.edu/R/library/vegan/html/capscale.html

# remove any samples with missing metadata from phyloseq object
fungi.rfy_noNA <- fungi.rfy %>%
  subset_samples(
      !is.na(SMC_perc_GravWaterCont_OutlierRmv) & 
      !is.na(Avg_root_width_diam) &
      !is.na(Network_length) & 
      !is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv  ) & 
      !is.na(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv ) 
  )

# OTU table - OTUs should be columns
fungi.rfy.otu <-  t(as.data.frame(as(otu_table(fungi.rfy_noNA),"matrix")))

# Distance matrix - should be as matrix for dbrda
set.seed(2)
fungi.rfy_bray <- vegdist(fungi.rfy.otu, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE)
fungi.rfy_bray.m <- as.matrix(fungi.rfy_bray)

# sampledata as dataframe 
sampledata_fungi.rfy <- as.data.frame(sample_data(fungi.rfy_noNA))

sampledata_fungi.rfy_sub <- select(sample_data(fungi.rfy_noNA), SMC_perc_GravWaterCont_OutlierRmv,Avg_root_width_diam,Network_length,NO3_ugN_g_drysoil_K2SO4_SMCOutRmv,NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)


# scale all variables by z score (#/mean of variable/sd variable)
#https://www.r-bloggers.com/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/

sampledata_fungi.rfy_sub_scale <- scale(sampledata_fungi.rfy_sub, center = TRUE, scale = TRUE)
sampledata_fungi.rfy_sub_scale <- merge(sampledata_fungi.rfy_sub_scale, sampledata_fungi.rfy[,c(7,16)], by = 0) # add block back to df

# run RDA (distance matrix ~ dataframe$factor1 + dataframe$factor2... + Condition(dataframe$Block)) 
  # Condition is a random factor that is controlle for first before evaluating effect of factors 
set.seed(2)
dbRDA <- dbrda(fungi.rfy_bray.m ~sampledata_fungi.rfy_sub_scale$NO3_ugN_g_drysoil_K2SO4 +sampledata_fungi.rfy_sub_scale$NH4_ugN_g_drysoil_K2SO4 + sampledata_fungi.rfy_sub_scale$SMC_perc_GravWaterCont  + sampledata_fungi.rfy_sub_scale$Avg_root_width_diam + sampledata_fungi.rfy_sub_scale$Network_length +  Condition(sampledata_fungi.rfy_sub_scale$Block))



dbRDA
#Inertia Proportion Rank RealDims
#Total         30.46122    1.00000              
#Conditional    1.89693    0.06227    3         
#Constrained    1.53055    0.05025    5        5
#Unconstrained 27.03374    0.88748  124      121
#Inertia is squared Unknown distance 


plot(dbRDA)

anova(dbRDA) # overal significance 
#  Df SumOfSqs      F Pr(>F)    
#Model      5   1.5305 1.4041  0.001 ***
#Residual 124  27.0337       

set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 

# Calculate R2 
x$R2 <- (x$SumOfSqs/sum(x$SumOfSqs)) *100

```
## 4 cultivars Root 
```{r}
#http://finzi.psych.upenn.edu/R/library/vegan/html/capscale.html

# remove any samples with missing metadata from phyloseq object
# Subset for just roots 
rhizoendo_roots <- subset_samples(rhizoendo, SampleSite == "Endosphere")
nsamples(rhizoendo_roots) # 42 

# remove any data with missing metadata 
rhizoendo_roots_noNA <- rhizoendo_roots %>%
  subset_samples(
      !is.na(SMC_perc_GravWaterCont_OutlierRmv) & 
      !is.na(Avg_root_width_diam) &
      !is.na(Network_length) & 
      !is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv  ) & 
      !is.na(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv ) 
      )

set.seed(2)
rhizoendo_roots_noNA_wuni <- phyloseq::distance(rhizoendo_roots_noNA, method="wunifrac")

# Distance matrix - should be as matrix for dbrda
rhizoendo_roots_noNA_w.uni.m <- as.matrix(rhizoendo_roots_noNA_wuni)

# sampledata as dataframe 
sampledata_rhizoendo_roots_noNA <- as.data.frame(sample_data(rhizoendo_roots_noNA))

sampledata_rhizoendo_roots_noNA_sub <- select(sampledata_rhizoendo_roots_noNA, SMC_perc_GravWaterCont_OutlierRmv,Avg_root_width_diam,Network_length,NO3_ugN_g_drysoil_K2SO4_SMCOutRmv,NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)

# scale all variables by z score (#/mean of variable/sd variable)
#https://www.r-bloggers.com/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/

sampledata_rhizoendo_roots_noNA_sub_scale <- scale(sampledata_rhizoendo_roots_noNA_sub, center = TRUE, scale = TRUE)
sampledata_rhizoendo_roots_noNA_sub_scale <- merge(sampledata_rhizoendo_roots_noNA_sub_scale, sampledata_rhizoendo_roots_noNA[,c(7,16)], by = 0) # add block back to df

# run RDA (distance matrix ~ dataframe$factor1 + dataframe$factor2... + Condition(dataframe$Block)) 
  # Condition is a random factor that is controlle for first before evaluating effect of factors 

set.seed(2)
dbRDA <- dbrda(rhizoendo_roots_noNA_w.uni.m ~sampledata_rhizoendo_roots_noNA_sub_scale$NO3_ugN_g_drysoil_K2SO4+ +sampledata_rhizoendo_roots_noNA_sub_scale$NH4_ugN_g_drysoil_K2SO4 + sampledata_rhizoendo_roots_noNA_sub_scale$SMC_perc_GravWaterCont  + sampledata_rhizoendo_roots_noNA_sub_scale$Avg_root_width_diam + sampledata_rhizoendo_roots_noNA_sub_scale$Network_length +  Condition(sampledata_rhizoendo_roots_noNA_sub_scale$Block))

dbRDA
# Inertia Proportion Rank RealDims
#Total         1.97606    1.00000              
#Conditional   0.12773    0.06464    3         
#Constrained   0.25251    0.12778    5        5
#Unconstrained 1.59583    0.80758   33       29
#Inertia is squared Unknown distance 

plot(dbRDA)

anova(dbRDA) # overal significance 
#     Df SumOfSqs      F Pr(>F)
#Model     5  0.25251 1.0443  0.374
#Residual 33  1.59583 

set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 
#  Df SumOfSqs      F Pr(>F)  
#sampledata_rhizoendo_roots_noNA_sub_scale$SMC_perc_GravWaterCont   1  0.02710 0.5604  0.918  
#sampledata_rhizoendo_roots_noNA_sub_scale$Avg_root_width_diam      1  0.03503 0.7243  0.729  
#sampledata_rhizoendo_roots_noNA_sub_scale$Network_length           1  0.06265 1.2956  0.213  
#sampledata_rhizoendo_roots_noNA_sub_scale$NO3_ugN_g_drysoil_K2SO4  1  0.09990 2.0658  0.026 *
#sampledata_rhizoendo_roots_noNA_sub_scale$NH4_ugN_g_drysoil_K2SO4  1  0.02783 0.5755  0.900  
#Residual                                                          33  1.59583                            


# Calculate R2 
x$R2 <- (x$SumOfSqs/sum(x$SumOfSqs)) *100

```

## 4 cultivars Soil
```{r}
#http://finzi.psych.upenn.edu/R/library/vegan/html/capscale.html

# remove any samples with missing metadata from phyloseq object
# Subset for just soil 
rhizoendo_soil <- subset_samples(rhizoendo, SampleSite == "Rhizosphere")
nsamples(rhizoendo_soil) # 46 

# remove any data with missing metadata 
rhizoendo_soil_noNA <- rhizoendo_soil %>%
  subset_samples(
      !is.na(SMC_perc_GravWaterCont_OutlierRmv) & 
      !is.na(Avg_root_width_diam) &
      !is.na(Network_length) & 
      !is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv  ) & 
      !is.na(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv ) 
      )

rhizoendo_soil_noNA_wuni <- phyloseq::distance(rhizoendo_soil_noNA, method="wunifrac")

# Distance matrix - should be as matrix for dbrda
rhizoendo_soil_noNA_w.uni.m <- as.matrix(rhizoendo_soil_noNA_wuni)

# sampledata as dataframe 
sampledata_rhizoendo_soil_noNA <- as.data.frame(sample_data(rhizoendo_soil_noNA))

sampledata_rhizoendo_soil_noNA_sub <- select(sampledata_rhizoendo_soil_noNA, SMC_perc_GravWaterCont_OutlierRmv,Avg_root_width_diam,Network_length,NO3_ugN_g_drysoil_K2SO4_SMCOutRmv,NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)

# scale all variables by z score (#/mean of variable/sd variable)
#https://www.r-bloggers.com/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/

sampledata_rhizoendo_soil_noNA_sub_scale <- scale(sampledata_rhizoendo_soil_noNA_sub, center = TRUE, scale = TRUE)
sampledata_rhizoendo_soil_noNA_sub_scale <- merge(sampledata_rhizoendo_soil_noNA_sub_scale, sampledata_rhizoendo_soil_noNA[,c(7,16)], by = 0) # add block back to df

# run RDA (distance matrix ~ dataframe$factor1 + dataframe$factor2... + Condition(dataframe$Block)) 
  # Condition is a random factor that is controlle for first before evaluating effect of factors 

set.seed(2)
dbRDA <- dbrda(rhizoendo_soil_noNA_w.uni.m ~sampledata_rhizoendo_soil_noNA_sub_scale$NO3_ugN_g_drysoil_K2SO4+ +sampledata_rhizoendo_soil_noNA_sub_scale$NH4_ugN_g_drysoil_K2SO4 + sampledata_rhizoendo_soil_noNA_sub_scale$SMC_perc_GravWaterCont  + sampledata_rhizoendo_soil_noNA_sub_scale$Avg_root_width_diam + sampledata_rhizoendo_soil_noNA_sub_scale$Network_length +  Condition(sampledata_rhizoendo_soil_noNA_sub_scale$Block))

dbRDA
#    Inertia Proportion Rank
#Total         0.37831    1.00000     
#Conditional   0.03718    0.09829    3
##Constrained   0.05791    0.15307    5
#Unconstrained 0.28322    0.74864   37  

plot(dbRDA)

anova(dbRDA) # overal significance 
##  Df SumOfSqs     F Pr(>F)   
#Model     5 0.057908 1.513  0.002 **
#Residual 37 0.283218 

set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 
# Df SumOfSqs      F Pr(>F)   
#sampledata_rhizoendo_soil_noNA_sub_scale$SMC_perc_GravWaterCont   1 0.010635 1.3894  0.104   
#sampledata_rhizoendo_soil_noNA_sub_scale$Avg_root_width_diam      1 0.007690 1.0046  0.386   
#sampledata_rhizoendo_soil_noNA_sub_scale$Network_length           1 0.011707 1.5294  0.067 . 
#sampledata_rhizoendo_soil_noNA_sub_scale$NO3_ugN_g_drysoil_K2SO4  1 0.020245 2.6449  0.002 **
#sampledata_rhizoendo_soil_noNA_sub_scale$NH4_ugN_g_drysoil_K2SO4  1 0.007631 0.9969  0.420   
#Residual                                                         37 0.283218                 

# Calculate R2 
x$R2 <- (x$SumOfSqs/sum(x$SumOfSqs)) *100

```
---
# Identify variable Phyla with MVabund
  *What bacterial phyla (or classes of proteobacteria) differ among cultivars?*

## Soil Bacteria 
```{r}
#http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)
library(dplyr)
library(tibble)

# What bacterial phyla (or classes of proteobacteria) differ among cultivars?

# Populate Proteobacteria classes in phylum level 
#https://github.com/joey711/phyloseq/issues/659

rhizo.rfy2 = rhizo.rfy #create a new phyloseq object, so the old isn't disrupted
plot.tax = as.data.frame(tax_table(rhizo.rfy2))
plot.tax = data.frame(lapply(plot.tax, as.character), stringsAsFactors = F)
plot.tax$Phylum[plot.tax$Phylum=="Proteobacteria"] = plot.tax$Class[plot.tax$Phylum=="Proteobacteria"] 
plot.tax[] = lapply(plot.tax, factor)
plot.tax.table = tax_table(plot.tax)
rownames(plot.tax.table) = rownames(tax_table(rhizo.rfy2))
colnames(plot.tax.table) = colnames(tax_table(rhizo.rfy2))
identical(plot.tax.table[1:14590,3:7], tax_table(rhizo.rfy2)[1:14590,3:7]) #final check; should be true. Ensure the number of rows is adjusted to the dimensions of your tax table. 
tax_table(rhizo.rfy2) = plot.tax.table


# glom at phyla level (with proteo classes)
rhizo.rfy_phyla <- tax_glom(rhizo.rfy2, "Phylum")
ntaxa(rhizo.rfy_phyla) #43 

#save(rhizo.rfy_phyla, file = "R_output/VariableTaxa/2019.08.07_rhizo.rfy_PhylaProteoClass_phyloseq.Rdata")

rhizo.rfy_phyla_otu <- as.data.frame(as(otu_table(rhizo.rfy_phyla),"matrix"))
rhizo.rfy_phyla_otu <- rownames_to_column(rhizo.rfy_phyla_otu, "OTU")

rhizo.rfy_phyla_tax <- as.data.frame(as(tax_table(rhizo.rfy_phyla),"matrix"))
rhizo.rfy_phyla_tax <- rownames_to_column(rhizo.rfy_phyla_tax, "OTU")
rhizo.rfy_phyla_names <- merge(rhizo.rfy_phyla_otu, rhizo.rfy_phyla_tax[,1:3], by = "OTU")

#write.csv(rhizo.rfy_phyla_names, "R_output/VariableTaxa/2019.06.24_MVabund_PhylaProteoClass_OTUtable.csv")
  # I exported this file, transposed it so that the samples are rows and columns are phyla; column 1 should be labeled with 'X.SampleID' to merge with metadata all other columns should be phyla names. 

###################################
#read in file from above
rhizo.rfy_phyla <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/MVabund_Analysis/Bacteria/2019.06.24_MVabund_rhizo.rfy_PhylaProteoClass_OTUtable.csv", header = TRUE)
nrow(rhizo.rfy_phyla) #138

# merge with metadata (need to have variety in the dataframe for mvabund analysis)
rhizo.rfy_phyla_meta <- merge(metadata[,c(2,3,4,5,7,10,14)], rhizo.rfy_phyla, by = "X.SampleID")

# sample in rows and variables in column 
# convert df to mvabund object
rhizo.rfy_phyla_MVA <- mvabund(rhizo.rfy_phyla_meta[,c(8:50)])

# visualize the spread of the data
par(mar=c(3,10,3,3)) # adjusts the margins
boxplot(rhizo.rfy_phyla_meta[,c(8:50)],horizontal = TRUE,las=2, main="Abundance")

# Look at the mean-variance relationship
meanvar.plot(rhizo.rfy_phyla_MVA) # very linear relationship - those with a high means (x-axis) have a high-variance (y-axis)


########### Create models to look at variability in phyla by cultivar or root traits ###################

# Variety only 
phylum_Var <- manyglm(rhizo.rfy_phyla_MVA ~ rhizo.rfy_phyla_meta$Variety, family = "negative_binomial")
plot(phylum_Var)

# Analysis 
#This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.
# IF OTUS differ among cultivars - To examine this further, and see which  otus are more likely to be found on which cultivars, we can run univariate tests for each species separately. This is done by using the p.uni="adjusted" argument in the anova function. The "adjusted" part of the argument refers to the resampling method used to compute the p values, taking into account the correlation between the response variables; e.g. after adjusting for multiple testing, is there still an effect of cultivar on the otus? 

# NOTE: I can't do block because block needs to be balanced - and it's not.
phylum_Var_results <- anova(phylum_Var, p.uni="adjusted", nBoot = 999)
phylum_Var_results$table #anova table 


## Save results into a dataframe 
t <- as.data.frame(phylum_Var_results$table)
t <- rownames_to_column(t)
x <- as.data.frame(phylum_Var_results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(phylum_Var_results$uni.p)
y <- rownames_to_column(y)

library(plyr)
phylum_Var_results.df <- rbind.fill(t,x,y)

#write.csv(phylum_Var_results.df, "R_output/VariableTaxa/2019.06.24_Mvabund_rhizo.rfy_phylum_proteoClass_Variety_Padjust_results.csv")

```
### Statistical tests for variable phyla 
```{r}
# Populate Proteobacteria classes in phylum level 
#https://github.com/joey711/phyloseq/issues/659

rhizo.rfy2 = rhizo.rfy #create a new phyloseq object, so the old isn't disrupted
plot.tax = as.data.frame(tax_table(rhizo.rfy2))
plot.tax = data.frame(lapply(plot.tax, as.character), stringsAsFactors = F)
plot.tax$Phylum[plot.tax$Phylum=="Proteobacteria"] = plot.tax$Class[plot.tax$Phylum=="Proteobacteria"] 
plot.tax[] = lapply(plot.tax, factor)
plot.tax.table = tax_table(plot.tax)
rownames(plot.tax.table) = rownames(tax_table(rhizo.rfy2))
colnames(plot.tax.table) = colnames(tax_table(rhizo.rfy2))
identical(plot.tax.table[1:14590,3:7], tax_table(rhizo.rfy2)[1:14590,3:7]) #final check; should be true. Ensure the number of rows is adjusted to the dimensions of your tax table. 
tax_table(rhizo.rfy2) = plot.tax.table


# glom at phyla level (with proteo classes)
rhizo.rfy_phyla <- tax_glom(rhizo.rfy2, "Phylum")
ntaxa(rhizo.rfy_phyla) #43 

############
# Subset dataset for phyla that significantly differ among cultivars (determined with MVabunda bove)
sig_phyla <- as.vector(c("Acidobacteria","Actinobacteria","Verrucomicrobia","Planctomycetes","Firmicutes","Bacteroidetes","Gemmatimonadetes", "Deltaproteobacteria", "Betaproteobacteria")) 

# change otu names to phylumn names to subset by sig phylums
taxnames <- as.data.frame(as(tax_table(rhizo.rfy_phyla), "matrix"))
taxa_names(rhizo.rfy_phyla) <- taxnames$Phylum

rhizo.rfy_phy_sig <- prune_taxa(sig_phyla, rhizo.rfy_phyla)
ntaxa(rhizo.rfy_phy_sig) # 9

# convert to dataframe 
rhizo.rfy_phy_sig.df <- psmelt(rhizo.rfy_phy_sig)

# Rename Cave-in-Rock with CinR 
library(plyr)
rhizo.rfy_phy_sig.df$Variety <- revalue(rhizo.rfy_phy_sig.df$Variety, c("Cave-in-Rock" = "CinR"))

# Sum each phyla by sample 
rhizo.rfy_phy_sig.df_sums <- rhizo.rfy_phy_sig.df %>% group_by(X.SampleID, Variety, Ecotype, Phylum, Date,SMC_perc_GravWaterCont_OutlierRmv ) %>%  summarize_at("Abundance",sum)

# use this dataset for all phyla-specific analyses below

```

#### Actinobacteria
```{r}
# Actinobacteria
rhizo.rfy_phy_sig.df_sums_Actino <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Actinobacteria")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Actino$Abundance)


Actino_Var <- lm(log(Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Actino )
Actino_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Actino )

res = Actino_VarSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.7

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Actino) 
# p <0.001
# log p = 0.16


AIC(Actino_VarSMC, Actino_Var)
# df      AIC
#Actino_VarSMC 14 53.27872
#Actino_Var    13 55.35742

Anova(Actino_VarSMC, type = "III")
# Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       123.010   1 1577.8144 < 2.2e-16 ***
#Variety                             8.926  11   10.4082 2.828e-13 ***
#SMC_perc_GravWaterCont_OutlierRmv   0.320   1    4.0986   0.04508 *  
#Residuals                           9.589 123                        

marginal = lsmeans(Actino_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean     SE  df lower.CL upper.CL .group
# EG1102        5.11 0.0827 123     4.87     5.35  a    
# NE28          5.17 0.0810 123     4.94     5.41  a    
# Alamo         5.24 0.0814 123     5.00     5.48  a    
# Kanlow        5.28 0.0851 123     5.03     5.53  ab   
# EG1101        5.53 0.0811 123     5.29     5.76   bc  
# EG2101        5.63 0.0850 123     5.39     5.88    c  
# Trailblazer   5.64 0.0893 123     5.37     5.90    cd 
# Southlow      5.69 0.0813 123     5.45     5.92    cd 
# CinR          5.69 0.0885 123     5.43     5.95    cd 
# Shelter       5.70 0.0894 123     5.44     5.97    cd 
# Dacotah       5.92 0.0925 123     5.65     6.19     de
# Blackwell     6.03 0.0821 123     5.79     6.27      e

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Actino$Abundance)

Actino_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Actino )
Actino_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Actino )

res = Actino_EcoSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.4

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Actino) 
# p <0.001
# log p = 0.06

AIC(Actino_EcoSMC, Actino_Eco)
#                     df       AIC
#Actino_EcoSMC  4  95.57694
#Actino_Eco     3 116.35475


Anova(Actino_EcoSMC, type = "III")
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       250.119   1 2194.175 < 2.2e-16 ***
#Ecotype                             3.354   1   29.425 2.657e-07 ***
#SMC_perc_GravWaterCont_OutlierRmv   2.459   1   21.570 8.093e-06 ***
#Residuals                          15.161 133                                            

marginal = lsmeans(Actino_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Ecotype lsmean     SE  df lower.CL upper.CL .group
# Lowland   5.32 0.0502 133     5.21     5.44  a    
# Upland    5.67 0.0366 133     5.58     5.75   b   

```
#### Acidobacteria
```{r}
rhizo.rfy_phy_sig.df_sums_Acido <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Acidobacteria")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Acido$Abundance)


Acido_Var <- lm((Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Acido )
Acido_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Acido )

res = Acido_VarSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.2

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Acido) 
# log p = 0.11


AIC(Acido_VarSMC, Acido_Var)
# df       AIC
#Acido_VarSMC 14 -156.6605
#Acido_Var    13 -156.2354

Anova(Acido_VarSMC, type = "III")
#  Sum Sq  Df    F value    Pr(>F)    
#(Intercept)                       212.486   1 12760.1283 < 2.2e-16 ***
#Variety                             0.826  11     4.5083  1.03e-05 ***
#SMC_perc_GravWaterCont_OutlierRmv   0.054   1     3.2701     0.073 .  
#Residuals                           2.048 123            

marginal = lsmeans(Acido_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean     SE  df lower.CL upper.CL .group
# Blackwell     7.05 0.0380 123     6.94     7.16  a    
 #Dacotah       7.09 0.0427 123     6.97     7.22  ab   
 #Southlow      7.21 0.0376 123     7.10     7.32   bc  
 #Kanlow        7.23 0.0393 123     7.11     7.34   bcd 
 #Alamo         7.24 0.0376 123     7.13     7.35   bcd 
 #Trailblazer   7.24 0.0413 123     7.12     7.36   bcd 
 #EG2101        7.26 0.0393 123     7.14     7.37    cd 
 #CinR          7.26 0.0409 123     7.14     7.38    cd 
 #EG1101        7.27 0.0375 123     7.16     7.38    cd 
 #NE28          7.29 0.0375 123     7.18     7.39    cd 
 #Shelter       7.31 0.0413 123     7.19     7.43    cd 
 #EG1102        7.36 0.0382 123     7.24     7.47     d 

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Acido$Abundance)

Acido_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Acido )
Acido_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Acido )

res = Acido_EcoSMC$residuals
shapiro.test(res) 
# log p = 0.84

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Acido) 
# log p = 0.003

AIC(Acido_EcoSMC, Acido_Eco)
#    df       AIC
#Acido_EcoSMC  4 -133.9932
#Acido_Eco     3 -122.3018


Anova(Acido_EcoSMC, type = "III")
#Sum Sq  Df    F value    Pr(>F)    
#(Intercept)                       377.55   1 17914.2826 < 2.2e-16 ***
#Ecotype                             0.07   1     3.3683  0.068698 .  
#SMC_perc_GravWaterCont_OutlierRmv   0.24   1    11.2135  0.001057 ** 
#Residuals                           2.80 133                         

marginal = lsmeans(Acido_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


```
#### Bacteroidetes
```{r}
rhizo.rfy_phy_sig.df_sums_Bacto <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Bacteroidetes")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Bacto$Abundance)


Bacto_Var <- lm((Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Bacto )
Bacto_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Bacto )

res = Bacto_VarSMC$residuals
shapiro.test(res) 
# log 0.2  

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Bacto) 
# log 0.97


AIC(Bacto_VarSMC, Bacto_Var)
#    df        AIC
#Bacto_VarSMC 14   50.25287
#Bacto_Var    13 1420.21821

Anova(Bacto_VarSMC, type = "III")
# Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       113.027   1 1482.3922 < 2.2e-16 ***
#Variety                             5.322  11    6.3449 3.043e-08 ***
#SMC_perc_GravWaterCont_OutlierRmv   0.003   1    0.0433    0.8355    
#Residuals                           9.378 123                           9.589 123                        

marginal = lsmeans(Bacto_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#Variety     lsmean     SE  df lower.CL upper.CL .group
# CinR          4.61 0.0875 123     4.36     4.87  a    
# EG2101        4.66 0.0840 123     4.42     4.91  ab   
# EG1101        4.71 0.0802 123     4.47     4.94  ab   
# Shelter       4.72 0.0884 123     4.46     4.97  ab   
# Southlow      4.72 0.0804 123     4.49     4.96  ab   
# Blackwell     4.77 0.0812 123     4.53     5.00  ab   
# Dacotah       4.84 0.0914 123     4.57     5.11  ab   
# Trailblazer   4.84 0.0883 123     4.58     5.10  ab   
# EG1102        4.92 0.0817 123     4.68     5.16   b   
# NE28          5.19 0.0802 123     4.96     5.42    c  
# Alamo         5.20 0.0805 123     4.96     5.43    c  
# Kanlow        5.20 0.0841 123     4.96     5.45    c  


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Bacto$Abundance)

Bacto_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Bacto )
Bacto_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Bacto )

res = Bacto_EcoSMC$residuals
shapiro.test(res) 
# log p = 0.4

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Bacto) 
# log = 0.36

AIC(Bacto_EcoSMC, Bacto_Eco)
#  df      AIC
#Bacto_EcoSMC  4 82.55517
#Bacto_Eco     3 85.25760


Anova(Bacto_EcoSMC, type = "III")
#Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       175.821   1 1697.3730 < 2.2e-16 ***
#Ecotype                             0.923   1    8.9122  0.003373 ** 
#SMC_perc_GravWaterCont_OutlierRmv   0.191   1    1.8419  0.177023    
#Residuals                          13.777 133                                                

marginal = lsmeans(Bacto_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#Ecotype lsmean     SE  df lower.CL upper.CL .group
# Upland    4.81 0.0349 133     4.73     4.89  a    
# Lowland   4.99 0.0479 133     4.88     5.10   b   


```
#### Betaproteobacteria
```{r}
rhizo.rfy_phy_sig.df_sums_Beta <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Betaproteobacteria")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Beta$Abundance)


Beta_Var <- lm((Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Beta )
Beta_VarSMC <- lm((Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Beta )

res = Beta_VarSMC$residuals
shapiro.test(res) 
# 0.09 
#log = 0.007

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Beta) 
# p = 0.9

AIC(Beta_VarSMC, Beta_Var)
#df      AIC
#Beta_VarSMC 14 1425.070
#Beta_Var    13 1449.602

Anova(Beta_VarSMC, type = "III")
#um Sq  Df  F value    Pr(>F)    
#(Intercept)                       220626   1 117.8085 < 2.2e-16 ***
#Variety                           107323  11   5.2098 1.072e-06 ***
#SMC_perc_GravWaterCont_OutlierRmv   4748   1   2.5351    0.1139    
#Residuals                         230349 123                   

marginal = lsmeans(Beta_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean   SE  df lower.CL upper.CL .group 
 #EG1102         228 12.8 123      191      265  a     
 #NE28           235 12.6 123      199      272  ab    
 #Alamo          255 12.6 123      218      292  abc   
 #Kanlow         258 13.2 123      219      296  abcd  
 #Blackwell      273 12.7 123      236      311   bcde 
 #EG2101         281 13.2 123      242      319    cde 
 #Trailblazer    295 13.8 123      255      336    cdef
 #Dacotah        298 14.3 123      256      339    cdef
 #CinR           303 13.7 123      263      343     def
 #Shelter        304 13.9 123      264      345     def
 #Southlow       305 12.6 123      268      342      ef
 #EG1101         325 12.6 123      288      361       f

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Beta$Abundance)

Beta_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Beta )
Beta_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Beta )

res = Beta_EcoSMC$residuals
shapiro.test(res) 
# 0.18

# homoegeneity of Ecoiances
bartlett.test((Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Beta) 
# p = 0.08

AIC(Beta_EcoSMC, Beta_Eco)
# df      AIC
#Beta_EcoSMC  4 1454.452
#Beta_Eco     3 1475.509


Anova(Beta_EcoSMC, type = "III")
#Sum Sq  Df  F value Pr(>F)    
#(Intercept)                       536103   1 215.2883 <2e-16 ***
#Ecotype                             6480   1   2.6022 0.1091    
#SMC_perc_GravWaterCont_OutlierRmv     11   1   0.0045 0.9466    
#Residuals                         331192 133                



```

#### DeltaProteobacteira

```{r}
rhizo.rfy_phy_sig.df_sums_Delta <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Deltaproteobacteria")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Delta$Abundance)


Delta_Var <- lm((Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Delta )
Delta_VarSMC <- lm((Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Delta )

res = Delta_VarSMC$residuals
shapiro.test(res) 
# 0.5

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Delta) 
# 0.2


AIC(Delta_VarSMC, Delta_Var)
#  df      AIC
#Delta_VarSMC 14 1391.534
#Delta_Var    13 1420.836

Anova(Delta_VarSMC, type = "III")
#Sum Sq  Df F value    Pr(>F)    
#(Intercept)                       191871   1 131.105 < 2.2e-16 ***
#Variety                           185201  11  11.504 1.676e-14 ***
#SMC_perc_GravWaterCont_OutlierRmv   9399   1   6.422   0.01253 *  
#Residuals                         180009 123                      

marginal = lsmeans(Delta_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean   SE  df lower.CL upper.CL .group
 #EG1102         204 11.3 123      171      237  a    
 #NE28           222 11.1 123      190      254  ab   
 #Kanlow         225 11.7 123      191      259  ab   
 #Alamo          250 11.1 123      217      282   bc  
 #CinR           271 12.1 123      236      307    cd 
 #Southlow       277 11.1 123      244      309    cd 
 #Shelter        287 12.3 123      251      323     d 
 #EG2101         293 11.6 123      259      327     d 
 #Trailblazer    293 12.2 123      257      329     d 
 #EG1101         307 11.1 123      274      339     de
 #Blackwell      334 11.3 123      302      367      e
 #Dacotah        342 12.7 123      305      379      e


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Delta$Abundance)

Delta_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Delta )
Delta_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Delta )

res = Delta_EcoSMC$residuals
shapiro.test(res) 
# p = 0.4

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Delta) 
# 0.04

AIC(Delta_EcoSMC, Delta_Eco)
#   df      AIC
#Delta_EcoSMC  4 1455.362
#Delta_Eco     3 1481.180


Anova(Delta_EcoSMC, type = "III")
# Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       554661   1 221.2552 < 2.2e-16 ***
#Ecotype                            31795   1  12.6832 0.0005126 ***
#SMC_perc_GravWaterCont_OutlierRmv   4481   1   1.7873 0.1835328    
#Residuals                         333415 133                                                      

marginal = lsmeans(Delta_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD



```
#### Firmicutes
```{r}
rhizo.rfy_phy_sig.df_sums_Firmi <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Firmicutes")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Firmi$Abundance)


Firmi_Var <- lm(log(Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Firmi )
Firmi_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Firmi )

res = Firmi_VarSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.14

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Firmi) 
# p = 0.009


AIC(Firmi_VarSMC, Firmi_Var)
#  df      AIC
#Firmi_VarSMC 14 165.9924
#Firmi_Var    13 175.0133

Anova(Firmi_VarSMC, type = "III")
# Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       84.725   1 474.4552 < 2.2e-16 ***
#Variety                            9.265  11   4.7165 5.239e-06 ***
#SMC_perc_GravWaterCont_OutlierRmv  1.787   1  10.0096  0.001962 ** 
#Residuals                         21.965 123                                

marginal = lsmeans(Firmi_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean    SE  df lower.CL upper.CL .group
# EG1101        3.67 0.123 123     3.31     4.03  a    
# Shelter       3.70 0.135 123     3.31     4.10  a    
# EG2101        3.72 0.129 123     3.34     4.09  a    
# CinR          3.74 0.134 123     3.35     4.13  a    
# Southlow      3.75 0.123 123     3.39     4.10  a    
# Trailblazer   4.00 0.135 123     3.61     4.40  ab   
# NE28          4.02 0.123 123     3.66     4.38  ab   
# Alamo         4.04 0.123 123     3.68     4.40  ab   
# EG1102        4.12 0.125 123     3.76     4.49  ab   
# Kanlow        4.24 0.129 123     3.87     4.62   b   
# Dacotah       4.35 0.140 123     3.94     4.76   b   
# Blackwell     4.45 0.124 123     4.09     4.82   b   

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Firmi$Abundance)

Firmi_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Firmi )
Firmi_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Firmi )

res = Firmi_EcoSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.3

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Firmi) 
# log p 0.005

AIC(Firmi_EcoSMC, Firmi_Eco)
#    df      AIC
#Firmi_EcoSMC  4 193.5379
#Firmi_Eco     3 207.0636


Anova(Firmi_EcoSMC, type = "III")
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       154.019   1 657.4704 < 2.2e-16 ***
#Ecotype                             0.073   1   0.3097 0.5788072    
#SMC_perc_GravWaterCont_OutlierRmv   3.021   1  12.8978 0.0004617 ***
#Residuals                          31.157 133                                                                 

marginal = lsmeans(Firmi_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


```
#### Gemmatimonadetes
```{r}
rhizo.rfy_phy_sig.df_sums_Gemma <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Gemmatimonadetes")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Gemma$Abundance)


Gemma_Var <- lm(log(Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Gemma )
Gemma_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Gemma )

res = Gemma_VarSMC$residuals
shapiro.test(res) 
# log p = 0.7

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Gemma) 
# log p = 0.6


AIC(Gemma_VarSMC, Gemma_Var)
# df      AIC
#Gemma_VarSMC 14 51.94668
#Gemma_Var    13 51.40848

Anova(Gemma_VarSMC, type = "III")
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       75.656   1 979.9778 < 2.2e-16 ***
##Variety                            4.833  11   5.6916 2.326e-07 ***
#SMC_perc_GravWaterCont_OutlierRmv  0.000   1   0.0000    0.9981    
#Residuals                          9.496 123                                             

marginal = lsmeans(Gemma_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#  Variety     lsmean     SE  df lower.CL upper.CL .group
# NE28          4.11 0.0807 123     3.87     4.34  a    
# EG1102        4.16 0.0823 123     3.92     4.40  a    
# Kanlow        4.17 0.0846 123     3.92     4.42  a    
# Alamo         4.27 0.0810 123     4.03     4.50  ab   
# Southlow      4.47 0.0809 123     4.24     4.71   bc  
# EG2101        4.48 0.0846 123     4.23     4.72   bc  
# Blackwell     4.55 0.0817 123     4.31     4.79    c  
# Dacotah       4.56 0.0920 123     4.29     4.83    c  
# EG1101        4.58 0.0807 123     4.34     4.81    c  
# Shelter       4.62 0.0890 123     4.36     4.88    c  
# CinR          4.65 0.0880 123     4.39     4.91    c  
# Trailblazer   4.71 0.0889 123     4.45     4.97    c  


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Gemma$Abundance)

Gemma_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Gemma )
Gemma_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Gemma )

res = Gemma_EcoSMC$residuals
shapiro.test(res) 
# log p = 0.4

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Gemma) 
# log p = 0.8

AIC(Gemma_EcoSMC, Gemma_Eco)
#    df      AIC
#Gemma_EcoSMC  4 78.15673
#Gemma_Eco     3 79.98720

Anova(Gemma_EcoSMC, type = "III")
#Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       148.775   1 1483.4818 < 2.2e-16 ***
#Ecotype                             0.991   1    9.8821  0.002058 ** 
#SMC_perc_GravWaterCont_OutlierRmv   0.269   1    2.6784  0.104081    
#Residuals                          13.338 133                                                   

marginal = lsmeans(Gemma_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#Ecotype lsmean     SE  df lower.CL upper.CL .group
# Lowland   4.31 0.0471 133     4.21     4.42  a    
# Upland    4.50 0.0343 133     4.42     4.58   b   
```
#### Planctomycetes
```{r}
rhizo.rfy_phy_sig.df_sums_Plancto <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Planctomycetes")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Plancto$Abundance)


Plancto_Var <- lm((Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Plancto )
Plancto_VarSMC <- lm((Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Plancto )

res = Plancto_VarSMC$residuals
shapiro.test(res) 
# 0.4

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Plancto) 
# p = 0.5

AIC(Plancto_VarSMC, Plancto_Var)
# df      AIC
#Plancto_VarSMC 14 1569.187
#Plancto_Var    13 1589.685

Anova(Plancto_VarSMC, type = "III")
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       1205590   1 223.1011 < 2.2e-16 ***
#Variety                            185001  11   3.1123 0.0009953 ***
#SMC_perc_GravWaterCont_OutlierRmv    4704   1   0.8704 0.3526659    
#Residuals                          664665 123                                             

marginal = lsmeans(Plancto_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Variety     lsmean   SE  df lower.CL upper.CL .group
# Blackwell      399 21.6 123      336      462  a    
# Shelter        411 23.5 123      342      479  ab   
# EG2101         421 22.4 123      356      487  abc  
# Dacotah        426 24.3 123      355      497  abc  
# CinR           428 23.3 123      360      496  abc  
# EG1101         429 21.4 123      367      491  abc  
# Trailblazer    434 23.5 123      366      503  abc  
# Southlow       453 21.4 123      391      516  abcd 
# Kanlow         490 22.4 123      425      556   bcd 
# NE28           499 21.3 123      436      561    cd 
# Alamo          514 21.4 123      452      577     d 
# EG1102         520 21.8 123      456      583     d 
############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Plancto$Abundance)

Plancto_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Plancto )
Plancto_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Plancto )

res = Plancto_EcoSMC$residuals
shapiro.test(res) 
# p = 0.06

# homoegeneity of Ecoiances
bartlett.test((Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Plancto) 
# p = 0.2

AIC(Plancto_EcoSMC, Plancto_Eco)
#  df      AIC
#Plancto_EcoSMC  4 1572.379
#Plancto_Eco     3 1592.732


Anova(Plancto_EcoSMC, type = "III")
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)                       1644677   1 277.5040 < 2.2e-16 ***
#Ecotype                             61418   1  10.3629  0.001616 ** 
#SMC_perc_GravWaterCont_OutlierRmv    2173   1   0.3666  0.545897    
#Residuals                          788248 133                                                                  

marginal = lsmeans(Plancto_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#Ecotype lsmean    SE  df lower.CL upper.CL .group
 #Upland     437  8.34 133      418      456  a    
 #Lowland    484 11.45 133      458      510   b    

```
#### Verrucomicrobia
```{r}
rhizo.rfy_phy_sig.df_sums_Verruco <- filter(rhizo.rfy_phy_sig.df_sums, Phylum == "Verrucomicrobia")

### By Variety 
ggqqplot(rhizo.rfy_phy_sig.df_sums_Verruco$Abundance)


Verruco_Var <- lm(log(Abundance) ~ Variety , data =rhizo.rfy_phy_sig.df_sums_Verruco )
Verruco_VarSMC <- lm(log(Abundance) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Verruco )

res = Verruco_VarSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.3

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Verruco) 
# p <0.001
# log p = 0.86


AIC(Verruco_VarSMC, Verruco_Var)
#  df       AIC
#Verruco_VarSMC 14 -31.15767
#Verruco_Var    13 -22.44540

Anova(Verruco_VarSMC, type = "III")
# S Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       181.283   1 4326.2068 < 2.2e-16 ***
#Variety                             2.422  11    5.2542 9.303e-07 ***
#SMC_perc_GravWaterCont_OutlierRmv   0.441   1   10.5256  0.001516 ** 
#Residuals                           5.154 123                        

marginal = lsmeans(Verruco_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

# Dacotah       5.92 0.0678 123     5.72     6.12  a    
 #Trailblazer   6.05 0.0655 123     5.86     6.24  ab   
 #Blackwell     6.07 0.0602 123     5.89     6.25  ab   
 #EG1101        6.13 0.0595 123     5.95     6.30   b   
 #Shelter       6.14 0.0656 123     5.94     6.33   b   
 #CinR          6.18 0.0649 123     6.00     6.37   bc  
 #EG2101        6.22 0.0623 123     6.03     6.40   bcd 
 #Southlow      6.24 0.0596 123     6.06     6.41   bcd 
 #EG1102        6.36 0.0606 123     6.19     6.54    cd 
 #Alamo         6.37 0.0596 123     6.20     6.54    cd 
 #NE28          6.41 0.0594 123     6.24     6.58     d 
 #Kanlow        6.41 0.0624 123     6.23     6.60     d 


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Verruco$Abundance)

Verruco_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Verruco )
Verruco_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Verruco )

res = Verruco_EcoSMC$residuals
shapiro.test(res) 
# 0.001
# log p = 0.7

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Verruco) 
# p <0.001
# log p = 0.7

AIC(Verruco_EcoSMC, Verruco_Eco)
#  df       AIC
#Verruco_EcoSMC  4 -7.461436
#Verruco_Eco     3 -7.451777


Anova(Verruco_EcoSMC, type = "III")
# Sum Sq  Df   F value    Pr(>F)    
#(Intercept)                       300.629   1 5625.8993 < 2.2e-16 ***
#Ecotype                             0.469   1    8.7747  0.003619 ** 
#SMC_perc_GravWaterCont_OutlierRmv   0.033   1    0.6110  0.435789    
#Residuals                           7.107 133                                                                  

marginal = lsmeans(Verruco_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

#Ecotype lsmean     SE  df lower.CL upper.CL .group
 #Upland    6.17 0.0251 133     6.11     6.23  a    
 #Lowland   6.30 0.0344 133     6.22     6.38   b  

```

## Soil Fungi 
```{r}
library(mvabund)
library(dplyr)
library(tibble)

# glom at phyla level 
fungi.rfy_phyla <- tax_glom(fungi.rfy, "Phylum")
ntaxa(fungi.rfy_phyla)  #21

fungi.rfy_phyla_otu <- as.data.frame(t(as(otu_table(fungi.rfy_phyla),"matrix")))
fungi.rfy_phyla_otu <- rownames_to_column(fungi.rfy_phyla_otu, "X.SampleID")

# merge with metadata 
fungi.rfy_phyla_all <- merge(sampledata_fungi.rfy[,c(2:5,10)], fungi.rfy_phyla_otu, by= "X.SampleID")

# sample in rows and variables in column 
# convert df to mvabund object
fungi.rfy_phyla_all.MVA <- mvabund(fungi.rfy_phyla_all[,c(6:26)])

# Look at the mean-variance relationship
meanvar.plot(fungi.rfy_phyla_all.MVA) # very linear relationship - those with a high means (x-axis) have a high-variance (y-axis)

# Variety only 
fungi.phyla_var <- manyglm(fungi.rfy_phyla_all.MVA ~ fungi.rfy_phyla_all$Variety, family = "negative_binomial")
plot(fungi.phyla_var)

# Ecotype only 
fungi.phyla_eco <- manyglm(fungi.rfy_phyla_all.MVA ~ fungi.rfy_phyla_all$Ecotype, family = "negative_binomial")
plot(fungi.phyla_eco)


# Analysis 
#This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.

fungi.phyla_var_results <-anova(fungi.phyla_var, p.uni = "adjusted", nBoot =999 )

fungi.phyla_eco_results <-anova(fungi.phyla_eco, p.uni = "adjusted", nBoot =999 )

results$table
           #  Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                    134      NA       NA       NA
#fungi.rfy_phyla_all$Ecotype    133       1 21.90787      0.7

t <- as.data.frame(results$table)
t <- rownames_to_column(t)
x <- as.data.frame(results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(results$uni.p)
y <- rownames_to_column(y)

library(plyr)
fungi.phyla_var_results.df <- rbind.fill(t,x,y)

#########################
#### Rozellomycota######

# Rozellomycota rel.abundance differed among fungi - what % of reads did it make up ? 

fungi.phyla <- tax_glom(fungi.rfy,"Phylum")
ntaxa(fungi.phyla) #21

fungi.Rozello <- subset_taxa(fungi.phyla, Phylum == "Rozellomycota")
ntaxa(fungi.Rozello) #1

fungi.Rozello.df <- psmelt(fungi.Rozello)

# Relative % of reads that Rozello makes up? 
sum(sample_sums(fungi.Rozello))/ sum(sample_sums(fungi.phyla)) *100 # 0.73% 


### Variation in abundance among cultivars? 
fungi.Rozello.df_Var <- lm(log(Abundance + 1) ~ Variety , data =fungi.Rozello.df)
fungi.Rozello.df_VarSMC <- lm(log(Abundance + 1) ~ Variety +SMC_perc_GravWaterCont_OutlierRmv, data =fungi.Rozello.df )

res = fungi.Rozello.df_Var$residuals
shapiro.test(res) 
# 0.001
# log + 1 p > 0.05

# homoegeneity of variances
bartlett.test(log(Abundance + 1) ~ Variety, data = fungi.Rozello.df) 
# p <0.001
# log p >0.05


AIC(fungi.Rozello.df_Var, fungi.Rozello.df_VarSMC)
# Df      AIC
#fungi.Rozello.df_Var    13 458.3770
#fungi.Rozello.df_VarSMC 14 444.2556

Anova(fungi.Rozello.df_VarSMC, type = "III")
#Sum Sq  Df F value    Pr(>F)    
#(Intercept)                        54.079   1 36.4418 1.805e-08 ***
#Variety                            48.576  11  2.9758  0.001585 ** 
###SMC_perc_GravWaterCont_OutlierRmv   2.560   1  1.7253  0.191516    
#Residuals                         178.077 120                        

```
---

# Orders or Core OTUS that vary with root traits 
## Soil Bacteria 
```{r}
#http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)
library(dplyr)
library(tibble)

# glom at order level 
rhizo.rfy_order <- tax_glom(rhizo.rfy, "Order")
ntaxa(rhizo.rfy_order) #236

rhizo.rfy_order_otu <- as.data.frame((as(otu_table(rhizo.rfy_order),"matrix")))
rhizo.rfy_order_otu <- rownames_to_column(rhizo.rfy_order_otu, "OTU")

rhizo.rfy_order_tax <- as.data.frame(as(tax_table(rhizo.rfy_order),"matrix"))
rhizo.rfy_order_tax <- rownames_to_column(rhizo.rfy_order_tax, "OTU")
rhizo.rfy_order_names <- merge(rhizo.rfy_order_otu, rhizo.rfy_order_tax[,1:5], by = "OTU")

#write.csv(rhizo.rfy_order_names, "R_output/VariableTaxa/2019.05.08_MVabund_OrderOTUtable.csv")
# I exported this file, transposed it so that the samples are rows and columns are phyla; column 1 should be labeled with 'X.SampleID' to merge with metadata; all other columns should be phyla names. 

################
# Order ~ Root Diam
################

rhizo.rfy_order <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/MVabund_Analysis/Bacteria/2019.05.08_MVabund_OrderOTUtable_XsampleID.csv", header = TRUE)

  # merge to add variety, soil moisture content, and root traits to dataframe
sampledata_rhizo.rfy <- as(sample_data(rhizo.rfy), "data.frame")
rhizo.rfy_order_meta <- merge(sampledata_rhizo.rfy[,c(2,3,10,18,23)], rhizo.rfy_order, by = "X.SampleID")

library(mvabund)
rhizo.rfy_order.MVA <- mvabund(rhizo.rfy_order_meta[,c(9:243)])

# Do any orders differ with root diameter?
Order_Rootdiam <- manyglm(rhizo.rfy_order.MVA ~ rhizo.rfy_order_meta$Avg_root_width_diam, family = "negative_binomial")
plot(Order_Rootdiam) # residuals look normal

Order_Rootdiam_results <-anova(Order_Rootdiam, p.uni = "adjusted", nBoot =999)

Order_Rootdiam_results$table
 #                                        Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                                 137      NA       NA       NA
#rhizo.rfy_order_meta$Avg_root_width_diam    136       1 369.4605    0.012

# make results into a dataframe 
t <- as.data.frame(Order_Rootdiam_results$table)
t <- rownames_to_column(t)
x <- as.data.frame(Order_Rootdiam_results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(Order_Rootdiam_results$uni.p)
y <- rownames_to_column(y)

library(plyr)
Order_Rootdiam_results.df <- rbind.fill(t,x,y)

#write.csv(Order_Rootdiam_results.df, "R_output/VariableTaxa/Final/2019.08.07_MVabund_rhizo.rfyOrder_RootDiameter_results_p999_Padjust.csv")
# export to see if any orders significnatly vary by root diameter.

###########################################
### Variation in core taxa (present in at least 80% of samples) by root traits

# convert all abundance data to rel.abund.
rhizo.rfy_ra = transform_sample_counts(rhizo.rfy, function(x) 100 * x/sum(x))


# Filter for taxa present in 80% of samples
# e.g. detection = 0 and prevalence = .8; or taxa abundant in over 80% of the samples and with relative abundance in their sample at > 0 (no threshold)
rhizo80perc <- core(rhizo.rfy_ra, prevalence = .8, detection = 0)
rhizo80perc_names <- taxa(rhizo80perc) # list of otus 
length(rhizo80perc_names) #316

rhizo.rfy_ra.80perc <- prune_taxa(rhizo80perc_names,rhizo.rfy)
ntaxa(rhizo.rfy_ra.80perc) #316 

rhizo.rfy_ra.80perc_otu <- as.data.frame(as(t(otu_table(rhizo.rfy_ra.80perc)),"matrix"))
rhizo.rfy_ra.80perc_otu <- rownames_to_column(rhizo.rfy_ra.80perc_otu, "X.SampleID")

rhizo.rfy_ra.80perc_samp <- as.data.frame(as(sample_data(rhizo.rfy),"matrix"))
rhizo.rfy_ra.80perc_samp.p <- rhizo.rfy_ra.80perc_samp[,c(1,2,3,7,18,21,23)] # select columns
rhizo.rfy_ra.80perc_otu_samp <- merge(rhizo.rfy_ra.80perc_samp.p, rhizo.rfy_ra.80perc_otu,by = "X.SampleID")
#View(rhizo.rfy_ra.80perc_otu_samp)
# for some reason the root traits are factors....convert to numeric
 rhizo.rfy_ra.80perc_otu_samp$Avg_root_width_diam <- as.numeric(as.character(rhizo.rfy_ra.80perc_otu_samp$Avg_root_width_diam ))
 rhizo.rfy_ra.80perc_otu_samp$SRL_giaroots <- as.numeric(as.character(rhizo.rfy_ra.80perc_otu_samp$SRL_giaroots ))
  rhizo.rfy_ra.80perc_otu_samp$Network_length <- as.numeric(as.character(rhizo.rfy_ra.80perc_otu_samp$Network_length ))
str(rhizo.rfy_ra.80perc_otu_samp)

###################### MVAbund 

# sample in rows and variables in column 
# convert df to mvabund object
library(mvabund)
rhizo.rfy_core.MVA <- mvabund(rhizo.rfy_ra.80perc_otu_samp[,c(8:323)])

# Look at the mean-variance relationship
meanvar.plot(rhizo.rfy_core.MVA) # very linear relationship - those with a high means (x-axis) have a high-variance (y-axis)

# Are some phyla more abundant in particular cultivars? 
plot(rhizo.rfy_core.MVA~rhizo.rfy_ra.80perc_otu_samp$Variety, cex.axis=0.8, cex=0.8)


#do any otus differ with volume-weighted SRL?
rhizo.rfy_core_SRL <- manyglm(rhizo.rfy_core.MVA ~  rhizo.rfy_ra.80perc_otu_samp$SRL_giaroots, family = "poisson")

# do any otus differ with root diameter?
rhizo.rfy_core_diam <- manyglm(rhizo.rfy_core.MVA ~  rhizo.rfy_ra.80perc_otu_samp$Avg_root_width_diam, family = "poisson")


plot(rhizo.rfy_core_diam) # residuals look slightly fan shaped...

# Analysis 
#This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.

rhizo.rfy_core_SRL_results <-anova(rhizo.rfy_core_SRL, p.uni = "adjusted", nBoot =9999)

rhizo.rfy_core_SRL_results$table
#                                           
 #                                    Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                                  137      NA       NA       NA
#rhizo.rfy_ra.80perc_otu_samp$SRL_giaroots    136       1 1223.051    1e-04

rhizo.rfy_core_diam_results <-anova(rhizo.rfy_core_diam, p.uni = "adjusted", nBoot =9999)
rhizo.rfy_core_diam_results$table
#                                                 Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                                         137      NA       NA       NA
#rhizo.rfy_ra.80perc_otu_samp$Avg_root_width_diam    136       1 1495.062    1e-04


##############################
# make results into a dataframe 
t <- as.data.frame(rhizo.rfy_core_diam_results$table)
t <- rownames_to_column(t)
x <- as.data.frame(rhizo.rfy_core_diam_results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(rhizo.rfy_core_diam_results$uni.p)
y <- rownames_to_column(y)

library(plyr)
rhizo.rfy_core_diam_results.df <- rbind.fill(t,x,y)

#write.csv(rhizo.rfy_core_SRL_results.df,"R_output/VariableTaxa/    2019.07.28_MVabund_rhizo.rfy_80perCORE_SRL_p9999_results_Padjust.csv")

#write.csv(rhizo.rfy_core_diam_results.df, "R_output/VariableTaxa/2019.07.28_MVabund_rhizo.rfy_80perCORE_RootDiam_p9999_results_Padjust_readme.csv")


# looked at these files to see if any OTUS significantly differ with volume-weighted SRL (giaroots) or root diameter -- none did. 


```
## Soil Fungi 
```{r}
#http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)
library(dplyr)
library(tibble)

# glom at order level 
fungi.rfy_order <- tax_glom(fungi.rfy, "Order")
ntaxa(fungi.rfy_order)  #129

fungi.rfy_order_otu <- as.data.frame(as(otu_table(fungi.rfy_order),"matrix"))
fungi.rfy_order_otu <- rownames_to_column(fungi.rfy_order_otu, "OTU")

fungi.rfy_order_tax <- as.data.frame(as(tax_table(fungi.rfy_order),"matrix"))
fungi.rfy_order_tax <- rownames_to_column(fungi.rfy_order_tax, "OTU")
fungi.rfy_order_names <- merge(fungi.rfy_order_otu, fungi.rfy_order_tax[,1:5], by = "OTU")

#write.csv(fungi.rfy_order_names, "R_output/VariableTaxa/2019.05.10_MVabund_fungi.rfy_OrderOTUtable.csv")
# I exported this file, transposed it so that the samples are rows and columns are phyla; column 1 should be labeled with 'X.SampleID' to merge with metadata; all other columns should be phyla names. 

################
# Order ~ Root Diam
################

# save as csv and then add variety to table
fungi.rfy_order <- read.csv("C:/Users/Tayler/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/ITS/R_output/VariableTaxa/Final/2019.05.10_MVabund_fungi.rfy_OrderOTUtable_clean.csv")

# merge with metadata to add variety and root traits to dataframe 
fungi.rfy_order_meta <- merge(sampledata_fungi.rfy[,c(2,3,7,8,10,18,23)], fungi.rfy_order, by = "X.SampleID")

library(mvabund)
fungi.rfy_order_all.MVA <- mvabund(fungi.rfy_order_meta[,c(8:136)])

# Do any orders differ with root diameter? 
Order_Rootdiam <- manyglm(fungi.rfy_order_all.MVA ~ fungi.rfy_order_meta$Avg_root_width_diam, family = "negative_binomial")
plot(Order_Rootdiam) # residuals look normal


Order_Rootdiam_results <-anova(Order_Rootdiam, p.uni = "adjusted", nBoot =999)

Order_Rootdiam_results$table
#                                         Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                                 134      NA       NA       NA
#fungi.rfy_order_meta$Avg_root_width_diam    133       1 238.8551    0.016

# put results into dataframe 
t <- as.data.frame(Order_Rootdiam_results$table)
t <- rownames_to_column(t)
x <- as.data.frame(Order_Rootdiam_results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(Order_Rootdiam_results$uni.p)
y <- rownames_to_column(y)

library(plyr)
Order_Rootdiam_results.df <- rbind.fill(t,x,y)

#write.csv(Order_Rootdiam_results.df, "C:/Users/Tayler/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/ITS/R_output/VariableTaxa/Final/2019.08.07_MVabund_fungi.rfyOrder_RootDiameter_results_p999_Padjust.csv")

################
# Order ~ Root Length
################

fungi.rfy_order_meta <- merge(sampledata_fungi.rfy[,c(2,3,7,8,10,18,23)], fungi.rfy_order, by = "X.SampleID")

library(mvabund)
fungi.rfy_order_all.MVA <- mvabund(fungi.rfy_order_meta[,c(8:136)])

# Do any Orders differ with root length?
Order_RootLength <- manyglm(fungi.rfy_order_all.MVA ~ fungi.rfy_order_meta$Network_length, family = "negative_binomial")
plot(Order_RootLength) # residuals look normal



Order_RootLength_results <-anova(Order_RootLength, p.uni = "adjusted", nBoot =999)

Order_RootLength_results$table
#                                    Res.Df Df.diff      Dev Pr(>Dev)
#(Intercept)                            134      NA       NA       NA
#fungi.rfy_order_meta$Network_length    133       1 313.2811    0.001
t <- as.data.frame(Order_RootLength_results$table)
t <- rownames_to_column(t)
x <- as.data.frame(Order_RootLength_results$uni.test)
x <- rownames_to_column(x)
y <- as.data.frame(Order_RootLength_results$uni.p)
y <- rownames_to_column(y)

library(plyr)
Order_RootLength_results.df <- rbind.fill(t,x,y)

#write.csv(Order_RootLength_results.df, "C:/Users/Tayler/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/ITS/R_output/VariableTaxa/Final/2019.08.07_MVabund_fungi.rfyOrder_RootLength_results_p999_Padjust.csv")

###################################################################################
####  Mortierellales correlated with root traits; confirm this with correlation

# glom by order 
fungi.rfy_order <- tax_glom(fungi.rfy, taxrank = "Order")

# relative abundance
fungi.rfy_order.ra = transform_sample_counts(fungi.rfy_order, function(x) 100 * x/sum(x))

# Put into a dataframe
fungi.rfy_order.ra.df<- psmelt(fungi.rfy_order.ra) #create dataframe from phyloseq


# Calculate mean realtive abundance by sample
fungi.rfy_order.ra.df1 <- fungi.rfy_order.ra.df %>% group_by(Order,X.SampleID) %>%  summarize_at("Abundance",mean)

fungi.rfy_order.ra.df1_meta <- merge(fungi.rfy_order.ra.df1, sampledata_fungi.rfy, by = "X.SampleID")

fungi.rfy_order.ra.df1_meta_Morti <- filter(fungi.rfy_order.ra.df1_meta, Order == "Mortierellales" )

# Mortierellales correlates - lets look at a linear regression
ggqqplot(fungi.rfy_order.ra.df1_meta_Morti$Abundance)
ggqqplot(sqrt(fungi.rfy_order.ra.df1_meta_Morti$Network_length))


cor.test(fungi.rfy_order.ra.df1_meta_Morti$Abundance, sqrt(fungi.rfy_order.ra.df1_meta_Morti$Network_length), method = "pearson")

#	Pearson's product-moment correlation
#data:  fungi.rfy_order.ra.df1_meta_Morti$Abundance and sqrt(fungi.rfy_order.ra.df1_meta_Morti$Network_length)
#t = -5.147, df = 133, p-value = 9.292e-07
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.5393730 -0.2562446
#sample estimates:
 #      cor -0.4075569 


#############################
### Variation in core taxa (present in at least 80% of samples) by root traits
##############################

library(microbiome)

# convert all abundance data to rel.abund.
fungi.rfy_ra = transform_sample_counts(rhizo.rfy, function(x) 100 * x/sum(x))


# phyloseq object with core micobiota (those in 80% of the samples)
# e.g. detection = .00 and prevalence = .8; or taxa abundant in over 80% of the samples and with relative abundance in their sample at > 0% (no threshold) 
fungi80perc <- core(fungi.rfy_ra, prevalence = .8, detection = 0)
fungi80perc_names <- taxa(fungi80perc) # list of otus 
length(fungi80perc_names) #70

fungi.rfy.80perc <- prune_taxa(fungi80perc_names,fungi.rfy)
ntaxa(fungi.rfy.80perc) #70

fungi.rfy.80perc_otu <- as.data.frame(as(t(otu_table(fungi.rfy.80perc)),"matrix"))
fungi.rfy.80perc_otu <- rownames_to_column(fungi.rfy.80perc_otu, "X.SampleID")

fungi.rfy.80perc_samp <- as.data.frame(as(sample_data(fungi.rfy),"matrix"))
fungi.rfy.80perc_samp.p <- fungi.rfy.80perc_samp[,c(2,3,18,23)]
fungi.rfy.80perc_otu_samp <- merge(fungi.rfy.80perc_samp.p, fungi.rfy.80perc_otu,by = "X.SampleID")

###################### MVAbund 

## Do core otus vary by Root Traits 

# sample in rows and variables in column 
# convert df to mvabund object
fungi.rfy_core.MVA <- mvabund(fungi.rfy.80perc_otu_samp[,c(5:74)])

# Look at the mean-variance relationship
meanvar.plot(fungi.rfy_core.MVA) # very linear relationship - those with a high means (x-axis) have a high-variance (y-axis)

# Are some phyla more abundant in particular cultivars? 
plot(fungi.rfy_core.MVA~fungi.rfy.80perc_otu_samp$Variety, cex.axis=0.8, cex=0.8)


# Do any Otus vary with root diameter? 
fungi.rfy_core_rootDiam <- manyglm(fungi.rfy_core.MVA ~ fungi.rfy.80perc_otu_samp$Avg_root_width_diam, family = "negative_binomial")

plot(fungi.rfy_core_rootDiam) # residuals look slightly fan shaped...

# Analysis 
#This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.

fungi.rfy_core_rootDiam_results <-anova(fungi.rfy_core_rootDiam, p.uni = "adjusted", nBoot =9)

fungi.rfy_core_VarSRLSMC_results$table
#                                              Res.Df Df.diff          Dev Pr(>Dev)

#fungi.rfy_ra.80perc_otu_samp_noNA$Variety         121      11 7.247022e+03      0.1
#fungi.rfy_ra.80perc_otu_$SRL_giaroots_OutlierRmv -12     133 7.969945e+04      0.1
#fungi.rfy_ra.80per$SMC_GravWaterCont_OutlierRmv      0      21 9.057587e-03      1.0

```
---
# Shared Taxa across all cultivars

*We defined shared taxa as those present in at least 75% of the samples within each cultivar (e.g.,9/12 sample units, 3 cores x 4 blocks = 12 samples/cultivar)*


## Soil Bacteria 
```{r}

# put dataframe into presnce absence to identify OTUS present in all cultivars 
pa <- function(x)(ifelse(x>0,1,0))
rhizo.rfy_pa <- transform_sample_counts(otu_table(rhizo.rfy),pa)
  # this is an otu table and it needs remerged with sample data 
rhizo.rfy_pa <- merge_phyloseq(rhizo.rfy_pa,sample_data(rhizo.rfy), tax_table(rhizo.rfy))
ntaxa(rhizo.rfy_pa)
nsamples(rhizo.rfy_pa) #12

###########################

# MELT Presence/abs phyloseq object to a dataframe
rhizo.rfy_pa.df<- psmelt(rhizo.rfy_pa)

rhizo.rfy_pa.df2 <- rhizo.rfy_pa.df[,c(1:6)]

library(dplyr)
# what is the mean abundance for each otu by variety (this will = the % reps the otus is found in per variety)
rhizo.rfy_pa.mean <- rhizo.rfy_pa.df2 %>%
  group_by(Variety,OTU) %>%
  summarise(PercAbund_Reps = mean(Abundance))
head(rhizo.rfy_pa.mean)

# Dataframe of the abundance of OTU/# of replicates for each otu by variety
library(tidyr)
rhizo.rfy_pa.mean2 <- spread(rhizo.rfy_pa.mean, Variety, PercAbund_Reps)

# Add taxononmy info to df
rhizo.rfy_pa.df_tax <- rhizo.rfy_pa.df[,c(1,76:82)]
rhizo.rfy_pa.df_tax <- rhizo.rfy_pa.df_tax %>% distinct(OTU, .keep_all = TRUE)
nrow(rhizo.rfy_pa.df_tax)

rhizo.rfy_pa.mean_tax<- merge(rhizo.rfy_pa.mean2, rhizo.rfy_pa.df_tax, by = "OTU")
head(rhizo.rfy_pa.mean_tax)
nrow(rhizo.rfy_pa.mean_tax) #14590 

#save this file 
#write.csv(rhizo.rfy_pa.mean_tax, "R_output/VariableTaxa/2019.05.18_rhizo.rfy_percRepsVar_taxa.csv")

### Subset all OTUS that are present in at least 75% of replicates for each cultivar 
rhizo.rfy_pa.mean_core <- rhizo.rfy_pa.mean2[apply(rhizo.rfy_pa.mean2 > 0.75,1, all),]
nrow(rhizo.rfy_pa.mean_core) #160 
head(rhizo.rfy_pa.mean_core)
# Add taxonomy to these core 
rhizo.rfy_pa.mean_core_tax <- merge(rhizo.rfy_pa.mean_core, rhizo.rfy_pa.df_tax, by = "OTU")
nrow(rhizo.rfy_pa.mean_core_tax) # 160 

#write.csv(rhizo.rfy_pa.mean_core_tax, "R_output/VariableTaxa/2019.05.18_rhizo.rfy_percRepsVar_75cutoff_AllVarshared.csv")


## What % of total sequences are comprised of the shared OTUS present in 75% of all cultivars replicates? 
shared_tax <- rhizo.rfy_pa.mean_core_tax$OTU
length(shared_tax) #160 
rhizo.rfy_shared <- prune_taxa(shared_tax, rhizo.rfy)
ntaxa(rhizo.rfy_shared)

sum(sample_sums(rhizo.rfy_shared))/sum(sample_sums(rhizo.rfy))*100 #44.59% 


#### What classes are the most abundant of the core shared taxa? (load summary stats function in previous code)
# Create summary data on the dominate taxa, (load function above) from:https://github.com/joey711/phyloseq/issues/418

# Relative abundance of bacteria in rhizosphere 
plot_taxa_summary(rhizo.rfy_shared, "Class") # mean relative abundance 
rhizo.rfy_shared_summary <- summarize_taxa(rhizo.rfy_shared, "Class")

#write.csv(rhizo.rfy_shared_summary, "R_output/VariableTaxa/2019.05.18_rhizo.rfy_percRepsVar_75cutoff_AllVarshared_ClassSummary.csv")


```
## Soil Fungi 
```{r}
# Convert to p/a - now we know which otus are found at least once in a cultivar 
pa <- function(x)(ifelse(x>0,1,0))
fungi.rfy_pa <- transform_sample_counts(otu_table(fungi.rfy),pa)
  # this is an otu table and it needs remerged with sample data 
fungi.rfy_pa <- merge_phyloseq(fungi.rfy_pa,sample_data(fungi.rfy), tax_table(fungi.rfy))
ntaxa(fungi.rfy_pa)
nsamples(fungi.rfy_pa) #135

###########################

# MELT Presence/abs phyloseq object to a dataframe
fungi.rfy_pa.df<- psmelt(fungi.rfy_pa)
colnames(fungi.rfy_pa.df)

fungi.rfy_pa.df2 <- fungi.rfy_pa.df[,c(1:6)] #keep only otu, sample, abundance, variety

library(dplyr)
# what is the mean abundance for each otu by variety (this will = the % reps the otus is found in per variety)
fungi.rfy_pa.mean <- fungi.rfy_pa.df2 %>%
  group_by(Variety,OTU) %>%
  summarise(PercAbund_Reps = mean(Abundance))
head(fungi.rfy_pa.mean)

# Dataframe of the abundance of OTU/# of replicates for each otu by variety
library(tidyr)
fungi.rfy_pa.mean2 <- spread(fungi.rfy_pa.mean, Variety, PercAbund_Reps)
head(fungi.rfy_pa.mean2)

# Add taxononmy info to df
fungi.rfy_pa.df_tax <- fungi.rfy_pa.df[,c(1,77:83)]
fungi.rfy_pa.df_tax <- fungi.rfy_pa.df_tax %>% distinct(OTU, .keep_all = TRUE)
nrow(fungi.rfy_pa.df_tax) #4064

fungi.rfy_pa.mean_tax<- merge(fungi.rfy_pa.mean2, fungi.rfy_pa.df_tax, by = "OTU")
head(fungi.rfy_pa.mean_tax)
nrow(fungi.rfy_pa.mean_tax) #4064 

#save this file 
#write.csv(fungi.rfy_pa.mean_tax, "R_output/VariableTaxa/2019.05.19_fungi.rfy_percRepsVar_taxa.csv")
############

### Subset all OTUS that are present in at least 75% of replicates for each cultivar 
fungi.rfy_pa.mean_core <- fungi.rfy_pa.mean2[apply(fungi.rfy_pa.mean2 > 0.75,1, all),]
nrow(fungi.rfy_pa.mean_core) #37
head(fungi.rfy_pa.mean_core)
# Add taxonomy to these core 
fungi.rfy_pa.mean_core_tax <- merge(fungi.rfy_pa.mean_core, fungi.rfy_pa.df_tax, by = "OTU")
nrow(fungi.rfy_pa.mean_core_tax) #37

#write.csv(fungi.rfy_pa.mean_core_tax, "R_output/VariableTaxa/2019.05.18_fungi.rfy_percRepsVar_75cutoff_AllVarshared.csv")


## What % of total sequences are comprised of the shared OTUS present in 75% of all cultivars replicates? 
shared_tax <- fungi.rfy_pa.mean_core_tax$OTU
fungi.rfy_shared <- prune_taxa(shared_tax, fungi.rfy)
ntaxa(fungi.rfy_shared)

sum(sample_sums(fungi.rfy_shared))/sum(sample_sums(fungi.rfy))*100 #35.29


#### What classes are the most abundant of the core shared taxa? (load summary stats function in previous code)

# Relative abundance of bacteria in rhizosphere 
plot_taxa_summary(fungi.rfy_shared, "Class") # mean relative abundance 
fungi.rfy_shared_summary <- summarize_taxa(fungi.rfy_shared, "Phylum")

#write.csv(fungi.rfy_shared_summary, "R_output/VariableTaxa/2019.05.18_fungi.rfy_percRepsVar_75cutoff_AllVarshared_ClassSummary.csv")

```
---
# Indicator Species Analysis 
  *We defined Indicator taxa (after removing singleton OTUs only) using the ‘multiplatt’ function in the indicspecies R package (Caceres and Legendre 2009) and defined as OTUs present in at least 25% of the samples (3/12 sample units, or indicspecies specificity parameter= 0.25)*

## Soil Bacteria 
```{r}
# Remove anything less than something present in 1 sample (singletons)
rhizo.rfy.removesing <- prune_taxa(taxa_sums(rhizo.rfy) > 1, rhizo.rfy)
ntaxa(rhizo.rfy) - ntaxa(11086) #3504
ntaxa(rhizo.rfy.removesing) # 6801

# data needs to have samples in rows and species in columns (opposite of phyloseq)
rhizo.rfy.removesing_otu_t <- as.data.frame(t(as(otu_table(rhizo.rfy.removesing), "matrix")))
dim(rhizo.rfy.removesing_otu_t) # 138 11086
sampledata_rhizo.rfy.removesing<- as.data.frame(as(sample_data(rhizo.rfy.removesing),"matrix"))

#############
# multiplatt is the function used to conduct indicator species analysis
  # duleg = T makes it so we avoid considering grouping combinations
  # func = IndVal.g sets our test statistic as the INdVal index (returns sq.rt of IndVal)
  # indvalcomp=TRUE will give us the separate A and B parameters for the indicator species. A indicates the probablity that the Sample belongs to the identified sample given the species in the sample (specificity); B is the probablity of finding the species in the sample belonging to that group (sensitivity) 

library(indicspecies)
# by variety 
Var_IndSp <- multipatt(rhizo.rfy.removesing_otu_t,
                       sampledata_rhizo.rfy.removesing$Variety, duleg = TRUE, # in future make duleg = TRUE so it doesn't consider other site combinations
                       control = how(nperm = 9999),func = "IndVal.g")



# send ISA summary results to a file 
 sink("R_output/IndicatorTaxa/2019.07.18_rhizo.rfy_nosingles_AllOTUs_perm9999.csv")
print(summary(Var_IndSp, indvalcomp=TRUE))
sink()
# exported this table and organized into columns 

###############################
### Indicator Padjust#########

# read in indicator taxa organized csv file 
IndSp <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/SharedIndicator_taxa/2019.07.18_rhizo.rfy_nosingles_AllOTUs_perm9999_clean.csv", header = TRUE)

# read in list of pvalues from primer pairwise analysis 
IndSp_p <-IndSp$p_value

# confirm that read list is a vector
is.vector(IndSp_p)
IndSp_p <- as.numeric(IndSp_p)

library(stats)
IndSp_p_fdr <- as.matrix(p.adjust(IndSp_p, method = "fdr", length(IndSp_p)))
colnames(IndSp_p_fdr) <- "p_value_fdr"


Tax <- as.data.frame(as(tax_table(rhizo.rfy),"matrix")) # taxonomy table with OTU column to merge with indsp results 
Tax <- rownames_to_column(Tax, "OTU")
IndSp_tax <- merge(Tax, IndSp, by = "OTU")


IndSp_p_fdr_merge <- cbind(IndSp_tax, IndSp_p_fdr)

#write.csv(IndSp_p_fdr_merge, "R_output/IndicatorTaxa/2019.07.19_rhizo.rfy_nosingles_allOTUS_indsp_Var_perm999_Cultivar_tax_pvalAdjust.csv")


### Now filter for those OTUS that are present in at least 25% of the replicates 
IndSp_p_fdr_merge_25 <- filter(IndSp_p_fdr_merge, IndSp_p_fdr_merge$B >= 0.25)
dim(IndSp_p_fdr_merge_25) #687 

# write.csv(IndSp_p_fdr_merge_25, "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/SharedIndicator_taxa/2019.07.19_rhizo.rfy_nosing_allOTUS_indsp_Var_perm9999_tax_pvalAdj_25percRep.csv")

# get taxa table to merge with indicator species table 
rhizo.rfy_tax <- as.data.frame(as(tax_table(rhizo.rfy),"matrix"))
rhizo.rfy_tax <- rownames_to_column(rhizo.rfy_tax,"OTU")
IndSp_Var_Tax <- merge(IndSp_p_fdr_merge_25, rhizo.rfy_tax, by = "OTU")
dim(IndSp_Var_Tax) # 687 28 

#write.csv(IndSp_Var_Tax,"R_output/IndicatorTaxa/2019.07.19_rhizo.rfy_nosing_allOTUS_indsp_Var_perm9999_tax_pvalAdj_25percRep_Tax.csv")

# Get rid of indicator taxa that aren't sig @ alpha = 0.05 
IndSp_Var_Tax_sig <- filter(IndSp_Var_Tax, IndSp_Var_Tax$p_value_fdr < 0.05)
dim(IndSp_Var_Tax_sig) #685 29

#write.csv(IndSp_Var_Tax_sig,"R_output/IndicatorTaxa/2019.07.19_rhizo.rfy_nosing_allOTUS_indsp_Var_perm9999_tax_pvalAdj_25percRep_Tax_Pfilter.csv")


#Extract these OTUs as a vector 
IndSp_OTU <- as.vector(IndSp_Var_Tax_sig$OTU)
length(IndSp_OTU) #685


## Prune rel.abund of rhizo.rfy by indicator taxa 
rhizo.rfy_ra  = transform_sample_counts(rhizo.rfy, function(x) x / sum(x))
sample_sums(otu_table(rhizo.rfy_ra)) #sanity check - sample sum should equal 1

rhizo.rfy_ra.indsp <- prune_taxa(IndSp_OTU, rhizo.rfy_ra) #subset rhizo.rfy.ra by indsp taxa 
ntaxa(rhizo.rfy_ra.indsp) #685

# What % of sequences are indicator otus? 
rhizo.rfy.indsp <- prune_taxa(IndSp_OTU, rhizo.rfy) #subset rhizo.rfy.ra by indsp taxa 
ntaxa(rhizo.rfy.indsp) #683
sum(sample_sums(rhizo.rfy.indsp))/sum(sample_sums(rhizo.rfy)) *100 #21.19


#### What classes are the most abundant of the indicator taxa? (load summary stats function in previous code)

# Relative abundance of bacteria in rhizosphere 
plot_taxa_summary(rhizo.rfy.indsp, "Class") # mean relative abundance 
rhizo.rfy.indsp_summary <- summarize_taxa(rhizo.rfy.indsp, "Class")
#View(rhizo.rfy.indsp_summary)
# Top Classes = Acidobacteria (33%), Alphaproteobacteria (10%) and Deltaproteobacteria (7%)

#write.csv(rhizo.rfy.indsp_summary, "R_output/IndicatorTaxa/2019.07.19_rhizo.rfy_nosing_allOTUS_indsp_Var_perm9999_tax_pvalAdj_25percRep_ClassSummary.csv")


```

## Soil Fungi 
```{r}
# Remove anything less than something present in 1 samples
fungi.rfy.removeSing <- prune_taxa(taxa_sums(fungi.rfy) > 1, fungi.rfy)
ntaxa(fungi.rfy) - ntaxa(fungi.rfy.removeSing) #633
ntaxa(fungi.rfy.removeSing) # 3431



# multiplatt is the function used to conduct indicator species analysis
  # duleg = T makes it so we avoid considering grouping combinations
  # func = IndVal.g sets our test statistic as the INdVal index (returns sq.rt of IndVal)
  # indvalcomp=TRUE will give us the separate A and B parameters for the indicator species. A indicates the probablity that the Sample belongs to the identified sample given the species in the sample (specificity); B is the probablity of finding the species in the sample belonging to that group (sensitivity) 

library(indicspecies)

# data needs to have samples in rows and species in columns (opposite of phyloseq)
fungi.rfy.removeSing_otu_t <- as.data.frame(t(as(otu_table(fungi.rfy.removeSing), "matrix")))
dim(fungi.rfy.removeSing_otu_t) # 135 3431

sampledata_fungi.rfy.removeSing<- as.data.frame(as(sample_data(fungi.rfy.removeSing),"matrix"))


# by variety 
Var_IndSp <- multipatt(fungi.rfy.removeSing_otu_t,
                       sampledata_fungi.rfy.removeSing$Variety, duleg = TRUE, # duleg = TRUE so it doesn't consider other site combinations
                       control = how(nperm = 9999),func = "IndVal.g")


# send ISA summary results to a file 
 sink("R_output/IndicatorTaxa/2019.07.19_fungi.rfy_rmvSing_indsp_Var_perm9999.csv")
print(summary(Var_IndSp, indvalcomp=TRUE))
sink()
# editted this file into a table (output is funny)

#######################
### Indicator Padjust##

# Merge indicator species table with taxonomy 
IndSp <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/SharedIndicator_taxa/2019.07.19_fungi.rfy_rmvSing_indsp_Var_perm9999_clean.csv", header = TRUE)

fungi.Tax <- as.data.frame(as(tax_table(fungi.rfy),"matrix"))
library(tibble)
fungi.Tax <- rownames_to_column(fungi.Tax, "OTU")
IndSp_Tax <- merge(IndSp, fungi.Tax, by = "OTU")
dim(IndSp_Tax) #288 13

# List of pvalues
IndSp_Tax_p <-IndSp_Tax$p_value
# confirm that read list is a vector
is.vector(IndSp_Tax_p)
IndSp_Tax_p <- as.numeric(IndSp_Tax_p)

# Padjust with FDR 
library(stats)
IndSp_p_fdr <- as.matrix(p.adjust(IndSp_Tax_p, method = "fdr", length(IndSp_Tax_p)))
colnames(IndSp_p_fdr) <- "p_value_fdr"


# merge padjusted values with datatable
IndSp_Tax_fdr <- cbind(IndSp_Tax, IndSp_p_fdr)


### Now filter for those OTUS that are present in at least 25% of the replicates 

IndSp_Tax_fdr_25 <- filter(IndSp_Tax_fdr, IndSp_Tax_fdr$B >= 0.25)
dim(IndSp_Tax_fdr_25) #213

#write.csv(IndSp_Tax_fdr_25, "R_output/IndicatorTaxa/2019.07.19_fungi.rfy_nosingles_indsp_Var_perm9999_tax_pvalAdjust_25percRep.csv")

##################################
### Indicator Speices Taxonomy###

# I exported Indicator Species analysis results and organized to have a column that lists variety, percentreps across cultivars for the otu and taxonomy 

IndSp_Var <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/SharedIndicator_taxa/2019.07.19_fungi.rfy_nosingles_indsp_Var_perm9999_tax_pvalAdjust_25percRep.csv", header = TRUE)

#Extract these OTUs as a vector 
IndSp_OTU <- as.vector(IndSp_Var$OTU)
length(IndSp_OTU) #213 

## Prune rel.abund of fungi.rfy by indicator taxa 
fungi.rfy_ra  = transform_sample_counts(fungi.rfy, function(x) x / sum(x))
sample_sums(otu_table(fungi.rfy_ra)) #sanity check - sample sum should equal 1

fungi.rfy_ra.indsp <- prune_taxa(IndSp_OTU, fungi.rfy_ra) #subset fungi.rfy.ra by indsp taxa 
ntaxa(fungi.rfy_ra.indsp) #213

# What % of sequences are indicator otus? 
fungi.rfy.indsp <- prune_taxa(IndSp_OTU, fungi.rfy) #subset fungi.rfy.ra by indsp taxa 
ntaxa(fungi.rfy.indsp) #213
sum(sample_sums(fungi.rfy.indsp))/sum(sample_sums(fungi.rfy)) *100 #24.9749


#### What classes are the most abundant of the indicator taxa? (load summary stats function in previous code)
# Load summarize taxa function earlier in the code 

# Relative abundance of bacteria in rhizosphere 
plot_taxa_summary(fungi.rfy.indsp, "Class") # mean relative abundance 
fungi.rfy.indsp_summary <- summarize_taxa(fungi.rfy.indsp, "Class")
# Top classes:  Sordariomycetes (19%), Dothideomycetes (17%), and 27% were unclassified at class level

```
---

# Variation in N-fixation Potential
## Rel.Abund. of N-fixing Orders 
*Does the relative abundance of common N-fixing orders (Burkholderiales and Rhizobiales) differ among cultivars or sample types (roots vs. soils)?*
```{r}
##########################################################################
# Do the rel.abund. of soil N-fixers differ among cultivars/ecotype? (n = 12 cultivars )
#########################################################################

# transform data to relative abundance 
rhizo.rfy_ra <- transform_sample_counts(rhizo.rfy, function(x) x/sum(x))
# subset relative abundance for N-fixing orders 
Nfixers_ra<- subset_taxa(rhizo.rfy_ra, Order == "Rhizobiales" | Order == "Burkholderiales")
ntaxa(Nfixers_ra) #494 

Nfixers_ra.df <- psmelt(Nfixers_ra)

# Take mean relative abundance of Nfixer orders by sample
Nfixers_ra.df_mean <- Nfixers_ra.df %>% group_by(X.SampleID, Order) %>%  summarize_at("Abundance",mean)

Nfixers_ra.df_mean_meta <- merge(Nfixers_ra.df_mean, sampledata_rhizo.rfy[,c(1,2,3,7,4,10,19,32,39,40,51,72)], by = "X.SampleID")

# ANALYSIS
kruskal.test(Nfixers_ra.df_mean_meta$Abundance ~ Nfixers_ra.df_mean_meta$Variety)
	#Kruskal-Wallis rank sum test
#data:  Nfixers_ra.df_mean_meta$Abundance by Nfixers_ra.df_mean_meta$Variety
#Kruskal-Wallis chi-squared = 8.6892, df = 11, p-value = 0.6506

kruskal.test(Nfixers_ra.df_mean_meta$Abundance ~ Nfixers_ra.df_mean_meta$Ecotype)
#data:  Nfixers_ra.df_mean_meta$Abundance by Nfixers_ra.df_mean_meta$Ecotype
#Kruskal-Wallis chi-squared = 0.034911, df = 1, p-value = 0.8518

##########################################################################
# Do the rel.abund. of soil N-fixers differ among sample sites (root or soil for n = 4 cultivars)?
#########################################################################
rhizoendo_order <- tax_glom(rhizoendo, "Order")

rhizoendo_Nfixers <- subset_taxa(rhizoendo_order, Order == "Rhizobiales" | Order == "Burkholderiales")
ntaxa(rhizoendo_order) #230


rhizoendo_Nfixers.df <- psmelt(rhizoendo_Nfixers)

rhizoendo_Nfixers.df_sum <- rhizoendo_Nfixers.df %>% group_by(X.SampleID,SampleSite, Order) %>%  summarize_at("Abundance",sum)
dim(rhizoendo_Nfixers.df_sum) #176 4

metadata2 <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)

rhizoendo_Nfixers_ra.df_mean_meta <- merge(rhizoendo_Nfixers.df_sum, metadata2[,c(1,2,3,7,10,4,19,32,39,40,51,72)], by = "X.SampleID")


kruskal.test(rhizoendo_Nfixers_ra.df_mean_meta$Abundance ~ rhizoendo_Nfixers_ra.df_mean_meta$SampleSite)
#Kruskal-Wallis rank sum test
#rhizoendo_Nfixers_ra.df_mean_meta$SampleSite
#Kruskal-Wallis chi-squared = 43.201, df = 1, p-value = 4.938e-11

rhizoendo_Nfixers_ra.df_mean_meta %>% group_by(SampleSite) %>%    summarise( mean(Abundance)) 
# Endosphere              305. 
#Rhizosphere              92.8

```
## Prop. of putative N-fixers (PICRUST)
```{r}
# We used a closed refrence Green-genes referenced OTU table with PICRUSt (Langille et al. 2013) to predict the relative abundance of N-fixers. 
 # We first normalized all OTUs by their predicted 16S rRNA gene copy number, which provides a pseudo-abundance estimate for each OTU.
 # Then we used ‘metagenome_predictions’ to obtain OTU-specific gene counts for N-fixation using the following KEGG pathway orthologs K02588, K02586, K02591, K00531. 
 # We calculated each samples’ predicted proportion of N-fixation genes by dividing the number of OTUs with at least one predicted N-fixation pathway for each sample by the normalized abundance of OTUS (e.g., the total 16S-gene normalized OTU counts).


# Read in Nfix contributions from 'metagenome_predictions' output 
Rhizo_Nfix <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Nfix_PiCRUST/Rhizo_NitrogenFixation_contributions.csv", header = TRUE)

# Take the sum of 'CountContributedByOTU' for each sample
Rhizo_Nfix_CountCont_SampSum <- Rhizo_Nfix %>% group_by(X.SampleID,OTU)%>% summarise(Sum_CountContbyOTU = sum(CountContributedByOTU))

# Reshape so that each column is a sample with the contribution sum from all 4 gene pathways for each otu. 
# Then make anything >1 == 1, we are interested in the presence of at least 1 nifH gene and are not accounting for differnces in the # of nifH genes (e.g. 1-4)
library(reshape)
x2 <- cast(Rhizo_Nfix_CountCont_SampSum, OTU ~ X.SampleID)
x2[is.na(x2)] <- 0 # make nas = 0 
x2[,-c(1)][x2[,-c(1)] >1] <- 1 #make anything >1 = 1 

# rename columns 
Nfix_OTUcounts <- as.data.frame(colSums(x2))
names(Nfix_OTUcounts)[1] <- "Nfix_OTUcounts"
Nfix_OTUcounts <- rownames_to_column(Nfix_OTUcounts, "X.SampleID")

# Read in 16S copy #-normalized OTU table (result from 'normalize_by_copy_number')
Rhizo_NormalizedOTUS<- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGvariety/Nfix_PiCRUST/Rhizo_OTU_corrected.csv", header = TRUE)
dim(Rhizo_NormalizedOTUS) #7905 OTUs

# find the normalized OTU sequence sum for each sample
Rhizo_sequenceSum <- as.data.frame(colSums(Rhizo_NormalizedOTUS[2:ncol(Rhizo_NormalizedOTUS)]))

dim(Rhizo_sequenceSum) #138
Rhizo_sequenceSum <- rownames_to_column(Rhizo_sequenceSum, "X.SampleID")
names(Rhizo_sequenceSum)[2] <- "NormalizedSeqSum"

# Now merge the total # of normalized sequences (e.g. pseudo taxa abundance) for each sample with # presence/absence of Nfixers 
Rhizo_NfixOTUcounts_NormSeqSum <- merge(Rhizo_sequenceSum,Nfix_OTUcounts, by = "X.SampleID", all = TRUE)

# Divide the # of OTUs with Nfix pathway by the total number of normalized "OTU individuals"/sequences per sample 
Rhizo_NfixOTUcounts_NormSeqSum$Sum_NfixCount_Normalized <- Rhizo_NfixOTUcounts_NormSeqSum$Nfix_OTUcounts/Rhizo_NfixOTUcounts_NormSeqSum$NormalizedSeqSum

# merge with metadata 
Rhizo_NfixOTUcounts_NormSeqSum_meta <- merge(Rhizo_NfixOTUcounts_NormSeqSum, metadata, by = "X.SampleID")
dim(Rhizo_NfixOTUcounts_NormSeqSum_meta) # 138 

###################################################################
## Does the predicted rel.abund. of soil Nfixers differ among cultivars?
####################################################################

# ANOVA 
# normal residuals?
ggqqplot(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized)
hist(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized)

# shapiro test - null hypothesis: residuals are normally distributed 
Nfix_rhizo_model = lmer(Sum_NfixCount_Normalized ~ Variety +  (1|Block/Variety) , data = Rhizo_NfixOTUcounts_NormSeqSum_meta)
res = resid(Nfix_rhizo_model)
shapiro.test(res) 
# 0.63

# homoegeneity of variances
bartlett.test(Sum_NfixCount_Normalized ~ Variety, data = Rhizo_NfixOTUcounts_NormSeqSum_meta) 
# p = 0.16

plot(Nfix_rhizo_model)

require(lme4)
require(lmerTest)

# Anova
anova(Nfix_rhizo_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)    
#Variety 0.0072158 0.00065598    11 32.499  5.9558 3.247e-05 ***


# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(Nfix_rhizo_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

#  Variety      lsmean      SE   df lower.CL upper.CL .group
# EG1102        0.104 0.00495 11.7   0.0861    0.121  a    
# NE28          0.114 0.00495 11.7   0.0960    0.131  ab   
# Kanlow        0.118 0.00495 11.7   0.1000    0.135   bc  
# Alamo         0.119 0.00495 11.7   0.1010    0.136   bcd 
# EG2101        0.123 0.00495 11.7   0.1053    0.140   bcde
# Southlow      0.124 0.00495 11.7   0.1062    0.141   bcde
# Shelter       0.128 0.00515 13.6   0.1104    0.146    cde
# Cave-in-Rock  0.129 0.00515 13.6   0.1109    0.146    cde
# EG1101        0.129 0.00495 11.7   0.1117    0.147    cde
# Blackwell     0.131 0.00495 11.7   0.1138    0.149     de
# Trailblazer   0.133 0.00515 13.6   0.1149    0.150      e
# Dacotah       0.136 0.00495 11.7   0.1183    0.153      e


########################################
# Does Nfix differ by Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Nfix_Rhizo_model_eco <- lmer(Sum_NfixCount_Normalized ~ Ecotype + (1|Block/Plot), data = Rhizo_NfixOTUcounts_NormSeqSum_meta)
plot(Nfix_Rhizo_model_eco)

# shapiro test 
res = resid(Nfix_Rhizo_model_eco)
shapiro.test(res) # p= 0.7

# homoegeneity of variances
bartlett.test(Sum_NfixCount_Normalized ~ Ecotype, data = Rhizo_NfixOTUcounts_NormSeqSum_meta) 
# p = 0.1


anova(Nfix_Rhizo_model_eco, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF  DenDF F value   Pr(>F)   
#Ecotype 0.0010218 0.0010218     1 42.637  9.3103 0.003911 **


emmeans(Nfix_Rhizo_model_eco, pairwise~Ecotype, adjust = "fdr")
# Ecotype emmean      SE   df lower.CL upper.CL
# Lowland  0.117 0.00398 5.70    0.107    0.127
# Upland   0.127 0.00355 3.64    0.117    0.137#

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 

#$contrasts
# contrast         estimate      SE   df t.ratio p.value
# Lowland - Upland -0.00971 0.00318 42.2 -3.051  0.0039 



#### IS predicted # of Nfixers correlated with soil Nfixation and soil N measured by paired study
#Roley et al. 2020 Phytobiomes https://doi.org/10.1094/PBIOMES-11-19-0064-FI 
  
cor.test(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized, Rhizo_NfixOTUcounts_NormSeqSum_meta$sfix_rate_gN_d, method = "pearson")

#t = -0.78528, df = 136, p-value = 0.4337
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval: -0.2316884  0.1010547
#sample estimates:   cor -0.06718471 


### Does Nfixer abundance correlate with nitrate availability?

cor.test(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized, Rhizo_NfixOTUcounts_NormSeqSum_meta$NO3_ugN_g_drysoil_K2SO4, method = "pearson")

#data:  Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized and Rhizo_NfixOTUcounts_NormSeqSum_meta$NO3_ugN_g_drysoil_K2SO4
#t = -4.0891, df = 136, p-value = 7.373e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.4718971 -0.1733656
#sample estimates:
#       cor -0.3308846 


```
---
# Microbial Biomass Carbon  
## Normality & Stats
```{r}
library("ggpubr")
# test for normal residuals and equal variances (ANOVA assumptions)

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
hist(metadata$ugC_MicBiomass_g_dry_soil) 


# Anova - Does MBC differ by Variety? should SMC be in the model?  
require(lmerTest)
MBC_Variety_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Variety + (1|Block/Variety), data = metadata)
MBC_Variety_SMC_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)
MBC_Variety_NO3_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Variety + NO3_ugN_g_drysoil_K2SO4_SMCOutRmv + (1|Block/Variety), data = metadata)
#boundary (singular) fit: see ?isSingular

isSingular(MBC_Variety_SMC_model) # TRUE

plot(MBC_Variety_SMC_model) # residuals look okay 

# normal residuals?
res = resid(MBC_Variety_SMC_model)
shapiro.test(res) 
# MBC - p = 0.04


#  homoegeneity of variances?
bartlett.test(ugC_MicBiomass_g_dry_soil~ Variety, data = metadata)
# MBC, p = 0.88

MBC_Variety_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Variety +  (1|Block/Variety), data = metadata)

# normal residuals?
res = resid(MBC_Variety_model)
shapiro.test(res) 
# MBC - p = 0.3


## Is the model improved with SMC as covariate?

AIC(MBC_Variety_NO3_model, MBC_Variety_model)
#         df      AIC
##MBC_Variety_SMC_model 16 1434.961 # YES
#MBC_Variety_model     15 1473.867
#MBC_Variety_NO3_model 16 1440.766

anova(MBC_Variety_SMC_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Variety                    251333   22848    11 123.22  7.8478 3.349e-10 ***
#SMC_perc_GravWaterCont      45910   45910     1 113.30 15.7689 0.0001261 ***


# tukeys comparison will give letter differentiation
library(multcompView)
library(lsmeans)
marginal = lsmeans(MBC_Variety_SMC_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

#Variety      lsmean   SE   df lower.CL upper.CL .group
#  Cave-in-Rock   96.8 18.3 18.3     36.9      157  a    
# Trailblazer    99.9 19.0 20.7     38.7      161  a    
# Southlow      108.0 18.5 18.7     47.7      168  a    
# EG1102        129.7 18.8 19.5     68.7      191  ab   
# Dacotah       146.5 20.7 25.8     81.5      212  abc  
# Blackwell     151.4 20.0 24.3     88.1      215  abcd 
# Kanlow        180.1 19.3 20.9    117.9      242   bcde
# NE28          190.3 18.4 18.5    130.1      250    cde
# Shelter       194.8 18.6 18.9    134.3      255    cde
# Alamo         204.3 18.5 18.7    144.0      265    cde
# EG2101        205.1 19.2 21.1    143.5      267     de
# EG1101        219.9 18.4 18.6    159.7      280      e

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 
#Conf-level adjustment: bonferroni method for 12 estimates 
#P value adjustment: fdr method for 66 tests 
#significance level used: alpha = 0.05

###########################################################################
# does MBC difer by ecotype? 
###########################################################################
MBC_Ecotype_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot),data = metadata)

# check normality of residuals
res = MBC_Ecotype_model$resdiduals   
shapiro.test(res) # p = 0.04

anova(MBC_Ecotype_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#                                  Sum Sq Mean Sq NumDF   DenDF F value    Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)
# Ecotype                     18328   18328     1  43.492  6.0636 0.017834 * 
#SMC_perc_GravWaterCo         30304   30304     1 125.527 10.0259 0.001938 **
emmeans(MBC_Ecotype_model, pairwise~Ecotype)
##cotype emmean   SE   df lower.CL upper.CL
# Lowland    185 15.6 8.22      150      221
# Upland     148 13.0 4.08      112      183
 

```

###  MBC correlation with SMC or SRL or root biomass?
```{r}
# Pearson Correlation 

# Specific Root Length 
ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(log(metadata$SRL_giaroots_OutlierRmv))


cor.test(metadata$ugC_MicBiomass_g_dry_soil, log(metadata$SRL_giaroots_OutlierRmv), method = "pearson")

#t = 0.91581, df = 136, p-value = 0.3614
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:-0.08999299  0.24222511
#sample estimates:      cor 0.07828918 

######################
# avg. root diam 

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(metadata$Avg_root_width_diam)


cor.test(metadata$ugC_MicBiomass_g_dry_soil, metadata$Avg_root_width_diam, method = "pearson")

#t = -0.77428, df = 139, p-value = 0.4401
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:-0.2283705  0.1008726
#sample estimates:       cor -0.06553232 

plot(Avg_root_width_diam ~ ugN_MicBiomass_g_dry_soil, data = metadata)

##############
# Root Biomass 

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(log(metadata$DryRootWt_total_g))


cor.test(metadata$ugC_MicBiomass_g_dry_soil, log(metadata$DryRootWt_total_g), method = "pearson")

#data:  metadata$ugC_MicBiomass_g_dry_soil and log(metadata$DryRootWt_total_g)
#t = 2.4101, df = 139, p-value = 0.01726
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.03616484 0.35387542 
#sample estimates: cor 0.2002796 

plot( DryRootWt_total_g~ ugC_MicBiomass_g_dry_soil, data = metadata)

################
# Root length  

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(sqrt(metadata$Network_length))


cor.test(metadata$ugC_MicBiomass_g_dry_soil, sqrt(metadata$Network_length), method = "pearson")


#data:  metadata$ugC_MicBiomass_g_dry_soil and sqrt(metadata$Network_length)
#t = 2.7226, df = 139, p-value = 0.007308
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval: 0.06199804 0.37631873
#sample estimates:   cor 0.2250041 

plot( Network_length~ ugC_MicBiomass_g_dry_soil, data = metadata)

#############

#############
# MBC correlate with SMC? 

cor.test(metadata$ugC_MicBiomass_g_dry_soil, (metadata$SMC_perc_GravWaterCont_OutlierRmv), method = "pearson")

##data:  metadata$ugC_MicBiomass_g_dry_soil and (metadata$SMC_perc_GravWaterCont_OutlierRmv)
#t = 4.6448, df = 137, p-value = 7.894e-06
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2155932 0.5043782
#sample estimates:
#      cor 0.3688535

###########################
# MBC correlate with NO3? 

ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4)
cor.test(metadata$ugC_MicBiomass_g_dry_soil, (metadata$NO3_ugN_g_drysoil_K2SO4), method = "pearson")

#data:  metadata$ugC_MicBiomass_g_dry_soil and (metadata$NO3_ugN_g_drysoil_K2SO4)
#t = 3.1025, df = 139, p-value = 0.002324
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.09309429 0.40285390
#sample estimates:
 #     cor 0.2544899 


```
# Microbial Biomass Nitrogen 
## Normality & Stats
```{r}
library("ggpubr")
# test for normal residuals and equal variances (ANOVA assumptions)

ggqqplot(metadata$ugN_MicBiomass_g_dry_soil)
hist(metadata$ugN_MicBiomass_g_dry_soil) 

# Anova - Does MBN differ by Variety? should SMC be in the model?  
require(lme4)
require(lmerTest)
MBN_Variety_SMC_model <- lmer(ugN_MicBiomass_g_dry_soil ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)


# normal residuals?
res = resid(MBN_Variety_SMC_model)
shapiro.test(res) 
# MBN - p = 0.4

#  homoegeneity of variances?
bartlett.test(ugN_MicBiomass_g_dry_soil~ Variety, data = metadata)
# MBN, p = 0.63

MBN_Variety_model <- lmer(ugN_MicBiomass_g_dry_soil ~ Variety +  (1|Block/Variety), data = metadata)

# normal residuals?
res = resid(MBN_Variety_model)
shapiro.test(res) 
# MBN - p = 0.38



AIC(MBN_Variety_SMC_model, MBN_Variety_model)
#   df      AIC
#MBN_Variety_SMC_model 16 934.3634 # YES 
#MBN_Variety_model     15 951.8143

anova(MBN_Variety_SMC_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Variety                     9562.3  869.30    11 120.19 13.6107 <2e-16 ***
#SMC_perc_GravWaterCont_OutlierRmv  144.8  144.84     1 122.09  2.2677 0.1347 

# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(MBN_Variety_SMC_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

#Variety      lsmean   SE   df lower.CL upper.CL .group
# Trailblazer    13.4 3.43 14.4     1.69     25.0  a    
# Southlow       13.6 3.34 13.1     2.05     25.1  ab   
# Cave-in-Rock   16.0 3.21 11.5     4.51     27.4  ab   
# Dacotah        19.2 3.44 14.3     7.44     30.9  ab   
# NE28           21.1 3.15 10.6     9.65     32.5  ab   
# Blackwell      21.6 3.18 10.9    10.09     33.0   b   
# EG1101         29.4 3.15 10.6    17.99     40.9    c  
# Shelter        30.4 3.16 10.8    18.94     41.8    c  
# EG2101         31.1 3.24 11.8    19.65     42.6    c  
# EG1102         32.4 3.19 11.1    20.94     43.9    c  
# Alamo          36.5 3.15 10.7    25.05     47.9    cd 
# Kanlow         40.0 3.27 11.9    28.48     51.6     d 


#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 
#Conf-level adjustment: bonferroni method for 12 estimates 
#P value adjustment: fdr method for 66 tests 
#significance level used: alpha = 0.05 

###########################################################################
# does MBN difer by ecotype? 
###########################################################################
MBN_Ecotype_model <- lmer(ugN_MicBiomass_g_dry_soil ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot),data = metadata)

# check normality of residuals
res = resid(MBN_Ecotype_model)
shapiro.test(res) # p = 0.5

anova(MBN_Ecotype_model, type = "III") 

#Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#Ecotype                   2648.9  2648.9     1  44.447 37.9547 1.886e-07 ***
#SMC_perc_GravWaterCont     100.8   100.8     1 127.486  1.4443    0.2317   



emmeans(MBN_Ecotype_model, pairwise~Ecotype)
#  Ecotype emmean   SE   df lower.CL upper.CL
# Lowland   34.7 2.66 5.85     28.2     41.2
# Upland    21.0 2.35 3.65     14.3     27.8

  #contrast         estimate   SE   df t.ratio p.value
 #Lowland - Upland     14.9 2.18 42.6 6.828   <.0001 

```
---
# Edaphic Variables by Variety, Block and Date 
## Block
```{r}

##### soil moisture content (gravimetric, SMC) ###########
##########################################################
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)
hist(log(metadata$SMC_perc_GravWaterCont_OutlierRmv))

model = lm(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Block, data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = model$residuals
shapiro.test(res)
# SMC- no outlier p = 0.5
plot(model)

# 2) homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Block , data = metadata)
# SMC- no outliers - p <0.001
# log SMC = 0.19 


# Does SMC differ by BLock
SMC_Block_model <- lm(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Block, data = metadata)
SMC_Block_model
Anova(SMC_Block_model, type = "III") 
#Anova Table (Type III tests)
#Response: SMC_perc_GravWaterCont_OutlierRmv
#   Sum Sq  Df  F value    Pr(>F)    
#(Intercept) 170.117   1 2379.357 < 2.2e-16 ***
#Block         2.203   3   10.273 3.785e-06 ***
#Residuals     9.867 138                         

# pairwise test (gives you p.values for every combination)
emmeans(SMC_Block_model, pairwise ~ Block, adjust = "fdr")
#$contrasts
# contrast     estimate    SE  df t.ratio p.value
#Four - One     0.2429 0.0635 138  3.827  0.0004 
# Four - Three   0.3172 0.0639 138  4.963  <.0001 
# Four - Two     0.2881 0.0635 138  4.540  <.0001 
# One - Three    0.0743 0.0635 138  1.171  0.3655 
# One - Two      0.0453 0.0630 138  0.718  0.5687 
# Three - Two   -0.0291 0.0635 138 -0.458  0.6478 
plot(metadata$SMC_perc_GravWaterCont_OutlierRmv ~ metadata$Block)

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(SMC_perc_GravWaterCont_OutlierRmv))


SMC_blockmean <- metadata2 %>%  group_by(Block) %>% summarise_at("SMC_perc_GravWaterCont_OutlierRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
SMC_blockmean
#Block  mean   var    sd   min   max
#  <fct> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 Four   9.45  7.71  2.78  4.82 17.6 
#2 One    7.34  3.44  1.86  3.94 11.3 
#3 Three  6.87  3.59  1.89  3.48 10.0 
#4 Two    6.94  2.08  1.44  4.59  9.93

######################################
##### Soil nitrate (NO3) #############

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv)
hist(log(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1))

model = lm(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = model$residuals
shapiro.test(res)
# NO3 p = 0<0.001
# log+1 NO3 p = 0.12
plot(model)

# 2) homoegeneity of variances
bartlett.test(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block , data = metadata)
# NO3, p = 0.96

# Does NO3 differ by BLock
NO3_Block_model <- lm(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
NO3_Block_model
Anova(NO3_Block_model, type ="III") 
 # Anova Table (Type III tests)
#Response: log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv + 1)
#             Sum Sq  Df   F value Pr(>F)    
#(Intercept) 24.3327   1 1129.6337 <2e-16 ***
#Block        0.0335   3    0.5182 0.6704    
#Residuals    2.9726 138  

# pairwise test (gives you p.values for every combination)
emmeans(NO3_Block_model, pairwise ~ Block, adjust = "fdr")
# contrast     estimate     SE  df t.ratio p.value
# Four - One    0.00929 0.0348 138  0.267  0.8144 
# Four - Three  0.03830 0.0351 138  1.092  0.8129 
# Four - Two    0.03011 0.0348 138  0.864  0.8129 
# One - Three   0.02901 0.0348 138  0.833  0.8129 
# One - Two     0.02082 0.0346 138  0.602  0.8144 
# Three - Two  -0.00819 0.0348 138 -0.235  0.8144 

plot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv ~ metadata$Block)

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv))


NO3_blockmean <- metadata2 %>%  group_by(Block) %>% summarise_at("NO3_ugN_g_drysoil_K2SO4_SMCOutRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
NO3_blockmean
#  Block  mean   var    sd   min   max
#  <fct> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 Four   1.32 0.109 0.330  0.74  2.21
#2 One    1.31 0.119 0.345  0.77  2.12
#3 Three  1.24 0.116 0.341  0.68  1.97
#4 Two    1.26 0.115 0.340  0.73  2.09


##############################
##### Soil Ammonium (NH4) ####

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)
hist(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1))

model = lm(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = model$residuals
shapiro.test(res)
# log + 1 NH4 p <0.001

# 2) homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block , data = metadata)
# NH4 p<0.001
# log + 1 p = 0.09 

plot(model)

# Does NH4 differ by BLock
NH4_Block_model <- lm(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
NH4_Block_model
Anova(NH4_Block_model, type = "III") 
  
#Response: log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv + 1)
#            Sum Sq  Df  F value Pr(>F)    
#(Intercept) 43.172   1 863.9327 <2e-16 ***
#Block        0.345   3   2.3015 0.0799 .  
#Residuals    6.896 138    

# pairwise test (gives you p.values for every combination)
emmeans(NH4_Block_model, pairwise ~ Block, adjust = "fdr")
#$contrasts
# Four - One     0.0565 0.0531 138  1.065  0.3467 
# Four - Three   0.1386 0.0534 138  2.593  0.0632 
# Four - Two     0.0787 0.0531 138  1.484  0.2803 
# One - Three    0.0821 0.0531 138  1.547  0.2803 
# One - Two      0.0222 0.0527 138  0.422  0.6736 
# Three - Two   -0.0598 0.0531 138 -1.127  0.3467 


plot(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv ~ metadata$Block)


```
## Variety
```{r}

##############################
##### Soil Moisture Content####
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)
hist(log(metadata$SMC_perc_GravWaterCont_OutlierRmv))

model = lm(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Variety , data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# SMC- no outlier p = 0.5
plot(model)


# 2) homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Variety , data = metadata)
# log SMC by Var p >0.05 

# Does SMC differ by Variety
SMC_Var_Block_model <- lmer(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Variety + (1|Block/Variety), data = metadata)
SMC_Var_Block_model
anova(SMC_Var_Block_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Variety 3.3854 0.30777    11 33.31  7.8595 1.672e-06 ***                     


# pairwise test (gives you p.values for every combination)
emmeans(SMC_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(SMC_Var_Block_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD
#Variety      lsmean     SE    df lower.CL upper.CL .group 
# Dacotah        1.55 0.0957 10.06     1.20     1.90  a     
# Blackwell      1.84 0.0939  9.38     1.49     2.20   b    
# EG2101         1.88 0.0957 10.06     1.52     2.23   bc   
# Shelter        1.89 0.0939  9.38     1.53     2.24   bcd  
# EG1101         1.91 0.0939  9.38     1.56     2.26   bcde 
# Cave-in-Rock   2.01 0.0939  9.38     1.66     2.37   bcdef
# Trailblazer    2.05 0.0939  9.38     1.70     2.40   bcdef
# NE28           2.09 0.0939  9.38     1.73     2.44    cdef
# Southlow       2.11 0.0939  9.38     1.76     2.47     def
# Alamo          2.12 0.0939  9.38     1.76     2.47      ef
# EG1102         2.19 0.0939  9.38     1.84     2.54       f
# Kanlow         2.23 0.0939  9.38     1.88     2.58       f


#### BY ecotype? 
SMC_Ecotype_model <- lmer(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Ecotype + (1|Block/Plot),data = metadata)

# check normality of residuals
res = resid(SMC_Ecotype_model)
shapiro.test(res) # p = 0.9

anova(SMC_Ecotype_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#Ecotype 0.34868 0.34868     1 42.515  8.8729 0.004767 **

emmeans(SMC_Ecotype_model, pairwise~Ecotype)
#Ecotype emmean     SE   df lower.CL upper.CL
# Lowland   2.11 0.0815 5.32     1.91     2.32
# Upland    1.93 0.0735 3.53     1.71     2.14

#$contrasts
# contrast         estimate     SE   df t.ratio p.value
# Lowland - Upland    0.183 0.0616 42.7 2.979   0.0048 


##############################
##### Soil Nitrate (NO3) #####

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1)
hist(log(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1))

model = lmer(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# NO3 p = 0.004

# 2) homoegeneity of variances
bartlett.test(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety , data = metadata)


# Does NO3 differ by Variety
NO3_Var_Block_model <- lmer(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
 # singular fit
anova(NO3_Var_Block_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF  DenDF F value   Pr(>F)    
#Variety 0.84872 0.077156    11 36.175  11.625 7.57e-09 ***

# pairwise test (gives you p.values for every combination)
emmeans(NO3_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(NO3_Var_Block_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD
#Variety      lsmean     SE    df lower.CL upper.CL .group 
# Dacotah       0.653 0.0359 38.4    0.543    0.762  a     
# EG2101        0.666 0.0359 38.4    0.557    0.776  ab    
# EG1101        0.717 0.0351 35.5    0.610    0.825  ab    
#Shelter       0.743 0.0351 35.5    0.636    0.851  abc   
# Blackwell     0.760 0.0351 35.5    0.652    0.867  abc   
# Southlow      0.761 0.0351 35.5    0.653    0.868  abc   
# Trailblazer   0.778 0.0351 35.5    0.670    0.885   bcd  
# Cave-in-Rock  0.836 0.0351 35.5    0.729    0.944    cde 
# EG1102        0.891 0.0351 35.5    0.783    0.998     def
# Alamo         0.947 0.0351 35.5    0.839    1.054      ef
# NE28          0.992 0.0351 35.5    0.885    1.100       f
# Kanlow        0.999 0.0351 35.5    0.892    1.107       f


#### BY ecotype? 
NO3_Ecotype_model <- lmer(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Ecotype + (1|Block/Plot),data = metadata)
# is singular 

# check normality of residuals
res = resid(NO3_Ecotype_model)
shapiro.test(res) # p = 0.2

anova(NO3_Ecotype_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF  DenDF F value   Pr(>F)   
#Ecotype 0.063818 0.063818     1 45.701  9.5923 0.003335 **

emmeans(NO3_Ecotype_model, pairwise~Ecotype)
#Ecotype emmean     SE    df lower.CL upper.CL
# Lowland  0.889 0.0302 20.99    0.826    0.951
# Upland   0.774 0.0214  6.65    0.723    0.825
##$contrasts
## contrast         estimate    SE   df t.ratio p.value
# Lowland - Upland    0.115 0.037 42.9 3.097   0.0034



##############################
##### Soil Ammonium (NH4) ####

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv_SMCOutRmv)
hist(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv))

model = lmer(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
# shapiro test - null hypothesis: residuals are normally distributed 
res = resid(model)
shapiro.test(res)
# NH4 p = 0.002

# 2) homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety , data = metadata)


plot(model)


# Does NH4 differ by Variety
NH4_Var_Block_model <- lmer(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
 # singular fit
anova(NH4_Var_Block_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Variety 0.25777 0.023434    11 33.117  0.5289 0.8695

# pairwise test (gives you p.values for every combination)
emmeans(NH4_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(NH4_Var_Block_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD


#### BY ecotype? 
NH4_Ecotype_model <- lmer(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Ecotype + (1|Block/Plot),data = metadata)

# check normality of residuals
res = resid(NH4_Ecotype_model)
shapiro.test(res) # p = 0.2

plot(NH4_Ecotype_model)

anova(NH4_Ecotype_model, type = "III") 
#Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Ecotype 0.044921 0.044921     1 42.772   1.014 0.3196#

emmeans(NH4_Ecotype_model, pairwise~Ecotype)
# Ecotype emmean     SE    df lower.CL upper.CL
# Lowland   1.01 0.0408 11.85    0.923     1.10
# Upland    1.06 0.0321  4.86    0.974     1.14
##$contrasts
# contrast         estimate     SE   df t.ratio p.value
# Lowland - Upland  -0.0446 0.0443 42.5 -1.007  0.3197 

```

## Sampling Date  
```{r}
##############################
##### Soil Moisture Content###

# test for normal residuals and equal variances (ANOVA assumptions)
# normal residuals?
ggqqplot(log(metadata$SMC_perc_GravWaterCont_OutlierRmv))
hist(log(metadata$SMC_perc_GravWaterCont_OutlierRmv))

# shapiro test - null hypothesis: residuals are normally distributed 
model = lm(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Date, data =   metadata)
res = model$residuals
shapiro.test(res) 
# p = 0.9 

# homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Date, data = metadata) 
# p = 0.7 

Anova(model, type = "III")
#Response: log(SMC_perc_GravWaterCont_OutlierRmv)
#            Sum Sq  Df  F value    Pr(>F)    
#(Intercept) 125.50   1 2122.506 < 2.2e-16 ***
#Date          3.91   3   22.043 1.009e-11 ***

emmeans(model, pairwise~Date, adjust = "fdr")
 #contrast        estimate     SE  df t.ratio p.value
 #13-Jul - 20-Jul   -0.171 0.0540 138 -3.172  0.0028 
 #13-Jul - 27-Jul   -0.286 0.0577 138 -4.962  <.0001 
 #13-Jul - 28-Jun    0.200 0.0653 138  3.058  0.0032 
 #20-Jul - 27-Jul   -0.115 0.0536 138 -2.145  0.0337 
 #20-Jul - 28-Jun    0.371 0.0617 138  6.017  <.0001 
 #27-Jul - 28-Jun    0.486 0.0649 138  7.488  <.0001 

plot(SMC_perc_GravWaterCont_OutlierRmv ~Date_Julian, metadata)

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(SMC_perc_GravWaterCont_OutlierRmv))


SMC_blockmean <- metadata2 %>%  group_by(VarRep) %>% summarise_at("SMC_perc_GravWaterCont_OutlierRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
#Date    mean   var    sd   min   max
#  <fct>  <dbl> <dbl> <dbl> <dbl> <dbl>
#1 13-Jul  6.80  2.59  1.61  4.59 12.2 
#2 20-Jul  8.13  4.20  2.05  4.19 12.7 
#3 27-Jul  9.11  5.52  2.35  5.49 17.6 
#4 28-Jun  5.62  2.14  1.46  3.48  7.97


metadata2  %>%  group_by(VarRep) %>% summarise_at("pct_water_soil", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
#<fct>  <dbl> <dbl> <dbl> <dbl> <dbl>
#1 13-Jul  15.8 0.533 0.730  13.9  16.7
#2 20-Jul  16.1 0.678 0.824  14.3  17.4
#3 27-Jul  15.4 0.962 0.981  12.9  16.5
#4 28-Jun  14.8 0.283 0.532  14.0  15.7


# as correlation with julian date 
## SMC 
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)

cor.test(metadata$SMC_perc_GravWaterCont_OutlierRmv, metadata$Date_Julian, method = "pearson")
cor.test(metadata$pct_water_soil, metadata$Date_Julian, method = "pearson")

#t = 7.2232, df = 140, p-value = 2.975e-11
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval: 0.3897837 0.6315732
#sample estimates:     cor 0.5210551 

plot(SMC_perc_GravWaterCont_OutlierRmv ~ Date_Julian, data = metadata)

########################
### Soil Nitrate (NO3)##

# normal residuals?
ggqqplot((metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv))
hist((metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv))

# shapiro test - null hypothesis: residuals are normally distributed 
model = lm(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv ~ Date, data =   metadata)
res = model$residuals
shapiro.test(res) 
# p <0.001 
plot(model)

# homoegeneity of variances
bartlett.test((NO3_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data = metadata) 
# p = 0.23

Anova(model, type = "III")
#
#Response: NO3_ugN_g_drysoil_K2SO4_SMCOutRmv
#            Sum Sq  Df F value    Pr(>F)    
#(Intercept) 38.231   1 561.420 < 2.2e-16 ***
#Date         6.610   3  32.356 6.695e-16 ***
#Residuals    9.397 138 

emmeans(model, pairwise~Date, adjust = "fdr")
# 13-Jul - 20-Jul  -0.2894 0.0580 138 -4.990  <.0001 
# 13-Jul - 27-Jul  -0.5412 0.0619 138 -8.737  <.0001 
# 13-Jul - 28-Jun  -0.0105 0.0700 138 -0.150  0.8810 
# 20-Jul - 27-Jul  -0.2518 0.0575 138 -4.377  <.0001 
# 20-Jul - 28-Jun   0.2789 0.0662 138  4.215  0.0001 
# 27-Jul - 28-Jun   0.5307 0.0697 138  7.619  <.0001 

# Summarize means by date 
metadata2  %>%  group_by(Date) %>% summarise_at("NO3_ugN_g_drysoil_K2SO4_SMCOutRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
#Date    mean    var    sd   min   max
#  <fct>  <dbl>  <dbl> <dbl> <dbl> <dbl>
#1 13-Jul  1.05 0.0458 0.214  0.73  1.65
#2 20-Jul  1.33 0.0733 0.271  1.01  2.21
#3 27-Jul  1.59 0.0634 0.252  1.13  2.12
#4 28-Jun  1.06 0.0988 0.314  0.68  1.67

## NO3
ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv)

cor.test(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv, metadata$Date_Julian, method = "pearson")

#data:  metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv and metadata$Date_Julian
#t = 7.7787, df = 140, p-value = 1.45e-12
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:0.4228780 0.6548113
#sample estimates:     cor 0.5493386

##########################
### Soil ammonium (NH4)###

# normal residuals?
ggqqplot((metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv))
hist(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv))

# shapiro test - null hypothesis: residuals are normally distributed 
model = lm(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data =   metadata)
res = model$residuals
shapiro.test(res) 
# p <0.0001 
# log p = 0.02 
plot(model)

# homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data = metadata) 
# log p = 0.03

Anova(model, type = "III")
#
#Response: log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)
#             Sum Sq  Df  F value Pr(>F)    
#(Intercept) 13.7638   1 118.0062 <2e-16 ***
#Date         0.4769   3   1.3629 0.2568    
#Residuals   16.0958 138 

emmeans(model, pairwise~Date, adjust = "fdr")

## NH4
ggqqplot(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)

cor.test(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv), metadata$Date_Julian, method = "pearson")

#data:  log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv) and metadata$Date_Julian
#t = -1.505, df = 140, p-value = 0.1346
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:-0.28498306  0.03936631
#sample estimates:  cor -0.1261792

```
---
# Root Trait Analysis 
### Check Correlations among traits 
```{r}
# SRLgia = volume-weighted SRL 
# SRL = mass-weighted SRL 

##########################
###SRLgia ~ SRL ##########
# Specific Root Length 
ggqqplot(metadata$SRL_giaroots_log)
ggqqplot(log(metadata$SRL_log))


cor.test(metadata$SRL_log, log(metadata$SRL_giaroots_log), method = "pearson")

#data:  metadata$SRL_log and log(metadata$SRL_giaroots_log)
#t = 11.621, df = 142, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6035419 0.7734336
#sample estimates:
#      cor 0.6981891 


plot(metadata$SRL_giaroots,metadata$SRL) 

################################
####Root biomass and length ####

cor.test(sqrt(metadata$Network_length), sqrt(metadata$DryRootWt_total_g), method = "pearson")

#data:  sqrt(metadata$Network_length) and sqrt(metadata$DryRootWt_total_g)
#t = 13.653, df = 142, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6727250 0.8163652
#sample estimates:
#  cor 0.753396 

################################
# Average root diameter and root weight

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$DryRootWt_total_g), method = "pearson")

#data:  metadata$Avg_root_width_diam and sqrt(metadata$DryRootWt_total_g)
#t = 6.7582, df = 142, p-value = 3.342e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3586880 0.6078453
#sample estimates:
#cor = 0.4933195 

######################################
# Average root diameter and root length

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$Network_length), method = "pearson")

#data:  metadata$Avg_root_width_diam and sqrt(metadata$DryRootWt_total_g)
#t = 6.7582, df = 142, p-value = 3.342e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3586880 0.6078453
#sample estimates:
#cor = 0.4933195 


##################################
### Vol-weighted SRL & diameter ##

ggqqplot(log(metadata$SRL_giaroots))
ggqqplot(metadata$Avg_root_width_diam)


cor.test(log(metadata$SRL_giaroots), metadata$Avg_root_width_diam, method = "pearson")


#data:  log(metadata$SRL_giaroots) and metadata$Avg_root_width_diam
#t = -40.859, df = 142, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:-0.9710883 -0.9447951
#sample estimates: cor -0.9600061 


#####################################
#### mass-weighted SRL & diameter ###

ggqqplot(log(metadata$SRL))
ggqqplot(metadata$Avg_root_width_diam)


cor.test(log(metadata$SRL), metadata$Avg_root_width_diam, method = "pearson")

#	Pearson's product-moment correlation

#data:  log(metadata$SRL) and metadata$Avg_root_width_diam
#t = -9.929, df = 142, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:-0.7275307 -0.5322962
#sample estimates: cor -0.6401351 

######################################
# Average root diameter and root weight

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$DryRootWt_total_g), method = "pearson")

#data:  metadata$Avg_root_width_diam and sqrt(metadata$DryRootWt_total_g)
#t = 6.7582, df = 142, p-value = 3.342e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3586880 0.6078453
#sample estimates:
#cor = 0.4933195 

######################################
# Average root diameter and root length

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$Network_length), method = "pearson")

#data:  metadata$Avg_root_width_diam and sqrt(metadata$DryRootWt_total_g)
#t = 6.7582, df = 142, p-value = 3.342e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3586880 0.6078453
#sample estimates:
#cor = 0.4933195 

```
## Volume-weighted SRL 
```{r}

# SRL giaroots = volume-weigthed SRL 
# Summarize 
metadata  %>%  group_by(Variety) %>% summarise_at("SRL_giaroots", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))


# Is SRL normally distributed? 
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(log(metadata$SRL_giaroots)) # decent
hist(log(metadata$SRL_giaroots)) # log improves


## ANALYSIS

# Does SRL differ by Variety? 
require(lme4)
require(lmerTest)

SRL_VarietyBl_model <- lmer(log(SRL_giaroots) ~ Variety + (1|Block/Variety), data = metadata)
plot(SRL_VarietyBl_model)

SRL_VarietyBl_model_SMC <- lmer(log(SRL_giaroots) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(SRL_VarietyBl_model_SMC, SRL_VarietyBl_model)
#df      AIC
#SRL_VarietyBl_model_SMC 16 205.1903
#SRL_VarietyBl_model     15 197.8613 # NO 

# check normality 
plot(SRL_VarietyBl_model)

res = resid(SRL_VarietyBl_model)
shapiro.test(res) 
# log(SRL_giaroots) p = 0.06 

# homoegeneity of variances
bartlett.test(log(SRL_giaroots) ~ Variety, data = metadata) 
# log(SRL_giaroots_OutlierRmv) p = 0.04

anova(SRL_VarietyBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
 #       Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
#Variety 6.1844 0.56222    11    36  3.6184 0.001667 **


# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(SRL_VarietyBl_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD
#Variety      lsmean    SE df lower.CL upper.CL .group
# EG1101         5.84 0.129 36     5.45     6.24  a    
# Cave-in-Rock   5.87 0.129 36     5.47     6.26  a    
# EG2101         5.88 0.129 36     5.48     6.27  a    
# Alamo          6.01 0.129 36     5.62     6.41  ab   
# EG1102         6.05 0.129 36     5.65     6.44  ab   
# Blackwell      6.07 0.129 36     5.67     6.47  ab   
# Dacotah        6.09 0.129 36     5.69     6.49  ab   
# Shelter        6.13 0.129 36     5.74     6.53  ab   
# Southlow       6.33 0.129 36     5.93     6.72  ab   
# Trailblazer    6.39 0.129 36     5.99     6.78  ab   
# Kanlow         6.52 0.129 36     6.13     6.92   b   
# NE28           6.54 0.129 36     6.14     6.93   b   

#Degrees-of-freedom method: kenward-roger 
#Results are given on the log (not the response) scale. 
#Confidence level used: 0.95 
#Conf-level adjustment: bonferroni method for 12 estimates 
#P value adjustment: fdr method for 66 tests 
#significance level used: alpha = 0.05

library(emmeans)
emmeans(SRL_VarietyBl_model, list(pairwise ~ Variety), adjust = "fdr")


############### By Ecotype 
SRL_EcotypeBl_model <- lmer(log(SRL_giaroots) ~ Ecotype + (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
isSingular(SRL_EcotypeBl_model) #FALSE


# shapiro test 
res = resid(SRL_EcotypeBl_model)
shapiro.test(res) 
# log p = 0.02

bartlett.test(log(SRL_giaroots) ~ Ecotype, data = metadata) 
# log(SRL_giaroots) p = 0.3


plot(SRL_EcotypeBl_model)
# log residuals are at slight angle -- maybe this is causing issues 


anova(SRL_EcotypeBl_model, type = "III") # not perfectly normal 
#        Sum Sq  Mean Sq NumDF DenDF F value Pr(>F)
#Ecotype 0.044672 0.044672     1    46  0.2875 0.5944


emmeans(SRL_EcotypeBl_model, pairwise~Ecotype)
# $contrasts
 #contrast         estimate    SE df t.ratio p.value
 #Lowland - Upland  -0.0543 0.101 43 -0.536  0.5946 


```
#### Volume-weighted SRL correlated with SMC or NO3? 
```{r}
# SMC
reg_SRL_SMC<-lm(log(SRL_giaroots)~SMC_perc_GravWaterCont_OutlierRmv, data=metadata)

#Look at the residuals to make sure they are normal
hist(stdres(reg_SRL_SMC))
qqPlot(stdres(reg_SRL_SMC))
#It mostly falls within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(reg_SRL_SMC)
#Since the confidence interval isn't within 1 - need to transform, log is better

summary(reg_SRL_SMC) #this will give you the statistics for the regression
#esidual standard error: 0.4037 on 185 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.0312,	Adjusted R-squared:  0.02596 
#F-statistic: 5.958 on 1 and 185 DF,  p-value: 0.01559

#Run the next 2 lines together to plot the results
plot(metadata$SRL_giaroots_OutlierRmv,metadata$SMC_perc_GravWaterCont_OutlierRmv) 
abline(reg_SRL_SMC)



# NO3
reg_SRL_NO3<-lm(NO3_ugN_g_drysoil_K2SO4~SRL_giaroots_OutlierRmv, data=metadata)
reg_SRL_NO32<-lm(NO3_ugperG.T0~SRL_giaroots_OutlierRmv, data=metadata)


#Look at the residuals to make sure they are normal
hist(stdres(reg_SRL_NO32))
qqPlot(stdres(reg_SRL_NO3))
#It mostly falls within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(reg_SRL_NO3)

summary(reg_SRL_NO3) #this will give you the statistics for the regression
# NO3_K2SO4
#Residual standard error: 0.3282 on 139 degrees of freedom
 # (3 observations deleted due to missingness)
#Multiple R-squared:  0.06524,	Adjusted R-squared:  0.05851 
#F-statistic: 9.701 on 1 and 139 DF,  p-value: 0.002238

summary(reg_SRL_NO32) #this will give you the statistics for the regression
# Residual standard error: 0.7473 on 139 degrees of freedom
 # (3 observations deleted due to missingness)
#Multiple R-squared:  0.0009668,	Adjusted R-squared:  -0.006221 
#F-statistic: 0.1345 on 1 and 139 DF,  p-value: 0.7144

# PLOT
plot(metadata$NO3_ugN_g_drysoil_K2SO4,metadata$SRL_giaroots_OutlierRmv) 
abline(reg_SRL_NO3)




```
## Mass-weighted SRL   
```{r}
# summarize
metadata %>% group_by(Variety)  %>% summarise_at("SRL", list(~mean(.),~var(.), ~min(.),~max(.)))


# Is SRL normally distributed? 
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$SRL) # outliers
ggqqplot(log(metadata$SRL)) # slightly better with log transformation
hist(metadata$SRL)
hist(log(metadata$SRL)) #better with log 


## ANALYSIS
## inspect data
metadata %>% group_by(Variety)  %>% summarise_at("SRL", list(~mean(.),~var(.)))
                                               

# Does SRL differ by Variety? 
require(lme4)
require(lmerTest)

SRL_VarietyBl_model <- lmer(log(SRL) ~ Variety + (1|Block/Variety), data = metadata)
#boundary (singular) fit: see ?isSingular

isSingular(SRL_VarietyBl_model) # TRUE 
# Look at variance - is block 0? 
library(lme4)
VarCorr(SRL_VarietyBl_model)
#Groups        Name        Std.Dev.
# Variety:Block (Intercept) 0.1062  
# Block         (Intercept) 0.0000  
# Residual                  0.6177



plot(SRL_VarietyBl_model)

# shapiro test
res = resid(SRL_VarietyBl_model)
shapiro.test(res) 
#log(SRL), p = 0.2


# homoegeneity of variances
bartlett.test(SRL_log ~ Variety, data = metadata) 
#SRL_log, p = 0.1

# Should SMC be included to improve model fit?

SRL_VarietyBl_model <- lmer(log(SRL) ~ Variety + (1|Block/Variety), data = metadata)

SRL_VarietyBl_model_SMC <- lmer(log(SRL) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(SRL_VarietyBl_model_SMC, SRL_VarietyBl_model)
# df      AIC
#SRL_VarietyBl_model_SMC 16 313.0181 # NO 
#SRL_VarietyBl_model     15 310.2930


anova(SRL_VarietyBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
 #       Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Variety 6.7905 0.61732    11 36.033  1.6179 0.1353
  

############### By Ecotype 
SRL_EcotypeBl_model <- lmer(log(SRL) ~ Ecotype + (1|Block/Plot), data = metadata)

plot(SRL_EcotypeBl_model)
res = resid(SRL_EcotypeBl_model)
shapiro.test(res)
# p = 0.4


anova(SRL_EcotypeBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Ecotype 0.92581 0.92581     1 46.055  2.4279  0.126


```

## Total Root Biomass 
```{r}
# Values back-calculated from dry:wet ratios of roots subset for endophyte DNA extractions 

# Is DryRootWt_total_g normally distributed? 
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$DryRootWt_total_g) # outliers
ggqqplot(log(metadata$DryRootWt_total_g)) # slightly better with log transformation
hist(metadata$DryRootWt_total_g)
hist(log(metadata$DryRootWt_total_g)) #better with log &sqrt


## ANALYSIS

require(lme4)
require(lmerTest)

BGdry_VarietyBl_model <- lmer(sqrt(DryRootWt_total_g) ~ Variety + (1|Block/Variety), data = metadata)

plot(BGdry_VarietyBl_model)

res = resid(BGdry_VarietyBl_model)
shapiro.test(res) 
#log(DryRootWt_total_g), p <0.001
#sqrt - p = 0.69

# homoegeneity of variances
bartlett.test(sqrt(DryRootWt_total_g) ~ Variety, data = metadata) 
#DryRootWt_total_g, p <0.001 
# sqrt p <0.01

# Should SMC be included to improve model fit?

BGdry_VarietyBl_model <- lmer(sqrt(DryRootWt_total_g) ~ Variety + (1|Block/Variety), data = metadata)

BGdry_VarietyBl_model_SMC <- lmer(sqrt(DryRootWt_total_g) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(BGdry_VarietyBl_model, BGdry_VarietyBl_model_SMC)
# df      AIC
#BGdry_VarietyBl_model     15 27.58935
#BGdry_VarietyBl_model_SMC 16 37.53703 # NO 


anova(BGdry_VarietyBl_model, type = "III")
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF DenDF F value Pr(>F)
#Variety 0.79601 0.072365    11    36  1.6071 0.1385

################ ecotype

BGdry_EcoBl_model <- lmer(sqrt(DryRootWt_total_g) ~ Ecotype + (1|Block/Plot), data = metadata)

plot(BGdry_EcoBl_model)

res = resid(BGdry_EcoBl_model)
shapiro.test(res) 
# sqrt p =0.012

# homoegeneity of variances
bartlett.test(sqrt(DryRootWt_total_g) ~ Ecotype, data = metadata) 
#sqrt p = 0.8


anova(BGdry_EcoBl_model, type = "III")
 # Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
#Ecotype 0.11115 0.11115     1    46  2.4685  0.123


```

## Average root diameter 
```{r}

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$Avg_root_width_diam) # great
hist(metadata$Avg_root_width_diam)

## ANALYSIS

require(lme4)
require(lmerTest)

rootdiam_VarietyBl_model <- lm(Avg_root_width_diam ~ Variety + (1|Block/Variety), data = metadata)
#boundary (singular) fit: see ?isSingular

plot(rootdiam_VarietyBl_model) 

res = resid(rootdiam_VarietyBl_model)
shapiro.test(res) 
# p = 0.9


# homoegeneity of variances
bartlett.test(Avg_root_width_diam ~ Variety, data = metadata) 
# p <0.07 

# Should SMC be included to improve model fit?

rootdiam_VarietyBl_model <- lmer(Avg_root_width_diam ~ Variety + (1|Block/Variety), data = metadata)

rootdiam_VarietyBl_model_SMC <- lmer(Avg_root_width_diam ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootdiam_VarietyBl_model, rootdiam_VarietyBl_model_SMC)
#   df       AIC
#rootdiam_VarietyBl_model     15 -923.9182
#rootdiam_VarietyBl_model_SMC 16 -890.7907 #NO

anova(rootdiam_VarietyBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq    Mean Sq NumDF DenDF F value    Pr(>F)    
#Variety 0.0016101 0.00014637    11    36  4.4326 0.0003226 ***

# tukeys comparison will give letter differentiation
library(multcomp)
library(lsmeans)
marginal = lsmeans(rootdiam_VarietyBl_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

# Variety      lsmean      SE   df lower.CL upper.CL .group
# Kanlow       0.0333 0.00144 21.9   0.0287   0.0379  a    
# NE28         0.0337 0.00179 52.0   0.0283   0.0390  ab   
# Trailblazer  0.0364 0.00179 52.0   0.0310   0.0418  abc  
# Southlow     0.0367 0.00144 21.9   0.0321   0.0413  abc  
# Shelter      0.0393 0.00179 52.0   0.0340   0.0447   bcd 
# EG1102       0.0398 0.00179 52.0   0.0344   0.0452   bcd 
# Dacotah      0.0400 0.00179 52.0   0.0346   0.0453   bcd 
# Blackwell    0.0403 0.00179 52.0   0.0349   0.0457    cd 
# Alamo        0.0405 0.00144 21.9   0.0359   0.0451    cd 
#EG1101       0.0435 0.00179 52.0   0.0382   0.0489     d 
# Cave-in-Rock 0.0440 0.00144 21.9   0.0394   0.0486     d 
# EG2101       0.0441 0.00179 52.0   0.0387   0.0494     d 

#Degrees-of-freedom method: kenward-roger 
#Confidence level used: 0.95 
#Conf-level adjustment: bonferroni method for 12 estimates 
#P value adjustment: fdr method for 66 tests 
#significance level used: alpha = 0.05  


### BY Ecotype? 
rootdiam_EcotypeBl_model <- lmer(Avg_root_width_diam ~ Ecotype + (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
res = resid(rootdiam_EcotypeBl_model)
shapiro.test(res) 
# p = 0.95

plot(rootdiam_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(Avg_root_width_diam ~ Ecotype, data = metadata) 
# p =0.5

anova(rootdiam_EcotypeBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#Ecotype 2.1529e-09 2.1529e-09     1    46   1e-04 0.99367

#######################
# Does root diameter correlate with soil moisture content
ggqqplot(metadata$Avg_root_width_diam)
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)


cor.test(metadata$SMC_perc_GravWaterCont_OutlierRmv, metadata$Avg_root_width_diam, method = "pearson")

#data:  metadata$SMC_perc_GravWaterCont_OutlierRmv and metadata$Avg_root_width_diam
#t = -2.2161, df = 140, p-value = 0.0283
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.33855297 -0.01997065
#sample estimates:
#       cor 
#-0.1840923 



```

## Root Length 
```{r}

# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(sqrt(metadata$Network_length)) # great
hist(sqrt(metadata$Network_length))

## ANALYSIS

require(lme4)
require(lmerTest)

rootlength_VarietyBl_model <- lmer(sqrt(Network_length) ~ Variety + (1|Block/Variety), data = metadata)
res = resid(rootlength_VarietyBl_model)
shapiro.test(res) 
# p = <0.001
# log(length), p = 0.04 
# sqrt(length),  p = 0.6

plot(rootlength_VarietyBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_length) ~ Variety, data = metadata) 
#log -  p = 0.02 
# sqrt - p <0.001

#####
# Should SMC be included to improve model fit?

rootlength_VarietyBl_model <- lmer(sqrt(Network_length) ~ Variety + (1|Block/Variety), data = metadata)

rootlength_VarietyBl_model_SMC <- lmer(sqrt(Network_length) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootlength_VarietyBl_model, rootlength_VarietyBl_model_SMC)
#         df     AIC
#rootlength_VarietyBl_model     15 836.643
#rootlength_VarietyBl_model_SMC 16 827.529 # YES 

anova(rootlength_VarietyBl_model_SMC, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Variety                       277.501  25.227    11 34.518  1.2158 0.3137
#SMC_perc_GravWaterCont         15.113  15.113     1 47.173  0.7283 0.3977





### BY Ecotype? 
length_EcotypeBl_model <- lmer(sqrt(Network_length) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
res = resid(length_EcotypeBl_model)
shapiro.test(res) 
# p = 0.06

plot(length_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_length) ~ Ecotype, data = metadata) 
# p =0.7



anova(length_EcotypeBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#                            Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#Ecotype                     66.892  66.892     1 46.891  3.2315 0.07867 .
#SMC_perc_GravWaterCont 20.616  20.616     1 61.639  0.9959 0.32220 


######################################
### CORRELATE ROOT LENGTH WITH SMC 

reg_rootlength_SMC<-lm(sqrt(Network_length)~SMC_perc_GravWaterCont_OutlierRmv, data=metadata)

#Look at the residuals to make sure they are normal
hist(stdres(reg_rootlength_SMC))
qqPlot(stdres(reg_rootlength_SMC))

summary(reg_rootlength_SMC) #this will give you the statistics for the regression
#Residual standard error: 4.645 on 140 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.001674,	Adjusted R-squared:  -0.005457 
#F-statistic: 0.2348 on 1 and 140 DF,  p-value: 0.6288

# PLOT 
plot(metadata$SRL_giaroots,metadata$SRL) 
abline(reg_SRL_SRLgia)


```

## Root network volume 
```{r}
# test for normal residuals and equal variances (ANOVA assumptions)
# 1) normal residuals
ggqqplot(metadata$Network_volume) # great
hist(metadata$Network_volume) # right skew


## ANALYSIS

require(lme4)
require(lmerTest)

vol_VarietyBl_model <- lmer(sqrt(Network_volume) ~ Variety + (1|Block/Variety), data = metadata)
res = resid(vol_VarietyBl_model)
shapiro.test(res) 
# p = <00.001
# sqrt p = 0.6

# homoegeneity of variances
bartlett.test(sqrt(Network_volume) ~ Variety, data = metadata) 
# sqrt - p= 0.7

#####
# Should SMC be included to improve model fit?

rootvol_VarietyBl_model <- lmer(sqrt(Network_volume) ~ Variety + (1|Block/Variety), data = metadata)

rootvol_VarietyBl_model_SMC <- lmer(sqrt(Network_volume) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootvol_VarietyBl_model, rootvol_VarietyBl_model_SMC)
# df      AIC
#rootvol_VarietyBl_model     15 126.0228 
#rootvol_VarietyBl_model_SMC 16 133.7509 #NO 


anova(vol_VarietyBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
 #       Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#Variety 2.0636  0.1876    11 36.001  1.9964 0.05854 .


### BY Ecotype? 
Nvol_EcotypeBl_model <- lmer(sqrt(Network_volume) ~ Ecotype + (1|Block/Plot), data = metadata)
res = resid(Nvol_EcotypeBl_model)
shapiro.test(res) 
# p = 0.5

plot(Nvol_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_volume) ~ Ecotype, data = metadata) 
# p =0.95



anova(Nvol_EcotypeBl_model, type = "III")
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
#Ecotype 0.14023 0.14023     1    46  1.4923 0.2281


```
