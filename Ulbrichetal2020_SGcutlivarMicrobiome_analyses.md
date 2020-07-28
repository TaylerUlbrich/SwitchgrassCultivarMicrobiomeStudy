---
title: "Ulbrich et al. 2020 SG cultivar analyses"
author: "Tayler Ulbrich"
date: "7/1/2020"
---


*Resources Used:*
EDAMAME Tutorial: https://github.com/raleva/Edamame_phyloseq/blob/master/Edamame_phyloseq.Rmd
Phyloseq Tutorial: http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
https://rpubs.com/maddieSC/R_SOP_UCR_Jan_2018
Microbiome course tutorial: https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html#subset-by-taxonomy

*metadata file:*  PhyloseqObjects_forAnalysis/2019.04.26_SGvariety_Metadata_OutliersRmv.csv

*Bacterial data files for analysis:*
- PhyloseqObjects_forAnalysis/SGvar_16S_rhizo.rfy.Rdata - rarified soil-associated samples (12 cultivars)
- PhyloseqObjects_forAnalysis/SGvar_16S_RhizoEndo_rfy.Rdata - combined and rarified root and soil samples (4 cultivars); subset into 4 cultivars = rhizoendo; or rhizoendo_root for the 4 cultivars' root samples; or rhizoendo_soil for the 4 cultivars' soil samples 
- vst.Rdata files are files transformed with DESEQ2 variance stabalizing transformation; used to confirm that results were robust to rarefaction (protest analysis)	

*Fungal data files for analysis:* 
- PhyloseqObjects_forAnalysis/SGvar_ITS_fungi.rfy.Rdata - rarified soil-associated samples (12 cultivars)
---


# Load Metadata	
```{r}

metadata <- read.csv("PhyloseqObjects_forAnalysis/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)
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
load("PhyloseqObjects_forAnalysis/SGvar_16S_Rhizo_rfy.Rdata")
# extract metadata for future analysis 
sampledata_rhizo.rfy <- as(sample_data(rhizo.rfy), "data.frame")


### Combined rhizosphere and endosphere  
#  rarified data (bact.rfy)
load("PhyloseqObjects_forAnalysis/SGvar_16S_RhizoEndo_rfy.Rdata")
# extract metadata for future analysis 
sampledata_bact.rfy <- as(sample_data(bact.rfy), "data.frame")


### Rhizosphere Fungi 
#  rarified data (fungi.rfy)
load("PhyloseqObjects_forAnalysis/SGvar_ITS_Rhizo_rfy.Rdata")

# extract metadata for future analysis 
sampledata_fungi.rfy <- as(sample_data(fungi.rfy), "data.frame")
# reorder cultivars so that they always appear in the same order 
#sampledata_fungi.rfy$Variety <- factor(sampledata_fungi.rfy$Variety, c("Alamo", "Kanlow", "EG1101", "EG1102", "Blackwell", "Cave-in-Rock", "Dacotah", "NE28", "EG2101", "Shelter", "Southlow","Trailblazer"))

###########################################
# DESEQ Variance Stablizing Transformation data 

### Add the lowest minimum value to the matrix so you don't have neg values. 

#Rhizo.vst
load("PhyloseqObjects_forAnalysis/SGvar_16S_Rhizo_vst.Rdata")
sampledata_rhizo.vst <- as(sample_data(rhizo.vst), "data.frame")

# Bact (Rhizoendo, bact.vst)
load("PhyloseqObjects_forAnalysis/SGvar_16S_RhizoEndo_vst.Rdata")
sampledata_bact.vst <- as(sample_data(bact.vst), "data.frame")


# Fung.vst (soil fungi)
load("PhyloseqObjects_forAnalysis/SGvar_ITS_Rhizo_vst.Rdata")
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
#Total                            137   1.04248                  1.00000                                  

#VST 
set.seed(2)
adonis(rhizo.vst.wunifrac~sampledata_rhizo.rfy$Date_Julian + sampledata_rhizo.rfy$Variety, permutations = 999)
#Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)    
#sampledata_rhizo.rfy$Date_Julian   1 0.00008258 8.2584e-05 15.2577 0.09358  0.001 ***
#sampledata_rhizo.rfy$Variety      10 0.00011797 1.1797e-05  2.1795 0.13367  0.001 ***
#Residuals                        126 0.00068199 5.4130e-06         0.77276           
#Total                            137 0.00088254 1.00000      


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
# Create summary data for taxonomic groups (https://github.com/joey711/phyloseq/issues/418)


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
rhizo.rfy_diversity <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Figures/rhizo.rfy_diversity.csv", header = TRUE)

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
 
 # editted metafile in powerpoint 

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


# 2) homoegeneity of variances
bartlett.test(log(Shannon) ~ Variety , data = fungi.rfy_diversity) 

# check if anova is sig. with mixed effects despite normality
fungi_shannon = lmer(Shannon ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)

fungi_shannon_SMC = lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = fungi.rfy_diversity)


AIC(fungi_shannon,fungi_shannon_SMC)


# ANOVA without SMC as covariate
Anova(fungi_shannon_SMC, type = "III")

### Non-parametric stats for Shannon diversity (cannot add Block or SMC as covariate)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
fungi.shannon_Kruskal <- kruskal.test(Shannon ~ Variety , data = fungi.rfy_diversity)
fungi.shannon_Kruskal

### Non-parametric t-test (Wilcoxon Signed-ranked test)
wilcox.test(Shannon ~ Ecotype, data = fungi.rfy_diversity) 

# check if results are similar when block is in model with parametric test
model = lmer(Shannon ~ Ecotype + (1|Block/Plot), data = fungi.rfy_diversity)
emmeans(model,pairwise~Ecotype, adjust = "fdr")

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


# 2) homoegeneity of variances
bartlett.test(Observed ~ Variety , data = fungi.rfy_diversity)

# Check if the model is stronger with SMC as covariate 
Obs_Variety_model <- lmer(Observed ~ Variety +(1|Block/Variety), data = fungi.rfy_diversity)


# model with SMC as covariate
Obs_Variety_model_SMC = lmer(Observed ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = fungi.rfy_diversity)

# check AIC values 
AIC(Obs_Variety_model,Obs_Variety_model_SMC)
	# lower with SMC 

anova(Obs_Variety_model_SMC, type = "III")

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
	# lower with SMC

anova(Observed_EcoBl_model_SMC, type = "III")

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")

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

# mixed effects model despite non-normal data 
Fungi_Even = lmer(log(Evenness_Pielou) ~ Variety + (1|Block/Variety), data = fungi.rfy_diversity)

Anova(Fungi_Even, type = "III")


### Non-parametric stats for Evenness diversity (can't add block or soil moisture content)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
fungi.even_Kruskal <- kruskal.test(Evenness_Pielou ~ Variety , data = fungi.rfy_diversity)
fungi.even_Kruskal


library(FSA)
fungi.even_dunn <- dunnTest(Evenness_Pielou ~ Variety, data = fungi.rfy_diversity, method= "bh") # bh is the same as fdr 
fungi.even_dunn


### Non-parametric t-test (Wilcoxon Signed-ranked test)
wilcox.test(Evenness_Pielou ~ Ecotype, data = fungi.rfy_diversity) 

# see if similar results with parametric test (controlling for block)
 model = lmer(Evenness_Pielou ~ Ecotype + (1|Block/Plot), data = fungi.rfy_diversity)
 emmeans(model,pairwise~Ecotype, adjust = "fdr")

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
# Diversity seems to be driven by Dacotah Cultivar; Remove Dacotah Cultivar and see if trends hold
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

# homoegeneity of variances
bartlett.test(Shannon ~ Variety, data = rhizo.rfy_diversity) 


# Does Shannon differ by Variety? (Shannon is normally distributed)
require(lme4)
require(lmerTest)

# first check if SMC as covariate improves model fit 
Shannon_VarietyBl_model <- lmer(Shannon ~ Variety + (1|Block/Variety), data = rhizo.rfy_diversity)

Shannon_VarietyBl_model_SMC <- lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data =   rhizo.rfy_diversity)

AIC(Shannon_VarietyBl_model,  Shannon_VarietyBl_model_NO3)
	# lower with SMC


anova(Shannon_VarietyBl_model_date, type = "III")

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

# homoegeneity of variances
bartlett.test(Shannon ~ Variety, data = rhizo.rfy_NoDac_diversity) 

# Does Shannon differ by Variety? (Shannon is normally distributed)
require(lme4)
require(lmerTest)
Shannon_VarietyBl_model_noDac <- lmer(Shannon ~ Variety +  (1|Block/Variety), data = rhizo.rfy_NoDac_diversity)

Shannon_VarietyBl_model_noDac_SMC <- lmer(Shannon ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = rhizo.rfy_NoDac_diversity)

AIC(Shannon_VarietyBl_model_noDac, Shannon_VarietyBl_model_noDac_SMC)
	# lower with SMC

anova(Shannon_VarietyBl_model_noDac_SMC, type = "III")

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

########################################
# Does Shannon Diversity differ by Ecotype? Upland > Lowland (but driven by Dacotah)
#########################################
require(lme4)
require(lmerTest)
Shannon_EcoBl_model <- lmer(Shannon ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_diversity)

Shannon_EcoBl_model_SMC <- lmer(Shannon ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Shannon_EcoBl_model, Shannon_EcoBl_model_SMC)
 # lower with SMC

plot(Shannon_EcoBl_model)

anova(Shannon_EcoBl_model, type = "III")

# pairise test
library(emmeans)
emmeans(Shannon_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")


################ Is this driven by Dacotah? 
require(lme4)
require(lmerTest)
Shannon_EcoBl_noDac_model <- lmer(Shannon ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

anova(Shannon_EcoBl_noDac_model, type = "III")

emmeans(Shannon_EcoBl_noDac_model, pairwise~Ecotype, p.adjust = "fdr")

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

## Is model fit improved with SMC as covariate?
Observed_VarietyBl_model <- lmer(Observed ~ Variety + (1|Block/Variety), data = rhizo.rfy_diversity)

Observed_VarietyBl_SMC_model <- lmer(Observed ~ SMC_perc_GravWaterCont_OutlierRmv + Variety  +  (1|Block/Variety), data = rhizo.rfy_diversity)

AIC(Observed_VarietyBl_date, Observed_VarietyBl_SMC_model)
 #lower with SMC


anova(Observed_VarietyBl_SMC_model, type = "III")



emmeans(Observed_VarietyBl_SMC_model, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Observed_VarietyBl_SMC_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

####################################
# Is it driven by dakota? 

## By Variety? 
Obs_Variety_noDak_SMC <- lmer(Observed ~ Variety +SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizo.rfy_NoDac_diversity)


anova(Obs_Variety_noDak_SMC, type = "III")

emmeans(Obs_Variety_noDak_SMC, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Obs_Variety_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD


########################################
# By Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Observed_EcoBl_model <- lmer(Observed ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_diversity)

Observed_EcoBl_model_SMC <- lmer(Observed ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Observed_EcoBl_model_SMC, Observed_EcoBl_model)
#lower with SMC

plot(Observed_EcoBl_model_SMC)

anova(Observed_EcoBl_model_SMC, type = "III")

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")


############### 
# Without dacotah 
###############
Observed_EcoBl_model <- lmer(Observed ~ Ecotype + (1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

plot(Observed_EcoBl_model)

anova(Observed_EcoBl_model, type = "III")

# pairise test
library(emmeans)
emmeans(Observed_EcoBl_model, pairwise ~ Ecotype, adjust = "fdr")

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

# 2) homoegeneity of variances
bartlett.test(Evenness_Pielou ~ Variety, data = rhizo.rfy_diversity) 

## Analysis

## Will SMC as covariate improve model fit? 
Even_Variety_model <- lmer(Evenness_Pielou ~ Variety +(1|Block/Variety), data = rhizo.rfy_diversity)

Even_Variety_model_SMC <- lmer(Evenness_Pielou ~ Variety + SMC_perc_GravWaterCont_OutlierRmv+(1|Block/Variety), data = rhizo.rfy_diversity)

AIC(Even_Variety_model, Even_Variety_model_SMC)
# df       AIC
 # not lower with SMC

anova(Even_Variety_model, type = "III")

emmeans(Even_Variety_model, pairwise ~ Variety, adjust = "fdr")

# OR DO THIS (gives you letter assignments)
marginal = lsmeans(Even_Variety_model,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons

CLD

########### Is it driven by dakota? 

## By Variety? 
Even_Variety_model <- lmer(Evenness_Pielou ~ Variety +(1|Block/Variety), data = rhizo.rfy_NoDac_diversity)
Even_Variety_model

anova(Even_Variety_model, type = "III")


########################################
# By Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Even_Eco_model <- lmer(Evenness_Pielou ~ Ecotype +(1|Block/Plot), data = rhizo.rfy_diversity)

Even_Eco_model_SMC <- lmer(Evenness_Pielou ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Plot), data = rhizo.rfy_diversity)

AIC(Even_Eco_model, Even_Eco_model_SMC)
	# not lower with SMC-

plot(Even_Eco_model)

anova(Even_Eco_model, type = "III")

# pairise test
library(emmeans)
emmeans(Even_Eco_model, pairwise ~ Ecotype, adjust = "fdr")

###################
## without dacotah 
Even_Eco_model <- lmer(Evenness_Pielou ~ Ecotype +(1|Block/Plot), data = rhizo.rfy_NoDac_diversity)

plot(Even_Eco_model)

anova(Even_Eco_model, type = "III")

emmeans(Even_Eco_model, pairwise ~ Ecotype, adjust = "fdr")

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


# 2) homoegeneity of variances
bartlett.test(Shannon^2 ~ Variety , data = rhizoendo_diversity)

# check if SMC as covariate improves model fit 
rhizoendo_shannon = lmer((Shannon)^2 ~ Variety + (1|Block/Variety), data = rhizoendo_diversity)
# singular fit 

rhizoendo_shannon_SMC = lmer((Shannon)^2 ~ Variety +  SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

AIC(rhizoendo_shannon,rhizoendo_shannon_SMC)
 # not lower with SMC

# Anova - Does Shannon differ by Variety? 
require(lme4)
require(lmerTest)
Shannon_Variety_model <- lmer(Shannon^2 ~ Variety*SampleSite + (1|Block/Variety), data = rhizoendo_diversity)
Shannon_Variety_model
anova(Shannon_Variety_model, type = "III")

emmeans(Shannon_Variety_model, pairwise ~ SampleSite, adjust = "fdr")

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

# 2) homoegeneity of variances
bartlett.test(Observed ~ Variety , data = rhizoendo_diversity) 


# Anova - Does Richness differ by Variety?  AND does SMC include model fit?
rhizoendo_obs = lmer(Observed ~ Variety*SampleSite + (1|Block/Variety), data = rhizoendo_diversity)

rhizoendo_obs_SMC = lmer(Observed ~ Variety*SampleSite +  SMC_perc_GravWaterCont_OutlierRmv +(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

AIC(rhizoendo_obs,rhizoendo_obs_SMC)
 # Lower with SMC 
 
require(lme4)
require(lmerTest)

anova(rhizoendo_obs_SMC, type = "III")

emmeans(rhizoendo_obs_SMC, pairwise ~ SampleSite, adjust = "fdr")

####################################
# does Observed difer by ecotype? 
################################
Obs_EcoSampleSite_model <- lmer(Observed ~ Ecotype*SampleSite + (1|Block/Plot), data = rhizoendo_diversity)
anova(Obs_EcoSampleSite_model)


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

# 2) homoegeneity of variances
bartlett.test(Evenness_Pielou^2 ~ Variety , data = rhizoendo_diversity)
# Evenness_Pielou p =0.6256 - equal variance
# sq Evenness_Pielou p = 0.71407



# Test with proper mixed-effects model with block despite non-normal data 
rhizoendo_even = lmer(Evenness_Pielou^2 ~ Variety*SampleSite+(1|Block/Variety), data = rhizoendo_diversity) # singular fit 

Anova(rhizoendo_even, Type = "III")


### Non-parametric stats for Shannon diversity (cannot add block)
# Note this compares medians/not means so can be affected by unequal distributions 
#http://influentialpoints.com/Training/Kruskal-Wallis_ANOVA_use_and_misuse.htm
rhizoendo.even <- kruskal.test(Evenness_Pielou~ SampleSite, rhizoendo_diversity)
rhizoendo.even

wilcox.test(Evenness_Pielou ~ SampleSite, data = rhizoendo_diversity) 
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
## Final pairwise p-adjust (output from Primer pairwise permanova)
```{r}
### BACTERIA 

# read in list of pvalues from primer pairwise analysis 
Variety_pairwise_Bact_wuni <-read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/PrimerAnalysis/Bacteria/2019.12.17_SMCoutliersRmvd_primer_pairwise_BACT.csv", header = TRUE)

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

Variety_pairwise_Fungi_Bray<-read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/PrimerAnalysis/Fungi/2019.12.17_SMCoutliersRmvd_primer_pairwise_FUNGI.csv", header = TRUE)


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

summary(dbRDA)

plot(dbRDA)

anova(dbRDA) # overall significance 
           
           
anova(dbRDA, by = "axis", perm.max = 500)

# significance & % contribution by each predictor variable
set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 
x

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
plot(dbRDA)

anova(dbRDA) # overal significance 

# significance & % contribution by each predictor variable
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
plot(dbRDA)

anova(dbRDA) # overall significance 
 # not significant 


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

plot(dbRDA)

anova(dbRDA) # overall significance 

# significance & % contribution by each predictor variable
set.seed(2)
x <- anova(dbRDA, by="terms", permu=999) # signifcance by terms 

# Calculate R2 
x$R2 <- (x$SumOfSqs/sum(x$SumOfSqs)) *100

```
---
# Identify variable phyla with MVabund (Figure 4)
  *What bacterial phyla (or classes of proteobacteria) differ among cultivars?*

## Soil Bacteria 
```{r}
#http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)
library(dplyr)
library(tibble)

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
rhizo.rfy_phyla <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/MVabund_Analysis/Bacteria/2019.06.24_MVabund_rhizo.rfy_PhylaProteoClass_OTUtable.csv", header = TRUE)
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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Actino) 


AIC(Actino_VarSMC, Actino_Var)
	# lower with SMC 
	
Anova(Actino_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Actino_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Actino$Abundance)

Actino_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Actino )
Actino_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Actino )

res = Actino_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Actino) 

AIC(Actino_EcoSMC, Actino_Eco)
	# lower with SMC

Anova(Actino_EcoSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Actino_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Acido) 


AIC(Acido_VarSMC, Acido_Var)
	 # slightly lower with SMC 

Anova(Acido_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Acido_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Acido$Abundance)

Acido_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Acido )
Acido_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Acido )

res = Acido_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Acido) 

AIC(Acido_EcoSMC, Acido_Eco)
	# lower with SMC


Anova(Acido_EcoSMC, type = "III")

# get p-adjusted significance letters
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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Bacto) 


AIC(Bacto_VarSMC, Bacto_Var)
	# lower with SMC
	
Anova(Bacto_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Bacto_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Bacto$Abundance)

Bacto_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Bacto )
Bacto_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Bacto )

res = Bacto_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Bacto) 

AIC(Bacto_EcoSMC, Bacto_Eco)
	# lower with SMC

Anova(Bacto_EcoSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Bacto_EcoSMC,
                   ~ Ecotype)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


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

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Beta) 

AIC(Beta_VarSMC, Beta_Var)
	# lower with SMC 
	
Anova(Beta_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Beta_VarSMC,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Beta$Abundance)

Beta_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Beta )
Beta_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Beta )

res = Beta_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test((Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Beta) 

AIC(Beta_EcoSMC, Beta_Eco)
	 #lower with SMC

Anova(Beta_EcoSMC, type = "III")


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

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Delta) 

AIC(Delta_VarSMC, Delta_Var)
	# lower with SMC 	

Anova(Delta_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Delta_VarSMC,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Delta$Abundance)

Delta_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Delta )
Delta_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Delta )

res = Delta_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Delta) 

AIC(Delta_EcoSMC, Delta_Eco)
	# lower with SMC	

Anova(Delta_EcoSMC, type = "III")

# get p-adjusted significance letters
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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Firmi) 

AIC(Firmi_VarSMC, Firmi_Var)
	# lower with SMC
	
Anova(Firmi_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Firmi_VarSMC,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Firmi$Abundance)

Firmi_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Firmi )
Firmi_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Firmi )

res = Firmi_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Firmi) 

AIC(Firmi_EcoSMC, Firmi_Eco)
	# lower with SMC


Anova(Firmi_EcoSMC, type = "III")

 # get p-adjusted significance letters
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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Gemma) 


AIC(Gemma_VarSMC, Gemma_Var)
	# lower with SMC 
	
Anova(Gemma_VarSMC, type = "III")


 # get p-adjusted significance letters
marginal = lsmeans(Gemma_VarSMC,
                   ~ Variety)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Gemma$Abundance)

Gemma_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Gemma )
Gemma_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Gemma )

res = Gemma_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Gemma) 

AIC(Gemma_EcoSMC, Gemma_Eco)
	 # lower with SMC 
	 
Anova(Gemma_EcoSMC, type = "III")


 # get p-adjusted significance letters
marginal = lsmeans(Gemma_EcoSMC,
                   ~ Ecotype)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

  
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

# homoegeneity of variances
bartlett.test((Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Plancto) 

AIC(Plancto_VarSMC, Plancto_Var)
	# lower with SMC 

Anova(Plancto_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Plancto_VarSMC,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD
############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Plancto$Abundance)

Plancto_Eco <- lm((Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Plancto )
Plancto_EcoSMC <- lm((Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Plancto )

res = Plancto_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test((Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Plancto) 

AIC(Plancto_EcoSMC, Plancto_Eco)
	# lower with SMC

Anova(Plancto_EcoSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Plancto_EcoSMC,
                   ~ Ecotype)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

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

# homoegeneity of variances
bartlett.test(log(Abundance) ~ Variety, data = rhizo.rfy_phy_sig.df_sums_Verruco) 

AIC(Verruco_VarSMC, Verruco_Var)
	# lower with SMC
	
Anova(Verruco_VarSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Verruco_VarSMC,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD

############################################################
# ECOTYPE 

ggqqplot(rhizo.rfy_phy_sig.df_sums_Verruco$Abundance)

Verruco_Eco <- lm(log(Abundance) ~ Ecotype , data =rhizo.rfy_phy_sig.df_sums_Verruco )
Verruco_EcoSMC <- lm(log(Abundance) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv, data =rhizo.rfy_phy_sig.df_sums_Verruco )

res = Verruco_EcoSMC$residuals
shapiro.test(res) 

# homoegeneity of Ecoiances
bartlett.test(log(Abundance) ~ Ecotype, data = rhizo.rfy_phy_sig.df_sums_Verruco) 

AIC(Verruco_EcoSMC, Verruco_Eco)
# lower with SMC

Anova(Verruco_EcoSMC, type = "III")

# get p-adjusted significance letters
marginal = lsmeans(Verruco_EcoSMC,
                   ~ Ecotype)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="fdr")         ###  Tukey-adjusted comparisons
CLD


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
    
	# export as table
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

# homoegeneity of variances
bartlett.test(log(Abundance + 1) ~ Variety, data = fungi.Rozello.df) 

AIC(fungi.Rozello.df_Var, fungi.Rozello.df_VarSMC)
	# lower with SMC 
	
Anova(fungi.Rozello.df_VarSMC, type = "III")

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

rhizo.rfy_order <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/MVabund_Analysis/Bacteria/2019.05.08_MVabund_OrderOTUtable_XsampleID.csv", header = TRUE)

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

##################### Analysis  with SRL 
#This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.

rhizo.rfy_core_SRL_results <-anova(rhizo.rfy_core_SRL, p.uni = "adjusted", nBoot =9999)

rhizo.rfy_core_SRL_results$table

############## Analysis with Root Diameter	                                     

rhizo.rfy_core_diam_results <-anova(rhizo.rfy_core_diam, p.uni = "adjusted", nBoot =9999)
rhizo.rfy_core_diam_results$table

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

################ Fungal order ~ Root diameter

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

################ Fungal Order ~ Root Length

fungi.rfy_order_meta <- merge(sampledata_fungi.rfy[,c(2,3,7,8,10,18,23)], fungi.rfy_order, by = "X.SampleID")

library(mvabund)
fungi.rfy_order_all.MVA <- mvabund(fungi.rfy_order_meta[,c(8:136)])

# Do any Orders differ with root length?
Order_RootLength <- manyglm(fungi.rfy_order_all.MVA ~ fungi.rfy_order_meta$Network_length, family = "negative_binomial")
plot(Order_RootLength) # residuals look normal

Order_RootLength_results <-anova(Order_RootLength, p.uni = "adjusted", nBoot =999)

Order_RootLength_results$table

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

# read in output from multiplatt function (organized by cultivar in excel) 	
IndSp <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/SharedIndicator_taxa/2019.07.18_rhizo.rfy_nosingles_AllOTUs_perm9999_clean.csv", header = TRUE)

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

# write.csv(IndSp_p_fdr_merge_25, "C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/SharedIndicator_taxa/2019.07.19_rhizo.rfy_nosing_allOTUS_indsp_Var_perm9999_tax_pvalAdj_25percRep.csv")

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
IndSp <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/SharedIndicator_taxa/2019.07.19_fungi.rfy_rmvSing_indsp_Var_perm9999_clean.csv", header = TRUE)

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

IndSp_Var <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/SharedIndicator_taxa/2019.07.19_fungi.rfy_nosingles_indsp_Var_perm9999_tax_pvalAdjust_25percRep.csv", header = TRUE)

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
	
kruskal.test(Nfixers_ra.df_mean_meta$Abundance ~ Nfixers_ra.df_mean_meta$Ecotype)

##########################################################################
# Do the rel.abund. of soil N-fixers differ among sample sites (root or soil for n = 4 cultivars)?
#########################################################################
rhizoendo_order <- tax_glom(rhizoendo, "Order")

rhizoendo_Nfixers <- subset_taxa(rhizoendo_order, Order == "Rhizobiales" | Order == "Burkholderiales")
ntaxa(rhizoendo_order) #230


rhizoendo_Nfixers.df <- psmelt(rhizoendo_Nfixers)

rhizoendo_Nfixers.df_sum <- rhizoendo_Nfixers.df %>% group_by(X.SampleID,SampleSite, Order) %>%  summarize_at("Abundance",sum)
dim(rhizoendo_Nfixers.df_sum) #176 4

metadata2 <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/2019.04.26_SGvariety_Metadata_OutliersRmv.csv", header = T)

rhizoendo_Nfixers_ra.df_mean_meta <- merge(rhizoendo_Nfixers.df_sum, metadata2[,c(1,2,3,7,10,4,19,32,39,40,51,72)], by = "X.SampleID")


kruskal.test(rhizoendo_Nfixers_ra.df_mean_meta$Abundance ~ rhizoendo_Nfixers_ra.df_mean_meta$SampleSite)

# summarize mean by sample site 
rhizoendo_Nfixers_ra.df_mean_meta %>% group_by(SampleSite) %>%    summarise( mean(Abundance)) 


```
## Prop. of putative N-fixers (PICRUST)
```{r}
# We used a closed refrence Green-genes referenced OTU table with PICRUSt (Langille et al. 2013) to predict the relative abundance of N-fixers. 
 # We first normalized all OTUs by their predicted 16S rRNA gene copy number, which provides a pseudo-abundance estimate for each OTU.
 # Then we used ‘metagenome_predictions’ to obtain OTU-specific gene counts for N-fixation using the following KEGG pathway orthologs K02588, K02586, K02591, K00531. 
 # We calculated each samples’ predicted proportion of N-fixation genes by dividing the number of OTUs with at least one predicted N-fixation pathway for each sample by the normalized abundance of OTUS (e.g., the total 16S-gene normalized OTU counts).


# Read in Nfix contributions from 'metagenome_predictions' output 
Rhizo_Nfix <- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Nfix_PiCRUST/Rhizo_NitrogenFixation_contributions.csv", header = TRUE)

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
Rhizo_NormalizedOTUS<- read.csv("C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/Git_SGcultivar/Nfix_PiCRUST/Rhizo_OTU_corrected.csv", header = TRUE)
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

# homoegeneity of variances
bartlett.test(Sum_NfixCount_Normalized ~ Variety, data = Rhizo_NfixOTUcounts_NormSeqSum_meta) 


plot(Nfix_rhizo_model)

require(lme4)
require(lmerTest)

# Anova
anova(Nfix_rhizo_model, type = "III")


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


########################################
# Does Nfix differ by Ecotype? 
#########################################
require(lme4)
require(lmerTest)
Nfix_Rhizo_model_eco <- lmer(Sum_NfixCount_Normalized ~ Ecotype + (1|Block/Plot), data = Rhizo_NfixOTUcounts_NormSeqSum_meta)
plot(Nfix_Rhizo_model_eco)

# shapiro test 
res = resid(Nfix_Rhizo_model_eco)
shapiro.test(res) 

# homoegeneity of variances
bartlett.test(Sum_NfixCount_Normalized ~ Ecotype, data = Rhizo_NfixOTUcounts_NormSeqSum_meta) 

anova(Nfix_Rhizo_model_eco, type = "III")


emmeans(Nfix_Rhizo_model_eco, pairwise~Ecotype, adjust = "fdr")


#### IS predicted # of Nfixers correlated with soil Nfixation and soil N measured by paired study
#Roley et al. 2020 Phytobiomes https://doi.org/10.1094/PBIOMES-11-19-0064-FI 
  
cor.test(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized, Rhizo_NfixOTUcounts_NormSeqSum_meta$sfix_rate_gN_d, method = "pearson")


### Does Nfixer abundance correlate with nitrate availability?

cor.test(Rhizo_NfixOTUcounts_NormSeqSum_meta$Sum_NfixCount_Normalized, Rhizo_NfixOTUcounts_NormSeqSum_meta$NO3_ugN_g_drysoil_K2SO4, method = "pearson")


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

plot(MBC_Variety_SMC_model) 

# normal residuals?
res = resid(MBC_Variety_SMC_model)
shapiro.test(res) 

#  homoegeneity of variances?
bartlett.test(ugC_MicBiomass_g_dry_soil~ Variety, data = metadata)


MBC_Variety_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Variety +  (1|Block/Variety), data = metadata)

# normal residuals?
res = resid(MBC_Variety_model)
shapiro.test(res) 


## Is the model improved with SMC as covariate?

AIC(MBC_Variety_NO3_model, MBC_Variety_model)
 # lower with SMC 

anova(MBC_Variety_SMC_model, type = "III")

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

###########################################################################
# does MBC difer by ecotype? 
###########################################################################
MBC_Ecotype_model <- lmer(ugC_MicBiomass_g_dry_soil ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot),data = metadata)

# check normality of residuals
res = MBC_Ecotype_model$resdiduals   
shapiro.test(res) 

anova(MBC_Ecotype_model, type = "III") 

emmeans(MBC_Ecotype_model, pairwise~Ecotype)

 
```

###  MBC correlation with edaphic conditions and root traits?
```{r}
* use a pearson correlation to determine if microbial biomass carbon correlates with edaphic cconditions and root traits*

# Specific Root Length 
ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(log(metadata$SRL_giaroots_OutlierRmv))


cor.test(metadata$ugC_MicBiomass_g_dry_soil, log(metadata$SRL_giaroots_OutlierRmv), method = "pearson")

######################
# avg. root diam 

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(metadata$Avg_root_width_diam)

cor.test(metadata$ugC_MicBiomass_g_dry_soil, metadata$Avg_root_width_diam, method = "pearson")

##############
# Root Biomass 

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(log(metadata$DryRootWt_total_g))

cor.test(metadata$ugC_MicBiomass_g_dry_soil, log(metadata$DryRootWt_total_g), method = "pearson")

plot( DryRootWt_total_g~ ugC_MicBiomass_g_dry_soil, data = metadata)

################
# Root length  

ggqqplot(metadata$ugC_MicBiomass_g_dry_soil)
ggqqplot(sqrt(metadata$Network_length))

cor.test(metadata$ugC_MicBiomass_g_dry_soil, sqrt(metadata$Network_length), method = "pearson")

plot( Network_length~ ugC_MicBiomass_g_dry_soil, data = metadata)

#############
# MBC correlate with SMC? 

cor.test(metadata$ugC_MicBiomass_g_dry_soil, (metadata$SMC_perc_GravWaterCont_OutlierRmv), method = "pearson")

###########################
# MBC correlate with NO3? 

ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4)
cor.test(metadata$ugC_MicBiomass_g_dry_soil, (metadata$NO3_ugN_g_drysoil_K2SO4), method = "pearson")

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


#  homoegeneity of variances?
bartlett.test(ugN_MicBiomass_g_dry_soil~ Variety, data = metadata)

MBN_Variety_model <- lmer(ugN_MicBiomass_g_dry_soil ~ Variety +  (1|Block/Variety), data = metadata)

# normal residuals?
res = resid(MBN_Variety_model)
shapiro.test(res) 


AIC(MBN_Variety_SMC_model, MBN_Variety_model)
#   df      AIC
	# lower with SMC 

anova(MBN_Variety_SMC_model, type = "III")


# p-adjust with significant letters
library(multcomp)
library(lsmeans)
marginal = lsmeans(MBN_Variety_SMC_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

###########################################################################
# does MBN difer by ecotype? 
###########################################################################
MBN_Ecotype_model <- lmer(ugN_MicBiomass_g_dry_soil ~ Ecotype + SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot),data = metadata)

# check normality of residuals
res = resid(MBN_Ecotype_model)
shapiro.test(res) # p = 0.5

anova(MBN_Ecotype_model, type = "III") 

emmeans(MBN_Ecotype_model, pairwise~Ecotype)

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


# 2) homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Block , data = metadata)


# Does SMC differ by BLock
SMC_Block_model <- lm(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Block, data = metadata)
SMC_Block_model
Anova(SMC_Block_model, type = "III") 

# pairwise test (gives you p.values for every combination)
emmeans(SMC_Block_model, pairwise ~ Block, adjust = "fdr")

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(SMC_perc_GravWaterCont_OutlierRmv))

SMC_blockmean <- metadata2 %>%  group_by(Block) %>% summarise_at("SMC_perc_GravWaterCont_OutlierRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))


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


# 2) homoegeneity of variances
bartlett.test(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block , data = metadata)


# Does NO3 differ by BLock
NO3_Block_model <- lm(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
NO3_Block_model
Anova(NO3_Block_model, type ="III") 


# pairwise test (gives you p.values for every combination)
emmeans(NO3_Block_model, pairwise ~ Block, adjust = "fdr")

plot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv ~ metadata$Block)

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv))

NO3_blockmean <- metadata2 %>%  group_by(Block) %>% summarise_at("NO3_ugN_g_drysoil_K2SO4_SMCOutRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))
NO3_blockmean


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

# 2) homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block , data = metadata)


plot(model)

# Does NH4 differ by BLock
NH4_Block_model <- lm(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Block, data = metadata)
NH4_Block_model
Anova(NH4_Block_model, type = "III") 
  

# pairwise test (gives you p.values for every combination)
emmeans(NH4_Block_model, pairwise ~ Block, adjust = "fdr")

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


# 2) homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Variety , data = metadata)


# Does SMC differ by Variety
SMC_Var_Block_model <- lmer(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Variety + (1|Block/Variety), data = metadata)
SMC_Var_Block_model
anova(SMC_Var_Block_model, type = "III") 

# pairwise test (gives you p.values for every combination)
emmeans(SMC_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# p-adjust with significant letters
library(multcomp)
library(lsmeans)
marginal = lsmeans(SMC_Var_Block_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

#### BY ecotype? 
SMC_Ecotype_model <- lmer(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Ecotype + (1|Block/Plot),data = metadata)

# check normality of residuals
res = resid(SMC_Ecotype_model)


anova(SMC_Ecotype_model, type = "III") 


emmeans(SMC_Ecotype_model, pairwise~Ecotype)


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

# 2) homoegeneity of variances
bartlett.test(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety , data = metadata)


# Does NO3 differ by Variety
NO3_Var_Block_model <- lmer(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
 # singular fit
anova(NO3_Var_Block_model, type = "III") 

# pairwise test (gives you p.values for every combination)
emmeans(NO3_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# p-adjust with significant letters
library(multcomp)
library(lsmeans)
marginal = lsmeans(NO3_Var_Block_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

#### BY ecotype? 
NO3_Ecotype_model <- lmer(log(NO3_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Ecotype + (1|Block/Plot),data = metadata)


# check normality of residuals
res = resid(NO3_Ecotype_model)
shapiro.test(res) 

anova(NO3_Ecotype_model, type = "III") 

emmeans(NO3_Ecotype_model, pairwise~Ecotype)


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


# 2) homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety , data = metadata)

plot(model)

# Does NH4 differ by Variety
NH4_Var_Block_model <- lmer(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv+1) ~ Variety + (1|Block/Variety), data = metadata)
 # singular fit
anova(NH4_Var_Block_model, type = "III") 

# pairwise test (gives you p.values for every combination)
emmeans(NH4_Var_Block_model, pairwise ~ Variety, adjust = "fdr")

# p-adjust with significant letters
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

emmeans(NH4_Ecotype_model, pairwise~Ecotype)

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

# homoegeneity of variances
bartlett.test(log(SMC_perc_GravWaterCont_OutlierRmv) ~ Date, data = metadata) 


Anova(model, type = "III")

emmeans(model, pairwise~Date, adjust = "fdr")


plot(SMC_perc_GravWaterCont_OutlierRmv ~Date_Julian, metadata)

# Get means by date
# remove any data with missing metadata 
metadata2 <- metadata %>%
  filter(!is.na(SMC_perc_GravWaterCont_OutlierRmv))

SMC_blockmean <- metadata2 %>%  group_by(VarRep) %>% summarise_at("SMC_perc_GravWaterCont_OutlierRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))


metadata2  %>%  group_by(VarRep) %>% summarise_at("pct_water_soil", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))


# as correlation with julian date 
## SMC 
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)

cor.test(metadata$SMC_perc_GravWaterCont_OutlierRmv, metadata$Date_Julian, method = "pearson")
cor.test(metadata$pct_water_soil, metadata$Date_Julian, method = "pearson")

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

plot(model)

# homoegeneity of variances
bartlett.test((NO3_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data = metadata) 


Anova(model, type = "III")


emmeans(model, pairwise~Date, adjust = "fdr")

# Summarize means by date 
metadata2  %>%  group_by(Date) %>% summarise_at("NO3_ugN_g_drysoil_K2SO4_SMCOutRmv", list(~mean(.),~var(.), ~sd(.), ~min(.),~max(.)))


## NO3
ggqqplot(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv)

cor.test(metadata$NO3_ugN_g_drysoil_K2SO4_SMCOutRmv, metadata$Date_Julian, method = "pearson")


##########################
### Soil ammonium (NH4)###

# normal residuals?
ggqqplot((metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv))
hist(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv))

# shapiro test - null hypothesis: residuals are normally distributed 
model = lm(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data =   metadata)
res = model$residuals
shapiro.test(res) 


# homoegeneity of variances
bartlett.test(log(NH4_ugN_g_drysoil_K2SO4_SMCOutRmv) ~ Date, data = metadata) 


Anova(model, type = "III")


emmeans(model, pairwise~Date, adjust = "fdr")

## NH4
ggqqplot(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv)

cor.test(log(metadata$NH4_ugN_g_drysoil_K2SO4_SMCOutRmv), metadata$Date_Julian, method = "pearson")


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


plot(metadata$SRL_giaroots,metadata$SRL) 

################################
####Root biomass and length ####

cor.test(sqrt(metadata$Network_length), sqrt(metadata$DryRootWt_total_g), method = "pearson")

################################
# Average root diameter and root weight

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$DryRootWt_total_g), method = "pearson")

######################################
# Average root diameter and root length

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$Network_length), method = "pearson")


##################################
### Vol-weighted SRL & diameter ##

ggqqplot(log(metadata$SRL_giaroots))
ggqqplot(metadata$Avg_root_width_diam)


cor.test(log(metadata$SRL_giaroots), metadata$Avg_root_width_diam, method = "pearson")

#####################################
#### mass-weighted SRL & diameter ###

ggqqplot(log(metadata$SRL))
ggqqplot(metadata$Avg_root_width_diam)


cor.test(log(metadata$SRL), metadata$Avg_root_width_diam, method = "pearson")

######################################
# Average root diameter and root weight

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$DryRootWt_total_g), method = "pearson")

######################################
# Average root diameter and root length

cor.test(metadata$Avg_root_width_diam, sqrt(metadata$Network_length), method = "pearson")


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
	 # not lower with SMC

# check normality 
plot(SRL_VarietyBl_model)

res = resid(SRL_VarietyBl_model)
shapiro.test(res) 


# homoegeneity of variances
bartlett.test(log(SRL_giaroots) ~ Variety, data = metadata) 


anova(SRL_VarietyBl_model, type = "III")


# p-adjust with significant letters
library(multcomp)
library(lsmeans)
marginal = lsmeans(SRL_VarietyBl_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD

library(emmeans)
emmeans(SRL_VarietyBl_model, list(pairwise ~ Variety), adjust = "fdr")


############### By Ecotype 
SRL_EcotypeBl_model <- lmer(log(SRL_giaroots) ~ Ecotype + (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
isSingular(SRL_EcotypeBl_model) #FALSE


# shapiro test 
res = resid(SRL_EcotypeBl_model)
shapiro.test(res) 


bartlett.test(log(SRL_giaroots) ~ Ecotype, data = metadata) 



plot(SRL_EcotypeBl_model)


anova(SRL_EcotypeBl_model, type = "III") # not perfectly normal 

emmeans(SRL_EcotypeBl_model, pairwise~Ecotype)


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


plot(SRL_VarietyBl_model)

# shapiro test
res = resid(SRL_VarietyBl_model)
shapiro.test(res) 



# homoegeneity of variances
bartlett.test(SRL_log ~ Variety, data = metadata) 


# Should SMC be included to improve model fit?
SRL_VarietyBl_model <- lmer(log(SRL) ~ Variety + (1|Block/Variety), data = metadata)

SRL_VarietyBl_model_SMC <- lmer(log(SRL) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(SRL_VarietyBl_model_SMC, SRL_VarietyBl_model)
	 # not lower with SMC

anova(SRL_VarietyBl_model, type = "III")
  

############### By Ecotype 
SRL_EcotypeBl_model <- lmer(log(SRL) ~ Ecotype + (1|Block/Plot), data = metadata)

plot(SRL_EcotypeBl_model)
res = resid(SRL_EcotypeBl_model)
shapiro.test(res)



anova(SRL_EcotypeBl_model, type = "III")

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

# homoegeneity of variances
bartlett.test(sqrt(DryRootWt_total_g) ~ Variety, data = metadata) 


# Should SMC be included to improve model fit?
BGdry_VarietyBl_model <- lmer(sqrt(DryRootWt_total_g) ~ Variety + (1|Block/Variety), data = metadata)

BGdry_VarietyBl_model_SMC <- lmer(sqrt(DryRootWt_total_g) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(BGdry_VarietyBl_model, BGdry_VarietyBl_model_SMC)
	# not lower with SMC

anova(BGdry_VarietyBl_model, type = "III")


################ ecotype

BGdry_EcoBl_model <- lmer(sqrt(DryRootWt_total_g) ~ Ecotype + (1|Block/Plot), data = metadata)

plot(BGdry_EcoBl_model)

res = resid(BGdry_EcoBl_model)
shapiro.test(res) 

# homoegeneity of variances
bartlett.test(sqrt(DryRootWt_total_g) ~ Ecotype, data = metadata) 

anova(BGdry_EcoBl_model, type = "III")


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

# homoegeneity of variances
bartlett.test(Avg_root_width_diam ~ Variety, data = metadata) 


# Should SMC be included to improve model fit?

rootdiam_VarietyBl_model <- lmer(Avg_root_width_diam ~ Variety + (1|Block/Variety), data = metadata)

rootdiam_VarietyBl_model_SMC <- lmer(Avg_root_width_diam ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootdiam_VarietyBl_model, rootdiam_VarietyBl_model_SMC)
	# not lower with SMC

anova(rootdiam_VarietyBl_model, type = "III")


# p-adjust with significant letters
library(multcomp)
library(lsmeans)
marginal = lsmeans(rootdiam_VarietyBl_model,
                   ~ Variety)
CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,        ### Use lower-case letters for .group
          adjust="FDR")         ###  Tukey-adjusted comparisons
CLD


### BY Ecotype? 
rootdiam_EcotypeBl_model <- lmer(Avg_root_width_diam ~ Ecotype + (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
res = resid(rootdiam_EcotypeBl_model)
shapiro.test(res) 

plot(rootdiam_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(Avg_root_width_diam ~ Ecotype, data = metadata) 

anova(rootdiam_EcotypeBl_model, type = "III")

#######################
# Does root diameter correlate with soil moisture content
ggqqplot(metadata$Avg_root_width_diam)
ggqqplot(metadata$SMC_perc_GravWaterCont_OutlierRmv)


cor.test(metadata$SMC_perc_GravWaterCont_OutlierRmv, metadata$Avg_root_width_diam, method = "pearson")


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

plot(rootlength_VarietyBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_length) ~ Variety, data = metadata) 


#####
# Should SMC be included to improve model fit?

rootlength_VarietyBl_model <- lmer(sqrt(Network_length) ~ Variety + (1|Block/Variety), data = metadata)

rootlength_VarietyBl_model_SMC <- lmer(sqrt(Network_length) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootlength_VarietyBl_model, rootlength_VarietyBl_model_SMC)
 # lower with SMC

anova(rootlength_VarietyBl_model_SMC, type = "III")


### BY Ecotype? 
length_EcotypeBl_model <- lmer(sqrt(Network_length) ~ Ecotype +SMC_perc_GravWaterCont_OutlierRmv+ (1|Block/Plot), data = metadata)
#boundary (singular) fit: see ?isSingular
res = resid(length_EcotypeBl_model)
shapiro.test(res) 


plot(length_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_length) ~ Ecotype, data = metadata) 


anova(length_EcotypeBl_model, type = "III")

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

# homoegeneity of variances
bartlett.test(sqrt(Network_volume) ~ Variety, data = metadata) 

#####
# Should SMC be included to improve model fit?

rootvol_VarietyBl_model <- lmer(sqrt(Network_volume) ~ Variety + (1|Block/Variety), data = metadata)

rootvol_VarietyBl_model_SMC <- lmer(sqrt(Network_volume) ~ Variety + SMC_perc_GravWaterCont_OutlierRmv + (1|Block/Variety), data = metadata)

# AIC check with SMC
AIC(rootvol_VarietyBl_model, rootvol_VarietyBl_model_SMC)
	# not lower with SMC	

anova(vol_VarietyBl_model, type = "III")

### BY Ecotype? 
Nvol_EcotypeBl_model <- lmer(sqrt(Network_volume) ~ Ecotype + (1|Block/Plot), data = metadata)
res = resid(Nvol_EcotypeBl_model)
shapiro.test(res) 

plot(Nvol_EcotypeBl_model)

# homoegeneity of variances
bartlett.test(sqrt(Network_volume) ~ Ecotype, data = metadata) 

anova(Nvol_EcotypeBl_model, type = "III")

```
