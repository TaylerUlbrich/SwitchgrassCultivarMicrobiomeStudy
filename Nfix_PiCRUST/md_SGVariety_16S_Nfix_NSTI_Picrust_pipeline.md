# Tayler Ulbrich
# Picrust Metagenome Prediction Pipeline for SG variety 16S data
#2019.06.03
# OTU table used for this analysis was run through usearch make_contigs and closed_ref with greengenes_97_otu database
#OTU table filtered for singletons, non-bacterial samples, and rarified at 4117 sequences, n =138 samples"

**Resources**
*https://github.com/LangilleLab/microbiome_helper/wiki/PICRUSt-tutorial*
*https://picrust.github.io/picrust/tutorials/metagenome_prediction.html*

# Soil Samples 
## 1) Convert txt file to .biom
		biom convert -i 2019.06.03_SGvar_GGbactrhizo_nosing_rfy4117_OTUtable.txt -o 2019.06.03_SGvar_GGbactrhizo_nosing_rfy4117_otu.biom --table-type="OTU table" --to-json
		## OR: look at the head of the biom file
		biom head -i 2019.06.03_SGvar_GGbactrhizo_nosing_rfy4117_otu.biom

## 2) Correct the OTU table based on predicted 16s copy number for each organism in the otu table
		python /home/ruben/picrust/scripts/normalize_by_copy_number.py -i 2019.06.03_SGvar_GGbactrhizo_nosing_rfy4117_otu.biom -o otu_corrected.biom 
		# convert this into a text file 
		biom convert -i otu_corrected.biom -o otu_corrected.txt --to-tsv --header-key taxonomy

## 3) Make function predictions of KEGG Ortholog (KOs) using the otus_corrected file
		# NOTE: default is KO pathway, but can also look at COGs and Rfams."
		# use -a to get the NSTI value predictions for quality checking 
		python /home/ruben/picrust/scripts/predict_metagenomes.py -i otu_corrected.biom -o ko_predictions.biom -a nsti_per_sample.tab
		# convert this to text (note the different header option)
		biom convert -i ko_predictions.biom -o ko_predictions.txt --to-tsv --header-key KEGG_Description


## 4) Map the KO pathway to KEGG pathways to get a better prediction of functions
		python /home/ruben/picrust/scripts/categorize_by_function.py -i ko_predictions.biom -c KEGG_Pathways -l 3 -o KEGGpathway_predictions.biom

		# convert this to text (note the different header option)
		biom convert -i KEGGpathway_predictions.biom -o KEGGpathway_predictions.txt --to-tsv --header-key KEGG_Pathways

## 5) Use metagenome_contributions.py to connect the OTUS that are contributing to each KO
      # Choose by the KO that you are specifically interested in by looking at ko_map_descsriptions file or  http://www.genome.jp/kegg-bin/get_htext?query=00910&htext=br08901.keg&option=-a
      # This file relates how much a single OTU contributes to a KO pathway within a single sample.
        # The 5th column relates actual rel.abundance contributed by this otu and the percentage contribution
     # N-fixation pathways: K02588, K02586, K02591, K00531  http://www.kegg.jp/kegg-bin/show_module?M00175"
      python /home/ruben/picrust/scripts/metagenome_contributions.py -i otu_corrected.biom -l K02588,K02586,K02591,K00531 -o NitrogenFixation_contributions.txt



# Root and soil samples

## 1) Convert txt file to .biom

	cd C:/Users/Tayler Ulbrich/Documents/01_Grad School/02_Projects/01_SwitchgrassVariety/03_Data/Sequencing/16S/Picrust/Filtered_greengenes_PicrustInput/RhizoEndo

     ## Convert txt to biom
	biom convert -i 2019.09.23_SGvar_RhizoEndo_GGbact_nosing_rfy4573_OTUtable.txt -o 2019.09.23_SGvar_RhizoEndo_GGbact_nosing_rfy4573_OTUtable.biom --table-type="OTU table" --to-json
	## OR: look at the head of the biom file
	biom head -i 2019.09.23_SGvar_RhizoEndo_GGbact_nosing_rfy4573_OTUtable.biom

# 2) Correct the OTU table based on predicted 16s copy number for each organism in the otu table
      python /home/tayler/picrust/scripts/normalize_by_copy_number.py -i 2019.09.23_SGvar_RhizoEndo_GGbact_nosing_rfy4573_OTUtable.biom -o otu_corrected.biom 
        biom convert -i otu_corrected.biom -o otu_corrected.txt --to-tsv 


# 3) Make function predictions of KEGG Ortholog (KOs) using the otus_corrected file
	# NOTE: default is KO pathway, but can also look at COGs and Rfams."
     # use -a to get the NSTI value predictions for quality checking 
		python /home/tayler/picrust/scripts/predict_metagenomes.py -i otu_corrected.biom -o ko_predictions.biom -a nsti_per_sample.tab
		# convert this to text (note the different header option)
     biom convert -i ko_predictions.biom -o ko_predictions.txt --to-tsv --header-key KEGG_Description


# 4) Map the KO pathway to KEGG pathways to get a better prediction of functions
		python /home/tayler/picrust/scripts/categorize_by_function.py -i ko_predictions.biom -c KEGG_Pathways -l 3 -o KEGGpathway_predictions.biom

	  # convert this to text (note the different header option)
	   biom convert -i KEGGpathway_predictions.biom -o KEGGpathway_predictions.txt --to-tsv --header-key KEGG_Pathways

# 5) Use metagenome_contributions.py to connect the OTUS that are contributing to each KO
      # Choose by the KO that you are specifically interested in by looking at ko_map_descsriptions file or  http://www.genome.jp/kegg-bin/get_htext?query=00910&htext=br08901.keg&option=-a
      # This file relates how much a single OTU contributes to a KO pathway within a single sample.
        # The 5th column relates actual rel.abundance contributed by this otu and the percentage contribution
     # N-fixation pathways: K02588, K02586, K02591, K00531  http://www.kegg.jp/kegg-bin/show_module?M00175"
      python /home/tayler/picrust/scripts/metagenome_contributions.py -i otu_corrected.biom -l K02588,K02586,K02591,K00531 -o NitrogenFixation_contributions.txt

