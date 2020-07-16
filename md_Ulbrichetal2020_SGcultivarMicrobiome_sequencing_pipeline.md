---
title: "SGcultivar sequencing pipeline"
author: "Tayler Ulbrich"
date: "7/6/2020"
output: html_document
---

# Bacterial Sequencing Pipeline (primers 515f/806r)
" February 1, 2018 by Lukas Bell-Dereske and Tayler Ulbrich"
## 1) Check read quality & merge reads for all runs 
```{r warning=TRUE}
# Soil samples (n = 144) were in run 3, runs 1 and 2 have root samples (n = 48)
  # Run1: /mnt/research/EvansLab/MMPRNT_2017/20180113_16S_V4_PE (contains SGendo samples)
  # Run2: /mnt/research/EvansLab/MMPRNT_2017/20180119_B_16S_V4_PE (contains SGendo samples)
  # Run3: /mnt/research/EvansLab/Namibia_LTER/20170220_16S_V4_PE (contains all SGvarietyRhizo16s)

######### Repeat these steps for all runs

# Make sure you are in the folder with the extract forward and reverse reads from Run 1

# 1) unzip the fasta.gz files
  gunzip *.fastq.gz

# 2a) look at the raw unmerged seqs https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

  mkdir fastq_info


  nano fasta_info_fq.sh

  !#/bin/bash

  for fq in *.fastq

  do
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_info $fq -output fastq_info/$fq

  done

  chmod +x fasta_info_fq.sh

  #run the for loop created above

  ./fasta_info_fq.sh


  #move to the fastq_info and make a file for all the fastq quality scores (EE)
  cd fastq_info
  grep ""^EE"" * > run2_fastq_EE.txt"
  grep ""^File"" * > run2_fastq_length.txt"

      # Lukas did this for run3 (Namibia_LTER/20170220_1620170220_16S_V4_PE/)


  #look for any forward and reverse reads that look especially bad in terms of quality.
  #this info will also be really helpful for troubleshooting later on (e.g. why some samples have extermely low read numbers)

#Repeat the above for Run2 and Run3 as well.

# 2b) rename the .fastq files
  # We were dumb and added - and _ into our file names. This creates issue with downstream merge step.
  # We need to rename all of the .fastq files to not have these names.
  # make a readme file to explain this to future users.
  nano readme.txt
  # Some of these original fastq file names have underscores (_) which cannot be read with usearch mergepais step.
# Need to rename all of the ""TMC_"" file names by removing the first two underscores:"
# Original names: ""TMC_SGendo_"" will change to ""TMCSGendo"" and ""TMC_neighID_"" will change to ""TMCneighI$"
      rename ""TMC_SGendo_"" ""TMCSGendo"" TMC_SGendo*"
      rename ""TMC_neighID_"" ""TMCneighID"" TMC_neighID*"



# 3) MERGING Sequences and checking Quality (repeat steps for all runs)
 
    # run1 MERGE -  Merging the sequences from Run 1 - Make sure you are in the folder with the extracted forward and reverse reads from Run 1
    #https://www.drive5.com/usearch/manual/merge_options.html
    mkdir mergedfastq_run1

    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/research/EvansLab/MMPRNT_2017/20180113_16S-V4_PE/*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run1/combined_merged_run1.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run1/combined_merged_merge_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run1/combined_merged_merge_alnout.txt

    ## OUTPUT: WARNING: Max OMP threads 1
    #03:42 128Mb   100.0% 80.3% merged
    #Merged length distribution:
    #        16  Min
    #       253  Low quartile
    #       253  Median
    #       253  High quartile
    #       484  Max
    #Totals:
    #   9405243  Pairs (9.4M)
    #   7550191  Merged (7.6M, 80.28%)"
    #   4095510  Alignments with zero diffs (43.54%)
    #   1837774  Too many diffs (> 5) (19.54%)
    #     17278  No alignment found (0.18%)
    #         0  Alignment too short (< 16) (0.00%)
    #     37572  Staggered pairs (0.40%) merged & trimmed
    #    246.69  Mean alignment length
    #    253.08  Mean merged length
    #      0.31  Mean fwd expected errors
    #      0.92  Mean rev expected errors
    #      0.07  Mean merged expected errors


    #let's check sequence quality of the merged seqs - USEARCH has a great command for checking the expected error and lengths of the sequences.
    #https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html

    mkdir /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run1/fastq_info/
    cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run1
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 combined_merged_run1.fastq -output fastq_info/combined_merged_run1_eestats2.txt

    # Output:
     7550191 reads, max len 484, avg 253.1"
    Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
    ------   ----------------   ----------------   ----------------
        50    7524799( 99.7%)    7549278(100.0%)    7549899(100.0%)
       100    7474062( 99.0%)    7545430( 99.9%)    7549550(100.0%)
       150    7434766( 98.5%)    7538997( 99.9%)    7549155(100.0%)
       200    7404935( 98.1%)    7534666( 99.8%)    7548586(100.0%)
       250    7328140( 97.1%)    7497031( 99.3%)    7520211( 99.6%)
       300      10102(  0.1%)      12050(  0.2%)      12905(  0.2%)
       350       2059(  0.0%)       3717(  0.0%)       4881(  0.1%)
       400       1326(  0.0%)       2806(  0.0%)       4252(  0.1%)
       450          4(  0.0%)          9(  0.0%)         22(  0.0%)


    # Rename the merged file to get rid of hyphens in sample names ""MMPRNT-#"" should be ""MMPRNT#"""
         sed -e 's/MMPRNT-/MMPRNT/g' mergedfastq_run1/combined_merged_run1.fastq > mergedfastq_run1/renamed_combined_merged_run1.fastq
         # check sample names to make sure hyphen is gone:
         less mergedfastq_run1/renamed_combined_merged_run1.fastq


     # cut-adapt to get rid of adapters (adapters for 515F and reverse complement of 806R)
       # FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT (reverse compliment: ATTAGAWACCCBDGTAGTCC)
       cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
        module load cutadapt/1.8.1
       cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o mergedfastq_run1/cut_renamed_combined_merged_run1.fastq mergedfastq_run1/renamed_combined_merged_run1.fastq > mergedfastq_run1/cutadpt_results_cut_renamed_combined_merged_run1.txt

       # Summary
       #Total reads processed:               7,550,191"
       #Reads with adapters:                     5,822 (0.1%)"
       #Reads written (passing filters):     7,550,191 (100.0%)"
       #Total basepairs processed: 1,910,814,015 bp"
       #Total written (filtered):  1,910,788,452 bp (100.0%)"


    # Look at the quality stats of this file
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq_run1/cut_renamed_combined_merged_run1.fastq -output mergedfastq_run1/cut_renamed_combined_merged_run1_eestats.txt

    Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
    ------   ----------------   ----------------   ----------------
        50    7524796( 99.7%)    7549275(100.0%)    7549896(100.0%)
       100    7474051( 99.0%)    7545419( 99.9%)    7549539(100.0%)
       150    7434738( 98.5%)    7538967( 99.9%)    7549125(100.0%)
       200    7404908( 98.1%)    7534635( 99.8%)    7548555(100.0%)
       250    7328084( 97.1%)    7496965( 99.3%)    7520145( 99.6%)
       300      10081(  0.1%)      12023(  0.2%)      12877(  0.2%)
       350       2040(  0.0%)       3694(  0.0%)       4853(  0.1%)
       400       1326(  0.0%)       2806(  0.0%)       4252(  0.1%)
       450          4(  0.0%)          9(  0.0%)         22(  0.0%)



# run2 MERGE -  Merging the sequences from Run 2 - Make sure you are in the folder with the extracted forward and reverse reads from Run 2
    #https://www.drive5.com/usearch/manual/merge_options.html
    mkdir mergedfastq_run2
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/research/EvansLab/MMPRNT_2017/20180119_B_16S-V4_PE/*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/combined_merged_run2.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/combined_merged_merge_run2_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/combined_merged_merge_run2_alnout_report.txt

    #02:25 162Mb   100.0% 76.0% merged
    #Merged length distribution:
    #        16  Min
    #       253  Low quartile
    #       253  Median
    #       253  High quartile
    #       484  Max

    #Totals:
    #  11690596  Pairs (11.7M)
    #   8886999  Merged (8.9M, 76.02%)"
    #   4335488  Alignments with zero diffs (37.09%)
    #   2766904  Too many diffs (> 5) (23.67%)
    #     36693  No alignment found (0.31%)
    #         0  Alignment too short (< 16) (0.00%)
    #     69233  Staggered pairs (0.59%) merged & trimmed
    #    246.51  Mean alignment length
    #    253.04  Mean merged length
    #      0.37  Mean fwd expected errors
    #      1.15  Mean rev expected errors
    #      0.08  Mean merged expected errors


    #let's check sequence quality of the merged seqs - USEARCH has a great command for checking the expected error and lengths of the sequences.
    #https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html
    mkdir /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/fastq_info/
    cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 combined_merged_run2.fastq -output fastq_info/combined_merged_run2_eestats2.txt

        WARNING: Max OMP threads 1

    00:27 38Mb    100.0% Reading reads

    8886999 reads, max len 484, avg 253.0"

    Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
    ------   ----------------   ----------------   ----------------
        50    8852082( 99.6%)    8884817(100.0%)    8885695(100.0%)
       100    8786614( 98.9%)    8879154( 99.9%)    8884573(100.0%)
       150    8731933( 98.3%)    8869639( 99.8%)    8883200(100.0%)
       200    8680162( 97.7%)    8859378( 99.7%)    8879980( 99.9%)
       250    8567947( 96.4%)    8797815( 99.0%)    8832906( 99.4%)
       300      10623(  0.1%)      14200(  0.2%)      16016(  0.2%)
       350       3644(  0.0%)       7121(  0.1%)       9822(  0.1%)
       400       2002(  0.0%)       4684(  0.1%)       8001(  0.1%)
       450         24(  0.0%)         87(  0.0%)        239(  0.0%)

  # Rename the merged file to get rid of hyphens in sample names ""MMPRNT-#"" should be ""MMPRNT#"""
           sed -e 's/MMPRNT-/MMPRNT/g' mergedfastq_run2/combined_merged_run2.fastq > mergedfastq_run2/renamed_combined_merged_run2.fastq
           # check sample names to make sure hyphen is gone:
           less mergedfastq_run2/renamed_combined_merged_run2.fastq

     # cut-adapt to get rid of adapters (adapters for 515F and reverse complement of 806R)
       # FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
       cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
        module load cutadapt/1.8.1
        cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o mergedfastq_run2/cut_renamed_combined_merged_run2.fastq mergedfastq_run2/renamed_combined_merged_run2.fastq > mergedfastq_run2/cut_adpt_results_cut_renamed_combined_merged_run2.txt

        #Total reads processed:               8,886,999"
       #Reads with adapters:                     5,090 (0.1%)"
        #Reads written (passing filters):     8,886,999 (100.0%)"
        #Total basepairs processed: 2,248,777,762 bp"
        #Total written (filtered):  2,248,661,956 bp (100.0%)"



# run3 MERGE -     Merging the sequences from Run 3 - Make sure you are in the folder with the extracted forward and reverse reads from Run 3
    # NOTE: these samples had really low merge (30%) with defaults (maxdiffs 5) so troubleshooting is needed first
    # TROUBLESHOOTING:
      # merge a subsample of the sequences to invesitgate issues:
      /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/home/chicoin1/SGvarietyRhizo16S/raw/*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test_report.fastq -alnout  /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test_align.fastq

      # try to merge pairs with fastq_minqual option
      /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/*R1*.fastq -fastq_minqual 1 -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test_report.fastq -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/test/merged_test_align.fastq
        # didn't improve merge
      # https://drive5.com/usearch/manual/pipe_readprep_merge.html
       # maxdiffs (# of bp mismatches allowed, default is 5)
        # pctid (Minimum %id of alignment. Default 90
      # also, keep in mind that the fwd reads have eemax of .39 whereas reverse are 1.8 (fwd reads are in much better quality),"
        # so this may suggest that we complete the analysis on just the forward read for these sequences
        # or trim the reverse reads before merging (look at alignment with -alnout .txt in mergepairs step to see where to trim to), in this case it looks like at bp 80 the quality gets shitty"
        # use fastx_truncate to trim the Sequences (https://www.drive5.com/usearch/manual/cmd_fastx_truncate.html)
          # we are going to strip right 80 bp on our reverse read  (-stripright 80)
          cd to folder with the reverse reads:
          nano reverse_truncate_trim.sh
            ###
            for fq in *R2_001.fastq
            do
            /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_truncate $fq -stripright 80 -fastqout trim_$fq
            done
            ###
            chmod +x reverse_truncate_trim.sh
            ./reverse_truncate_trim.sh

       # Next, the samples were merged with usearch (Lukas did this)"
    # cut-adapt to get rid of adapters (adapters for 515F and reverse complement of 806R)
      # FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
      # note: this step was done on lukas's directory for the SGvariety Rhizo/blank namibia files
      cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o mergedfastq_run2/cut_combined_merged_run2.fastq mergedfastq_run2/combined_merged_run2.fastq > mergedfastq_run2/cut_adpt_results_cut_combined_merged_run2.txt



############## RUN 3 MERGE
       mkdir mergedfastq_run3

      /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/research/EvansLab/MMPRNT_2017/20180119_B_16S-V4_PE/*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/combined_merged_run2.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/combined_merged_merge_run2_report.txt

      #let's check sequence quality of the merged seqs - USEARCH has a great command for checking the expected error and lengths of the sequences.
      #https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html
      mkdir /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/mergedfastq_run2/fastq_info/

      /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 combined_merged_run2.fastq -output fastq_info/combined_merged_run2_eestats2.txt
```
## 2) Combine the three runs with merged reads 
```{r}
    cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
    cat mergedfastq_run1/cut_renamed_combined_merged_run1.fastq mergedfastq_run2/cut_renamed_combined_merged_run2.fastq mergedfastq_SGvarietyRhizo/cut_trimR_combined_merged_SGvarietyRhizo.fastq> combined_cut_merged_ALLruns.fastq

    # how big is this file?
    wc -l combined_cut_merged_ALLruns.fastq    # 91122628

# 4b) Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step."
    #Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs"
    #https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html

    cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_get_sample_names combined_cut_merged_ALLruns.fastq -output combined_cut_merged_ALLruns_samplenames.txt
      # 100.0% 794 samples found

    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 combined_cut_merged_ALLruns.fastq -output ccombined_cut_merged_ALLruns_eestats2.txt

    #22780657 reads, max len 484, avg 251.9"

    #Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
    ------   ----------------   ----------------   ----------------
        50   22715806( 99.7%)   22777251(100.0%)   22778831(100.0%)
       100   22530470( 98.9%)   22763021( 99.9%)   22777195(100.0%)
       150   22363282( 98.2%)   22736147( 99.8%)   22774991(100.0%)
       200   22198437( 97.4%)   22707461( 99.7%)   22770652(100.0%)
       250   20559621( 90.3%)   21199671( 93.1%)   21306200( 93.5%)
       300      81083(  0.4%)     107838(  0.5%)     121655(  0.5%)
       350      48530(  0.2%)      76189(  0.3%)      95502(  0.4%)
       400       3327(  0.0%)       7490(  0.0%)      12252(  0.1%)
       450         28(  0.0%)         95(  0.0%)        256(  0.0%)
       
```       

## 3) Filter the merged and truncated seqs to MaxEE and bp length 
```{r}
# https://www.drive5.com/usearch/manual/cmd_fastq_filter.html
    # set maxee to 1 and trunclength to 250 (which is how many bp overlaps we should have with these primers)
    # NOTE: this file should be used for the otu clustering, but the unfiltered file combined_cut_merged_ALLruns.fastq should be used for otumapping step"
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter  combined_cut_merged_ALLruns.fastq -fastq_maxee 1 -fastq_trunclen 250 -fastqout combined_cut_merged_ALLruns_trunc250_maxee1.fastq

    #02:50 4.1Mb   100.0% Filtering, 93.1% passed"
    #  22780657  Reads (22.8M)
    #   1472015  Discarded reads length < 250
    #    108971  Discarded reads with expected errs > 1.00
    #  21199671  Filtered reads (21.2M, 93.1%)"
```

## 4) Dereplication (Filter so we only have unique sequences)
```{r}
    # https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques combined_cut_merged_ALLruns_trunc250_maxee1.fastq  -fastqout uniques_combined_cut_merged_ALLruns_trunc250_maxee1.fastq -sizeout

    #01:41 13.8Gb  100.0% DF
    #01:42 14.1Gb 21199671 seqs, 5132869 uniques, 3295410 singletons (64.2%)"
    #01:42 14.1Gb Min size 1, median 1, max 541909, avg 4.13"
    #03:22 13.8Gb  100.0% Writing uniques_combined_cut_merged_ALLruns_trunc250_maxee1.fastq
```

## 5) Denoise (Cluster into OTUS (97%) and filter out singletons)
```{}
    # First lets cluster at 97% identity using uparse https://www.drive5.com/usearch/manual/cmd_cluster_otus.html
    # this also denovo chimera filters and filters singletons (OTUs with less than 2 reads) by default
    # I submitted a job for this:
      #!/bin/sh -login
      #PBS -o /mnt/home/chicoin1/SGvarietyRhizo16S/ # output for hpcc
      #PBS -j oe #yes, send me error messages"
      #PBS -l nodes=1:ppn=20,walltime=06:00:00,mem=128gb # time in the que"
      #PBS -M ""chicoinet@gmail.com"" # where to send notice of job"
      #PBS -m abe # sends email when job is aborted
      #PBS -N cluster_otus.sh # job name
      #PBS -r n
      cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
otu
     /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus uniques_combined_cut_merged_ALLruns_trunc250_maxee1.fastq -otus combined_merged_both_runs_otus.fa -uparseout combined_merged_both_runs_otus_uparse.txt -relabel OTU
     ####
     chmod +x cluster_otus.sh
     qsub cluster_otus.sh # 52158489.mgr-04.i
```     
     

## 6) Map the reads back to OTUS (97%)
```{r}
    #https://www.drive5.com/usearch/manual/cmd_otutab.html
    #-id 0.97 -strand plus are defaults
    #-otutab (original unfiltered merged file)
    #-otus (otus.fa from cluster_otus step)

    # I submitted a job for this:
      nano otutab_map.sh
      #!/bin/sh -login
      #PBS -o /mnt/home/chicoin1/SGvarietyRhizo16S/ # output for hpcc
      #PBS -j oe #yes, send me error messages"
      #PBS -l nodes=1:ppn=20,walltime=06:00:00,mem=128gb # time in the que"
      #PBS -M ""chicoinet@gmail.com"" # where to send notice of job"
      #PBS -m abe # sends email when job is aborted
      #PBS -N otutab_map.sh # job name
      #PBS -r n

      cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/

      /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab combined_cut_merged_ALLruns.fastq -otus combined_merged_both_runs_otus.fa -uc combined_merged_ALLruns_OTU_map_merg.uc -otutabout combined_merged_ALLruns_OTU_table.txt -biomout combined_merged_ALLruns_OTU_jsn.biom -notmatchedfq combined_merged_ALLruns_otu_unmapped.fq
      ####
      chomd +x otutab_map.sh
      qsub otutab_map.sh # 52169240.mgr-04.i
```      

## 7) Classify clustered otus against the reference database (silva vs. 123)
```{r}
  #https://www.drive5.com/usearch/manual/cmd_sintax.html
  # sintax: otus.fa output from cluster_otus
  # Usearch silva file, vs. 123 (""/mnt/research/EvansLab/Databases/Silva132_release/silva_16s_v123.fa"")
  # can specify -sintax_cutoff 0.8 to improve predictiability of taxonomy ranking, but can also filter this later

  # I submitted a job for this:
    nano taxonomy_otus_zotus.sh
    #!/bin/sh -login
    #PBS -o /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/ # output for hpcc
    #PBS -j oe #yes, send me error messages"
    #PBS -l nodes=1:ppn=20,walltime=40:00:00,mem=128gb # time in the que"
    #PBS -M ""chicoinet@gmail.com"" # where to send notice of job"
    #PBS -m abe # sends email when job is aborted
    #PBS -N taxonomy_otus_zotus.sh # job name
    #PBS -r n

  cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax combined_merged_both_runs_otus.fa -db /mnt/research/EvansLab/Databases/Silva132_release/silva_16s_v123.fa -tabbedout combined_merged_ALLruns_otus_taxonomy.sintax -strand both
  
```  
## 8) Make a phylogentic tree file with PASTA
```{r}
        # https://github.com/smirarab/pasta/blob/master/pasta-doc/pasta-tutorial.md

        # First a tree file for the OTUs

        module load Python/2.7.2
        # Lukas's example: python /mnt/research/EvansLab/Software/pasta-code/pasta/run_pasta.py -i representative_seq.fa -j job_name_here -o output_directory_because_many_files_produced/
        # practice from pasta sample.
        python /mnt/research/EvansLab/Software/pasta-code/pasta/run_pasta.py -i /mnt/research/EvansLab/Software/pasta-code/pasta/data/small.fasta -j pasta_practice -o /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/Pasta/Practice
          # output: the tree is saved in a file called [jobname].tre and the alignment file is named [jobname].marker001.small.aln
           # tree: ""/mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/Pasta/Practice/pasta_practice.tre"""
            # alignment file:  ""/mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/Pasta/Practice/pasta_practice.marker001.small.aln""
            # pastajob.score.txt - liklihood score produced by PASTA's internal ML tool
            # pastajob_temp_pasta_config.txt. This file contains all the PASTA configurations used in the current run. You should 1) inspect this file to make sure the configurations are what you intended, and 2) keep this file for future reference (so that you know what exact options where used). We will see in the next steps how configuration files can be used for running PASTA with a reproducible set of configurations."

        # I submitted a job to the hpcc to run this code
        #!/bin/sh -login
           #PBS -o /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/ # output for $
           #PBS -j oe #yes, send me error messages"
           #PBS -l nodes=1:ppn=20,walltime=40:00:00,mem=128gb # time in the que"
           #PBS -M ""chicoinet@gmail.com"" # where to send notice of job"
           #PBS -m abe # sends email when job is aborted
           #PBS -N pasta_tree.sh # job name
           #PBS -r n

           cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
           module load Python/2.7.2
           python /mnt/research/EvansLab/Software/pasta-code/pasta/run_pasta.py -i /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/combined_merged_both_runs_otus.fasta -j SGvarietyRhizoEndo16s_pasta -o /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/Pasta/

        
```

# Fungal Sequencing Pipeline (primers ITS1-F/ITS2)
  # SGvarietyRhizo samples and RL samples
  # info about run in document: C:\Users\Owner\Documents\01_Grad School\06_Sequencing\SGvarietyRhizo\MiSeqRuns_16S_ITSinfo.txt
  # https://github.com/GLBRC-TeamMicrobiome/ITS-amplicon-pipeline/blob/master/ITS-amplicon-pipeline.md

## 1) Check read quality & merge reads 
```{r}
# Only 1 run for the fungal soil sequences (cd /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/), n = 144 samples

    # read this to learn more about read quality: https://www.drive5.com/usearch/manual/exp_errs.html
    # remember that a low EEmax is high quality, and that a normal threshold is 1.0
  cd /mnt/home/chicoin1/SGvarietyRhizoITS/

  # 1a) make directory to store fastq info files
    mkdir fastq_info

    # make nano file with code to move all fastq info into this directory
    nano fasta_info_fq.sh

    # code to put in fasta_info_fq.sh file
    !#/bin/bash
    for fq in *.fastq
    do
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_info $fq -output /mnt/home/chicoin1/SGvarietyRhizoITS/fastq_info/$fq
     done

    # activate .sh file
    chmod +x fasta_info_fq.sh

    # use fasta_info_fq.sh to get info for every sample
    # first navigate to the folder with the .fastq files
    # this will take a little while to calculate the fastainfo for every sample
    cd /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/
    /mnt/home/chicoin1/SGvarietyRhizoITS/fasta_info_fq.sh

    # go back to the fastaq_info file you just generated and pull out information for every sample
    /mnt/home/chicoin1/SGvarietyRhizoITS/fastq_info/

    # summarize the expected errors and puts into a text file
    # remember: EE is expected errors; small E is high quality and low E is high quality; we can use fastq_filter to filter all reads with a low quality (set threshold);  A natural choice is to set E_max = 1 because the most probable read errors = 0, but to be more stringent, can set E_max = 0.5 or 0.25
    grep "^EE" * > ITS_allraw_EE.txt
    # can also summarize the sequence length and size into a text file
    grep "^File" * > ITS_allraw_seqlength.txt


  # 1b) You can also use fastqc to look at your reads quality
    module load FastQC
    # merge all of your .fastq files into a single R1 and R2 file
      # note: original .fastq files are in EvansLab folder and I created a directory in my folder to store the cat file
    cat /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/*R1_001.fastq > /mnt/home/chicoin1/SGvarietyRhizoITS/fastqc/raw_reads_R1.fastq
    cat /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/*R2_001.fastq > /mnt/home/chicoin1/SGvarietyRhizoITS/fastqc/raw_reads_R2.fastq

    # produce reads quality graph using fastqc
      # per_base_quality is a useful figure and will give you an indication of the quality of the reads by position
    fastqc raw_reads_R1.fastq
    fastqc raw_reads_R2.fastq
      # this creates a .zip file with files that summarize the quality of your Reads; unzip the file
    unzip *.zip
    # download the html files to your computer and you can open them in your web-browser to look at the quality of your reads.


#################### 
# Merge reads
####################
    # you need to be in the folder with the .fastq file #https://www.drive5.com/usearch/manual/merge_options.html
    mkdir /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/
    cd /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -fastq_maxdiffs 5 -fastq_merge_maxee 1 -fastqout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxxee1_merged.fq -relabel @

    ## RESULTS (maxdiffs: 5, eemax = 1)
    02:42 302Mb   100.0% 54.3% merged

      Totals:
         7063977  Pairs (7.1M)
         3835431  Merged (3.8M, 54.30%)
         1634088  Alignments with zero diffs (23.13%)
         3030109  Too many diffs (> 5) (42.90%)
          198437  No alignment found (2.81%)
               0  Alignment too short (< 16) (0.00%)
               0  Exp.errs. too high (max=1.0) (0.00%)
          516579  Staggered pairs (7.31%) merged & trimmed
          213.02  Mean alignment length
          282.21  Mean merged length
            0.81  Mean fwd expected errors
            0.89  Mean rev expected errors
            0.16  Mean merged expected errors

#We need to improve the merge (we only had 54.3% merge)
    # https://www.drive5.com/usearch/manual/merge_badrev.html


    # Reverse truncate the reverse reads to try and improve merging
        # I cannot do this with the raw reads in the EvansLab folder, so I'm going to copy all of the raw reads to my directory
          mkdir /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads
          cp  /mnt/research/EvansLab/Namibia_LTER/20170224_ITS1_PE/*R*.fastq /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads

          # looking at the fastqc report, it seems that we have low quality bp calls from bp 250 - 209 (Q score median <25)
          # use fastx_truncate to trim the Sequences (https://www.drive5.com/usearch/manual/cmd_fastx_truncate.html)
          # we are going to strip right 50 bp from the end of our reverse read  (-stripright 50)

          # make folder in my directory for truncated reverse reads
            mkdir /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads/R2trunc50/

          # cd to folder with the reverse reads:
          cd  /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads
          nano TMC_R2_reverse_truncate.sh
            ###
            for fq in *R2_001.fastq
            do
            /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_truncate $fq -stripright 50 -fastqout /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads/R2trunc50/trim_$fq
            done
            ###
            chmod +x TMC_R2_reverse_truncate.sh
            ./TMC_R2_reverse_truncate.sh

        # Merge forward-R1 with 50bp-truncated-R2
        cd  /mnt/home/chicoin1/SGvarietyRhizoITS/rawreads
        /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -reverse R2trunc50/*R2*.fastq -fastqout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged.fastq  -report /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged_align.txt

        ## Increased our merge to 62% - let's go forward with this...
        #02:42 335Mb   100.0% 62.1% merged
        Merged length distribution:
                19  Min
               266  Low quartile
               278  Median
               301  High quartile
               434  Max
        Totals:
           7063977  Pairs (7.1M)
           4388857  Merged (4.4M, 62.13%)
           2266303  Alignments with zero diffs (32.08%)
           2418930  Too many diffs (> 5) (34.24%)
            256190  No alignment found (3.63%)
                 0  Alignment too short (< 16) (0.00%)
             37402  Staggered pairs (0.53%) merged & trimmed
            164.45  Mean alignment length
            283.24  Mean merged length
              0.92  Mean fwd expected errors
              0.48  Mean rev expected errors
              0.26  Mean merged expected errors


# check the quality of merged pairs
  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged.fastq -output /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged_eestats2.txt

  4388857 reads, max len 434, avg 283.2

  Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
  ------   ----------------   ----------------   ----------------
      50    4383781( 99.9%)    4387752(100.0%)    4387784(100.0%)
     100    4282148( 97.6%)    4382786( 99.9%)    4385195( 99.9%)
     150    4144047( 94.4%)    4362225( 99.4%)    4384288( 99.9%)
     200    4054132( 92.4%)    4311464( 98.2%)    4356218( 99.3%)
     250    3645168( 83.1%)    3932229( 89.6%)    3992208( 91.0%)
     300     883185( 20.1%)    1045949( 23.8%)    1101926( 25.1%)
     350     106990(  2.4%)     142985(  3.3%)     167127(  3.8%)
     400       8759(  0.2%)      14139(  0.3%)      18907(  0.4%)
```

## 2) Remove adapters with cutadapt before merging (ITS1 and ITS2 primers)
```{r}
    # https://cutadapt.readthedocs.io/en/stable/guide.html

      > ITS1F adapter: AATGATACGGCGACCACCGAGATCTACAC (http://www.earthmicrobiome.org/protocols-and-standards/its/)
      > Reverse compliment of ITS1F adapter: CAAGCAGAAGACGGCATACGAGAT
      > Forward: ITS1F primer
        CTTGGTCATTTAGAGGAAGTAA
      > Reverse: ITS2 primer
        GCTGCGTTCTTCATCGATGC
      > ITS2 Reverse-Compliment:
        GCATCGATGAAGAACGCAGC

  # cut-adapt to get rid of adapters
    # cutadapt -a "whattocut" -o "output file" "inputfile" > results.txt
    # - g ^forwardread -a revcomp -f fastq -n 2 --discard-untrimmed --match-read-wildcards
    # -n 2 : how many times you expect to see the primers in your sequences
    # --discard-untrimmed will remove any sequences that didn't find the adapter/primers
         cd /mnt/home/chicoin1/SGvarietyRhizo1TS/
     module load cutadapt/1.8.1
     cutadapt  -a GCATCGATGAAGAACGCAGC -g ^CTTGGTCATTTAGAGGAAGTAA -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/cutadapt_diff5_maxee1_R2trunc50_merged.fastq /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/diff5_maxee1_R2trunc50_merged.fastq > /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/cutadaptSummary_diff5_maxee1_R2trunc50_merged.txt

     #SUMMARY
     #Total reads processed:               4,388,857
     #Reads with adapters:                 4,388,809 (100.0%)
     #Reads written (passing filters):     4,388,809 (100.0%)
     #Total basepairs processed: 1,243,112,117 bp
     #Total written (filtered):  1,059,711,383 bp (85.2%)
```

## 3) Filter the merged and truncated seqs to MaxEE and bp length 
```{r}
  # no need to trim since the ITS can be super variable in length https://www.drive5.com/usearch/manual/global_trimming_its.html
  # let's check the quality of the merged and truncated seqs
  # can play with the eemax value that you filter to, don't go lower than 1.0!

  # first, check eescores for merged and cut sequences
  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/cutadapt_diff5_maxee1_R2trunc50_merged.fastq -output  /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/cutadapt_diff5_maxee1_R2trunc50_merged_eestats2.txt

  #  4388809 reads, max len 413, avg 241.5

  #Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
  ------   ----------------   ----------------   ----------------
  #    50    4375276( 99.7%)    4385409( 99.9%)    4385459( 99.9%)
  #   100    4265118( 97.2%)    4379094( 99.8%)    4384564( 99.9%)
  #   150    4159937( 94.8%)    4347074( 99.0%)    4373700( 99.7%)
  #   200    3809430( 86.8%)    4022554( 91.7%)    4066068( 92.6%)
  #   250    1320976( 30.1%)    1469731( 33.5%)    1517809( 34.6%)
  #   300     144126(  3.3%)     182536(  4.2%)     207542(  4.7%)
  #   350      20707(  0.5%)      29827(  0.7%)      37454(  0.9%)
  #   400          8(  0.0%)         12(  0.0%)         18(  0.0%)


  #now filtering the merged and truncated seqs to MaxEE (1)

  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/cutadapt_diff5_maxee1_R2trunc50_merged.fastq -fastq_maxee 1 -fastaout  /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fasta -fastqout  /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fastq

  #SUMMARY
    #00:00 38Mb   CPU has 20 cores, defaulting to 10 threads
    #01:02 627Mb   100.0% Filtering, 98.3% passed
    #4388809  Reads (4.4M)
    #76570  Discarded reads with expected errs > 1.00
    #4312137  Filtered reads (4.3M, 98.3%)
```
## 4) Dereplication (Filter so we only have unique sequences)
```{r}
  #https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html

  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fasta -fastaout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fa -sizeout -relabel Uniq

```
## 5) Denoise (Cluster OTUs (97%) and filter out singletons)
```{r}
  # https://www.drive5.com/usearch/manual/cmd_cluster_otus.html
  # this also chimera filters...
  # to get more fine-scale otu separation, can use unoise command (99% cluster)
  # -minsize # will filter out any otus that are present <# times (e.g. # = 2 filters out singletons)

  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fa  -otus /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/OTUS_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fa  -uparseout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/uparseOTUS_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.txt  -relabel OTU
```

## 6) Classify clustered otus against the reference database (Unite vs. 7.2)
```{r}
  #https://www.drive5.com/usearch/manual/cmd_sintax.html
  # input: otus file from cluster_otus step
  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/OTUS_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.fa -db /mnt/research/EvansLab/Databases/UNITE_7.2/utax_reference_dataset_10.10.2017.fasta -tabbedout /mnt/home/chicoin1/SGvarietyRhizoITS/Oct2018_usearchOUT/Taxonomy_unite7.2_OTUS_uniques_maxee_1_cutadapt_diff5_maxee1_R2trunc50_merged.sintax -strand both

  #output
    #00:00 37Mb      0.1% Reading /mnt/research/EvansLab/Databases/UNITE_7.2/utax_reference00:01 37Mb      0.2% Reading /mnt/research/EvansLab/Databases/UNITE_7.2/utax_reference00:01 94Mb    100.0% Reading /mnt/research/EvansLab/Databases/UNITE_7.2/utax_reference_dataset_10.10.2017.fasta
    #00:01 60Mb      0.1% Masking (fastnucleo)                                             00:02 60Mb    100.0% Masking (fastnucleo)
    #00:03 61Mb    100.0% Word stats
    #00:03 61Mb    100.0% Alloc rows
    #00:06 214Mb   100.0% Build index
    #00:06 224Mb   100.0% Initialize taxonomy data
    #00:06 226Mb   100.0% Building name table
    #00:06 226Mb  30794 names, tax levels min 1, avg 5.7, max 7
    #WARNING: 2 taxonomy nodes have >1 parent
    #00:06 259Mb  CPU has 20 cores, defaulting to 10 threads
    #00:39 891Mb   100.0% Processing



############# NOTE #####################
#########################################################################
# Lukas Bell-Dereske filtered the unite7.2 database to make it usable for our pipelines
###THIS CODE WAS RUN TO PREPARE THE UNITE DATABASE AND DOES NOT NEED TO BE RUN AGAIN 
#first need to prep the database. I am using the UNITE 7.1 2016-11-20 release https://unite.ut.ee/repository.php
#path is now /mnt/research/EvansLab/Databases/UNITE_7.1

#First I am going to split the utax_reference_dataset_20.11.2016.fasta into the regions ITS1 and ITS2
#UNITE recomends using the ITSx tool http://microbiology.se/publ/itsx_users_guide.pdf

ITSx -i /mnt/research/EvansLab/Databases/UNITE_7.2/utax_reference_dataset_10.10.2017.fasta -o /mnt/research/EvansLab/Databases/UNITE_7.1/ITSx_split -t "f"
#-t option set the taxon you would like to compare againest "f" is fungi.

#ITSx -- Identifies ITS sequences and extracts the ITS region
#by Johan Bengtsson-Palme et al., University of Gothenburg
#Version: 1.0.11
-----------------------------------------------------------------
#Wed May 31 16:26:39 2017 : Preparing HMM database (should be quick)...
#Wed May 31 16:26:39 2017 : Checking and handling input sequence data (should not take long)...
#Wed May 31 16:26:43 2017 : Comparing sequences to HMM database (this may take a long while)...
#    Wed May 31 18:25:29 2017 : Fungi analysis of main strand finished.
#    Wed May 31 20:08:58 2017 : Fungi analysis of complementary strand finished.
#Wed May 31 20:08:58 2017 : Analysing results of HMM-scan (this might take quite some time)...
#Wed May 31 20:09:10 2017 : Extraction finished!
#-----------------------------------------------------------------
#Thank you for using ITSx!
#Please report bugs or unsupported lineages to itsx@microbiology.se
#NOTE!!! this is not necessary for RDP. See below from the README file on the public RDP database
#       Based on our experience, trimming the sequences to a specific region does not improve accuracy. The ranks are not required to be uniform neither, which means you can define any number of ranks as necessary. The speed of the## Classifier is proportional to the number of genera, not the number of training sequences

#I am going to remove all of the extraneous info that ITSx add at the end of each sequence name (i.e. |F|ITS1....)

sed s/'|F|ITS1.*'// ITSx_split.ITS1.fasta > reference_dataset_20.11.2016.ITS1.fasta
#reg exp '|F|ITS1.*' finds occurence of |F|ITS1 followed by any number of occurences and then I am replacing them with nothin //

#now removing all of the empty taxon levels
sed -i s/o:,//  reference_dataset_20.11.2016.ITS1.fasta

#I am going to make a sintax database

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -makeudb_sintax /mnt/research/EvansLab/Databases/UNITE_7.1/reference_dataset_20.11.2016.ITS1.fasta -output /mnt/research/EvansLab/Databases/UNITE_7.1/sintax_reference_dataset_20.11.2016.ITS1.udb

#00:01 61Mb    100.0% Reading /mnt/research/EvansLab/Databases/UNITE_7.1/reference_dataset_20.11.2016.ITS1.fasta
#00:01 27Mb    100.0% Converting to upper case
#00:01 29Mb    100.0% Word stats
#00:01 29Mb    100.0% Alloc rows
#00:01 69Mb    100.0% Build index
#00:02 78Mb    100.0% Initialize taxonomy data
#00:02 80Mb    100.0% Building name table
#00:02 80Mb   28857 names, tax levels min 1, avg 5.5, max 7
#00:02 80Mb   Buffers (54234 seqs)
#00:02 97Mb    100.0% Seqs

```

## 7)  Map the reads back to OTUS (97%)
```{r}
  # https://www.drive5.com/usearch/manual/cmd_otutab.html
    # -otutab fastqoutput from -fastq_mergepairs
    # -otus otu.fa file from -cluster_otus

  # This step takes a long time so I am going to create a job
    #  https://wiki.hpcc.msu.edu/display/TEAC/Job+Scripts

  cd SGvarietyRhizoITS/Oct2018_usearchOUT/
  nano OTU_ITS_mapping_TMC.sb

  #!/bin/bash --login
  ########## Define Resources Needed with SBATCH Lines ##########

  #SBATCH --time=06:00:00             # limit of wall clock time - how long the job will takes
  #SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes)
  #SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
  #SBATCH --mem=2G                    # memory required per node - amount of memory
  #SBATCH --job-name OTU_ITS_mapping_TMC.sb      # you can give your job a name

  cd /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/

  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab Nash_merged.fq -otus /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/merged_NOTstripped_OTUs_nonchimera.fa -otutabout /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/TMC_NamibLTER_ITS_otutab.txt -biomout /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/TMC_NamibLTER_ITS_otutab.json -mapout /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/TMC_NamibLTER_ITS_OTU_map.txt -notmatchedfq /mnt/home/chicoin1/SGvarietyRhizoITS/usearchOUT/TMC_NamibLTER_ITS_otus_ref_nonchimera_unmapped.fq

  ####

  #start job
  sbatch OTU_ITS_mapping_TMC.sb
  #Submitted batch job 1246551
```




