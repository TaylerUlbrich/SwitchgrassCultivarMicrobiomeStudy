# Tayler Ulbrich
# Closedref greengenes reference for picrust analysis

# Resources:
  *http://picrust.github.io/picrust/tutorials/otu_picking.html#using-an-open-reference-otu-table-from-qiime*
  *https://drive5.com/usearch/manual/cmd_closed_ref.html http://picrust.github.io/picrust/tutorials/otu_picking.html#using-an-open-reference-otu-table-from-qiime*
  *https://drive5.com/usearch/manual/cmd_closed_ref.html*

# 1) Merge fastq files from each run
	mkdir TMCmergedONLY

	  #run1
		/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/research/EvansLab/MMPRNT_2017/20180113_16S-V4_PE/TMC*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_run1.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_merge_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_merge_alnout.txt

		00:43 76Mb    100.0% 76.9% merged

		Merged length distribution:
				16  Min
			   253  Low quartile
			   253  Median
			   253  High quartile
			   473  Max

		Totals:
			636432  Pairs (636.4k)
			489284  Merged (489.3k, 76.88%)
			250307  Alignments with zero diffs (39.33%)
			145334  Too many diffs (> 5) (22.84%)
			  1814  No alignment found (0.29%)
				 0  Alignment too short (< 16) (0.00%)
			  7486  Staggered pairs (1.18%) merged & trimmed
			244.49  Mean alignment length
			254.55  Mean merged length
			  0.47  Mean fwd expected errors
			  0.98  Mean rev expected errors
			  0.09  Mean merged expected errors


	  #run2
		/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/research/EvansLab/MMPRNT_2017/20180119_B_16S-V4_PE/TMC*R1*.fastq -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_run2.fastq  -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_merge_run2_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCcombined_merged_merge_run2_alnout_report.txt

			  03:37 98Mb    100.0% 75.4% merged

		  Merged length distribution:
				  16  Min
				 253  Low quartile
				 253  Median
				 253  High quartile
				 484  Max

		  Totals:
			 3795293  Pairs (3.8M)
			 2860291  Merged (2.9M, 75.36%)
			 1382145  Alignments with zero diffs (36.42%)
			  916968  Too many diffs (> 5) (24.16%)
			   18034  No alignment found (0.48%)
				   0  Alignment too short (< 16) (0.00%)
			   41080  Staggered pairs (1.08%) merged & trimmed
			  245.99  Mean alignment length
			  253.11  Mean merged length
				0.39  Mean fwd expected errors
				1.16  Mean rev expected errors
				0.08  Mean merged expected errors

	  #run3
		Need to truncate reverse reads first to improve merge
		mkdir /mnt/home/chicoin1/SGvarietyRhizo16S/raw/truncated_rev
		nano truncated_rev.sh
		  ###
		  for fq in *R2_001.fastq
		  do
		  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_truncate $fq -stripright 80 -fastqout truncated_rev/trim_$fq

		  done
		  ###
		  chmod +x reverse_truncate_trim.sh
		  ./reverse_truncate_trim.sh


		  # merge truncated reverse with forward read
		  /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs /mnt/home/chicoin1/SGvarietyRhizo16S/raw/*R1*.fastq -reverse /mnt/home/chicoin1/SGvarietyRhizo16S/raw/truncated_rev/*R2* -relabel @ -fastqout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCrhizo_merged.fastq -report /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCrhizo_merged_report.txt -alnout /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/TMCrhizo_merged_alnout_report.txt


		  02:55 101Mb   100.0% 69.4% merged

		  Merged length distribution:
				  16  Min
				 253  Low quartile
				 253  Median
				 253  High quartile
				 398  Max

		  Totals:
			 3837965  Pairs (3.8M)
			 2663829  Merged (2.7M, 69.41%)
			  981257  Alignments with zero diffs (25.57%)
			 1151127  Too many diffs (> 5) (29.99%)
			   23009  No alignment found (0.60%)
				   0  Alignment too short (< 16) (0.00%)
				 127  Staggered pairs (0.00%) merged & trimmed
			  166.89  Mean alignment length
			  253.09  Mean merged length
				0.60  Mean fwd expected errors
				0.85  Mean rev expected errors
				0.13  Mean merged expected errors


# 2) cutadapt all of these merged files to remove adapters
	  module load cutadapt/1.8.1
	  cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/
	  #run1
	  cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o TMCmergedONLY/cut_TMCcombined_merged_run1.fastq TMCmergedONLY/TMCcombined_merged_run1.fastq > TMCmergedONLY/cutadaptRESULTS_TMCcombined_merged_run1.txt

	  #run2
	  cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o TMCmergedONLY/cut_TMCcombined_merged_run2.fastq TMCmergedONLY/TMCcombined_merged_run2.fastq > TMCmergedONLY/cutadaptRESULTS_TMCcombined_merged_run2.txt

	  #run3
	  cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o TMCmergedONLY/cut_TMCrhizo_merged.fastq TMCmergedONLY/TMCrhizo_merged.fastq > TMCmergedONLY/cutadaptRESULTS_TMCrhizo_merged.txt


# 3) now combine all of these merged files
	  cd TMCmergedONLY
	  cat cut_TMCrhizo_merged.fastq cut_TMCcombined_merged_run1.fastq cut_TMCcombined_merged_run2.fastq > TMConly_cut_combined_merged_ALLruns.fastq

		# check to make sure samplenames are okay
    /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_get_sample_names TMConly_cut_combined_merged_ALLruns.fastq -output TMConly_cut_combined_merged_ALLruns_samplenames.txt


# 4)  Check naming convention of merged files
		  less
		 # Note: The output for this file will have every sequence in every row (ie: TMCA1.1; TMCA1.2) you need to change the naming convention to (TMCA1_, TMCA1_)
		  # example: sed -e 's/\(@TMC.*\)\..*/\1_/g' mergedpairs_original.fq > mergedpairs_reformatted.fq

		 sed -e 's/\(@TMCSGendo.*\)\..*/\1_/g' TMConly_cut_combined_merged_ALLruns.fastq > renamed1.fastq
		 sed -e 's/\(@TMC[A-L].*\)\..*/\1_/g' renamed1.fastq > renamed2.fastq
		 sed -e 's/\(@TMCneighID.*\)\..*/\1_/g' renamed2.fastq > combined_cut_merged_ALLruns_renamed.fastq



 # 5) Then, use closed_ref with the greengenes database (vs. 13_8, 97% otus)
	  a. I submitted a bash file to the hpcc for this:
	  nano closedref.sh

		  #!/bin/sh -login
		  #PBS -o /mnt/home/chicoin1/SGvarietyRhizo16S/ # output for hpcc
		  #PBS -j oe #yes, send me error messages
		  #PBS -l nodes=1:ppn=20,walltime=30:00:00,mem=128gb # time in the que
		  #PBS -M "chicoin1@msu.edu" # sends the notice of job start/end to your email
		  #PBS -m abe # sends email when job is aborted
		  #PBS -N closedref.sh # job name
		  #PBS -r n

		  # CODE:

		  cd /mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/

		   /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -closed_ref mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/combined_cut_merged_ALLruns_renamed.fastq -db /mnt/research/EvansLab/Databases/gg_otus_13_8_release/rep_set/97_otus.fasta -strand plus -otutabout mnt/home/chicoin1/SGvarietyRhizo16S/SGvarietyRhizoEndoALL/TMCmergedONLY/combined_cut_merged_ALLruns_renamed_closedrefOTUtable.fastq
		  ##

		  chmod +x closedref.sh
		  qsub closedref.sh #
		  52455729.mgr-04.i


