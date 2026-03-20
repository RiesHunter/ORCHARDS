This SNIFFLES_V3 folder should be uploaded to your personal CHTC environment. In each folder, 
there will always be .sh files (code) and .sub files (submit file). The code is referred to 
by the submit file, so do not change its name. The submit file contains instructions for 
your run, and is relayed to CHTC with the command: condor_submit <submit_file_name>. In the 
instance of a 1-preprocess run, it will submit 1 job containing all files. For all other runs 
(2-normalize, 3-call, etc.) it will submit a job for every file specified in states.txt. You 
will need to update states.txt with the file names you wish to be included in the run. 

These files must be in the same folder as your submit file (I use cyberduck for uploading/downloading), and these files must be in the proper input format:
Directory	file step	input file	states.txt	refseq specified?
1-preprocess	raw		.fastq.gz	NOT REQUIRED	NOT REQUIRED
2-normalize	filtered	.fastq.gz	REQUIRED	NOT REQUIRED
3-call		normalized	.fastq.gz	REQUIRED	REQUIRED
4-consensus	called		.vcf		REQUIRED	REQUIRED

######################################################

The architecture of this directory is as follows:
SNIFFLES_V3
|- README.txt

|- 1-preprocess
   |- sh_SNIFFLES-1-preprocess.sh
   |- sub-1-preprocess.sub
   |- _states.txt-not-required
   |- _specify-adaptor

|- 2-normalize
   |- sh_SNIFFLES-2-normalize.sh
   |- sub-2-foo.sub
   |- states.txt
   |- _states.txt-required

|- 3-call
   |- sh_SNIFFLES-3-call.sh
   |- sub-3-call.sub
   |- states.txt
   |- _states.txt-required
   |- _refseq-required-in-sub

|- 4-consensus
   |- sh_SNIFFLES-4-consensus.sh
   |- sub-4-consensus.sub
   |- states.txt
   |- _states.txt-required
   |- _refseq-required-in-sub

|- refseqs
   |- H1N1-A_Brisbane_02_2018_HA.fasta
   |- H1N1-A_Michigan_45_2015.fasta
   |- H3N2-A_Kansas_04_2017_H3N2_HA.fasta
   |- H3N2-A_Singapore_INFIMH-16-0019_2016.fasta
   |- Victoria-B_Brisbane_60_2008.fasta
   |- Victoria-B_Colorado_06_2017.fasta
