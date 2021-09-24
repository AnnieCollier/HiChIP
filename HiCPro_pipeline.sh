#!/bin/bash
#SBATCH --job-name=HiCPro --output=HiCPro.out --error=HiCPro.err --time=01:00:00 --qos=normal --nodes=1 --mem-per-cpu=4000 --ntasks-per-node=8    

ml biology hic-pro
ml samtools
ml bowtie2

## PATH/Main folder of fastq files
fastq="fastq/WT" 
## Name of output folder created
output="AHDC1_WT_HiChIP" 
##Config file name
config_file="config_test_latest.txt"

HiC-Pro -i $fastq -o $output -c $config_file  -p
