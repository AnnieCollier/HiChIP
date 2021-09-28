#!/bin/bash
#SBATCH --job-name=diff --output=diff.out --error=diff.err --time=10:00:00 --qos=normal --nodes=2 --mem-per-cpu=64G

ml R/3.6.1
ml biology bedtools samtools

# code executable for differential loop analysis 
currworkdir=/oak/stanford/groups/oro/anncoll/HiChIP_NovaSeq/differential_loops/
currscriptdir=/home/groups/oro/software/FitHiChIP/
cd $currscriptdir
echo 'current directory containing this script: '$currscriptdir

# Rscript version
RscriptExec=`which Rscript`
echo 'Rscript version installed in the system : '$RscriptExec

# source code (R) of the differential analysis
DiffAnalysisCodeExec=$currscriptdir'/Imp_Scripts/DiffAnalysisHiChIP.r'
echo 'R Code of differential analysis: '$DiffAnalysisCodeExec

BaseOutDir=$currworkdir'/Results_DiffLoops'
mkdir -p $BaseOutDir

$RscriptExec ${DiffAnalysisCodeExec} --AllLoopList $currworkdir'/AHDC1_WT_r1.interactions_FitHiC.bed.gz':$currworkdir'AHDC1_WT_r2.interactions_FitHiC.bed.gz':$currworkdir'/AHDC1_WT_r3.interactions_FitHiC.bed.gz':$currworkdir'/AHDC1_KO_r1.interactions_FitHiC.bed.gz':$currworkdir'/AHDC1_KO_r2.interactions_FitHiC.bed.gz':$currworkdir'/AHDC1_KO_r3.interactions_FitHiC.bed.gz' --ChrSizeFile /home/groups/oro/software/ChromHMM/CHROMSIZES/hg38.txt --FDRThr 0.05 --ChIPAlignFileList /oak/stanford/groups/oro/anncoll/NextSeqChIP/CTCF_WT_AHDC1KO/Mapping/WT1CTCF_srt_rmdup.bam:/oak/stanford/groups/oro/anncoll/NextSeqChIP/CTCF_WT_AHDC1KO/Mapping/KO1CTCF_srt_rmdup.bam --CovThr 25 --OutDir $BaseOutDir --CategoryList 'WT':'KO' --ReplicaCount 3:3 --ReplicaLabels1 "R1":"R2":"R3" --ReplicaLabels2 "R1":"R2":"R3" --FoldChangeThr 2 --DiffFDRThr 0.05 

# revert to the old directory
cd $currworkdir

