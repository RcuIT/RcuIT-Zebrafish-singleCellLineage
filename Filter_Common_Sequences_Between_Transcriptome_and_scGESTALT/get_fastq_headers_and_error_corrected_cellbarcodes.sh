#!/bin/bash
#BSUB -J scATAC_Sample_Data
#BSUB -o scATAC_Sample_Data.%J.out
#BSUB -e scATAC_Sample_Data.%J.error
#BSUB -n 8 
#BSUB -M 8000 # Set memory limit to 8 GB
#BSUB -R "span[hosts=1] rusage[mem=8000]"


module load module load  samtools/1.7
# example bam file - Brain_2_linD1_32_possorted_genome_bam.bam
samtools view Brain_1_Lin_All_Concat_possorted_genome_bam.bam | \
awk '{for (i=1; i<=NF; i++) {if ($i ~ /^CB:Z:/) print $1, $i}}' | \
sort -u > Brain.1.1597_1083_unommon.common.CB.Tr.Lin.txt


library(Seurat)


fb <-Read10X_h5("fb_filtered_feature_bc_matrix.h5")
mb <-Read10X_h5("mb_filtered_feature_bc_matrix.h5")
hb <-Read10X_h5("hb_filtered_feature_bc_matrix.h5")
wb <-Read10X_h5("wb_filtered_feature_bc_matrix.h5")
b1 <- Read10X_h5("brain_1_transcriptomic_filtered_feature_bc_matrix.h5")
b2 <- Read10X_h5("brain_2_transcriptomic_filtered_feature_bc_matrix.h5")


w3.merge <- merge(x=fb, y= c(mb,hb,wb,b1,b2), add.cell.ids = c("fb","mb","hb","b1","b2"), project = "week.3.brain"  )
sobj_w3 <- CreateSeuratObject(w3.merge, min.cells = 3, min.features = 300)

saveRDS(sobj_w3, "sobj_w3.rds")


#!/bin/bash
#BSUB -J scATAC_Sample_Data  # LSF job name
#BSUB -o Sample_Fastaq.%J.out     # Name of the job output file 
#BSUB -e Sample_Fastaq.%J.error   # Name of the job error file
#BSUB -n 8
#BSUB -M 8000 #10GB
#BSUB -R "span[hosts=1] rusage [mem=8000]" 

module load R/4.2
# Load Anaconda module (modify the module name if necessary)
module load anaconda

# Create a conda environment
conda create -n seurat_env

# Activate the conda environment
conda activate seurat_env

# Install Seurat package
conda install -c bioconda r-seurat

Rscript /home/roshanpe/Seurat/seurat_integration.R