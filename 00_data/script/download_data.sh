
############################################ Data set 1. Caspar Ohnmacht ############################################
# GEO study link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243483

# Sample: GSM7789315 control SI-LP steady state
# Illumina HiSeq 2500
wget -P ../ -O ../GSM7789315.barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789315&format=file&file=GSM7789315%5FMUC29369%5Fbarcodes%2Etsv%2Egz"
wget -P ../ -O ../GSM7789315.features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789315&format=file&file=GSM7789315%5FMUC29369%5Ffeatures%2Etsv%2Egz"
wget -P ../ -O ../GSM7789315.matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789315&format=file&file=GSM7789315%5FMUC29369%5Fmatrix%2Emtx%2Egz"

# Sample: GSM7789317 control mLN steady state
wget -P ../ -O ../GSM7789317.barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789317&format=file&file=GSM7789317%5FMUC29371%5Fbarcodes%2Etsv%2Egz"
wget -P ../ -O ../GSM7789317.features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789317&format=file&file=GSM7789317%5FMUC29371%5Ffeatures%2Etsv%2Egz"
wget -P ../ -O ../GSM7789317.matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7789317&format=file&file=GSM7789317%5FMUC29371%5Fmatrix%2Emtx%2Egz"

############################################ Data set 2. Vuk Cerovic ############################################
# GEO study link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283808

# Sample: GSM8672515	Small intestinal dendritic cells CCR7gfp SI LP DCs
# Illumina NextSeq 500
wget -P ../ -O ../GSM8672515.barcodes.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8672515&format=file&file=GSM8672515%5Fbarcodes%2Etsv%2Egz"
wget -P ../ -O ../GSM8672515.features.tsv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8672515&format=file&file=GSM8672515%5Ffeatures%2Etsv%2Egz"
wget -P ../ -O ../GSM8672515.matrix.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8672515&format=file&file=GSM8672515%5Fmatrix%2Emtx%2Egz"

# Sample: GSM9122899	Small intestinal dendritic cells from WT mice
# NextSeq 2000
wget -P ../ -O ../GSM9122899.feature.matrix.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9122899&format=file&file=GSM9122899%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

############################################ Data set 3. Vuk Cerovic ############################################
# GEO study link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160156

# Sample: GSM4861588 CD103+DC_gut_steadystate_R1
# Sample: GSM4861589 CD103+DC_gut_steadystate_R2
# Sample: GSM4861590 CD103+DC_gut_steadystate_R3
# Sample: GSM4861591 Double+DC_gut_steadystate_R1
# Sample: GSM4861592 Double+DC_gut_steadystate_R2
# Sample: GSM4861593 Double+DC_gut_steadystate_R3
# Sample: GSM4861585 CD11b+DC_gut_steadystate_R1
# Sample: GSM4861586 CD11b+DC_gut_steadystate_R2
# Sample: GSM4861587 CD11b+DC_gut_steadystate_R3

# Download all DC read counts
# Data should be filtered for the gut samples (exclude duct and node)
# gut:	Intestinal lamina propria	(LP)
# duct: Mesenteric lymphatic duct	(mL)
# node: Mesenteric lymph node	(mLN)
wget -P ../ -O ../GSE160156.read_counts.csv.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE160156&format=file&file=GSE160156%5Fread%5Fcounts%5FDC%2Ecsv%2Egz"

# ^^ THIS IS BULK DATA 

############################################ Data set 4. George Kollias ############################################
# GEO study link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255350

# GSE255350_myeloid_raw_counts.txt.gz 
wget -P ../ -O ../GSE255350_myeloid_raw_counts.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE255350&format=file&file=GSE255350%5Fmyeloid%5Fraw%5Fcounts%2Etxt%2Egz"

# Meta data
# 2 conditions. TNFDARE= Inflamed, WT= Wild type
wget -P ../ -O ../GSE255350_metadata.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE255350&format=file&file=GSE255350%5Fmetadata%5Ftable%2Etxt%2Egz"

############################################ Data set 5. Fiona Powrie ############################################
# Study link: https://treg-gut-niches.cellgeni.sanger.ac.uk/ 

# Myeloid. Contain multiple site such as Lamina Propria, MLN and others. 
# Two conditions: Helicobacter infection = Inflammation, and WT

# Myeloid - Processed data for scRNA-seq (Not raw counts)
# wget https://cellgeni.cog.sanger.ac.uk/treg-gut-niches/Myeloid_allgenes.h5ad

# Load h5ad (python scanpy object) in R:
# library(zellkonverter)
# adata <- readH5AD("Myeloid_allgenes.h5ad")
# adata

# Download fastq files. CellRanger would need to be done after. 
# Not sure which samples are what
# https://www.ebi.ac.uk/ena/browser/view/PRJEB57700?show=reads


# LP GEX files 
# wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR125/062/ERR12552062/ERR12552062_1.fastq.gz
# wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR125/061/ERR12552061/ERR12552061_2.fastq.gz
# wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR125/061/ERR12552061/ERR12552061_1.fastq.gz
# wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR125/062/ERR12552062/ERR12552062_2.fastq.gz

# Analyzed on computerome 
rsync -avzP helweg@transfer.computerome.dk:/home/projects/dtu_00062/people/helweg/LPcDC_cellranger/cellranger_analysis/outputs/LP_CRAM1/outs/raw_feature_bc_matrix ../CRAM1
rsync -avzP helweg@transfer.computerome.dk:/home/projects/dtu_00062/people/helweg/LPcDC_cellranger/cellranger_analysis/outputs/LP_CRAM2/outs/raw_feature_bc_matrix ../CRAM2

# 2025-10-15, email from Faidra: "The Fiona Powrie dataset will have to wait a bit since we might not be able to use it".

############################################ Data set 6. Rivera CA ############################################

# Study link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188379 

# GSM5678427 CD103+CD11b+ dendritic cells from lamina propria

wget -P ../ -O ../GSM5678427_LP_raw_gene_bc_matrices_h5.h5 "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5678427&format=file&file=GSM5678427%5FLP%5Fraw%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5" 

############################################ Colon data ############################################

# I am also going to attach some data from colon that we should keep in case we want to do an object for colon in the future.
# colon : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137927
# GSM4094548	SPF_scRNAseq
# GSM4094549	GF_scRNAseq


############################################ Data set 7. ############################################
# https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-9522?query=E-MTAB-9522
# wget -P ../ -O ../E_MTAB_9522_wt_S1_R1.fastq.gz https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/522/E-MTAB-9522/Files/wt_S1_R1.fastq.gz
# wget -P ../ -O ../E_MTAB_9522_wt_S1_R2.fastq.gz https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/522/E-MTAB-9522/Files/wt_S1_R2.fastq.gz
# wget -P ../ -O ../E_MTAB_9522_wt_S1_I1.fastq.gz https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/522/E-MTAB-9522/Files/wt_S1_I1.fastq.gz

# preprocessed in computerome 
rsync -avzP helweg@transfer.computerome.dk:/home/people/helweg/ciir/people/helweg/LPcDC_cellranger/cellranger_analysis_7/outputs/wt/outs/filtered_feature_bc_matrix ../E_MTAB_9522



