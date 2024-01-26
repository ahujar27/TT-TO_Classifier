library(maftools)
library(ggplot2)
library(dplyr)
require(reshape2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

indir = args[1]
sample_list = args[2]
geneList_path = args[3]
gene_bins = args[4]
cnv_info = args[5]

fullGeneList = read.csv(geneList_path, header = F)

somatic_table <- data.frame(matrix(ncol = 29, nrow = 0))

patient_list = readLines(sample_list)

for (patient_id in patient_list) {
  maf_path <- paste(indir, patient_id, "/",patient_id, ".maf", sep = '')
  tmp = read.maf(maf = maf_path)
  tmp@data <- subset(tmp@data, Hugo_Symbol %in% fullGeneList$V1)
  tmp@data <- subset(tmp@data, FILTER == "PASS")
  
  somatic_table <-rbind(somatic_table,tmp@data)
}

somatic_table <- somatic_table[somatic_table$t_alt_count>5,]
mut_mtx = dcast(somatic_table[ ,c('Tumor_Sample_Barcode', 'Hugo_Symbol')], Tumor_Sample_Barcode ~ Hugo_Symbol, value.var="Hugo_Symbol")
mut_mtx <- mut_mtx[!is.na(mut_mtx$Tumor_Sample_Barcode),]
rownames(mut_mtx) <- mut_mtx$Tumor_Sample_Barcode
mut_mtx$Tumor_Sample_Barcode <- NULL
mut_mtx <- mut_mtx[rowSums(mut_mtx[])>0,colSums(mut_mtx[])>0]
dim(mut_mtx)

#### Generate CNV feature matrix

CNV_table <- data.frame(matrix(ncol = 10, nrow = 0))

CNV_table1 <- read.csv(cnv_info, sep = '\t')
CNV_table1 <- CNV_table1 %>%
  pivot_longer(cols = -(1:4), names_to = "Column", values_to = "Value")

### Generate copy number alteration matrix using gistic numbers
cpy_mtx = dcast(CNV_table1[ ,c('Column', 'Gene.Symbol', 'Value')], Column ~ Gene.Symbol, value.var="Value")
rownames(cpy_mtx) <- cpy_mtx$sampleId
cpy_mtx$sampleId <- NULL

## remove rows and columns containing only 0
cpy_mtx<-cpy_mtx[!apply(cpy_mtx, 1, function(row) all(row == 0)),]
cpy_mtx<-cpy_mtx[,!apply(cpy_mtx, 2, function(col) all(col == 0))]
mut_mtx<-mut_mtx[!apply(mut_mtx, 1, function(row) all(row == 0)),]
mut_mtx<-mut_mtx[,!apply(mut_mtx, 2, function(col) all(col == 0))]

sample_use <- unique(intersect(row.names(cpy_mtx),row.names(mut_mtx))) #NEED TO FILL THIS
cpy_df <- cpy_mtx[sample_use,]
somatic_df <- mut_mtx[sample_use,]

## compute mutation RMD for each sample
somatic_table$Chromosome <- sub("^chr", "", somatic_table$Chromosome)
df_cols = c("Chromosome","Start_Position","Reference_Allele","Variant_Allele","Variant_Classification","Tumor_Sample_Barcode","Hugo_Symbol")
df <- somatic_table[,df_cols]
colnames(df) <- c("chrom","pos","ref","alt","mut_type","Sample","Gene")
sample_df <- df %>% group_by(Sample)

## import annotated bins genome_bin.txt
genome_bins <- read.csv("genome_bins.txt", sep='\t')

genome_bins$chrom[genome_bins$chrom=="X"] <- "23"
mt <- matrix(0,nrow = 3071, ncol = 0)
cnt = 0

## create count matrix to save RMDs
counts <- structure(
  rep(0, nrow(genome_bins)), 
  names=naturalsort::naturalsort(unique(genome_bins$bin_name))
)

## use shared samples in Mut and CNV table to generate RMD matrix
for (samp in sample_use){
  
  tmp_df <- sample_df[sample_df$Sample==samp,]
  tmp_df <- tmp_df[tmp_df$chrom %in% genome_bins$chrom,]
  uniq_chroms <- naturalsort::naturalsort(unique(genome_bins$chrom))
  genome_bins$chrom <- factor(genome_bins$chrom, uniq_chroms)
  tmp_df$chrom <- factor(tmp_df$chrom, uniq_chroms)
  
  max_pos_digits <- nchar(max(c(genome_bins$start, genome_bins$end, tmp_df$pos)))
  
  chromPosToInt <- function(chrom, position, n.digits=max_pos_digits){
    pos_char <- formatC(position, width=n.digits, format='d', flag='0')
    chrom_int <- as.integer(chrom)
    as.numeric(paste0(chrom_int, pos_char))
  }
  
  genome_bins$start_int <- chromPosToInt(genome_bins$chrom, genome_bins$start)
  genome_bins$end_int <- chromPosToInt(genome_bins$chrom, genome_bins$end)
  tmp_df$pos_int <- chromPosToInt(tmp_df$chrom, tmp_df$pos)
  
  matched_bins <- lapply(tmp_df$pos_int, function(i){
    genome_bins$bin_name[i>=genome_bins$start_int & i<=genome_bins$end_int ]
  })
  
  match_lengths <- sapply(matched_bins, length)
  n_warn_matches <- sum(match_lengths>1 | match_lengths==0)
  matched_bins[match_lengths==0] <- NA
  tmp_df$bin_name <- sapply(matched_bins,`[[`,1)
  tab <- table(tmp_df$bin_name[!is.na(tmp_df$bin_name)])
  counts[names(tab)] <- tab
  
  ### normalize
  counts <- counts/sum(counts)
  counts[is.na(counts)] <- 0
  
  out <- as.matrix(counts)
  colnames(out) <- samp
  
  mt <- cbind(out,mt)
  
}

mt = as.data.frame(mt)
mt <-mt[!apply(mt, 1, function(row) all(row == 0)),]
mt <-mt[,!apply(mt, 2, function(col) all(col == 0))]
write.csv(mt,"rmd_all_tcga.csv")

### generate SBS matrix
df_SBS_cols = c("Chromosome","Start_Position","Reference_Allele","Varinat_Allele","Tumor_Sample_Barcode")
df_SBS <- somatic_table[,df_SBS_cols]
colnames(df_SBS) <- c('chrom','pos','ref','alt','Sample')
df_SBS$chrom[df_SBS$chrom=="23"] <- "X"
sample_df_SBS <- df_SBS %>% group_by(Sample)

## create empty matrix to save SBS 
mt <- matrix(0,nrow = 96, ncol = 0)

for (samp in sample_use){
  
  tmp_df <- sample_df_SBS[sample_df_SBS$Sample==samp,]
  tmp_df <- tmp_df[,c('chrom','pos','ref','alt')]
  tmp_sbs <- mutSigExtractor::extractSigsSnv(df=tmp_df,  output='contexts', verbose=F)[,1]
  out_sbs <- as.matrix(tmp_sbs)
  colnames(out_sbs) <- samp
  mt <- cbind(out_sbs,mt)
}

mt_SBS = as.data.frame(mt_SBS)
mt_SBS <- mt_SBS[!apply(mt_SBS, 1, function(row) all(row == 0)),]
mt_SBS <- mt_SBS[,!apply(mt_SBS, 2, function(col) all(col == 0))]
write.csv(mt_SBS,"sbs_all_tcga.csv")
