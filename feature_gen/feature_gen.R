library(maftools)
library(ggplot2)
library(dplyr)
require(reshape2)
library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
#

args = commandArgs(trailingOnly=TRUE)

if ("-t" %in% args) {
  indir = args[2]
  sample_list = args[3]
  geneList_path = args[4]
  gene_bins = args[5]
  outdir = args[6]
} else {
  indir = args[1]
  sample_list = args[2]
  geneList_path = args[3]
  gene_bins = args[4]
  cnv_info = args[5]
  outdir = args[6]
}

fullGeneList = read.csv(geneList_path, header = F)
#fullGeneList = read.csv('/Volumes/jin.zhang/Active/rohil/projects/TT_TO/TT-TO_Classifier/feature_gen/fullGeneList.txt', header = F)
#indir = '/Volumes/jin.zhang/Active/rohil/projects/TT_TO/maf/output_to/'
somatic_table <- data.frame(matrix(ncol = 29, nrow = 0))

# Compile All MAF files together for all patients in Cohort
tmr_list = c('COAD', 'READ', 'LIHC', 'STAD', 'ESCA', 'PAAD')
file_list <- list()
i = 1
for (tumor_type in tmr_list) {
  subdirectories <- list.dirs(paste(indir, tumor_type, sep = ''), full.names = FALSE, recursive = FALSE)
  for (subdir in subdirectories) {
    maf_path <- paste(indir, tumor_type, "/", subdir, "/snv/", subdir, ".filtered.maf", sep = '')
    file_list[[i]] <- maf_path
    i = i + 1
  }
}

maf_list = lapply(file_list, function(x) data.table::fread(x))
remove_ind = list()
j = 1
for (i in 1:length(maf_list)) {
  if (!('Hugo_Symbol' %in% colnames(maf_list[[i]]))) {
    remove_ind[[j]] <- i
    j = j+1
  }
}
 
for (i in 1:length(remove_ind)) {
    removing = remove_ind[[i]] - (i-1)
    maf_list <- maf_list[-removing]
}

for (i in 1:length(maf_list)) {
  maf_list[[i]][, 'Tumor_Sample_Barcode' := gsub(".*/([^/]+)/snv/.*", "\\1", file_list[[i]])]
}

#Bind all MAF together and read final compiled maf file

maf = rbindlist(maf_list)
mafObj = read.maf(maf)

maf <- subset(maf, Variant_Classification %in% unique(mafObj@data$Variant_Classification))

maf <- subset(maf, Hugo_Symbol %in% fullGeneList$V1)
maf <- subset(maf, FILTER == "PASS")
class(as.data.frame(maf))
somatic_table <-rbind(somatic_table,maf)

#somatic_table <- somatic_table[somatic_table$t_alt_count>5,]
mut_mtx = dcast(somatic_table[ ,c('Tumor_Sample_Barcode', 'Hugo_Symbol')], Tumor_Sample_Barcode ~ Hugo_Symbol, value.var="Hugo_Symbol")
mut_mtx <- mut_mtx[!is.na(mut_mtx$Tumor_Sample_Barcode),]
rownames(mut_mtx) <- mut_mtx$Tumor_Sample_Barcode
mut_mtx$Tumor_Sample_Barcode <- NULL
#mut_mtx <- mut_mtx[rowSums(mut_mtx)>0,colSums(mut_mtx[])>0]
dim(mut_mtx)

#### Generate CNV feature matrix

# mut_mtx<-mut_mtx[!apply(mut_mtx, 1, function(row) all(row == 0)),]
# mut_mtx<-mut_mtx[,!apply(mut_mtx, 2, function(col) all(col == 0))]

if (!("-t" %in% args)) {
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
  
  sample_use <- unique(intersect(row.names(cpy_mtx),row.names(mut_mtx))) #NEED TO FILL THIS
  cpy_df <- cpy_mtx[sample_use,]
  somatic_df <- mut_mtx[sample_use,]
  write.csv(cpy_df, paste(outdir, "CNV_all_tcga.csv", sep = ""))
  write.csv(somatic_df, paste(outdir, "somatic_all_tcga.csv", sep = ""))
} else {
  sample_use <- row.names(mut_mtx)
  write.csv(mut_mtx, paste(outdir, "somatic_all_tcga.csv", sep = ""))
}
#Rscript feature_gen.R -t /opt/Active/rohil/projects/TT_TO/maf/output_to/ sample_names.txt fullGeneList.txt genome_bins.txt /opt/Active/rohil/projects/TT_TO/features/
# labels <- somatic_table[,c("sampleId","studyId")]
# labels <- labels[!duplicated(labels$sampleId), ]
# write.csv(labels, paste(outdir, "labels_all_tcga.csv", sep = ""))
                         
## compute mutation RMD for each sample
somatic_table$Chromosome <- sub("^chr", "", somatic_table$Chromosome)
df_cols = c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Tumor_Sample_Barcode","Hugo_Symbol")
df <- somatic_table[,..df_cols]
colnames(df) <- c("chrom","pos","ref","alt","mut_type","Sample","Gene")
sample_df <- df %>% group_by(Sample)

## import annotated bins genome_bin.txt
#genome_bins <- read.csv("/Volumes/jin.zhang/Active/rohil/projects/TT_TO/TT-TO_Classifier/feature_gen/genome_bins.txt", sep = '\t')
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
  if (dim(tmp_df)[1] == 0) {
    next
  }
  tab <- table(tmp_df$bin_name)
  counts[names(tab)] <- tab
  
  ### normalize
  counts <- counts/sum(counts)
  counts[is.na(counts)] <- 0
  
  out <- as.matrix(counts)
  colnames(out) <- samp
  
  mt <- cbind(out,mt)
  
}

mt = as.data.frame(mt)
# mt <-mt[!apply(mt, 1, function(row) all(row == 0)),]
# mt <-mt[,!apply(mt, 2, function(col) all(col == 0))]
write.csv(mt,paste(outdir, "rmd_all_tcga.csv", sep = ""))

### generate SBS matrix
df_SBS_cols = c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
df_SBS <- somatic_table[,..df_SBS_cols]
colnames(df_SBS) <- c('chrom','pos','ref','alt','Sample')
df_SBS$chrom[df_SBS$chrom=="23"] <- "X"
sample_df_SBS <- df_SBS %>% group_by(Sample)

## create empty matrix to save SBS 
mt <- matrix(0,nrow = 96, ncol = 0)

for (samp in sample_use){
  
  tmp_df <- sample_df_SBS[sample_df_SBS$Sample==samp,]
  tmp_df <- tmp_df[,c('chrom','pos','ref','alt')]
  tmp_sbs <- mutSigExtractor::extractSigsSnv(df=tmp_df,  output='contexts', verbose=F, ref.genome=BSgenome.Hsapiens.UCSC.hg38)[,1] # NEED TO DEBUG THIS
  out_sbs <- as.matrix(tmp_sbs)
  colnames(out_sbs) <- samp
  mt <- cbind(out_sbs,mt)
}

mt_SBS = as.data.frame(mt)
# mt_SBS <- mt_SBS[!apply(mt_SBS, 1, function(row) all(row == 0)),]
# mt_SBS <- mt_SBS[,!apply(mt_SBS, 2, function(col) all(col == 0))]
write.csv(mt_SBS,paste(outdir, "sbs_all_tcga.csv", sep = ""))

#RNA Section


