# TT/TO Classifier Instructions

There are 4 different steps required to run the Tumor Type/Tissue of Origin Classifier. Each one is contained in a separate folder with a separate Dockerfile that contains everything necessary to run the particular script. The SNV and CNV calls are created using the GATK best practices workflows. 

#NEED TO WRITE THIS SUCH THAT U RUN THE SCRIPTS TO PREPARE WDLs and THEN RUN ALL THE WDLs? dockers are only needed to actually run the WDL?


## 1. SNV
### 1a. Alignment
The alignment step is responsible for aligning the BAM/FASTQ file to the reference genome. The reference available in the docker container is GRCh38. The alignment step requires X arguments.

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. FASTQ or BAM (Input 1 if the file is a FASTQ, 2 if the file is a BAM)

