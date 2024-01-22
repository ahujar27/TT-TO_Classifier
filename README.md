# TT/TO Classifier Instructions

There are 4 different steps required to run the Tumor Type/Tissue of Origin Classifier. Each one is contained in a separate folder with a separate Dockerfile that contains everything necessary to run the particular script. The SNV and CNV calls are created using the GATK best practices workflows. 

#NEED TO WRITE THIS SUCH THAT U RUN THE SCRIPTS TO PREPARE WDLs and THEN RUN ALL THE WDLs? dockers are only needed to actually run the WDL?


## 1. SNV/CNV
*All files used in the SNV step can be found in the SNV/CNV folder*

### 1a. Pre-Alignment
The pre-alignment step converts either an aligned BAM or a paired-end FASTQ file and converts it to an unmapped BAM to be aligned. 

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. FASTQ or BAM (Input 1 if the file is a FASTQ, 2 if the file is a BAM)

If the file provided is a FASTQ file, it can be converted to a BAM file using the `paired-fastq-to-unmapped-bam.wdl` workflow. The JSON to run this particular WDL can be created using the python script `make.paired-fastq-to-unmapped-bam.input.json.py`.

The python script can be ran in the following way. `python `

If the file provided is an aligned BAM file, the BAM file needs to be unmapped so it can be realigned using the appropriate files. This can be done using the `bam-to-unmapped-bams.wdl`. The JSON to run this particular WDL can be created using the python script `make.bam-to-unmapped-bams.input.json.py`.

If the provided file is an unmapped BAM file this step can be skipped. 

### 1b. Alignment
The alignment step is responsible for aligning the BAM/FASTQ file to the reference genome. The reference available in the docker container is GRCh38. The alignment step requires X arguments.

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)

To align the BAM, first generate the JSON for the WDL using the python script `make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py`. Then, the WDL, `processing-for-variant-discovery-gatk4.wdl` can be run using the inputs. 

### 1c. SNV Calling
This step []

1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 1b)
3. TUMORBAMIND (aligned tumor index BAI from Step 1b)
4. NORMALBAM (aligned normal BAM from Step 1b)
5. NORMALBAMIND (aligned normal index BAI from Step 1b)

To generate SNV the calls, first run the python script `make.mutect2.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `mutect2.wdl` can be run using the JSON.

### 1d. CNV Calling

1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 1b)
3. TUMORBAMIND (aligned tumor index BAI from Step 1b)
4. NORMALBAM (aligned normal BAM from Step 1b)
5. NORMALBAMIND (aligned normal index BAI from Step 1b)

To generate CNV calls, first run the python script `make.cnv_somatic_pair_workflow.hg38.input.json.py` to create the JSON necessayr for the WDL run. Then the WDL, `cnv_somatic_pair_workflow.wdl` can be run using the JSON.

## 2. MAF Conversion

## 3. Feature Generation

## 4. Model

