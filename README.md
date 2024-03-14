# TT/TO Classifier Instructions

There are 4 different steps required to run the Tumor Type/Tissue of Origin Classifier. Each one is contained in a separate folder with a separate Dockerfile that contains everything necessary to run the particular script. The SNV and CNV calls are created using the GATK best practices workflows, for which documentation can be found at https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows

## 1. Download the References

The reference files for GRCh38 are required in each individual step of the pre-processing to run the classifier. To download the references, first, gsutil needs to be installed. This can be done from the following link: https://cloud.google.com/storage/docs/gsutil_install, selecting your specific operating system.

From there, run the commands:

`mkdir references` and then `cd references`. 

The references can then be downloaded as such:

```
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
```
Following the download of the references, download the files necessary for SNV calling:

```
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi
gsutil cp gs://gatk-test-data/mutect2/Homo_sapiens_assembly38.index_bundle
```

## 2. SNV/CNV
*All files used in the SNV step can be found in the SNV/CNV folder*

### 2a. Build the Docker

The Dockerfile is contained in the `/SNV_CNV/` folder. Build this docker container to run the SNV and the CNV portion of the workflow. Once the docker image is built, activate the conda environment to properly run the workflows: `conda activate gatk`.

### 2b. Pre-Alignment
The pre-alignment step converts either an aligned BAM or a paired-end FASTQ file and converts it to an unmapped BAM to be aligned. 

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. OUTDIR (Output Directory)

If the file provided is a FASTQ file, it can be converted to a BAM file using the `paired-fastq-to-unmapped-bam.wdl` workflow. The JSON to run this particular WDL can be created using the python script `make.paired-fastq-to-unmapped-bam.input.json.py`. The WDL can be found at this link: https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/paired-fastq-to-unmapped-bam.wdl

The python script can be ran in the following way. `python make.paired-fastq-to-unmapped-bam.input.json.py $SRCBAM $OUTDIR/$SAMPNMAME/ $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json`

If the file provided is an aligned BAM file, the BAM file needs to be unmapped so it can be realigned using the appropriate files. This can be done using the `bam-to-unmapped-bams.wdl`. The JSON to run this particular WDL can be created using the python script `make.bam-to-unmapped-bams.input.json.py`. THe WDL can be found at this link: https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/bam-to-unmapped-bams.wdl

The python script can be ran in the following way. `python make.bam-to-unmapped-bams.input.json.py $SRCBAM $OUTDIR/$SAMPNMAME/ $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json`

If the provided file is an unmapped BAM file this step can be skipped. 

### 2c. Alignment
The alignment step is responsible for aligning the BAM/FASTQ file to the reference genome. The reference available in the docker container is GRCh38. The alignment step requires X arguments.

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. OUTDIR (Output Directory)
4. BAMLIST (Path to list of BAM files)

The format of the BAM list can be found under /SNV_CNV/bam_list.txt

To align the BAM, first generate the JSON for the WDL using the python script `make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py`. Then, the WDL, `processing-for-variant-discovery-gatk4.wdl` can be run using the inputs. 

The python script can be ran in the following way. `python make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py $SRCBAM $BAMLIST $OUTDIR/$SAMPNAME.processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json`

The WDL for this particular step can be found at: https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/master/processing-for-variant-discovery-gatk4.wdl

### 2d. SNV Calling
This step calls the variants. If a normal and a tumor are both present follow the tumor-normal SNV calling. If a normal is not present, then follow the steps for tumor-only calling. 

#### Tumor-Normal
1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 2c)
3. TUMORBAMIND (aligned tumor index BAI from Step 2c)
4. NORMALBAM (aligned normal BAM from Step 2c)
5. NORMALBAMIND (aligned normal index BAI from Step 2c)
6. OUTDIR (Output Directory)
7. INTPATH (Path to TSO500 Intervals Bed)

To generate SNV the calls, first run the python script `make.mutect2.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `mutect2.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.mutect2.input.json.py $TUMORBAM $TUMORBAMIND $NORMALBAM $NORMALBAMIND $OUTDIR/$SAMPNAME.mutect2.inputs.json`

#### Tumor-only

1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 2c)
3. TUMORBAMIND (aligned tumor index BAI from Step 2c)
4. INTPATH (Path to TSO500 Intervals Bed)

To generate SNV the calls, first run the python script `make.mutect2.tumor_only.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `mutect2.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.mutect2.tumor_only.input.json.py $TUMORBAM $TUMORBAMIND $OUTDIR/$SAMPNAME.mutect2.tumor_only.inputs.json`

The WDL for the tumor-normal and tumor-only are the same, and can be found here: https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mutect2_wdl/mutect2.wdl

### 2e. CNV Panel of Normals Generation

**If you are running the pipeline in Tumor-Only mode, please skip both CNV steps (2e and 2f)**

1. NORMALBAMLIST (List of normal bams to generate the panel of normals)
2. OUTFILE (The name of the output JSON file)
3. INTERVALSPATH (location of TSO500 intervals file (contained in repository))

To generate CNV calls, a panel of normals must first be generated. To generate this, first run the python script `make.cnv_somatic_panel_workflow.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `cnv_somatic_panel_workflow.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.cnv_somatic_panel_workflow.input.json.py $NORMALBAMLIST $OUTFILE $REFPATH $INTERVALSPATH`

The WDL can be found at: https://github.com/gatk-workflows/gatk4-somatic-cnvs/raw/master/cnv_somatic_panel_workflow.wdl

### 2f. CNV Calling

1. TUMORBAM (aligned tumor BAM from Step 2c)
2. TUMORBAMIND (aligned tumor index BAI from Step 2c)
3. NORMALBAM (aligned normal BAM from Step 2c)
4. NORMALBAMIND (aligned normal index BAI from Step 2c)
5. OUTFILE (The name of the output JSON file)
6. INTERVALSPATH (location of TSO500 intervals file (contained in repository))
7. PONPATH (Path to PON from previous step)

Once the panel of normals is generated, the python script `make.cnv_somatic_pair_workflow.hg38.input.json.py` was ran to create the JSON necessary for the WDL run. Then the WDL, `cnv_somatic_pair_workflow.wdl` can be run using the JSON.

The python script can be ran in the following way: `python make.cnv_somatic_pair_workflow.hg38.input.json.py $TUMORBAM $TUMORBAMIND $NORMALBAM $NORMALBAMIND $OUTFILE $INTERVALSPATH $REFPATH $PONPATH`

The WDL can be found at: https://github.com/gatk-workflows/gatk4-somatic-cnvs/raw/master/cnv_somatic_pair_workflow.wdl

## 3. MAF/GISTIC Conversion

The Dockerfile for this portion is contained in the folder `/MAF/`. Build this docker container to run the SNV and the CNV portion of the workflow. Once the docker image is built, activate the conda environment to properly run the workflows: `conda activate gatk`.

### 3a. MAF Conversion

1. INVCF (Final Processed VCF)
2. SAMPNAME (Sample Name)
3. OUTDIR (Output Directory)
   
In the docker container, to convert a processed VCF to MAF, run the following command: `perl /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $INVCF --output-maf $OUTDIR/$SAMPNAME.maf --ncbi-build GRCh38 --ref-fasta /opt/references/GRCh38.fa --vep-path /opt/miniconda/envs/gatk/bin`

### 3b. GISTIC Conversion 

1. INFILE (Combined CNV calls from step 4a)

Once this is complete, GISTIC2 can be run using the following command: `./gistic2 -refgene /refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -cnv $INFILE`

## 4. Feature Generation

The Dockerfile for this portion is contained in the folder `/feature_gen/`. Additionally in the folder is the feature generation script, at `/feature_gen/feature_gen.R`. 

1. INDIR (Directory containing all SNV/CNV calls)
2. SAMPNAMES (List of sample names from previous steps)
3. GENEPATH (Gene List Path)
4. BINPATH (Genome Bins Path)
5. CNVPATH (all_thresholded.by_genes.txt Path)
6. OUTDIR (output directory where features should be stored)

Once the docker environment is run, the feature generation script can be run as such: `RSCRIPT feature_gen.R $INDIR $SAMPNAMES $GENEPATH $BINPATH $CNVPATH $OUTDIR`

## 5. Model

The Dockerfile for this portion is contained in the folder `/model/`. Additionally in the folder is the xgboost model, in a JSON format, under `/model/final_model.json`.

1. MODELPATH (Path to model)
2. GENELIST (List of genes from previous feature generation step)
3. RMDPATH (Path to RMD from feature generation (from step 4))
4. CNVPATH (Path to CNV from feature generation (from step 4))
5. SNVPATH (Path to SNV from feature generation (from step 4))
6. LABELPATH (Path to labels from feature generation (from step 4))
7. SBSPATH (Path to SBS from feature generation (from step 4))
8. PASTPATH (Path to the `tcga_500.csv` file)
9. OUTDIR (Path to the output directory for predictions to be stored)

Once the docker container is created, the model can be ran using the command: `python model.py $MODELPATH $GENELIST $RMDPATH $CNVPATH $SNVPATH $LABELPATH $SBSPATH $PASTPATH $OUTDIR`


