# TT/TO Classifier Instructions

There are 4 different steps required to run the Tumor Type/Tissue of Origin Classifier. Each one is contained in a separate folder with a separate Dockerfile that contains everything necessary to run the particular script. The SNV and CNV calls are created using the GATK best practices workflows, for which documentation can be found at https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows

## 1. SNV/CNV
*All files used in the SNV step can be found in the SNV/CNV folder*

### 1a. Build the Docker

The Dockerfile is contained in the `/SNV_CNV/` folder. Build this docker container to run the SNV and the CNV portion of the workflow. Once the docker image is built, activate the conda environment to properly run the workflows: `conda activate gatk`.

### 1b. Pre-Alignment
The pre-alignment step converts either an aligned BAM or a paired-end FASTQ file and converts it to an unmapped BAM to be aligned. 

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. OUTDIR (Output Directory)

If the file provided is a FASTQ file, it can be converted to a BAM file using the `paired-fastq-to-unmapped-bam.wdl` workflow. The JSON to run this particular WDL can be created using the python script `make.paired-fastq-to-unmapped-bam.input.json.py`. The WDL can be found at this link: https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/paired-fastq-to-unmapped-bam.wdl

The python script can be ran in the following way. `python make.paired-fastq-to-unmapped-bam.input.json.py $SRCBAM $OUTDIR/$SAMPNMAME/ $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json`

If the file provided is an aligned BAM file, the BAM file needs to be unmapped so it can be realigned using the appropriate files. This can be done using the `bam-to-unmapped-bams.wdl`. The JSON to run this particular WDL can be created using the python script `make.bam-to-unmapped-bams.input.json.py`. THe WDL can be found at this link: https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/bam-to-unmapped-bams.wdl

The python script can be ran in the following way. `python make.bam-to-unmapped-bams.input.json.py $SRCBAM $OUTDIR/$SAMPNMAME/ $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json`

If the provided file is an unmapped BAM file this step can be skipped. 

### 1c. Alignment
The alignment step is responsible for aligning the BAM/FASTQ file to the reference genome. The reference available in the docker container is GRCh38. The alignment step requires X arguments.

1. SAMPNAME (Sample Name)
2. SRCBAM (path to the input BAM file, e.g. a tumor or a normal WXS BAM downloaded from TCGA)
3. OUTDIR (Output Directory)
4. BAMLIST (Path to list of BAM files)

The format of the BAM list can be found under /SNV_CNV/bam_list.txt

To align the BAM, first generate the JSON for the WDL using the python script `make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py`. Then, the WDL, `processing-for-variant-discovery-gatk4.wdl` can be run using the inputs. 

The python script can be ran in the following way. `python make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py $SRCBAM $BAMLIST $OUTDIR/$SAMPNAME.processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json`

The WDL for this particular step can be found at: https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/master/processing-for-variant-discovery-gatk4.wdl

### 1d. SNV Calling
This step calls the variants. If a normal and a tumor are both present follow the tumor-normal SNV calling. If a normal is not present, then follow the steps for tumor-only calling. 

#### Tumor-Normal
1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 1b)
3. TUMORBAMIND (aligned tumor index BAI from Step 1b)
4. NORMALBAM (aligned normal BAM from Step 1b)
5. NORMALBAMIND (aligned normal index BAI from Step 1b)
6. OUTDIR (Output Directory)

To generate SNV the calls, first run the python script `make.mutect2.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `mutect2.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.mutect2.input.json.py $TUMORBAM $TUMORBAMIND $NORMALBAM $NORMALBAMIND $OUTDIR/$SAMPNAME.mutect2.inputs.json`

#### Tumor-only

1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 1b)
3. TUMORBAMIND (aligned tumor index BAI from Step 1b)

To generate SNV the calls, first run the python script `make.mutect2.tumor_only.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `mutect2.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.mutect2.tumor_only.input.json.py $TUMORBAM $TUMORBAMIND $OUTDIR/$SAMPNAME.mutect2.tumor_only.inputs.json`

The WDL for the tumor-normal and tumor-only are the same, and can be found here: https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mutect2_wdl/mutect2.wdl

### 1e. CNV Panel of Normals Generation

1. NORMALBAMLIST (List of normal bams to generate the panel of normals)
2. OUTFILE (The name of the output JSON file)

To generate CNV calls, a panel of normals must first be generated. To generate this, first run the python script `make.cnv_somatic_panel_workflow.input.json.py` to create the JSON necessary for the WDL run. Then the WDL, `cnv_somatic_panel_workflow.wdl` can be run using the JSON.

The python script can be ran in the following way. `python make.cnv_somatic_panel_workflow.input.json.py $NORMALBAMLIST $OUTFILE`

The WDL can be found at: https://github.com/gatk-workflows/gatk4-somatic-cnvs/raw/master/cnv_somatic_panel_workflow.wdl

### 1f. CNV Calling

1. SAMPNAME (sample name)
2. TUMORBAM (aligned tumor BAM from Step 1b)
3. TUMORBAMIND (aligned tumor index BAI from Step 1b)
4. NORMALBAM (aligned normal BAM from Step 1b)
5. NORMALBAMIND (aligned normal index BAI from Step 1b)
6. OUTFILE (The name of the output JSON file) 

Once the panel of normals is generated, the python script `make.cnv_somatic_pair_workflow.hg38.input.json.py` was ran to create the JSON necessary for the WDL run. Then the WDL, `cnv_somatic_pair_workflow.wdl` can be run using the JSON.

The python script can be ran in the following way: `python make.cnv_somatic_pair_workflow.hg38.input.json.py $TUMORBAM $TUMORBAMIND $NORMALBAM $NORMALBAMIND $OUTFILE`

The WDL can be found at: https://github.com/gatk-workflows/gatk4-somatic-cnvs/raw/master/cnv_somatic_pair_workflow.wdl

## 2. GISTIC2 + MAF Conversion

The Dockerfile for this portion is contained in the folder `/MAF/`. Build this docker container to run the SNV and the CNV portion of the workflow. Once the docker image is built, activate the conda environment to properly run the workflows: `conda activate gatk`.

1. INVCF (Final Processed VCF)
2. SAMPNAME (Sample Name)
3. OUTDIR (Output Directory)
   
In the docker container, to convert a processed VCF to MAF, run the following command: `perl /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $INVCF --output-maf $OUTDIR/$SAMPNAME.maf --ncbi-build GRCh38 --ref-fasta /opt/references/GRCh38.fa --vep-path /opt/miniconda/envs/gatk/bin`

## 3. Feature Generation

The Dockerfile for this portion is contained in the folder `/feature_gen/`. Additionally in the folder is the feature generation script, at `/feature_gen/feature_gen.R`. 

1. INDIR (Directory containing all SNV/CNV calls)
2. SAMPNAMES (List of sample names from previous steps)
3. GENEPATH (Gene List Path)
4. BINPATH (Genome Bins Path)
5. CNVPATH (all_thresholded.by_genes.txt Path)

Once the docker environment is run, the feature generation script can be run as such: `RSCRIPT feature_gen.R $INDIR $SAMPNAMES $GENEPATH $BINPATH $CNVPATH`

## 4. Model


