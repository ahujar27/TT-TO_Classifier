#!/bin/bash

fastq_bam_WDL=https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/paired-fastq-to-unmapped-bam.wdl
py_dir=$1
SRCBAM=$2
INDIR=$3
SAMPNAME=$4
OUTDIR=$5

#Fastq to BAM
python3 $py_dir/make.paired-fastq-to-unmapped-bam-input.json.py $SRCBAM $INDIR/$SAMPNAME $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json

/usr/bin/java -Dconfig.file=$INDIR/cromwell.config.$SAMPNAME -jar /opt/cromwell.jar run $fastq_bam_WDL -i $OUTDIR/$SAMPNAME.bam-to-unmapped-bams.inputs.json

#Once Fastq to BAM conversions are all complete, create the BAM list as outlined in the github 
#Alignment of BAM

BAMLIST=$6
alignment_WDL=https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/master/processing-for-variant-discovery-gatk4.wdl

python $py_dir/make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py $SRCBAM $BAMLIST $OUTDIR/$SAMPNAME.processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json

/usr/bin/java -Dconfig.file=$INDIR/cromwell.config.$SAMPNAME -jar /opt/cromwell.jar run $alignment_WDL -i $OUTDIR/$SAMPNAME.processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json

#SNV Calling Step
#Tumor Bam, and Bam Index are generated from the previous step.
TUMORBAM=$7
TUMORBAMIND=$8
mutect2_wdl=https://raw.githubusercontent.com/broadinstitute/gatk/141529b167fb73768314a3caccbbac3c5021741a/scripts/mutect2_wdl/mutect2.wdl

python $py_dir/make.mutect2.input.json.py $TUMORBAM $TUMORBAMIND $OUTDIR/$SAMPNAME.mutect2.inputs.json