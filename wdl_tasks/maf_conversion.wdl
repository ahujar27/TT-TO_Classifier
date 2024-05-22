version 1.0

# Step 3. [MAF Conversion](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L121-L129)
#
# The Dockerfile for this portion is contained in the folder `/MAF/`. Build this docker container to run the SNV and the CNV portion of the workflow. Once the docker image is built, activate the conda environment to properly run the workflows: `conda activate gatk`.

# 1. INVCF (Final Processed VCF)
# 2. SAMPNAME (Sample Name)
# 3. OUTDIR (Output Directory)

# In the docker container, to convert a processed VCF to MAF, run the following command: `perl /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $INVCF --output-maf $OUTDIR/$SAMPNAME.maf --ncbi-build GRCh38 --ref-fasta /opt/references/GRCh38.fa --vep-path /opt/miniconda/envs/gatk/bin`

task MafConversion {
  input {
    File inputVCF
    File ref_fasta

    String sampleName
    String vep_path
    String dockerPath
  }

  command <<<

    set -e

    perl /opt/mskcc-vcf2maf-754d68a/vcf2maf.pl \
      --input-vcf "${inputVCF}" \
      --output-maf "${dockerPath}.maf" \
      --ncbi-build GRCh38 \
      --ref-fasta "${ref_fasta}" \
      --vep-path "${vep_path}"
  >>>

  runtime {
    docker: "${dockerPath}"
  }

  output {
    File mafFile = "${sampleName}.maf"
  }
}