import json
import sys

tumorBam = sys.argv[1]
tumorBamInd = sys.argv[2]
normalBam = sys.argv[3]
normalBamInd = sys.argv[4]
outName = sys.argv[5]
intervalsPath = sys.argv[6]
ponPath = sys.argv[7]

jsonDict = {
  "CNVSomaticPairWorkflow.tumor_bam": tumorBam,
  "CNVSomaticPairWorkflow.tumor_bam_idx": tumorBamInd,
  "CNVSomaticPairWorkflow.normal_bam": normalBam,
  "CNVSomaticPairWorkflow.normal_bam_idx": normalBamInd,
  
  "CNVSomaticPairWorkflow.ref_fasta": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CNVSomaticPairWorkflow.ref_fasta_fai": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CNVSomaticPairWorkflow.ref_fasta_dict": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict",
  "CNVSomaticPairWorkflow.common_sites": "gs://gatk-test-data/cnv/somatic/common_snps.interval_list",
  "CNVSomaticPairWorkflow.read_count_pon": ponPath, 
  "CNVSomaticPairWorkflow.intervals": intervalsPath,

  "CNVSomaticPairWorkflow.gatk_docker": "broadinstitute/gatk:4.1.4.1",
  "CNVSomaticPairWorkflow.gatk4_jar_override": "/opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar",

  "CNVSomaticPairWorkflow.is_run_funcotator": "false"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
