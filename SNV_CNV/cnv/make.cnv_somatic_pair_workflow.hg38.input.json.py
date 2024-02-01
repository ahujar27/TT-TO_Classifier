import json
import sys

tumorBam = sys.argv[1]
tumorBamInd = sys.argv[2]
normalBam = sys.argv[3]
normalBamInd = sys.argv[4]
outName = sys.argv[5]
intervalsPath = sys.argv[6]
refPath = sys.argv[7]
ponPath = sys.argv[8]

appPath = /opt/gatk-4.4.0.0/

jsonDict = {
  "CNVSomaticPairWorkflow.tumor_bam": tumorBam,
  "CNVSomaticPairWorkflow.tumor_bam_idx": tumorBamInd,
  "CNVSomaticPairWorkflow.normal_bam": normalBam,
  "CNVSomaticPairWorkflow.normal_bam_idx": normalBamInd,
  
  "CNVSomaticPairWorkflow.ref_fasta": "{}/GRCh38.d1.vd1.fa".format(refPath),
  "CNVSomaticPairWorkflow.ref_fasta_fai": "{}/GRCh38.d1.vd1.fa.fai".format(refPath),
  "CNVSomaticPairWorkflow.ref_fasta_dict": "{}/GRCh38.d1.vd1.fa.dict".format(refPath),
  "CNVSomaticPairWorkflow.common_sites": "gs://gatk-test-data/cnv/somatic/common_snps.interval_list",
  "CNVSomaticPairWorkflow.read_count_pon": ponPath, 
  "CNVSomaticPairWorkflow.intervals": intervalsPath,

  "CNVSomaticPairWorkflow.gatk_docker": "broadinstitute/gatk:4.4.0.0",
  "CNVSomaticPairWorkflow.gatk4_jar_override": "{}/gatk-package-4.4.0.0-local.jar".format(appPath),

  "CNVSomaticPairWorkflow.is_run_funcotator": "false"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
