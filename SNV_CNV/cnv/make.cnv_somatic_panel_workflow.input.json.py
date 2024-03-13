import json
import sys

normalList = sys.argv[1]
outName = sys.argv[2]
intervalPath = sys.argv[3]

with open(normalList,'r') as fin:
    normBams = fin.readlines()
normBams = [fname[:-1] for fname in normBams]
normBais = [fname[:-1]+'i' for fname in normBams]

jsonDict = {
    "CNVSomaticPanelWorkflow.normal_bams": normBams,
    "CNVSomaticPanelWorkflow.normal_bais": normBais,
    "CNVSomaticPanelWorkflow.pon_entity_id": "tuo-cnv-pon",
  
    "CNVSomaticPanelWorkflow.ref_fasta": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
    "CNVSomaticPanelWorkflow.ref_fasta_fai": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    "CNVSomaticPanelWorkflow.ref_fasta_dict": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict",
    "CNVSomaticPanelWorkflow.intervals": intervalPath,

    "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "CNVSomaticPanelWorkflow.gatk4_jar_override": "/opt/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar".format(appPath),
    "CNVSomaticPanelWorkflow.preemptible_attempts": "3"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
