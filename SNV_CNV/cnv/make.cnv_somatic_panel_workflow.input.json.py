import json
import sys

normalList = sys.argv[1]
outName = sys.argv[2]
refPath = sys.argv[3]
intervalPath = sys.argv[4]
appPath = sys.argv[5]

with open(normalList,'r') as fin:
    normBams = fin.readlines()
normBams = [fname[:-1] for fname in normBams]
normBais = [fname[:-1]+'i' for fname in normBams]

jsonDict = {
    "CNVSomaticPanelWorkflow.normal_bams": normBams,
    "CNVSomaticPanelWorkflow.normal_bais": normBais,
    "CNVSomaticPanelWorkflow.pon_entity_id": "tuo-cnv-pon",
  
    "CNVSomaticPanelWorkflow.ref_fasta": "{}/GRCh38.d1.vd1.fa".format(refPath),
    "CNVSomaticPanelWorkflow.ref_fasta_fai": "{}/GRCh38.d1.vd1.fa.fai".format(refPath),
    "CNVSomaticPanelWorkflow.ref_fasta_dict": "{}/GRCh38.d1.vd1.fa.dict".format(refPath),
    "CNVSomaticPanelWorkflow.intervals": intervalPath,

    "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.1.4.0",
    "CNVSomaticPanelWorkflow.gatk4_jar_override": "{}/gatk-package-4.1.4.1-local.jar".format(appPath),
    "CNVSomaticPanelWorkflow.preemptible_attempts": "3"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
