import json
import sys

normalList = sys.argv[1]
outName = sys.argv[2]

with open(normalList,'r') as fin:
    normBams = fin.readlines()
normBams = [fname[:-1] for fname in normBams]
normBais = [fname[:-1]+'i' for fname in normBams]

jsonDict = {
    "CNVSomaticPanelWorkflow.normal_bams": normBams,
    "CNVSomaticPanelWorkflow.normal_bais": normBais,
    "CNVSomaticPanelWorkflow.pon_entity_id": "tuo-cnv-pon",
  
    "CNVSomaticPanelWorkflow.ref_fasta": "/storage1/fs1/jin.zhang/Active/tmp/rohil_docker_test/10_test/workflow_test.aligner/ref/Homo_sapiens_assembly38.fasta",
    "CNVSomaticPanelWorkflow.ref_fasta_fai": "/storage1/fs1/jin.zhang/Active/tmp/rohil_docker_test/10_test/workflow_test.aligner/ref/Homo_sapiens_assembly38.fasta.fai",
    "CNVSomaticPanelWorkflow.ref_fasta_dict": "/storage1/fs1/jin.zhang/Active/tmp/rohil_docker_test/10_test/workflow_test.aligner/ref/Homo_sapiens_assembly38.dict",
    "CNVSomaticPanelWorkflow.intervals": "/storage1/fs1/jin.zhang/Active/tmp/rohil_docker_test/tso500_intervals.bed",

    "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.1.4.0",
    "CNVSomaticPanelWorkflow.gatk4_jar_override": "/storage1/fs1/jin.zhang/Active/tmp/rohil_docker_test/apps/gatk-package-4.1.4.1-local.jar",
    "CNVSomaticPanelWorkflow.preemptible_attempts": "3"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
