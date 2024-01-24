import json
import sys

tumorBam = sys.argv[1]
tumorBamInd = sys.argv[2]
normalBam = sys.argv[3]
normalBamInd = sys.argv[4]
outName = sys.argv[5]

appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'

jsonDict = {
  "CNVSomaticPairWorkflow.tumor_bam": tumorBam,
  "CNVSomaticPairWorkflow.tumor_bam_idx": tumorBamInd,
  "CNVSomaticPairWorkflow.normal_bam": normalBam,
  "CNVSomaticPairWorkflow.normal_bam_idx": normalBamInd,
  
  "CNVSomaticPairWorkflow.ref_fasta": "{}/Homo_sapiens_assembly38.fasta".format(refPath),
  "CNVSomaticPairWorkflow.ref_fasta_fai": "{}/Homo_sapiens_assembly38.fasta.fai".format(refPath),
  "CNVSomaticPairWorkflow.ref_fasta_dict": "{}/Homo_sapiens_assembly38.dict".format(refPath),
  "CNVSomaticPairWorkflow.common_sites": "{}/cnv/snp151Common.bed".format(refPath),
  "CNVSomaticPairWorkflow.read_count_pon": "{}/cnv/tuo-cnv-pon.pon.hdf5".format(refPath), 
  "CNVSomaticPairWorkflow.intervals": intervalsPath,

  "CNVSomaticPairWorkflow.gatk_docker": "broadinstitute/gatk:4.1.4.0",
  "CNVSomaticPairWorkflow.gatk4_jar_override": "{}/gatk-package-4.1.4.1-local.jar".format(appPath),

  "CNVSomaticPairWorkflow.is_run_funcotator": "false"
}


with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
