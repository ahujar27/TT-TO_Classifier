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
  "Mutect2.gatk_docker": "broadinstitute/gatk:4.1.4.1",
  "Mutect2.gatk_override": "{}/gatk-package-4.1.4.1-local.jar".format(appPath),

  "Mutect2.intervals": intervalsPath,
  "Mutect2.scatter_count": 50,
  "Mutect2.m2_extra_args": "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

  "Mutect2.ref_fasta": "{}/Homo_sapiens_assembly38.fasta".format(refPath),
  "Mutect2.ref_dict": "{}/Homo_sapiens_assembly38.dict".format(refPath),
  "Mutect2.ref_fai": "{}/Homo_sapiens_assembly38.fasta.fai".format(refPath),
  "Mutect2.normal_reads": normalBam,
  "Mutect2.normal_reads_index": normalBamInd,
  "Mutect2.tumor_reads": tumorBam,
  "Mutect2.tumor_reads_index": tumorBamInd,

  "Mutect2.pon": "{}/somatic-hg38/1000g_pon.hg38.vcf".format(refPath),
  "Mutect2.pon_idx": "{}/somatic-hg38/1000g_pon.hg38.vcf.idx".format(refPath),
  "Mutect2.gnomad": "{}/somatic-hg38/af-only-gnomad.hg38.vcf".format(refPath),
  "Mutect2.gnomad_idx": "{}/somatic-hg38/af-only-gnomad.hg38.vcf.idx".format(refPath),
  "Mutect2.variants_for_contamination": "{}/somatic-hg38/small_exac_common_3.hg38.vcf".format(refPath),
  "Mutect2.variants_for_contamination_idx": "{}/somatic-hg38/small_exac_common_3.hg38.vcf.idx".format(refPath),
  "Mutect2.realignment_index_bundle": "{}/Homo_sapiens_assembly38.index_bundle".format(refPath)
}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
