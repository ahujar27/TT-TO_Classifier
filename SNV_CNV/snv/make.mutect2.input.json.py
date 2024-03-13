import json
import sys

tumorBam = sys.argv[1]
tumorBamInd = sys.argv[2]
normalBam = sys.argv[3]
normalBamInd = sys.argv[4]
outName = sys.argv[5]
intervalsPath = sys.argv[6]
refPath = sys.argv[7]

jsonDict = {
  "Mutect2.gatk_docker": "broadinstitute/gatk:4.1.4.1",
  "Mutect2.gatk_override": "/opt/gatk-4.4.0.0/gatk-package-4.1.4.1-local.jar",

  "Mutect2.intervals": intervalsPath,
  "Mutect2.scatter_count": 50,
  "Mutect2.m2_extra_args": "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

  "Mutect2.ref_fasta": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict",
  "Mutect2.ref_dict": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
  "Mutect2.ref_fai": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "Mutect2.normal_reads": normalBam,
  "Mutect2.normal_reads_index": normalBamInd,
  "Mutect2.tumor_reads": tumorBam,
  "Mutect2.tumor_reads_index": tumorBamInd,

  "Mutect2.pon": "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz",
  "Mutect2.pon_idx": "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.tbi",
  "Mutect2.gnomad": "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf",
  "Mutect2.gnomad_idx": "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.idx",
  "Mutect2.variants_for_contamination": "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf",
  "Mutect2.variants_for_contamination_idx": "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.idx",
  "Mutect2.realignment_index_bundle": "gs://gatk-test-data/mutect2/Homo_sapiens_assembly38.index_bundle"
}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
