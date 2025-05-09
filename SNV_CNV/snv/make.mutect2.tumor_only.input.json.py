import json
import sys

tumorBam = sys.argv[1]
tumorBamInd = sys.argv[2]
outName = sys.argv[3]
intervalsPath = sys.argv[4]
refPath = sys.argv[5]

jsonDict = {
  "Mutect2.gatk_docker": "broadinstitute/gatk:4.1.4.1",
  "Mutect2.gatk_override": "/opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar".format(appPath),

  "Mutect2.intervals": intervalsPath,
  "Mutect2.scatter_count": 50,
  "Mutect2.m2_extra_args": "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

  "Mutect2.ref_fasta": "/path-to/references/hg38_references/Homo_sapiens_assembly38.fasta",
  "Mutect2.ref_dict": "/path-to/references/hg38_references/Homo_sapiens_assembly38.dict",
  "Mutect2.ref_fai": "/path-to/references/hg38_references/Homo_sapiens_assembly38.fasta.fai",
#  "Mutect2.normal_reads": normalBam,
#  "Mutect2.normal_reads_index": normalBamInd,
  "Mutect2.tumor_reads": tumorBam,
  "Mutect2.tumor_reads_index": tumorBamInd,

  "Mutect2.pon": "/path-to/references/hg38_references/1000g_pon.hg38.vcf.gz",
  "Mutect2.pon_idx": "/path-to/references/hg38_references/1000g_pon.hg38.vcf.gz.tbi",
  "Mutect2.gnomad": "/path-to/references/hg38_references/af-only-gnomad.hg38.vcf.gz",
  "Mutect2.gnomad_idx": "/path-to/references/hg38_references/af-only-gnomad.hg38.vcf.gz.tbi",
  "Mutect2.variants_for_contamination": "/path-to/references/hg38_references/small_exac_common_3.hg38.vcf.gz",
  "Mutect2.variants_for_contamination_idx": "/path-to/references/hg38_references/small_exac_common_3.hg38.vcf.gz.tbi",
  "Mutect2.realignment_index_bundle": "/path-to/references/hg38_references/Homo_sapiens_assembly38.index_bundle"
}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
