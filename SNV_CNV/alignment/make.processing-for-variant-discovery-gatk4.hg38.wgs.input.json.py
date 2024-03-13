import pysam
import json
import sys

sourceBam = sys.argv[1]
bamsList = sys.argv[2]
outName = sys.argv[3]

bam = pysam.AlignmentFile(sourceBam, mode='rb')
bam_header = dict((k, v) for k, v in bam.header.items())
if type(bam_header['RG'])==list:
    dictRG = bam_header['RG'][0]
else:
    dictRG = bam_header['RG']

tagDict = {}
tags = ['SM','ID','LB','PU','DT','PL','CN']

for tag in tags:
    if tag in dictRG:
        tagDict[tag] = dictRG[tag]
    else:
        tagDict[tag] = 'NA'

jsonDict = {"PreProcessingForVariantDiscovery_GATK4.sample_name": tagDict['SM'],
            "PreProcessingForVariantDiscovery_GATK4.ref_name": "hg38",
            "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list": bamsList,
            "PreProcessingForVariantDiscovery_GATK4.unmapped_bam_suffix": ".bam",
            "PreProcessingForVariantDiscovery_GATK4.ref_dict": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict",
            "PreProcessingForVariantDiscovery_GATK4.ref_fasta": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
            "PreProcessingForVariantDiscovery_GATK4.ref_fasta_index": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
            "PreProcessingForVariantDiscovery_GATK4.ref_alt": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
            "PreProcessingForVariantDiscovery_GATK4.ref_sa": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa",
            "PreProcessingForVariantDiscovery_GATK4.ref_amb": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb",
            "PreProcessingForVariantDiscovery_GATK4.ref_bwt": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt",
            "PreProcessingForVariantDiscovery_GATK4.ref_ann": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann",
            "PreProcessingForVariantDiscovery_GATK4.ref_pac": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac",
            "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index": "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
            "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs": ["gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz","gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"],
            "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices": ["gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi","gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"]
}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
