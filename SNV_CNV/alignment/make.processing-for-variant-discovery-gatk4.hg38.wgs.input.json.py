import pysam
import json
import sys

sourceBam = sys.argv[1]
bamsList = sys.argv[2]
outName = sys.argv[3]

refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'

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
            "PreProcessingForVariantDiscovery_GATK4.ref_dict": "{}/Homo_sapiens_assembly38.dict".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_fasta": "{}/Homo_sapiens_assembly38.fasta".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_fasta_index": "{}/Homo_sapiens_assembly38.fasta.fai".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_alt": "{}/Homo_sapiens_assembly38.fasta.64.alt".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_sa": "{}/Homo_sapiens_assembly38.fasta.64.sa".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_amb": "{}/Homo_sapiens_assembly38.fasta.64.amb".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_bwt": "{}/Homo_sapiens_assembly38.fasta.64.bwt".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_ann": "{}/Homo_sapiens_assembly38.fasta.64.ann".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.ref_pac": "{}/Homo_sapiens_assembly38.fasta.64.pac".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf": "{}/Homo_sapiens_assembly38.dbsnp138.vcf".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index": "{}/Homo_sapiens_assembly38.dbsnp138.vcf.idx".format(refPath),
            "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs": ["{}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz".format(refPath),"{}/Homo_sapiens_assembly38.known_indels.vcf.gz".format(refPath)],
  "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices": ["{}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi".format(refPath),"{}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi".format(refPath)]
}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
