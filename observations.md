# Doubts for the researchers:

## 1. Regarding step: 2b. Pre-Alignment

There was a missing reference file: `sourceBAM`.
This file is required when starting with FASTQ files:
https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/SNV_CNV/alignment/make.paired-fastq-to-unmapped-bam.input.json.py#L5-L7

in the above mentioned script, and also required on `2c. Alignment` step:
https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L52

* Could you please send us the file or the proper path to this file?

---
## 2. Regarding step: 2c. Alignment

There were anothers missing references files for the task [processing-for-variant-discovery-gatk4](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L62):

    * ref_alt
    * dbSNP_vcf
    * dbSNP_vcf_index
    * known_indels_sites_VCFs
    * known_indels_sites_indices

The [downloaded files according to the README instructions](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L11-L22) are differently named than those obtained downloading those references:

| Task input required: | Reference file name on [make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py](SNV_CNV/alignment/make.processing-for-variant-discovery-gatk4.hg38.wgs.input.json.py) | [Reference files downloading according to README.md file](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L11-L22)
| --- | --- | --- |
| PreProcessingForVariantDiscovery_GATK4.ref_dict | Homo_sapiens_assembly38.dict | references/GRCh38.d1.vd1.dict
| PreProcessingForVariantDiscovery_GATK4.ref_fasta | Homo_sapiens_assembly38.fasta | references/GRCh38.d1.vd1.fa
| PreProcessingForVariantDiscovery_GATK4.ref_fasta_index | Homo_sapiens_assembly38.fasta.fai | references/GRCh38.d1.vd1.fa.fai
| PreProcessingForVariantDiscovery_GATK4.ref_alt | Homo_sapiens_assembly38.fasta.64.alt | `no equivalent file found`
| PreProcessingForVariantDiscovery_GATK4.ref_sa | Homo_sapiens_assembly38.fasta.64.sa | references/GRCh38.d1.vd1.fa.sa
| PreProcessingForVariantDiscovery_GATK4.ref_amb | Homo_sapiens_assembly38.fasta.64.amb | references/GRCh38.d1.vd1.fa.amb
| PreProcessingForVariantDiscovery_GATK4.ref_bwt | Homo_sapiens_assembly38.fasta.64.bwt | references/GRCh38.d1.vd1.fa.bwt
| PreProcessingForVariantDiscovery_GATK4.ref_ann | Homo_sapiens_assembly38.fasta.64.ann | references/GRCh38.d1.vd1.fa.ann
| PreProcessingForVariantDiscovery_GATK4.ref_pac | Homo_sapiens_assembly38.fasta.64.pac | references/GRCh38.d1.vd1.fa.pac
| PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf | Homo_sapiens_assembly38.dbsnp138.vcf | `no equivalent file found`
| PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index | Homo_sapiens_assembly38.dbsnp138.vcf.idx | `no equivalent file found`
| PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs | [Mills_and_1000G_gold_standard.indels.hg38.vcf.gz, Homo_sapiens_assembly38.known_indels.vcf.gz] | `no equivalent file found`
| PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices | [Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi, Homo_sapiens_assembly38.known_indels.vcf.gz.tbi] | `no equivalent file found`

* Could you please double-check the reference files and supply the missing ones in the `.tar` files?

---
## 3. Regarding '2d. SNV Calling Tumor Normal | Tumor-Only' steps:

[Currently we do not have access to this repos or files](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/SNV_CNV/snv/make.mutect2.tumor_only.input.json.py#L10-L12):
```
appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'
```

Therefore those files are missing for the Mutect task:

* gatk_override = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps/gatk-package-4.1.4.1-local.jar"
* pon = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf"
* pon_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf.idx"
* gnomad = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf"
* gnomad_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf.idx"
* variants_for_contamination = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf"
* variants_for_contamination_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf.idx"
* realignment_index_bundle = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/Homo_sapiens_assembly38.index_bundle"

---
## 4. Regarding step: 2e. CNV Panel of Normals Generation

### 4.1. Is the [**NORMALBAMLIST** (List of normal bams to generate the panel of normals)](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L93) a list of the [normal BAM files' that should have been created on item 2.c](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L71)

* ['normal_BAM_file_1.bam','normal_BAM_file_2.bam']
* ['normal_BAM_file_1.bam.bai','normal_BAM_file_2.bam.bai']
---

### 4.2. Is a PON going to be build for each sample processing?

---
### 4.3 If a sample is tumor-only, does it will go through the CNV Calling step?
#### 4.3.1 If so, which PON is supposed to be used?

---
## 5. WDL workflows from the paths described could not get validated through womtools:

### 5.1 `cnv_somatic_pair_workflow.wdl`
The workflow copied as described from:
https://github.com/ahujar27/TT-TO_Classifier/blob/0de949714a4e17ebdf74fd6b1ee6729ee6cdd1ed/README.md?plain=1#L119

In a validation attempt using `Womtools`:
```bash
(base) ➜  TT-TO_Classifier git:(vars_somatic_classifier) ✗ java -jar ../womtool-86.jar validate wdl_tasks/cnv_somatic_pair_workflow.wdl
```

Got the following error message:
```bash
ERROR: Unexpected symbol (line 56, col 5) when parsing 'setter'.

Expected equal, got "File".

    File? intervals
    ^

$setter = :equal $e -> $1
```

---
### 5.2 `cnv_somatic_panel_workflow.wdl`
The workflow copied as described from:
https://github.com/ahujar27/TT-TO_Classifier/blob/0de949714a4e17ebdf74fd6b1ee6729ee6cdd1ed/README.md?plain=1#L102

In a validation attempt using `Womtools`:
```bash
(base) ➜  TT-TO_Classifier git:(vars_somatic_classifier) ✗ java -jar ../womtool-86.jar validate wdl_tasks/cnv_somatic_panel_workflow.wdl
```

Got the following error message:
```bash
ERROR: Unexpected symbol (line 43, col 5) when parsing 'setter'.

Expected equal, got "File".

    File? blacklist_intervals
    ^

$setter = :equal $e -> $1
```

---
### 5.3 `paired-fastq-to-unmapped-bam.wdl`
The workflow copied as described from:
https://github.com/ahujar27/TT-TO_Classifier/blob/0de949714a4e17ebdf74fd6b1ee6729ee6cdd1ed/README.md?plain=1#L38

In a validation attempt using `Womtools`:
```bash
java -jar ../womtool-86.jar validate wdl_tasks/paired-fastq-to-unmapped-bam.wdl
```

Got the following error message:
```bash
Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 1 of 3): Cannot lookup value 'run_date', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 2 of 3): Cannot lookup value 'sequencing_center', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 3 of 3): Cannot lookup value 'platform_name', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
```

---
### 5.4 Global workflow: `somatic-classifier-workflow.wdl`

Thus, in a validation attempt using `Womtools`:
```bash
java -jar ../womtool-86.jar validate wdl_tasks/somatic-classifier-workflow.wdl
```

The workflow log errors are a combination of the dependencies errors:
```bash
Failed to import 'paired-fastq-to-unmapped-bam.wdl' (reason 1 of 3): Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 1 of 3): Cannot lookup value 'run_date', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
Failed to import 'paired-fastq-to-unmapped-bam.wdl' (reason 2 of 3): Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 2 of 3): Cannot lookup value 'sequencing_center', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
Failed to import 'paired-fastq-to-unmapped-bam.wdl' (reason 3 of 3): Failed to process workflow definition 'ConvertPairedFastQsToUnmappedBamWf' (reason 3 of 3): Cannot lookup value 'platform_name', it is never declared. Available values are: ['CreateFoFN.__after', 'sample_name', 'gatk_path', 'fastq_2', 'ubam_list_name', 'output_unmapped_bam', 'fastq_1', 'gatk_docker', 'PairedFastQsToUnmappedBAM.output_unmapped_bam', 'library_name', 'CreateFoFN.fofn_list', 'unmapped_bam_list', 'PairedFastQsToUnmappedBAM.__after', 'disk_multiplier', 'make_fofn', 'platform_unit', 'readgroup_name']
Failed to import 'cnv_somatic_panel_workflow.wdl' (reason 1 of 1): ERROR: Unexpected symbol (line 43, col 3) when parsing 'setter'.

Expected equal, got "File".

  File? blacklist_intervals
  ^

$setter = :equal $e -> $1

Failed to import 'cnv_somatic_pair_workflow.wdl' (reason 1 of 1): ERROR: Unexpected symbol (line 56, col 3) when parsing 'setter'.

Expected equal, got "File".

  File intervals
  ^

$setter = :equal $e -> $1

Failed to import 'maf_conversion.wdl' (reason 1 of 1): ERROR: Unexpected symbol (line 16, col 3) when parsing 'setter'.

Expected equal, got "String".

  String sampleName
  ^

$setter = :equal $e -> $1
```