version 1.0

import "bam-to-unmapped-bams.wdl" as bam_to_unmappedbam
import "paired-fastq-to-unmapped-bam.wdl" as convert_paired_fastqs_to_unmapped_bam
import "processing-for-variant-discovery-gatk4.wdl" as pre_processing_for_variant_discovery
import "mutect2.wdl" as mutect2
import "cnv_somatic_panel_workflow.wdl" as build_cnv_pon
import "cnv_somatic_pair_workflow.wdl" as CNVcalling
import "maf_conversion.wdl" as MAFconversion


workflow somaticClassifier {

  input {

    SampleInputs sampleInputs

  }

  # Step: 2b. Pre-Alignment starting with aligned BAM files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L42
  if defined(sampleInputs.aligned_bam) {
    call bam_to_unmappedbam.BamToUnmappedBam {
      input:
        input_bam = sampleInputs.aligned_bam,
    }
  }

  # Step: 2b. Pre-Alignment starting with Fastq Files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L38
  if defined(sampleInputs.fastq_1)
  && defined(sampleInputs.fastq_2)
  && defined(sampleInputs.readgroup_name)
  && defined(sampleInputs.library_name)
  && defined(sampleInputs.platform_unit) {
    call convert_paired_fastqs_to_unmapped_bam.ConvertPairedFastQsToUnmappedBamWf {
      input:
        sample_name = sampleInputs.sampleName,
        fastq_1 = sampleInputs.fastq_1,
        fastq_2 = sampleInputs.fastq_2,
        readgroup_name = sampleInputs.readgroupName,
        library_name = sampleInputs.libraryName,
        platform_unit = sampleInputs.platformUnit
    }
  }

  Array[File] flowcell_unmapped_bams_list = if (defined(BamToUnmappedBam.output_bams))
  then BamToUnmappedBam.output_bams
  else [select_first([sampleInputs.unmappedBam, ConvertPairedFastQsToUnmappedBamWf.output_unmapped_bam])],

  # Step: 2c. Alignment;
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L48-L62
  call pre_processing_for_variant_discovery.PreProcessingForVariantDiscovery_GATK4 {
    input:
      sample_name = sampleInputs.sampleName,
      flowcell_unmapped_bams_list = flowcell_unmapped_bams_list,
      ref_name = "hg38",
      unmapped_bam_suffix = ".bam",
      ref_fasta = "../references/GRCh38.d1.vd1.fa",
      ref_fasta_index = "../references/GRCh38.d1.vd1.fa.fai",
      ref_dict = "../references/GRCh38.d1.vd1.dict",
      ref_alt = "../references/",
      ref_sa = "../references/GRCh38.d1.vd1.fa.sa",
      ref_ann = "../references/GRCh38.d1.vd1.fa.ann",
      ref_bwt = "../references/GRCh38.d1.vd1.fa.bwt",
      ref_pac = "../references/GRCh38.d1.vd1.fa.pac",
      ref_amb = "../references/GRCh38.d1.vd1.fa.amb",
      dbSNP_vcf = "../references/",
      dbSNP_vcf_index = "../references/",
      known_indels_sites_VCFs = ["../references/"],
      known_indels_sites_indices = ["../references/"]
    }
  }

  ## Regarding step: 2d. SNV Calling
  ##
  ## This step is described as:
  ## [ This step calls the variants. If a normal and a tumor are both present follow the tumor-normal SNV calling. If a normal is not present, then follow the steps for tumor-only calling.](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L65)
  ##
  ## However, the outputs of the previous task [processing-for-variant-discovery-gatk4](https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L62) are:
  ##
  ## ```
  ## File duplication_metrics = MarkDuplicates.duplicate_metrics
  ## File bqsr_report = GatherBqsrReports.output_bqsr_report
  ## File analysis_ready_bam = GatherBamFiles.output_bam
  ## File analysis_ready_bam_index = GatherBamFiles.output_bam_index
  ## File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
  ## ```
  ##
  ## Therefore are no specifications for `tumour BAM files` or `normal BAM files`, only `analysis_ready_bam` thus, neither of the conditional steps in sequence will have their requirements fullfilled.

  String AnalysisMode = if (defined(normalBam) && defined(normalBamInd))
  then "Tumor-Normal"
  else "Tumor-Only"

  # Step: 2d. SNV Calling Tumor-Normal
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L67-L77

  if (AnalysisMode == "Tumor-Normal") {

    call mutect2.Mutect2 {
      input:
        normal_reads = normalBam,
        normal_reads_index = normalBamInd,
        tumor_reads = tumorBam,
        tumor_reads_index = tumorBamInd,

        scatter_count = 50,
        m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

        ref_fasta = "../references/GRCh38.d1.vd1.fa",
        ref_fai = "../references/GRCh38.d1.vd1.fa.fai",
        ref_dict = "../references/GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.1.4.1"

        # Currently we do not have access to this repos or files
        #
        # appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
        # intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
        # refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'
        #
        # therefore those are missing files for this task:
        gatk_override = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps/gatk-package-4.1.4.1-local.jar"

        intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
        pon = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf"
        pon_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf.idx"
        gnomad = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf"
        gnomad_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf.idx"
        variants_for_contamination = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf"
        variants_for_contamination_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf.idx"
        realignment_index_bundle = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/Homo_sapiens_assembly38.index_bundle"
    }

    # Considering that normal_bams and normal_bais on CNV Panel of Normals Generation
    # are lists composed of the normalBam, normalBamInd,
    # that should have been created on Step: 2c. Alignment

    Array[String] normal_bams = glob("*normal.bam")
    Array[String] normal_bais = glob("*normal.bam.bai")

    # Will a Panel of Normals be created for each sample processing?
    call build_cnv_pon.CNVSomaticPanelWorkflow {
      input:

        normal_bams = normal_bams
        normal_bais = normal_bais

        pon_entity_id = "tuo-cnv-pon",

        ref_fasta = "../references/GRCh38.d1.vd1.fa",
        ref_fai = "../references/GRCh38.d1.vd1.fa.fai",
        ref_dict = "../references/GRCh38.d1.vd1.dict",
        gatk_docker = "broadinstitute/gatk:4.1.4.1",
        preemptible_attempts = 3

        # Currently configuring only TSO500 intervals
        intervals = "tso500_intervals.bed"
    }

    call CNVcalling.CNVSomaticPairWorkflow {
      input:
        common_sites = "gs://gatk-test-data/cnv/somatic/common_snps.interval_list"
        intervals = "tso500_intervals.bed"

        tumor_bam = tumorBam,
        tumor_bam_idx = tumorBamInd,
        normal_bam = normalBam,
        normal_bam_idx = normalBamInd,

        read_count_pon = CNVSomaticPanelWorkflow.read_count_pon

        ref_fasta = "../references/GRCh38.d1.vd1.fa",
        ref_fai = "../references/GRCh38.d1.vd1.fa.fai",
        ref_dict = "../references/GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.4.0.0",
        gatk4_jar_override = "gatk-package-4.4.0.0-local.jar"
    }

    call MAFconversion.MafConversion {
      input:
        inputVCF = filtered_vcf,
        sampleName = sampleInputs.sampleName
    }

  }

  # Step: 2d. SNV Calling Tumor-Only
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L79-L89

  if (AnalysisMode == "Tumor-Only") {
    call mutect2.Mutect2 {
      input:
        tumor_reads = tumorBam,
        tumor_reads_index = tumorBamInd,

        scatter_count = 50,
        m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

        ref_fasta = "../references/GRCh38.d1.vd1.fa",
        ref_fai = "../references/GRCh38.d1.vd1.fa.fai",
        ref_dict = "../references/GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.1.4.1"

        # Currently we do not have access to this repos or files
        #
        # appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
        # intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
        # refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'
        #
        # therefore those are missing files for this task:
        gatk_override = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps/gatk-package-4.1.4.1-local.jar",

        intervalsPath = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed",
        pon = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf",
        pon_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf.idx",
        gnomad = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf",
        gnomad_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf.idx",
        variants_for_contamination = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf",
        variants_for_contamination_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf.idx",
        realignment_index_bundle = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/Homo_sapiens_assembly38.index_bundle"
    }

    # Will a tumor-only sample go through CNVCalling?
    # Which PON is supposed to be used?
    #
    # call CNVcalling.CNVSomaticPairWorkflow {
    #   input:
    #
    #     common_sites = "gs://gatk-test-data/cnv/somatic/common_snps.interval_list"
    #     intervals = "tso500_intervals.bed"
    #
    #     tumor_bam = tumorBam,
    #     tumor_bam_idx = tumorBamInd,
    #
    #     read_count_pon = ?
    #
    #     ref_fasta = "../references/GRCh38.d1.vd1.fa",
    #     ref_fai = "../references/GRCh38.d1.vd1.fa.fai",
    #     ref_dict = "../references/GRCh38.d1.vd1.dict",
    #
    #     gatk_docker = "broadinstitute/gatk:4.4.0.0",
    #     gatk4_jar_override = "gatk-package-4.4.0.0-local.jar"
    # }

    call MAFconversion.MafConversion {
      input:
        inputVCF = filtered_vcf,
        sampleName = sampleInputs.sampleName
    }


  }

  output {
    Array[File] output_bams = SortSam.sorted_bam
  }
}
