version 1.0

import "bam-to-unmapped-bams.wdl" as bam_to_unmappedbam
import "paired-fastq-to-unmapped-bam.wdl" as convert_paired_fastqs_to_unmapped_bam
import "processing-for-variant-discovery-gatk4.wdl" as pre_processing_for_variant_discovery
import "mutect2.wdl" as mutect2
import "cnv_somatic_panel_workflow.wdl" as build_cnv_pon
import "cnv_somatic_pair_workflow.wdl" as CNVcalling
import "maf_conversion.wdl" as MafConversion

workflow somaticClassifier {

  String sampleName = sampleName
  File? tumour_unmapped_bam = tumourSample.unmappedBam

  File? tumourAlignedBam = tumourSample.aligned_bam

  File? tumourFastq1 = tumourSample.fastqs1
  File? tumourFastq2 = tumourSample.fastqs2
  String? tumourLibraryName = tumourSample.libraryName
  String? tumourPlatformUnit = tumourSample.platformUnit

  File? normalUnmappedBam = normalSample.unmappedBam

  File? normalAlignedBam = normalSample.aligned_bam

  File? normalFastq1 = normalSample.fastqs1
  File? normalFastq2 = normalSample.fastqs2
  String? normalLibraryName = normalSample.libraryName
  String? normalPlatformUnit = normalSample.platformUnit

  String referencesPath = "../references/"

  File tso_bed_intervals_file = "../tso500_intervals.bed"

  # Tumour Sample PreProcessing
  # Step: 2b. Pre-Alignment starting with aligned BAM files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L42
  if ((!defined(tumour_unmapped_bam)) && (defined(tumourAlignedBam))) {
    call bam_to_unmappedbam.BamToUnmappedBam as TumourBam2UnmappedBam {
      input:
        input_bam = tumourAlignedBam
    }
  }

  # Step: 2b. Pre-Alignment starting with Fastq Files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L38
  if ((!defined(tumour_unmapped_bam)) &&
  (defined(tumourFastq1)) &&
  (defined(tumourFastq2)) &&
  (defined(tumourLibraryName)) &&
  (defined(tumourPlatformUnit))) {
    call convert_paired_fastqs_to_unmapped_bam.ConvertPairedFastQsToUnmappedBamWf as TumourFastq2UnmappedBam {
      input:
        sample_name = sampleName,
        fastq_1 = tumourSample.fastq_1,
        fastq_2 = tumourSample.fastq_2,
        readgroup_name = "tumour" + sampleName,
        library_name = tumourSample.libraryName,
        platform_unit = platformUnit,
        run_date = "2020-01-01"
    }
  }

  Array[File] tumour_unmapped_bams_list = if (defined(tumour_unmapped_bam))
  then [tumour_unmapped_bam]
  else select_first([TumourBam2UnmappedBam.output_bams,[TumourFastq2UnmappedBam.output_unmapped_bam]])

  # Step: 2c. Alignment;
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L48-L62
  call pre_processing_for_variant_discovery.PreProcessingForVariantDiscovery_GATK4 as TumourSamplePreProcessing {
    input:
      sample_name = sampleName,
      flowcell_unmapped_bams_list = tumour_unmapped_bams_list,
      ref_name = "hg38",
      unmapped_bam_suffix = ".bam",
      ref_fasta = path + "GRCh38.d1.vd1.fa",
      ref_fasta_index = path + "GRCh38.d1.vd1.fa.fai",
      ref_dict = path + "GRCh38.d1.vd1.dict",
      ref_alt = path,
      ref_sa = path + "GRCh38.d1.vd1.fa.sa",
      ref_ann = path + "GRCh38.d1.vd1.fa.ann",
      ref_bwt = path + "GRCh38.d1.vd1.fa.bwt",
      ref_pac = path + "GRCh38.d1.vd1.fa.pac",
      ref_amb = path + "GRCh38.d1.vd1.fa.amb",
      dbSNP_vcf = path,
      dbSNP_vcf_index = path,
      known_indels_sites_VCFs = [path],
      known_indels_sites_indices = [path]
  }

  # Normal Sample PreProcessing
    # Step: 2b. Pre-Alignment starting with aligned BAM files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L42
  if ((!defined(normalUnmappedBam)) &&
  (defined(normalAlignedBam))) {
    call bam_to_unmappedbam.BamToUnmappedBam as Bam2UnmappedBam {
      input:
        input_bam = normalAlignedBam
    }
  }

  # Step: 2b. Pre-Alignment starting with Fastq Files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L38
  if ((!defined(normalUnmappedBam)) &&
  (defined(normalFastq1)) &&
  (defined(normalFastq2)) &&
  (defined(normalLibraryName)) &&
  (defined(normalPlatformUnit))) {
    call convert_paired_fastqs_to_unmapped_bam.ConvertPairedFastQsToUnmappedBamWf as Fastq2UnmappedBam {
      input:
        sample_name = sampleName,
        fastq_1 = normalFastq1,
        fastq_2 = normalFastq2,
        readgroup_name = "normal" + sampleName,
        library_name = normalLibraryName,
        platform_unit = normalPlatformUnit
    }
  }

  Array[File] unmapped_bams_list = if (defined(normalUnmappedBam))
  then [normalUnmappedBam]
  else select_first([Bam2UnmappedBam.output_bams,[Fastq2UnmappedBam.output_unmapped_bam]])

  # Step: 2c. Alignment;
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L48-L62
  call pre_processing_for_variant_discovery.PreProcessingForVariantDiscovery_GATK4 as SamplePreProcessing {
    input:
      sample_name = sampleName,
      flowcell_unmapped_bams_list = unmapped_bams_list,
      ref_name = "hg38",
      unmapped_bam_suffix = ".bam",
      ref_fasta = path + "GRCh38.d1.vd1.fa",
      ref_fasta_index = path + "GRCh38.d1.vd1.fa.fai",
      ref_dict = path + "GRCh38.d1.vd1.dict",
      ref_alt = path,
      ref_sa = path + "GRCh38.d1.vd1.fa.sa",
      ref_ann = path + "GRCh38.d1.vd1.fa.ann",
      ref_bwt = path + "GRCh38.d1.vd1.fa.bwt",
      ref_pac = path + "GRCh38.d1.vd1.fa.pac",
      ref_amb = path + "GRCh38.d1.vd1.fa.amb",
      dbSNP_vcf = path,
      dbSNP_vcf_index = path,
      known_indels_sites_VCFs = [path],
      known_indels_sites_indices = [path]
  }

  ## 2d. SNV Calling
  String analysis_mode = if (defined(SamplePreProcessing.analysis_ready_bam)) then "Tumor-Normal" else "Tumor-Only"

  # Step: 2d. SNV Calling Tumor-Normal
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L67-L77

  if (analysis_mode == "Tumor-Normal") {

    call mutect2.Mutect2 {
      input:
        intervalsPath = tso_bed_intervals_file,
        normal_reads = SamplePreProcessing.analysis_ready_bam,
        normal_reads_index = SamplePreProcessing.analysis_ready_bam_index,
        tumor_reads = TumourSamplePreProcessing.analysis_ready_bam,
        tumor_reads_index = TumourSamplePreProcessing.analysis_ready_bam_index,

        scatter_count = 50,
        m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

        ref_fasta = path + "GRCh38.d1.vd1.fa",
        ref_fai = path + "GRCh38.d1.vd1.fa.fai",
        ref_dict = path + "GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.1.4.1",

        # Currently we do not have access to this repos or files
        #
        # appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
        # intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
        # refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'
        #
        # therefore those are missing files for this task:
        gatk_override = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps/gatk-package-4.1.4.1-local.jar",
        pon = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf",
        pon_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf.idx",
        gnomad = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf",
        gnomad_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf.idx",
        variants_for_contamination = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf",
        variants_for_contamination_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf.idx",
        realignment_index_bundle = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/Homo_sapiens_assembly38.index_bundle"
    }

    # call build_cnv_pon.CNVSomaticPanelWorkflow {
    #   input:
    #     normal_bams = SamplePreProcessing.analysis_ready_bam,
    #     normal_bais = SamplePreProcessing.analysis_ready_bam_index,

    #     pon_entity_id = "tuo-cnv-pon",

    #     ref_fasta = path + "GRCh38.d1.vd1.fa",
    #     ref_fai = path + "GRCh38.d1.vd1.fa.fai",
    #     ref_dict = path + "GRCh38.d1.vd1.dict",
    #     gatk_docker = "broadinstitute/gatk:4.1.4.1",
    #     preemptible_attempts = 3,

    #     # Currently configuring only TSO500 intervals
    #     intervals = tso_bed_intervals_file
    # }

    call CNVcalling.CNVSomaticPairWorkflow {
      input:
        common_sites = "gs://gatk-test-data/cnv/somatic/common_snps.interval_list",
        intervals = tso_bed_intervals_file,
        tumor_bam = TumourSamplePreProcessing.analysis_ready_bam,
        tumor_bam_idx = TumourSamplePreProcessing.analysis_ready_bam_index,
        normal_bam = SamplePreProcessing.analysis_ready_bam,
        normal_bam_idx = SamplePreProcessing.analysis_ready_bam_index,

        read_count_pon = CNVSomaticPanelWorkflow.read_count_pon,

        ref_fasta = path + "GRCh38.d1.vd1.fa",
        ref_fai = path + "GRCh38.d1.vd1.fa.fai",
        ref_dict = path + "GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.4.0.0",
        gatk4_jar_override = "gatk-package-4.4.0.0-local.jar"
    }

    call MAFconversion.MafConversion {
      input:
        inputVCF = Mutect2.unfiltered_vcf,
        sampleName = sampleName
    }
  }

  # Step: 2d. SNV Calling Tumor-Only
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L79-L89

  if (analysis_mode == "Tumor-Only") {
    call mutect2.Mutect2 {
      input:
        tumor_reads = TumourSamplePreProcessing.analysis_ready_bam,
        tumor_reads_index = TumourSamplePreProcessing.analysis_ready_bam_index,

        scatter_count = 50,
        m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

        ref_fasta = path + "GRCh38.d1.vd1.fa",
        ref_fai = path + "GRCh38.d1.vd1.fa.fai",
        ref_dict = path + "GRCh38.d1.vd1.dict",

        gatk_docker = "broadinstitute/gatk:4.1.4.1",

        # Currently we do not have access to this repos or files
        #
        # appPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps'
        # intervalsPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/tso500_intervals.bed'
        # refPath = '/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref'
        #
        # therefore those are missing files for this task:
        gatk_override = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/apps/gatk-package-4.1.4.1-local.jar",

        intervalsPath = tso_bed_intervals_file,
        pon = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf",
        pon_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/1000g_pon.hg38.vcf.idx",
        gnomad = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf",
        gnomad_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/af-only-gnomad.hg38.vcf.idx",
        variants_for_contamination = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf",
        variants_for_contamination_idx = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/somatic-hg38/small_exac_common_3.hg38.vcf.idx",
        realignment_index_bundle = "/storage1/fs1/jin.zhang/Active/rohil/gatk_best_practices_test/ref/Homo_sapiens_assembly38.index_bundle"
    }

    call MafConversion.MafConversion {
      input:
        inputVCF = Mutect2.unfiltered_vcf,
        sampleName = sampleInputs.sampleName
    }
  }

  output {
    Array[File] output_bams = SortSam.sorted_bam
  }
}
