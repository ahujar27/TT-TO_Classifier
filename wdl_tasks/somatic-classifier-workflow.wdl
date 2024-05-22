version 1.0

import "paired-fastq-to-unmapped-bam.wdl" as convert_paired_fastqs_to_unmapped_bam
import "processing-for-variant-discovery-gatk4.wdl" as pre_processing_for_variant_discovery
import "mutect2.wdl" as mutect2
import "maf_conversion.wdl" as MafConversion
import "input_struct.wdl"

workflow somaticClassifier {
  input {
    SampleInputs sampleInputs
    References references
    Mutect2inputs mutect2inputs
    MAFinputs mafInputs
  }

  # Step: 2b. Pre-Alignment starting with Fastq Files
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L38
  call convert_paired_fastqs_to_unmapped_bam.ConvertPairedFastQsToUnmappedBamWf as TumourFastq2UnmappedBam {
    input:
      sample_name = sampleInputs.sample_name,
      fastq_1 = sampleInputs.fastqs1,
      fastq_2 = sampleInputs.fastqs2,
      readgroup_name = sampleInputs.sample_name,
      library_name = sampleInputs.library_name,
      platform_unit = sampleInputs.platform_unit,
      run_date = sampleInputs.run_date
  }

  File flowcell_unmapped_bams_list = select_first([TumourFastq2UnmappedBam.unmapped_bam_list])

  # Step: 2c. Alignment;
  # as seen in https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L48-L62
  call pre_processing_for_variant_discovery.PreProcessingForVariantDiscovery_GATK4 as TumourSamplePreProcessing {
    input:
      sample_name = sampleInputs.sample_name,
      flowcell_unmapped_bams_list = flowcell_unmapped_bams_list,
      ref_name = references.ref_name,
      unmapped_bam_suffix = ".bam",
      ref_fasta = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      ref_dict = references.ref_dict,
      ref_alt = references.ref_alt,
      ref_sa = references.ref_sa,
      ref_ann = references.ref_ann,
      ref_bwt = references.ref_bwt,
      ref_pac = references.ref_pac,
      ref_amb = references.ref_amb,
      dbSNP_vcf = references.dbSNP_vcf,
      dbSNP_vcf_index = references.dbSNP_vcf_index,
      known_indels_sites_VCFs = references.known_indels_sites_VCFs,
      known_indels_sites_indices = references.known_indels_sites_indices
  }

  call mutect2.Mutect2 {
    input:
      tumor_reads = TumourSamplePreProcessing.analysis_ready_bam,
      tumor_reads_index = TumourSamplePreProcessing.analysis_ready_bam_index,

      scatter_count = 50,
      m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

      ref_fasta = references.ref_fasta,
      ref_fai = references.ref_fasta_index,
      ref_dict = references.ref_dict,

      realignment_index_bundle = mutect2inputs.realignment_index_bundle,
      gatk_docker = mutect2inputs.gatk_docker_version,
      gatk_override = mutect2inputs.gatk_docker_version
  }

  call MafConversion.MafConversion {
    input:
      inputVCF = Mutect2.filtered_vcf,
      sampleName = sampleInputs.sample_name,
      ref_fasta = references.ref_fasta,
      vep_path = mafInputs.vep_path,
      dockerPath = mafInputs.docker_path,
  }

  output {
    File mafFile = MafConversion.mafFile
  }
}
