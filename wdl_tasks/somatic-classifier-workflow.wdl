version 1.0

import "paired-fastq-to-unmapped-bam.wdl" as convert_paired_fastqs_to_unmapped_bam
import "processing-for-variant-discovery-gatk4.wdl" as pre_processing_for_variant_discovery
import "mutect2.wdl" as mutect2
import "maf_conversion.wdl" as MafConversion
import "input_struct.wdl"

workflow somaticClassifier {
  input {
    SampleInputs SampleInputs
    References References
    Mutect2inputs Mutect2inputs
    MAFinputs MAFinputs
  }

  # **Pre-Alignment starting with Fastq Files**
  # as seen in https://github.com/ahujar27/TT-TO_Classifier?tab=readme-ov-file#2b-pre-alignment
  call convert_paired_fastqs_to_unmapped_bam.ConvertPairedFastQsToUnmappedBamWf as TumourFastq2UnmappedBam {
    input:
      sample_name = SampleInputs.sample_name,
      fastq_1 = SampleInputs.fastqs1,
      fastq_2 = SampleInputs.fastqs2,
      readgroup_name = SampleInputs.sample_name,
      library_name = SampleInputs.library_name,
      platform_unit = SampleInputs.platform_unit,
      platform_name = SampleInputs.platform_name,
      sequencing_center = SampleInputs.sequencing_center,
      run_date = SampleInputs.run_date
  }

  File flowcell_unmapped_bams_list = select_first([TumourFastq2UnmappedBam.unmapped_bam_list])

  # **Alignment**
  # as seen in https://github.com/ahujar27/TT-TO_Classifier?tab=readme-ov-file#2c-alignment
  call pre_processing_for_variant_discovery.PreProcessingForVariantDiscovery_GATK4 as TumourSamplePreProcessing {
    input:
      sample_name = SampleInputs.sample_name,
      flowcell_unmapped_bams_list = flowcell_unmapped_bams_list,
      ref_name = References.ref_name,
      unmapped_bam_suffix = ".bam",
      ref_fasta = References.ref_fasta,
      ref_fasta_index = References.ref_fasta_index,
      ref_dict = References.ref_dict,
      ref_alt = References.ref_alt,
      ref_sa = References.ref_sa,
      ref_ann = References.ref_ann,
      ref_bwt = References.ref_bwt,
      ref_pac = References.ref_pac,
      ref_amb = References.ref_amb,
      dbSNP_vcf = References.dbSNP_vcf,
      dbSNP_vcf_index = References.dbSNP_vcf_index,
      known_indels_sites_VCFs = References.known_indels_sites_VCFs,
      known_indels_sites_indices = References.known_indels_sites_indices
  }

  # ** Tumour-Only SNV Calling **
  # as seen in https://github.com/ahujar27/TT-TO_Classifier?tab=readme-ov-file#2d-snv-calling
  call mutect2.Mutect2 {
    input:
      tumor_reads = TumourSamplePreProcessing.analysis_ready_bam,
      tumor_reads_index = TumourSamplePreProcessing.analysis_ready_bam_index,

      scatter_count = 50,
      m2_extra_args = "--downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6",

      ref_fasta = References.ref_fasta,
      ref_fai = References.ref_fasta_index,
      ref_dict = References.ref_dict,

      realignment_index_bundle = Mutect2inputs.realignment_index_bundle,
      gatk_docker = Mutect2inputs.gatk_docker_version,
      gatk_override = Mutect2inputs.gatk_docker_version
  }

  # ** MAF Conversion **
  # as seen in https://github.com/ahujar27/TT-TO_Classifier?tab=readme-ov-file#3a-maf-conversion
  call MafConversion.MafConversion {
    input:
      inputVCF = Mutect2.filtered_vcf,
      sampleName = SampleInputs.sample_name,
      ref_fasta = References.ref_fasta,
      vep_path = MAFinputs.vep_path,
      dockerPath = MAFinputs.docker_path,
  }

  # ** Tumour-Only Feature Generation **
  # as seen in https://github.com/ahujar27/TT-TO_Classifier?tab=readme-ov-file#4-feature-generation

  output {
    File mafFile = MafConversion.mafFile
  }
}
