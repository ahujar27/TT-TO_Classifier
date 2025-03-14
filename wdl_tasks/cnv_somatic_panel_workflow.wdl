version 1.0

# As recommended in: https://github.com/ahujar27/TT-TO_Classifier/blob/19667e90d4d1474b83c19865c8c5c033f7c3d28d/README.md?plain=1#L102C26-L102C120

# Workflow for creating a GATK CNV Panel of Normals given a list of normal samples. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
#   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
#   with panels containing only individuals of the same sex as the case samples).
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
#  A reasonable blacklist for excluded intervals (-XL) can be found at:
#   hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
#   hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.list (untested)
#
# - Example invocation:
#
#     java -jar cromwell.jar run cnv_somatic_panel_workflow.wdl -i my_parameters.json
#
#############

#import "../cnv_common_tasks.wdl" as CNVTasks

import "https://raw.githubusercontent.com/gatk-workflows/gatk4-somatic-cnvs/1.4.0/tasks/cnv_common_tasks.wdl" as CNVTasks

workflow CNVSomaticPanelWorkflow {

  ##################################
  #### required basic arguments ####
  ##################################
  File intervals
  File? blacklist_intervals
  Array[String] normal_bams
  Array[String] normal_bais
  String pon_entity_id
  File ref_fasta_dict
  File ref_fasta_fai
  File ref_fasta
  String gatk_docker

  ##################################
  #### optional basic arguments ####
  ##################################
  # If true, AnnotateIntervals will be run to create GC annotations and explicit GC correction
  # will be performed by the PoN generated by CreateReadCountPanelOfNormals before PCA is performed on subsequent cases
  Boolean? do_explicit_gc_correction
  File? gatk4_jar_override
  Int? preemptible_attempts

  ####################################################
  #### optional arguments for PreprocessIntervals ####
  ####################################################
  Int? padding
  Int? bin_length
  Int? mem_gb_for_preprocess_intervals

  ##################################################
  #### optional arguments for AnnotateIntervals ####
  ##################################################
  File? mappability_track_bed
  File? mappability_track_bed_idx
  File? segmental_duplication_track_bed
  File? segmental_duplication_track_bed_idx
  Int? feature_query_lookahead
  Int? mem_gb_for_annotate_intervals

  ##############################################
  #### optional arguments for CollectCounts ####
  ##############################################
  String? collect_counts_format
  Int? mem_gb_for_collect_counts

  ##############################################################
  #### optional arguments for CreateReadCountPanelOfNormals ####
  ##############################################################
  Float? minimum_interval_median_percentile
  Float? maximum_zeros_in_sample_percentage
  Float? maximum_zeros_in_interval_percentage
  Float? extreme_sample_median_percentile
  Boolean? do_impute_zeros
  Float? extreme_outlier_truncation_percentile
  Int? number_of_eigensamples
  Int? maximum_chunk_size
  Int? mem_gb_for_create_read_count_pon

  Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

  call CNVTasks.PreprocessIntervals {
    input:
      intervals = intervals,
      blacklist_intervals = blacklist_intervals,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_fasta_dict = ref_fasta_dict,
      padding = padding,
      bin_length = bin_length,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_preprocess_intervals,
      preemptible_attempts = preemptible_attempts
  }

  if (select_first([do_explicit_gc_correction, true])) {
    call CNVTasks.AnnotateIntervals {
      input:
        intervals = PreprocessIntervals.preprocessed_intervals,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        mappability_track_bed = mappability_track_bed,
        mappability_track_bed_idx = mappability_track_bed_idx,
        segmental_duplication_track_bed = segmental_duplication_track_bed,
        segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
        feature_query_lookahead = feature_query_lookahead,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_annotate_intervals,
        preemptible_attempts = preemptible_attempts
    }
  }

  scatter (normal_bam_and_bai in normal_bams_and_bais) {
    call CNVTasks.CollectCounts {
      input:
        intervals = PreprocessIntervals.preprocessed_intervals,
        bam = normal_bam_and_bai.left,
        bam_idx = normal_bam_and_bai.right,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        format = collect_counts_format,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_collect_counts,
        preemptible_attempts = preemptible_attempts
    }
  }

  call CreateReadCountPanelOfNormals {
    input:
      pon_entity_id = pon_entity_id,
      read_count_files = CollectCounts.counts,
      minimum_interval_median_percentile = minimum_interval_median_percentile,
      maximum_zeros_in_sample_percentage = maximum_zeros_in_sample_percentage,
      maximum_zeros_in_interval_percentage = maximum_zeros_in_interval_percentage,
      extreme_sample_median_percentile = extreme_sample_median_percentile,
      do_impute_zeros = do_impute_zeros,
      extreme_outlier_truncation_percentile = extreme_outlier_truncation_percentile,
      number_of_eigensamples = number_of_eigensamples,
      maximum_chunk_size = maximum_chunk_size,
      annotated_intervals = AnnotateIntervals.annotated_intervals,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_create_read_count_pon,
      preemptible_attempts = preemptible_attempts
  }

  output {
    File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals
    Array[String] read_counts_entity_ids = CollectCounts.entity_id
    Array[File] read_counts = CollectCounts.counts
    File read_count_pon = CreateReadCountPanelOfNormals.read_count_pon
  }
}

task CreateReadCountPanelOfNormals {
  String pon_entity_id
  Array[File] read_count_files
  Float? minimum_interval_median_percentile
  Float? maximum_zeros_in_sample_percentage
  Float? maximum_zeros_in_interval_percentage
  Float? extreme_sample_median_percentile
  Boolean? do_impute_zeros
  Float? extreme_outlier_truncation_percentile
  Int? number_of_eigensamples
  Int? maximum_chunk_size
  File? annotated_intervals   #do not perform explicit GC correction by default
  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? cpu
  Int? preemptible_attempts

  Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
  Int command_mem_mb = machine_mem_mb - 500

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

    gatk --java-options "-Xmx${command_mem_mb}m" CreateReadCountPanelOfNormals \
      --input ${sep=" --input " read_count_files} \
      --minimum-interval-median-percentile ${default="10.0" minimum_interval_median_percentile} \
      --maximum-zeros-in-sample-percentage ${default="5.0" maximum_zeros_in_sample_percentage} \
      --maximum-zeros-in-interval-percentage ${default="5.0" maximum_zeros_in_interval_percentage} \
      --extreme-sample-median-percentile ${default="2.5" extreme_sample_median_percentile} \
      --do-impute-zeros ${default="true" do_impute_zeros} \
      --extreme-outlier-truncation-percentile ${default="0.1" extreme_outlier_truncation_percentile} \
      --number-of-eigensamples ${default="20" number_of_eigensamples} \
      --maximum-chunk-size ${default="16777216" maximum_chunk_size} \
      ${"--annotated-intervals " + annotated_intervals} \
      --output ${pon_entity_id}.pon.hdf5
  >>>

  runtime {
    docker: "${gatk_docker}"
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
    cpu: select_first([cpu, 1])
    preemptible: select_first([preemptible_attempts, 2])
  }

  output {
    File read_count_pon = "${pon_entity_id}.pon.hdf5"
  }
}