version 1.0

struct SampleInputs {
    String sample_name
    String routine
    String run_date
    String fastqs1
    String fastqs2
    String library_name
    String platform_unit
    String sequencing_center
    String platform_name
    String output_dir
}

struct References {
    File intervals_path
    String ref_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_alt
    File ref_sa
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_amb
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
}

struct Mutect2inputs {
    String gatk_docker_version
    File realignment_index_bundle
}

struct MAFinputs {
    String vep_path
    String docker_path
}
