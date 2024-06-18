version 1.0

# 4. **Feature Generation**

# The Dockerfile for this portion is contained in the folder `/feature_gen/`. Additionally in the folder is the feature generation script, at `/feature_gen/feature_gen.R`.

# Once the docker environment is run, there is an option to run the tool in **Tumor-Only** mode.

# 1. INDIR (Directory containing all SNV/CNV calls)
# 2. SAMPNAMES (List of sample names from previous steps)
# 3. GENEPATH (Gene List Path)
# 4. BINPATH (Genome Bins Path)
# 5. CNVPATH (all_thresholded.by_genes.txt Path)
# 6. OUTDIR (output directory where features should be stored)

# **Tumor-Normal**

# To run the feature generation script in tumor-normal mode: `RSCRIPT feature_gen.R $INDIR $SAMPNAMES $GENEPATH $BINPATH $CNVPATH $OUTDIR`

# **Tumor-Only**

# To run the feature generation script in tumor-only mode: `RSCRIPT feature_gen.R -t $INDIR $SAMPNAMES $GENEPATH $BINPATH $OUTDIR`

task FeatureGen {
  input {
    String input_dir = "../feature_gen/indir/"
    File sample_names_files = "../feature_gen/sample_names.txt"
    File gene_list = "../feature_gen/fullGeneList.txt"
    File bins_list = "../feature_gen/genome_bins.txt"
    String output_dir = "../feature_gen/outdir/"
    String dockerPath = "../feature_gen/"
  }

  command <<<

    set -e

    Rscript feature_gen.R -t ${input_dir} ${sample_names_files} ${gene_list} ${bins_list} ${output_dir}

  >>>

  runtime {
    docker: "${dockerPath}"
  }

  output {
    File model_path = "${output_dir}/model.RData"
    File rmd_path = "${output_dir}/feature_gen_report.Rmd"
    File snv_path = "${output_dir}/snv_features.csv"
    File label_path = "${output_dir}/labels.csv"
    File sbs_path = "${output_dir}/sbs_features.csv"
  }
}