version 1.0

# Model
# The Dockerfile for this portion is contained in the folder /model/. Additionally in the folder is the xgboost model, in a JSON format, under /model/final_model.json.

# Once the docker environment is run, there is an option to run the tool in Tumor-Only mode.

# MODELPATH (Path to model)
# GENELIST (List of genes from previous feature generation step)
# RMDPATH (Path to RMD from feature generation (from step 4))
# CNVPATH (Path to CNV from feature generation (from step 4))
# SNVPATH (Path to SNV from feature generation (from step 4))
# LABELPATH (Path to labels from feature generation (from step 4))
# SBSPATH (Path to SBS from feature generation (from step 4))
# PASTPATH (Path to the tcga_500.csv file)
# OUTDIR (Path to the output directory for predictions to be stored)

# **Tumor-Normal**

# To run the model of in tumor-normal mode: python model.py $MODELPATH $GENELIST $RMDPATH $CNVPATH $SNVPATH $LABELPATH $SBSPATH $PASTPATH $OUTDIR

# **Tumor-Only**

# To run the model in tumor-only mode: python model.py -t $MODELPATH $GENELIST $RMDPATH $SNVPATH $LABELPATH $SBSPATH $PASTPATH $OUTDIR

# Once the docker container is created, the model can be ran using the command: python model.py $MODELPATH $GENELIST $RMDPATH $CNVPATH $SNVPATH $LABELPATH $SBSPATH $PASTPATH $OUTDIR

task RunModel {

  input {
    File model_path = "${output_dir}/model.RData"
    File rmd_path = "${output_dir}/feature_gen_report.Rmd"
    File snv_path = "${output_dir}/snv_features.csv"
    File label_path = "${output_dir}/labels.csv"
    File sbs_path = "${output_dir}/sbs_features.csv"
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