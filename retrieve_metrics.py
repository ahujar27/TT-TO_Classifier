# %%
import os
import json
import pandas as pd

metrics_list_of_samples_dicts = []
# %%
column_order = [
    "analysisDate",
    "sampleId",
    "analysisName",
    "PCT_PF_READS",
    "PCT_Q30_R1",
    "PCT_Q30_R2",
    "DNA_TOTAL_PF_READS",
    "DNA_PCT_PF_UQ_READS",
    "DNA_PCT_ALIGNED_READS",
    "DNA_PCT_READ_ENRICHMENT",
    "DNA_PCT_CHIMERIC_READS",
    "DNA_PCT_USABLE_UMI_READS",
    "DNA_CONTAMINATION_SCORE",
    "DNA_CONTAMINATION_P_VALUE",
    "DNA_PCT_CONTAMINATION_EST",
    "DNA_MEAN_FAMILY_SIZE",
    "DNA_MEDIAN_INSERT_SIZE",
    "DNA_MEDIAN_EXON_COVERAGE",
    "DNA_PCT_EXON_50X",
    "DNA_PCT_EXON_100X",
    "DNA_MEAN_TARGET_COVERAGE",
    "DNA_MEDIAN_TARGET_COVERAGE",
    "DNA_PCT_TARGET_04X_MEAN",
    "DNA_PCT_TARGET_100X",
    "DNA_PCT_TARGET_250X",
    "RNA_TOTAL_PF_READS",
    "RNA_MEDIAN_INSERT_SIZE",
    "RNA_MEDIAN_CV_GENE_500X",
    "RNA_SCALED_MEDIAN_GENE_COVERAGE",
    "RNA_TOTAL_ON_TARGET_READS",
    "RNA_PCT_ON_TARGET_READS",
    "RNA_PCT_CHIMERIC_READS",
    "COVERAGE_MAD",
    "MEDIAN_BIN_COUNT_CNV_TARGET",
    "USABLE_MSI_SITES",
    "TOTAL_MSI_SITES_UNSTABLE",
    "TMB_CODING_REGION_SIZE_MB",
    "TMB_SOMATIC_CODING_VARIANTS_COUNT",
    "analysisTime",
]

# %%
stgisilon_tso_path = r"/stgisilon_rotina/tso500/"


# %%
def add_tmb_metrics(tmb_metrics: list, sample_info_dict: dict):
    TMB_CODING_REGION_SIZE_MB = ""
    TMB_SOMATIC_CODING_VARIANTS_COUNT = ""

    for tmb_metric_dict in tmb_metrics:
        tmb_metric_name = tmb_metric_dict.get("name")

        if tmb_metric_name == "CodingRegionSizeMb":
            TMB_CODING_REGION_SIZE_MB = tmb_metric_dict.get("value")

        if tmb_metric_name == "SomaticCodingVariantsCount":
            TMB_SOMATIC_CODING_VARIANTS_COUNT = tmb_metric_dict.get("value")

    sample_info_dict["TMB_CODING_REGION_SIZE_MB"] = TMB_CODING_REGION_SIZE_MB
    sample_info_dict["TMB_SOMATIC_CODING_VARIANTS_COUNT"] = (
        TMB_SOMATIC_CODING_VARIANTS_COUNT
    )


# %%
def add_msi_metric(msi_metrics: list, sample_info_dict: dict):
    for msi_metrics_dict in msi_metrics:
        msi_metric_name = msi_metrics_dict.get("name")

        if msi_metric_name == "TotalMicrosatelliteSitesUnstable":
            TOTAL_MSI_SITES_UNSTABLE = msi_metrics_dict.get("value")
            sample_info_dict.update(
                {"TOTAL_MSI_SITES_UNSTABLE": TOTAL_MSI_SITES_UNSTABLE}
            )


# %%
def add_biomarker_metrics(biomarkers_dict: dict, sample_info_dict: dict):

    tmb_dict = biomarkers_dict.get("tumorMutationalBurden")

    if tmb_dict is not None:
        tmb_metrics = tmb_dict.get("additionalMetrics")
        add_tmb_metrics(tmb_metrics, sample_info_dict)

    msi_dict = biomarkers_dict.get("microsatelliteInstability")

    if msi_dict is not None:
        msi_metrics = msi_dict.get("additionalMetrics")
        add_msi_metric(msi_metrics, sample_info_dict)


# %%
def add_sample_general_metrics(general_metrics_list: list, sample_info_dict: dict):
    PCT_PF_READS = ""
    PCT_Q30_R1 = ""
    PCT_Q30_R2 = ""

    for metric_dict in general_metrics_list:
        metric_name = metric_dict.get("name")

        if metric_name == "PCT_PF_READS":
            PCT_PF_READS = metric_dict.get("value")

        if metric_name == "PCT_Q30_R1":
            PCT_Q30_R1 = metric_dict.get("value")

        if metric_name == "PCT_Q30_R2":
            PCT_Q30_R2 = metric_dict.get("value")

    sample_info_dict["PCT_PF_READS"] = PCT_PF_READS
    sample_info_dict["PCT_Q30_R1"] = PCT_Q30_R1
    sample_info_dict["PCT_Q30_R2"] = PCT_Q30_R2


# %%
def add_dna_contamination_metrics(quality_dict: dict, sample_info_dict: dict):
    DNA_CONTAMINATION_SCORE = ""
    DNA_CONTAMINATION_P_VALUE = ""

    metrics_list_of_dicts = quality_dict.get("metrics")

    for metric_dict in metrics_list_of_dicts:
        metric_name = metric_dict.get("name")

        if metric_name == "CONTAMINATION_SCORE":
            DNA_CONTAMINATION_SCORE = metric_dict.get("value")

        if metric_name == "CONTAMINATION_P_VALUE":
            DNA_CONTAMINATION_P_VALUE = metric_dict.get("value")

    sample_info_dict["DNA_CONTAMINATION_SCORE"] = DNA_CONTAMINATION_SCORE
    sample_info_dict["DNA_CONTAMINATION_P_VALUE"] = DNA_CONTAMINATION_P_VALUE


# %%
def add_exon_metrics(quality_dict: dict, sample_info_dict: dict):
    DNA_MEDIAN_INSERT_SIZE = ""
    DNA_MEDIAN_EXON_COVERAGE = ""
    DNA_PCT_EXON_50X = ""

    metrics_list_of_dicts = quality_dict.get("metrics")

    for metric_dict in metrics_list_of_dicts:
        metric_name = metric_dict.get("name")

        if metric_name == "MEDIAN_INSERT_SIZE":
            DNA_MEDIAN_INSERT_SIZE = metric_dict.get("value")

        if metric_name == "MEDIAN_EXON_COVERAGE":
            DNA_MEDIAN_EXON_COVERAGE = metric_dict.get("value")

        if metric_name == "PCT_EXON_50X":
            DNA_PCT_EXON_50X = metric_dict.get("value")

    sample_info_dict["DNA_MEDIAN_INSERT_SIZE"] = DNA_MEDIAN_INSERT_SIZE
    sample_info_dict["DNA_MEDIAN_EXON_COVERAGE"] = DNA_MEDIAN_EXON_COVERAGE
    sample_info_dict["DNA_PCT_EXON_50X"] = DNA_PCT_EXON_50X


# %%
def add_cnv_metrics(quality_dict: dict, sample_info_dict: dict):
    COVERAGE_MAD = ""
    MEDIAN_BIN_COUNT_CNV_TARGET = ""

    metrics_list_of_dicts = quality_dict.get("metrics")

    for metric_dict in metrics_list_of_dicts:
        metric_name = metric_dict.get("name")

        if metric_name == "COVERAGE_MAD":
            COVERAGE_MAD = metric_dict.get("value")

        if metric_name == "MEDIAN_BIN_COUNT_CNV_TARGET":
            MEDIAN_BIN_COUNT_CNV_TARGET = metric_dict.get("value")

    sample_info_dict["COVERAGE_MAD"] = COVERAGE_MAD
    sample_info_dict["MEDIAN_BIN_COUNT_CNV_TARGET"] = MEDIAN_BIN_COUNT_CNV_TARGET


# %%
def add_usable_msi_metric(quality_dict: dict, sample_info_dict: dict):
    USABLE_MSI_SITES = ""

    metrics_list_of_dicts = quality_dict.get("metrics")

    for metric_dict in metrics_list_of_dicts:
        metric_name = metric_dict.get("name")

        if metric_name == "USABLE_MSI_SITES":
            USABLE_MSI_SITES = metric_dict.get("value")

    sample_info_dict["USABLE_MSI_SITES"] = USABLE_MSI_SITES


# %%
def add_rna_qc_metrics(quality_dict: dict, sample_info_dict: dict):
    RNA_MEDIAN_CV_GENE_500X = ""
    RNA_TOTAL_ON_TARGET_READS = ""
    RNA_MEDIAN_INSERT_SIZE = ""

    metrics_list_of_dicts = quality_dict.get("metrics")

    for metric_dict in metrics_list_of_dicts:
        metric_name = metric_dict.get("name")

        if metric_name == "MEDIAN_CV_GENE_500X":
            RNA_MEDIAN_CV_GENE_500X = metric_dict.get("value")

        if metric_name == "TOTAL_ON_TARGET_READS":
            RNA_TOTAL_ON_TARGET_READS = metric_dict.get("value")

        if metric_name == "MEDIAN_INSERT_SIZE":
            RNA_MEDIAN_INSERT_SIZE = metric_dict.get("value")

    sample_info_dict["RNA_MEDIAN_CV_GENE_500X"] = RNA_MEDIAN_CV_GENE_500X
    sample_info_dict["RNA_TOTAL_ON_TARGET_READS"] = RNA_TOTAL_ON_TARGET_READS
    sample_info_dict["RNA_MEDIAN_INSERT_SIZE"] = RNA_MEDIAN_INSERT_SIZE


# %%
def add_quality_metrics(quality_control_metrics_list: list, sample_info_dict: dict):

    for quality_dict in quality_control_metrics_list:
        quality_dict_name = quality_dict.get("name")

        if quality_dict_name == "DNA Library QC Metrics":
            add_dna_contamination_metrics(quality_dict, sample_info_dict)

        if (
            quality_dict_name
            == "DNA Library QC Metrics for Small Variant Calling and TMB"
        ):
            add_exon_metrics(quality_dict, sample_info_dict)

        if (
            quality_dict_name
            == "DNA Library QC Metrics for Copy Number Variant Calling"
        ):
            add_cnv_metrics(quality_dict, sample_info_dict)

        if quality_dict_name == "DNA Library QC Metrics for MSI":
            add_usable_msi_metric(quality_dict, sample_info_dict)

        if quality_dict_name == "RNA Library QC Metrics":
            add_rna_qc_metrics(quality_dict, sample_info_dict)


# %%
def add_expanded_dna_metrics(dna_expanded_metrics_list: list, sample_info_dict: dict):
    DNA_TOTAL_PF_READS = ""
    DNA_MEAN_FAMILY_SIZE = ""
    DNA_MEDIAN_TARGET_COVERAGE = ""
    DNA_PCT_CHIMERIC_READS = ""
    DNA_PCT_EXON_100X = ""
    DNA_PCT_READ_ENRICHMENT = ""
    DNA_PCT_USABLE_UMI_READS = ""
    DNA_MEAN_TARGET_COVERAGE = ""
    DNA_PCT_ALIGNED_READS = ""
    DNA_PCT_CONTAMINATION_EST = ""
    DNA_PCT_PF_UQ_READS = ""
    DNA_PCT_TARGET_04X_MEAN = ""
    DNA_PCT_TARGET_100X = ""
    DNA_PCT_TARGET_250X = ""

    for dna_expanded_mectrics_dict in dna_expanded_metrics_list:
        dna_metric_name = dna_expanded_mectrics_dict.get("name")

        if dna_metric_name == "TOTAL_PF_READS":
            DNA_TOTAL_PF_READS = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "MEAN_FAMILY_SIZE":
            DNA_MEAN_FAMILY_SIZE = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "MEDIAN_TARGET_COVERAGE":
            DNA_MEDIAN_TARGET_COVERAGE = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_CHIMERIC_READS":
            DNA_PCT_CHIMERIC_READS = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_EXON_100X":
            DNA_PCT_EXON_100X = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_READ_ENRICHMENT":
            DNA_PCT_READ_ENRICHMENT = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_USABLE_UMI_READS":
            DNA_PCT_USABLE_UMI_READS = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "MEAN_TARGET_COVERAGE":
            DNA_MEAN_TARGET_COVERAGE = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_ALIGNED_READS":
            DNA_PCT_ALIGNED_READS = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_CONTAMINATION_EST":
            DNA_PCT_CONTAMINATION_EST = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_PF_UQ_READS":
            DNA_PCT_PF_UQ_READS = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_TARGET_0.4X_MEAN":
            DNA_PCT_TARGET_04X_MEAN = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_TARGET_100X":
            DNA_PCT_TARGET_100X = dna_expanded_mectrics_dict.get("value")

        if dna_metric_name == "PCT_TARGET_250X":
            DNA_PCT_TARGET_250X = dna_expanded_mectrics_dict.get("value")

    sample_info_dict["DNA_TOTAL_PF_READS"] = DNA_TOTAL_PF_READS
    sample_info_dict["DNA_MEAN_FAMILY_SIZE"] = DNA_MEAN_FAMILY_SIZE
    sample_info_dict["DNA_MEDIAN_TARGET_COVERAGE"] = DNA_MEDIAN_TARGET_COVERAGE
    sample_info_dict["DNA_PCT_CHIMERIC_READS"] = DNA_PCT_CHIMERIC_READS
    sample_info_dict["DNA_PCT_EXON_100X"] = DNA_PCT_EXON_100X
    sample_info_dict["DNA_PCT_READ_ENRICHMENT"] = DNA_PCT_READ_ENRICHMENT
    sample_info_dict["DNA_PCT_USABLE_UMI_READS"] = DNA_PCT_USABLE_UMI_READS
    sample_info_dict["DNA_MEAN_TARGET_COVERAGE"] = DNA_MEAN_TARGET_COVERAGE
    sample_info_dict["DNA_PCT_ALIGNED_READS"] = DNA_PCT_ALIGNED_READS
    sample_info_dict["DNA_PCT_CONTAMINATION_EST"] = DNA_PCT_CONTAMINATION_EST
    sample_info_dict["DNA_PCT_PF_UQ_READS"] = DNA_PCT_PF_UQ_READS
    sample_info_dict["DNA_PCT_TARGET_04X_MEAN"] = DNA_PCT_TARGET_04X_MEAN
    sample_info_dict["DNA_PCT_TARGET_100X"] = DNA_PCT_TARGET_100X
    sample_info_dict["DNA_PCT_TARGET_250X"] = DNA_PCT_TARGET_250X


# %%
def add_rna_expanded_metrics(rna_expanded_metrics_list: list, sample_info_dict: dict):
    RNA_PCT_CHIMERIC_READS = ""
    RNA_PCT_ON_TARGET_READS = ""
    RNA_SCALED_MEDIAN_GENE_COVERAGE = ""
    RNA_TOTAL_PF_READS = ""

    for rna_expanded_mectrics_dict in rna_expanded_metrics_list:
        rna_metric_name = rna_expanded_mectrics_dict.get("name")

        if rna_metric_name == "PCT_CHIMERIC_READS":
            RNA_PCT_CHIMERIC_READS = rna_expanded_mectrics_dict.get("value")

        if rna_metric_name == "PCT_ON_TARGET_READS":
            RNA_PCT_ON_TARGET_READS = rna_expanded_mectrics_dict.get("value")

        if rna_metric_name == "SCALED_MEDIAN_GENE_COVERAGE":
            RNA_SCALED_MEDIAN_GENE_COVERAGE = rna_expanded_mectrics_dict.get("value")

        if rna_metric_name == "TOTAL_PF_READS":
            RNA_TOTAL_PF_READS = rna_expanded_mectrics_dict.get("value")

    sample_info_dict["RNA_PCT_CHIMERIC_READS"] = RNA_PCT_CHIMERIC_READS
    sample_info_dict["RNA_PCT_ON_TARGET_READS"] = RNA_PCT_ON_TARGET_READS
    sample_info_dict["RNA_SCALED_MEDIAN_GENE_COVERAGE"] = (
        RNA_SCALED_MEDIAN_GENE_COVERAGE
    )
    sample_info_dict["RNA_TOTAL_PF_READS"] = RNA_TOTAL_PF_READS


# %%
def add_expanded_metrics(expanded_metrics_list: list, sample_info_dict: dict):

    for expanded_metrics_dict in expanded_metrics_list:
        dict_sample_type = expanded_metrics_dict.get("name")

        if dict_sample_type == "DNA Expanded Metrics":
            dna_expanded_metrics_list = expanded_metrics_dict.get("metrics")
            add_expanded_dna_metrics(dna_expanded_metrics_list, sample_info_dict)

        if dict_sample_type == "RNA Expanded Metrics":
            rna_expanded_metrics_list = expanded_metrics_dict.get("metrics")
            add_rna_expanded_metrics(rna_expanded_metrics_list, sample_info_dict)


# %%
def add_metrics_to_sample_info(metrics_data: dict, sample_info_dict: dict):

    biomarkers_dict = metrics_data.get("biomarkers")
    add_biomarker_metrics(biomarkers_dict, sample_info_dict)

    general_metrics_dict = metrics_data.get("sampleMetrics")
    run_metrics = general_metrics_dict.get("run")

    if run_metrics is not None:

        general_metrics_list = run_metrics.get("metrics")

        if general_metrics_list is not None:
            add_sample_general_metrics(general_metrics_list, sample_info_dict)

            quality_control_metrics_list = general_metrics_dict.get(
                "qualityControlMetrics"
            )
            add_quality_metrics(quality_control_metrics_list, sample_info_dict)

            expanded_metrics_list = general_metrics_dict.get("expandedMetrics")
            add_expanded_metrics(expanded_metrics_list, sample_info_dict)


# %%
def retrieve_samples_metrics(
    stgisilon_tso_path: str, metrics_list_of_samples_dicts: list
) -> list:

    stgisilon_content = os.listdir(stgisilon_tso_path)

    for routine in stgisilon_content:
        full_path = stgisilon_tso_path + routine

        if os.path.isdir(full_path):
            path_to_metrics = f"{full_path}/Logs_Intermediates/SampleAnalysisResults"

            if os.path.isdir(path_to_metrics):
                samples_analysis_files = os.listdir(path_to_metrics)

                for file in samples_analysis_files:
                    if not file.startswith("dsdm") and file.endswith(r".json"):
                        with open(f"{path_to_metrics}/{file}") as sample_metrics_json:
                            metrics_dict = json.load(sample_metrics_json)
                            metrics_data = metrics_dict.get("data")

                            sample_info_dict = metrics_data.get("sampleInformation")
                            add_metrics_to_sample_info(metrics_data, sample_info_dict)

                            metrics_list_of_samples_dicts.append(sample_info_dict)

    return metrics_list_of_samples_dicts


# %%
metrics_list_of_samples_dicts = retrieve_samples_metrics(
    stgisilon_tso_path, metrics_list_of_samples_dicts
)
metrics_df = pd.DataFrame(metrics_list_of_samples_dicts)
metrics_df = metrics_df[column_order]
metrics_df_sorted = metrics_df.sort_values(
    by=["analysisDate", "analysisName", "sampleId"]
)
metrics_df_sorted.to_csv("metrics_tso.tsv", index=False, sep="\t")
