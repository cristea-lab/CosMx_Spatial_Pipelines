import re

configfile: "config.yaml"

# Set output path- removing trailing "/"
output_path = config.get("output_path", "analysis").rstrip("/")
src_path = config.get("src_path", "CosMx_Spatial").rstrip("/")

def targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append(f"{output_path}/qc/{sample}/A_1_pre_QC_and_filtered_sample.rds")
        ls.append(f"{output_path}/umap/{sample}/{sample}_normScaledUMAP.rds")
	
        #ls.append(f"{output_path}/insitutype/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds")	
        ls.append(f"{output_path}/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds")
        ls.append(f"{output_path}/metrics/{sample}/metrics_summary.csv")
        ls.append(f"{output_path}/polarity/{sample}/4_1_plot_correlation_AUCell_scores.jpeg")

    #ADD report targets
    ls.extend(report_targets(wildcards))
    return ls

rule all:
    input:
        targets

def getAtoMx_path(wildcards):
    sample = wildcards.sample
    atomx_path = config['samples'][sample]
    return atomx_path

rule ST_QC:
    input:
        getAtoMx_path
    params:
        min_nCount_RNA = config['min_nCount_RNA'],
        max_nFeature_negprobes = config['max_nFeature_negprobes'],

        outdir=lambda wildcards,input,output: os.path.dirname(output[0]),
        script= f"{src_path}/scripts/1_QC_preprocessing_single_sample_pipeline.R",
    output:
        output_path + "/qc/{sample}/A_1_pre_QC_and_filtered_sample.rds",
        output_path + "/qc/{sample}/1_1_QC_probe_count_histogram.png",
        output_path + "/qc/{sample}/1_2_QC_feature_count_histogram.png",
        output_path + "/qc/{sample}/1_3_QC_negative_probe_count_histogram.png",
    shell:
        """Rscript {params.script} {wildcards.sample} {input[0]} {params.min_nCount_RNA} {params.max_nFeature_negprobes} {params.outdir}"""

rule ST_normalize_scale_umap:
    """Normalize, scale, and perform UMAP for one {sample}."""
    input:
        output_path + "/qc/{sample}/A_1_pre_QC_and_filtered_sample.rds"
    params:
        features = config['umap_features'],
        outdir = lambda wildcards, input, output: os.path.dirname(output[0]),
        script = f"{src_path}/scripts/CosMx_UMAP.R"
    output:
        seurat_obj = output_path + "/umap/{sample}/{sample}_normScaledUMAP.rds",
        umap       = output_path + "/umap/{sample}/umap.png",
        #LEN: cool trick to generate outputs based on user defined list in config
        feature_umaps = [
            output_path + "/umap/{sample}/" + feat.strip() + "_umap.png"
            for feat in config['umap_features'].split(",")
        ]
    shell:
        """Rscript {params.script} {input} {params.features} {output.seurat_obj} {params.outdir}"""

def ST_insitutype_inputFn(wildcards):
    sample = wildcards.sample
    atomx_path = config['samples'][sample]
    #rds_path = f"{output_path}/qc/{sample}/A_1_pre_QC_sample.rds"
    rds_path = f"{output_path}/umap/{sample}/{sample}_normScaledUMAP.rds"
    tmp = {'atomx_path': atomx_path, 'rds_path': rds_path}
    return tmp

rule ST_insitutype:
    input:
        unpack(ST_insitutype_inputFn)
    params:
        min_clusters = config['min_clusters'],
        max_clusters = config['max_clusters'],
        insitutype_profile_ref = config['insitutype_profile_ref'],
        #outdir=lambda wildcards,input,output: os.path.abspath(os.path.dirname(output[0])),
        outdir=lambda wildcards: f"{output_path}/insitutype/{wildcards.sample}/",
        script=f"{src_path}/scripts/2_insitutype_single_sample_processing_pipeline.R",
    output:
        #output_path + "/insitutype/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds",
        #LEN: Try remove hardcoded cluster number
        #HYP: it's inefficient but the script should just take one cluster num
        #as input and generate one instead of multiple clusters
        output_path + "/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds",
        output_path + "/insitutype/{sample}/12_clusters/2_1_flightpath_initial_cluster.png",
        output_path + "/insitutype/{sample}/12_clusters/2_2_flightpath_new_cluster.png",
        output_path + "/insitutype/{sample}/12_clusters/2_3_dotplot_marker_genes.png",
        output_path + "/insitutype/{sample}/12_clusters/2_4_barplot_level_1_clusters.png",
        output_path + "/insitutype/{sample}/12_clusters/2_5_barplot_level_2_clusters.png",
        output_path + "/insitutype/{sample}/12_clusters/2_6_barplot_level_3_clusters.png",
        output_path + "/insitutype/{sample}/12_clusters/3_1_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_1.png",
        output_path + "/insitutype/{sample}/12_clusters/3_2_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_2.png",
        output_path + "/insitutype/{sample}/12_clusters/4_1_spatial_plot.png",
        output_path + "/insitutype/{sample}/12_clusters/4_2_spatial_plot_level_1_clusters.png",
    shell:
        """Rscript {params.script} {wildcards.sample} {input.rds_path} {input.atomx_path} {params.outdir} {params.insitutype_profile_ref} {params.min_clusters} {params.max_clusters}"""

rule ST_meta:
    input:
        output_path + "/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds",
    params:
        script=f"{src_path}/scripts/3_meta_data_excel_sheet.R",
    output:
        output_path + "/metrics/{sample}/metrics_summary.csv"
    shell:
        """Rscript {params.script} {wildcards.sample} {input} {output}"""

rule ST_polarity:
    """Run tumor polarity script to probe tumor fibroblast status"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds",
    params:
        polarity_genesets= config['polarity_genesets'],
        outdir=lambda wildcards,input,output: os.path.dirname(output[0]),
        script= f"{src_path}/scripts/4_tumor_polarity_and_fibrolast_single_sample_processing_pipeline.R",
    output:
       output_path + "/polarity/{sample}/4_1_plot_correlation_AUCell_scores.jpeg",
       output_path + "/polarity/{sample}/4_2_plot_correlation_classical_basal_scores.jpeg",
       output_path + "/polarity/{sample}/4_3_plot_classical_scores_vs_subtype.jpeg",
       output_path + "/polarity/{sample}/4_4_plot_basal_scores_vs_subtype.jpeg",
       output_path + "/polarity/{sample}/4_5_plot_coexpressors_scores_vs_subtype.jpeg",
       output_path + "/polarity/{sample}/6_1_boxplots_distribution_basal_classical.jpeg",
       output_path + "/polarity/{sample}/8_1_boxplot_level_3_cluster_fibrablast_subtype_frequency.jpeg",
       output_path + "/polarity/{sample}/8_2_boxplot_level_3B_cluster_fibrablast_subtype_frequency.jpeg",
       output_path + "/polarity/{sample}/8_3_boxplot_level_2_cluster_fibrablast_subtype_frequency_cafs.jpeg",
       output_path + "/polarity/{sample}/8_4_boxplot_level_3B_cluster_fibrablast_subtype_frequency_cafs.jpeg",
       output_path + "/polarity/{sample}/plots/{sample}_FAP_CAF_AUC_Histogram.png",
       output_path + "/polarity/{sample}/plots/{sample}_Fibroblast_Subtype_Frequency.png",
    shell:
        """Rscript {params.script} {wildcards.sample} {input} {params.polarity_genesets} {params.outdir}"""

include: "./modules/report.snakefile"
