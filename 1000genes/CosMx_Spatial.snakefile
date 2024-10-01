configfile: "config.yaml"

def targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append(f"outputs/qc/{sample}/A_1_pre_QC_sample.rds")
        ls.append(f"outputs/insitutype/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds")
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
        outdir=lambda wildcards,input,output: os.path.dirname(output[0]),
    output:
        "outputs/qc/{sample}/A_1_pre_QC_sample.rds",
	#NOTE to Sam-- there should be more!
    shell:
        """Rscript 1_QC_preprocessing_single_sample_pipeline.R {wildcards.sample} {input[0]} {params.outdir}"""

def ST_insitutype_inputFn(wildcards):
    sample = wildcards.sample
    atomx_path = config['samples'][sample]
    rds_path = f"outputs/qc/{sample}/A_1_pre_QC_sample.rds"
    tmp = {'atomx_path': atomx_path, 'rds_path': rds_path}
    return tmp

rule ST_insitutype:
    input:
        unpack(ST_insitutype_inputFn)
    params:
        min_clusters = config['min_clusters'],
        max_clusters = config['max_clusters'],
        min_nCount_RNA = config['min_nCount_RNA'],
        max_nFeature_negprobes = config['max_nFeature_negprobes'],
        insitutype_profile_ref = config['insitutype_profile_ref'],
        #outdir=lambda wildcards,input,output: os.path.abspath(os.path.dirname(output[0])),
        outdir=lambda wildcards,input,output: os.path.dirname(output[0]),
    output:
        "outputs/insitutype/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds",
        #LEN: Try remove hardcoded cluster number
        "outputs/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds",
    	#MORE!
    shell:
        """Rscript 2_insitutype_single_sample_processing_pipeline.R {wildcards.sample} {input.rds_path} {input.atomx_path} {params.outdir} {params.insitutype_profile_ref} {params.min_clusters} {params.max_clusters} {params.min_nCount_RNA} {params.max_nFeature_negprobes}"""

rule ST_meta:
    input:
        "outputs/insitutype/{sample}/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds",
    output:
        "outputs/metrics/{sample}/metrics_summary.csv"
    shell:
        """Rscript 3_meta_data_excel_sheet.R {wildcards.sample} {input} {output}"""


#TODO: rule for tumor polarity and fibroblast
