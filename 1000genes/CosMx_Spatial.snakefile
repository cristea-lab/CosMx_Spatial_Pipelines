configfile: "config.yaml"

def targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append(f"outputs/qc/{sample}/A_1_pre_QC_sample.rds")
	ls.append(f"outputs/qc/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds")
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

rule ST_insitutype:
    input: 
        getAtoMx_path
    params:
        min_clusters = config['min_clusters'],
        max_clusters = config['max_clusters'],
	min_nCount_RNA = config['min_nCount_RNA'],
	max_nFeature_negprobes = config['max_nFeature_negprobes'],
	insitutype_profile_ref = config['insitutype_profile_ref'],
        #outdir=lambda wildcards,input,output: os.path.abspath(os.path.dirname(output[0])),
	outdir=lambda wildcards,input,output: os.path.dirname(output[0]),
    output: 
        "outputs/qc/{sample}/A_2_post_QC_and_normalization_20_counts_per_cell.rds",
        #"outputs/insitutype/{sample}/A_1_pre_QC_sample.rds",
    	#MORE!
    shell: 
        """Rscript 2_insitutype_single_sample_processing_pipeline.R {wildcards.sample} {input[0]} {params.outdir} {params.insitutype_profile_ref} {params.min_clusters} {params.max_clusters} {params.min_nCount_RNA} {params.max_nFeature_negprobes}"""


#TODO: rule for tumor polarity and fibroblast
