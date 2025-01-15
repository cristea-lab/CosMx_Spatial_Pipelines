#report module
from yaml import dump as yaml_dump

def report_targets_sansHTML(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append(f"{output_path}/report/{sample}/01_QC_probe_count.png")
        ls.append(f"{output_path}/report/{sample}/02_QC_feature_count.png")
        ls.append(f"{output_path}/report/{sample}/03_QC_neg_probe_count.png")
        ls.append(f"{output_path}/report/{sample}/04_umap.png")
        ls.append(f"{output_path}/report/{sample}/05_umap_features.tsv")
        ls.append(f"{output_path}/report/{sample}/06_initial_flightpath.png")
        ls.append(f"{output_path}/report/{sample}/07_new_flightpath.png")
        ls.append(f"{output_path}/report/{sample}/08_marker_genes.png")
        ls.append(f"{output_path}/report/{sample}/09_level1_clusters.png")
        ls.append(f"{output_path}/report/{sample}/10_level2_clusters.png")
        ls.append(f"{output_path}/report/{sample}/11_level3_clusters.png")
        ls.append(f"{output_path}/report/{sample}/12_violin_plot_1.png")
        ls.append(f"{output_path}/report/{sample}/13_violin_plot_2.png")
        ls.append(f"{output_path}/report/{sample}/14_spatial_plot.png")
        ls.append(f"{output_path}/report/{sample}/15_spatial_plot_level_1.png")
        
        ls.append(f"{output_path}/report/{sample}/16_cell_type_metrics.csv")
    return ls

def report_targets(wildcards):
    ls = report_targets_sansHTML(wildcards)
    #REPORT
    ls.append(f"{output_path}/report/report.html")
    ls.append(f"{output_path}/report.tar.gz")
    return ls

rule report_all:
    input:
        report_targets
    #benchmark: "benchmarks/report/report_all.txt"

###############################################################################
# QC plots
###############################################################################
rule report_probe_count:
    """Generate the probe count histogram for the report"""
    input:
        output_path + "/qc/{sample}/1_1_QC_probe_count_histogram.png"
    output:
        png= output_path + "/report/{sample}/01_QC_probe_count.png",
        details=output_path + "/report/{sample}/01_details.yaml",
    params:
        caption="""caption: 'This plot shows the probe counts per cell.'"""
    message:
        "REPORT: creating probe count histogram"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_feature_count:
    """Generate the feature count histogram for the report"""
    input:
        output_path + "/qc/{sample}/1_2_QC_feature_count_histogram.png"
    output:
        png= output_path + "/report/{sample}/02_QC_feature_count.png",
        details=output_path + "/report/{sample}/02_details.yaml",
    params:
        caption="""caption: 'This plot shows the feature counts per cell.'"""
    message:
        "REPORT: creating feature count histogram"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_neg_probe_count:
    """Generate the feature count histogram for the report"""
    input:
        output_path + "/qc/{sample}/1_3_QC_negative_probe_count_histogram.png"
    output:
        png= output_path + "/report/{sample}/03_QC_neg_probe_count.png",
        details=output_path + "/report/{sample}/03_details.yaml",
    params:
        caption="""caption: 'This plot shows the negative probe counts per cell.'"""
    message:
        "REPORT: creating negative probe count histogram"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

###############################################################################
# UMAP
###############################################################################
rule report_umap_plot:
    input:
        output_path + "/umap/{sample}/umap.png"
    output:
        png= output_path + "/report/{sample}/04_umap.png",
        details=output_path + "/report/{sample}/04_details.yaml",
    params:
        caption="""caption: 'UMAP plot.'"""
    message:
        "REPORT: creating UMAP plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

#TODO: integrate umap_feature plots
rule report_umap_feature:
    input:
        expand(output_path + "/umap/{{sample}}/{feature}_umap.png", feature=[f.strip() for f in config['umap_features'].split(",")])
    output:
        tsv=output_path + "/report/{sample}/05_umap_features.tsv",
        details=output_path + "/report/{sample}/05_details.yaml",
    params:
        caption="""caption: 'UMAP features plots.'""",
        #NOTE: the plots dir is relative to the reports dir
        report_dir = output_path + "/report/",
        plots_dir = lambda wildcards: f"{wildcards.sample}/plots",
        script = f"{src_path}/scripts/report_makeUMAPfeature.py",
        png_files = lambda wildcards, input: " -f ".join(input),
    message:
        "REPORT: creating UMAP feature plots"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && {params.script} -f {params.png_files} -r {params.report_dir} -p {params.plots_dir} -o {output.tsv}"""

###############################################################################
# Insitutype
###############################################################################
rule report_init_flighpath:
    """Generate the initial flightpath for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_1_flightpath_initial_cluster.png"
    output:
        png= output_path + "/report/{sample}/06_initial_flightpath.png",
        details=output_path + "/report/{sample}/06_details.yaml",
    params:
        caption="""caption: 'Initial flightpath plot assigning cells to cell types based on their posterior probabilities.'"""
    message:
        "REPORT: creating initial flightpath plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_new_flighpath:
    """Generate the initial flightpath for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_2_flightpath_new_cluster.png"
    output:
        png= output_path + "/report/{sample}/07_new_flightpath.png",
        details=output_path + "/report/{sample}/07_details.yaml",
    params:
        caption="""caption: 'New flightpath plot assigning cells to cell types based on their posterior probabilities.'"""
    message:
        "REPORT: creating new flightpath plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_marker_genes:
    """Generate the marker gene dotplots for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_3_dotplot_marker_genes.png"
    output:
        png= output_path + "/report/{sample}/08_marker_genes.png",
        details=output_path + "/report/{sample}/08_details.yaml",
    params:
        caption="""caption: 'Dot plot of marker genes.'"""
    message:
        "REPORT: creating marker genes plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_level1_clusters:
    """Generate the level1 clusters for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_4_barplot_level_1_clusters.png"
    output:
        png= output_path + "/report/{sample}/09_level1_clusters.png",
        details=output_path + "/report/{sample}/09_details.yaml",
    params:
        caption="""caption: 'Level 1 clusters.'"""
    message:
        "REPORT: creating level 1 cluster plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_level2_clusters:
    """Generate the level2 clusters for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_5_barplot_level_2_clusters.png"
    output:
        png= output_path + "/report/{sample}/10_level2_clusters.png",
        details=output_path + "/report/{sample}/10_details.yaml",
    params:
        caption="""caption: 'Level 2 clusters.'"""
    message:
        "REPORT: creating level 2 cluster plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_level3_clusters:
    """Generate the level3 clusters for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/2_6_barplot_level_3_clusters.png"
    output:
        png= output_path + "/report/{sample}/11_level3_clusters.png",
        details=output_path + "/report/{sample}/11_details.yaml",
    params:
        caption="""caption: 'Level 3 clusters.'"""
    message:
        "REPORT: creating level 3 cluster plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_violin_plot1:
    """Generate the violin plot1 for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/3_1_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_1.png"
    output:
        png= output_path + "/report/{sample}/12_violin_plot_1.png",
        details=output_path + "/report/{sample}/12_details.yaml",
    params:
        caption="""caption: 'Volin plot of marker scores stratified by collapsed cell type 1.'"""
    message:
        "REPORT: creating violin plot 1"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_violin_plot2:
    """Generate the violin plot2 for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/3_2_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_2.png"
    output:
        png= output_path + "/report/{sample}/13_violin_plot_2.png",
        details=output_path + "/report/{sample}/13_details.yaml",
    params:
        caption="""caption: 'Volin plot of marker scores stratified by collapsed cell type 2.'"""
    message:
        "REPORT: creating violin plot 2"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_spatial_plot:
    """Generate the spatial plot for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/4_1_spatial_plot.png"
    output:
        png= output_path + "/report/{sample}/14_spatial_plot.png",
        details=output_path + "/report/{sample}/14_details.yaml",
    params:
        caption="""caption: 'Spatial plot.'"""
    message:
        "REPORT: creating spatial plot"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

rule report_spatial_plot2:
    """Generate the spatial plot 2 for the report"""
    input:
        output_path + "/insitutype/{sample}/12_clusters/4_2_spatial_plot_level_1_clusters.png"
    output:
        png= output_path + "/report/{sample}/15_spatial_plot_level_1.png",
        details=output_path + "/report/{sample}/15_details.yaml",
    params:
        caption="""caption: 'Spatial plot of level 1 clusters.'"""
    message:
        "REPORT: creating spatial plot of level 1 clusters"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cp {input} {output.png}"""

###############################################################################
# META
###############################################################################
rule report_metrics:
    """Generate the cell type metrics for the report"""
    input:
        output_path + "/metrics/{sample}/metrics_summary.csv"
    output:
        csv= output_path + "/report/{sample}/16_cell_type_metrics.csv",
        details=output_path + "/report/{sample}/16_details.yaml",
    params:
        caption="""caption: 'Cell type metrics.'""",
        #AWK cmd to remove double-quotes from the file
        awk_cmd = "\'{gsub(/\"/, \"\"); print}\'"
    message:
        "REPORT: creating cell type metrics"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && awk {params.awk_cmd} {input} > {output.csv}"""

###############################################################################

def getSections():
    section_str = ",".join(config['samples'])
    #print(section_str)
    return section_str

rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_targets_sansHTML
    params:
        jinja2_template= src_path + "/report/index.sample.html",
        report_path = output_path + "/report",
        sections_list = getSections(),
        title="Summary Report",
    output:
        output_path + "/report/report.html"
    message:
        "REPORT: Generating report"
    group: "report"
    shell:
        src_path + """/scripts/report.py -d {params.report_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r """ + src_path +"""/report/static {params.report_path}"""

rule report_gzipReport:
    input:
        output_path + "/report/report.html"
    params:
        report_path = output_path +"/report"
    output:
        output_path + "/report.tar.gz"
    message: "REPORT: Zipping up report directory"
    group: "report"
    shell:
        "tar -c {params.report_path} | gzip > {output}"
