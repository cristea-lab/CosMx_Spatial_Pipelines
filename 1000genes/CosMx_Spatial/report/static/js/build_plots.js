// This script contains functions for building and displaying plotly.js plots
// The functions tend to follow a similar foundational structure:
// Loop through available samples, check if these samples aren't filtered out, populate data structures used in plotly plot
// The functions vary in complexity, some of which use danfo.js to transform and extract data from maf tables
// The build_mapping_plot() is an example of a more basic function which lays the groundwork for others and contains comments explaining basic steps

// check if elemennt exists in array
function checkAvailability(arr, val) {
    return arr.some(function (arrVal) {
        return val === arrVal;
    });
}

function build_data_quality(){
    build_mapping_plot();
    build_coverage_plot();
    build_mean_quality_plot();
    build_gc_content_plot();
    build_insert_size_plot();
}

function build_copy_number_variation(){
    build_clonality_plot();
    build_purity_plot();
    build_ploidy_plot();
}

function build_somatic_variants(){
    build_sv_summary_plots();
    build_ti_tv_plot();
    build_lego_plot();
    build_tmb_plot();
    build_oncoplot_plot();
}

function build_HLA(){
    build_hla_oncoplot_plot();
}

function build_MISC(){
    build_tcellextrect_plot();
    build_msisensor2_plot();
}


// build somatic variant summary plots
function build_sv_summary_plots(){
    build_variant_classification_plot();
    build_variant_type_plot();
    build_snv_class_plot();
    build_variant_per_sample_plot();
    build_variant_classification_summary_plot();
    build_top_10_genes_plot();
}
// sort by elements in array by frequency, remove duplicates
function sortByFrequency(array) {
    var frequency = {};

    array.forEach(function(value) { frequency[value] = 0; });

    var uniques = array.filter(function(value) {
        return ++frequency[value] == 1;
    });

    return uniques.sort(function(a, b) {
        return frequency[b] - frequency[a];
    });
}
// variant classification types
var vc_types = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Splice_Site','In_Frame_Del','Frame_Shift_Ins','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation'];
var bases = ['T','C','G','A']
// snv class conversion from maftools
var snv_conv = {
    'A>G':'T>C',
    'T>C':'T>C',
    'C>T':'C>T',
    'G>A':'C>T',
    'A>T':'T>A',
    'T>A':'T>A',
    'A>C':'T>G',
    'T>G':'T>G',
    'C>A':'C>A',
    'G>T':'C>A',
    'C>G':'C>G',
    'G>C':'C>G'
  };
// snv to ti-tv conversion from maftools
var ti_tv_conv = {
    "T>C":"Ti",
    "C>T":"Ti",
    "T>A":"Tv",
    "T>G":"Tv",
    "C>A":"Tv",
    "C>G":"Tv"
}

// find indices of top values
function findIndicesOfMax(inp, count) {
    var outp = new Array();
    for (var i = 0; i < inp.length; i++) {
        outp.push(i);
        if (outp.length > count) {
            outp.sort(function(a, b) { return inp[b] - inp[a]; });
            outp.pop();
        }
    }
    return outp;
}
// sum all values by key in an object
function sumObjectsByKey(...objs) {
    return objs.reduce((a, b) => {
      for (let k in b) {
        if (b.hasOwnProperty(k))
          a[k] = (a[k] || 0) + b[k];
      }
      return a;
    }, {});
  }

// Data quality plots

var sample_type = ['normal', 'tumor'];
var sample_type_suffix = ['.N', '.T'];

function build_mapping_plot() {

    let total_reads = [];
    let mapped_reads = [];
    let dedup_reads = [];
    let sample_ids = [];
// for each sample in wes_data...
    for (let i = 0; i < wes_data.length; ++i) {
        // check if sample hasn't been filtered out by user
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            // for normal and tumor type
            for (let j = 0; j < sample_type.length; ++j) {
                // check if type exists
                if (wes_data[i][sample_type[j]] != undefined) {
                    // push counts from different types of reads to arrays for plot
                    total_reads.push(wes_data[i][sample_type[j]]['alignment']['total_reads'])
                    mapped_reads.push(wes_data[i][sample_type[j]]['alignment']['mapped_reads'])
                    dedup_reads.push(wes_data[i][sample_type[j]]['alignment']['dedup_reads'])
                    // push sample (+type) name for x-axis of plot
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }
    // total reads trace
    let trace1 = {
        x: sample_ids,
        y: total_reads,
        name: 'Total Reads',
        type: 'bar'
    };
    // mapped reads trace
    let trace2 = {
        x: sample_ids,
        y: mapped_reads,
        name: 'Mapped Reads',
        type: 'bar'
    };
    // dedup reads trace
    let trace3 = {
        x: sample_ids,
        y: dedup_reads,
        name: 'Dedup Reads',
        type: 'bar'
    };
    // array of different traces for plot
    let data = [trace1, trace2, trace3];
    // layout object for defining appearance
    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Reads' },
            automargin: true
        },
    };
    // plot to div name "mapping_plot"
    Plotly.newPlot("mapping_plot", data, layout);

};

function build_coverage_plot() {

    let mean_depth = [];
    let sample_ids = [];

    let hover_text = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // read depth coverage
                    mean_depth.push(wes_data[i][sample_type[j]]['coverage']['mean_depth'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                    // extra information for hover
                    hover_text.push(
                        'Sample ID: ' + wes_data[i]['id'] + sample_type_suffix[j] + '<br>' +
                        'Mean Depth: ' + wes_data[i][sample_type[j]]['coverage']['mean_depth'] + '<br>' +
                        'Median Depth: ' + wes_data[i][sample_type[j]]['coverage']['median_depth'] + '<br>' +
                        'Q1 Depth: ' + wes_data[i][sample_type[j]]['coverage']['q1_depth'] + '<br>' +
                        'Q3 Depth: ' + wes_data[i][sample_type[j]]['coverage']['q3_depth'] + '<br>' +
                        '% Bases > 50: ' + wes_data[i][sample_type[j]]['coverage']['percent_bases_gt_50']
                    )
                }
            }
        }
    }

    let trace1 = {
        y: sample_ids,
        x: mean_depth,
        name: 'Mean Depth',
        type: 'bar',
        orientation: 'h',
        hovertemplate: '%{text}<extra></extra>',
        text: hover_text,
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Mean Depth' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
    };

    Plotly.newPlot("coverage_plot", data, layout);

};

function build_gc_content_plot() {

    let gc_content = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // gc content
                    gc_content.push({
                        y: wes_data[i][sample_type[j]]['alignment']['gc_content'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    let data = gc_content;

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        hovermode: 'closest',
        xaxis: {
            title: { text: '% GC bases' },
            automargin: true
        },
        yaxis: {
            title: { text: 'GC Content' },
            automargin: true
        },
        legend: { title: { "text": "Sample" } }
    };

    Plotly.newPlot("GC_Content_plot", data, layout);

};

function build_insert_size_plot() {

    let insert_size = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // insert size
                    insert_size.push({
                        y: wes_data[i][sample_type[j]]['alignment']['insert_size'],
                        mode: 'lines',
                        type: 'scatter',
                        name: wes_data[i]['id'] + sample_type_suffix[j]
                    })
                }
            }
        }
    }

    let data = insert_size;

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        hovermode: 'closest',
        xaxis: {
            title: { text: 'Insert Size' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Counts' },
            automargin: true
        },
        legend: { title: { "text": "Sample" } }
    };

    Plotly.newPlot("insert_size_plot", data, layout);

};



function build_mean_quality_plot() {

    let mean_quality = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (wes_data[i][sample_type[j]] != undefined) {
                    // mean quality score
                    mean_quality.push(wes_data[i][sample_type[j]]['alignment']['mean_quality_score'])
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }

    let trace1 = {
        x: sample_ids,
        y: mean_quality,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Mean Quality Score' },
            automargin: true
        },
    };

    Plotly.newPlot("mean_quality_plot", data, layout);

};

// CNV plots

function build_clonality_plot() {
    let traces = []; //cluster cell prevalence matrix, index 0 ~ cluster0 ccfs
    let sample_ids = [];
    let sample_indices = [];
    let max_clusters = 0;


    //first loop to figure out max clusters
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            // we don't differentiate between tumor and normal samples here
            if (wes_data[i]['copy_number'] != undefined) {
                // clonality
                sample_indices.push(i);
                sample_ids.push(wes_data[i]['id']);
                if (wes_data[i]['copy_number']['clonality']['cellular_prevalences'].length > max_clusters) {
                    max_clusters = wes_data[i]['copy_number']['clonality']['cellular_prevalences'].length;
                }
            }
        }
    }

    //Second loop to init traces
    for (let i = 0; i < max_clusters; i++) {
        traces.push(new Array(sample_ids.length));
    }
    //Add CCFs to build up traces
    for (let i = 0; i < sample_indices.length; i++) {
        for (let j = 0; j < wes_data[sample_indices[i]]['copy_number']['clonality']['cellular_prevalences'].length; j++) {
            traces[j][i] = wes_data[sample_indices[i]]['copy_number']['clonality']['cellular_prevalences'][j];
        }
    }

    let data = [];
    for (let i = 0; i < traces.length; i++) {
        data.push({name: 'cluster'+i, x: sample_ids, y: traces[i], type:'bar'})
    }
    let layout = {
        height: 450,
        width: 1044,
        barmode: 'stack',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Clonality' },
            automargin: true
        },
    };

    Plotly.newPlot("clonality_plot", data, layout);
};

function build_purity_plot() {

    let purity = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                // purity
                purity.push(wes_data[i]['copy_number']['purity'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    let trace1 = {
        x: sample_ids,
        y: purity,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Purity' },
            automargin: true
        },
    };

    Plotly.newPlot("purity_plot", data, layout);

};

function build_ploidy_plot() {

    let ploidy = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['copy_number'] != undefined) {
                // ploidy
                ploidy.push(wes_data[i]['copy_number']['ploidy'])
                sample_ids.push(wes_data[i]['id'])
            }
        }
    }
    let trace1 = {
        x: sample_ids,
        y: ploidy,
        mode: 'lines',
        type: 'scatter',
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Run' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Ploidy' },
            automargin: true
        },
    };

    Plotly.newPlot("ploidy_plot", data, layout);

};

// somatic variants plots

function build_variant_classification_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                // generate list of indexes in variant_classificaiton column which have values contained in vc_types array
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                // utilizing danfo.js to keep only values in vc_types
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Classification']);
            }
        }
    }
    // variant classification count
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
        transforms: [{
            type: 'sort',
            target: 'y',
            order: 'descending'
          }]
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        title: 'Variant Classificaiton',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Classification' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("variant_class_plot", data, layout);

}

function build_variant_type_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                // generate list of indexes in variant_classificaiton column which have values contained in vc_types array
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                // utilizing danfo.js to keep only values in vc_types
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Type']);
            }
        }
    }
    // variant type count
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
        transforms: [{
            type: 'sort',
            target: 'y',
            order: 'descending'
          }]
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        title: 'Variant Type',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Type' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("variant_type_plot", data, layout);

}

function build_snv_class_plot(){

    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                let ref = wes_data[i]['somatic']['maf']['Reference_Allele']['data'];
                let allele2 = wes_data[i]['somatic']['maf']['Tumor_Seq_Allele2']['data'];
                let class_list = [];
                // reference_allele > tumor_seq_allele2
                for (let j=0;j<ref.length; ++j) {
                    // keep bases only from bases array
                    if (checkAvailability(bases, ref[j]) && checkAvailability(bases, allele2[j])) {
                        // use snv_conv function to convert bases
                        class_list.push(snv_conv[ref[j]+'>'+allele2[j]])
                    }
                }
                df_list.push(new dfd.Series(class_list));
            }
        }
    }
    let counts = dfd.concat({ df_list: df_list, axis: 0 }).value_counts()

    let trace1 = {
        x: counts['index_arr'],
        y: counts['data'],
        type: 'bar',
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        title: 'SNV Class',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'SNV Class' },
            automargin: true,
            categoryorder: "array",
            categoryarray:  ['T>G','T>A','T>C','C>T','C>G','C>A']
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("SNV_Class_plot", data, layout);

}


function build_variant_per_sample_plot(){

    let sample_ids = [];
    let total_counts = [];
    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], x: sample_ids, type: 'bar',name: vc_types[v]}));

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                // generate variant counts
                vc = wes_data[i]['somatic']['maf']['Variant_Classification'].value_counts();
                let total_count = 0
                // for each element in vc_types array, push variant count from sample
                for (let j = 0; j < vc_types.length; ++j) {
                    data[j]['y'].push(vc['data'][vc['index_arr'].indexOf(vc_types[j])]);
                    // add to total count
                    if (typeof vc['data'][vc['index_arr'].indexOf(vc_types[j])] != 'undefined'){
                        total_count = total_count + vc['data'][vc['index_arr'].indexOf(vc_types[j])];
                    }
                }
                // sample by sample
                sample_ids.push(wes_data[i]['id']);
                total_counts.push(total_count);
            }
        }
    }
    // sort sample ids by total counts
    let indices = [...total_counts.keys()].sort((a, b) => total_counts[b] - total_counts[a]);
    let sample_ids_reordered = [sample_ids].map(a => indices.map(i => a[i]))[0];

    let layout = {
        height: 450,
        width: 1044,
        title: 'Variants Per Sample',
        barmode: 'stack',
        xaxis: {
            title: { text: 'Run' },
            automargin: true,
            categoryorder: 'array',
            categoryarray: sample_ids_reordered
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };
    Plotly.newPlot("variants_per_sample_plot", data, layout);
}


function build_variant_classification_summary_plot(){

    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], type: 'box',name: vc_types[v]}));

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                // vc count per sample
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                // indexes of variants from sample which are also in vc_types
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                // variant type counts
                let vc_counts = wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)})['Variant_Classification'].value_counts();
                // push counts to data structure for plotly
                for (let j = 0; j < vc_types.length; ++j) {
                    data[j]['y'].push(vc_counts['data'][vc_counts['index_arr'].indexOf(vc_types[j])]);
                }
            }
        }
    }

    let layout = {
        height: 450,
        width: 1044,
        title: 'Variant Classification Summary',
        barmode: 'overlay',
        xaxis: {
            title: { text: 'Variant Classification' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
        //legend: {traceorder: 'reversed'},
        showlegend: false
    };

    Plotly.newPlot("variant_class_summary_plot", data, layout);

}

function build_top_10_genes_plot(){
    // keep entries pertaining to vc_list
    let df_list = []
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                let vc_list = wes_data[i]['somatic']['maf']['Variant_Classification']['data'];
                let keep_idx = [];
                for (let j=0;j<vc_list.length; ++j) {
                    if (checkAvailability(vc_types, vc_list[j])) {
                        keep_idx.push(j)
                    }
                }
                df_list.push(wes_data[i]['somatic']['maf'].iloc({rows: Array.from(keep_idx)}));
            }
        }
    }
    // pool all samples with danfo, count top 10 gene occurances
    let combined_df = dfd.concat({ df_list: df_list, axis: 0 })
    let gene_counts = combined_df['Hugo_Symbol'].value_counts()
    let keep_idx = findIndicesOfMax(gene_counts['data'], 20);
    let top_10_genes = keep_idx.map(i => gene_counts['index_arr'][i]);

    let combined_df_hs = combined_df['Hugo_Symbol']['data'];
    df_list = [];
    // keep entries pertaining to top_10_genes
    for (let i = 0; i < top_10_genes.length; ++i) {
        let keep_idx = [];
        for (let j=0;j<combined_df_hs.length; ++j) {
            if (combined_df_hs[j] == top_10_genes[i]) {
                keep_idx.push(j)
            }
        }
        df_list.push(combined_df.iloc({rows: Array.from(keep_idx)}));
    }

    let data = Array(vc_types.length).fill().map((u,v) => ({y: [], x: top_10_genes, type: 'bar',name: vc_types[v]}));
    // vc counts by hugo_symbol
    for (let i = 0; i < top_10_genes.length; ++i) {
        vc = df_list[i]['Variant_Classification'].value_counts();
        for (let j = 0; j < vc_types.length; ++j) {
            data[j]['y'].push(vc['data'][vc['index_arr'].indexOf(vc_types[j])]);
        }
    }

    let layout = {
        height: 450,
        width: 1044,
        title: 'Top 10 Mutated Genes',
        barmode: 'stack',
        xaxis: {
            title: { text: 'Gene' },
            automargin: true
        },
        yaxis: {
            title: { text: 'Count' },
            automargin: true
        },
    };

    Plotly.newPlot("top_10_genes_plot", data, layout);

}

function build_ti_tv_plot(){

    let snv = {'T>G':[],'T>A':[],'T>C':[],'C>T':[],'C>G':[],'C>A':[]}
    let df_list = []
    let stack_colors = [
        '#1f77b4',  // muted blue
        '#ff7f0e',  // safety orange
        '#2ca02c',  // cooked asparagus green
        '#d62728',  // brick red
        '#9467bd',  // muted purple
        '#8c564b',  // chestnut brown
    ];
    let ti_tv_colors = [
        '#e377c2',  // raspberry yogurt pink
        '#7f7f7f',  // middle gray
    ];
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['maf'] != undefined) {
                let ref = wes_data[i]['somatic']['maf']['Reference_Allele']['data'];
                let allele2 = wes_data[i]['somatic']['maf']['Tumor_Seq_Allele2']['data'];
                let class_list = [];
                // reference_allele > tumor_seq_allele2
                for (let j=0;j<ref.length; ++j) {
                    // bases only
                    if (checkAvailability(bases, ref[j]) && checkAvailability(bases, allele2[j])) {
                        // use snv_conv function to convert bases
                        class_list.push(snv_conv[ref[j]+'>'+allele2[j]])
                    }
                }
                df_list.push(new dfd.Series(class_list));
            }
        }
    }
    // ti-tv data structure
    let ti_tv = {'Ti':Array(df_list.length).fill(0),'Tv':Array(df_list.length).fill(0)}

    for (let i = 0; i < df_list.length; ++i) {
        let counts = df_list[i].value_counts()
	//ensure that all SNV types are accounted for
        for (const snv_key in snv) {
            if (! counts['index_arr'].includes(snv_key)) {
                //ADD the missing key to the end of the list w/ 0 count
                counts['index_arr'].push(snv_key);
                counts['data'].push(0);
            }
        }
        let counts_sum = counts['data'].reduce((a, b) => a + b, 0)
        for (let j = 0; j < Object.keys(snv).length; ++j) {
            if (checkAvailability(counts['index_arr'], Object.keys(snv)[j])) {
                // normalize and add snv and ti-tv values
                snv[counts['index_arr'][j]].push((counts['data'][j]/counts_sum)*100);
                ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] = ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] + (counts['data'][j]/counts_sum)*100;
            } else {
                // if snv doesnt exist in sample, add 0 to snv and ti-tv data structures
                snv[counts['index_arr'][j]].push(0);
                ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] = ti_tv[ti_tv_conv[counts['index_arr'][j]]][i] + 0;
            }
        }
    }
    let data = [];

    for (let i = 0; i < Object.keys(snv).length; ++i) {
        // snv box plot
        data.push(
            {
                name: Object.keys(snv)[i],
                y: Object.values(snv)[i],
                type: 'box',
                marker: { color: stack_colors[i] },
            }
        )
        // stacked barplot
        data.push(
            {
                name: Object.keys(snv)[i],
                x: current_samples,
                y: Object.values(snv)[i],
                marker: { color: stack_colors[i] },
                type: 'bar',
                xaxis: 'x3',
                yaxis: 'y3',
            }
        )
    }
    // ti-tv boxplot
    for (let i = 0; i < Object.keys(ti_tv).length; ++i) {
        data.push(
            {
                name: Object.keys(ti_tv)[i],
                y: Object.values(ti_tv)[i],
                marker: { color: ti_tv_colors[i] },
                type: 'box',
                xaxis: 'x2',
                yaxis: 'y2',
            }
        )
    }

    let layout = {
        height: 450,
        width: 1044,
        title: 'Ti-Tv',
        barmode: 'stack',
        showlegend: false,
        // subplots
        xaxis: {domain: [0, 0.7], anchor: 'y'},
        yaxis: {domain: [0.55, 1], anchor: 'x', range: [0, 100], title: { text: '% Mutations' }},
        yaxis2: {domain: [0.55, 1], anchor: 'x2', range: [0, 100]},
        xaxis2: {domain: [0.8, 1], anchor: 'y2'},
        xaxis3: {domain: [0, 1], anchor: 'y3'},
        yaxis3: {domain: [0, 0.45], anchor: 'x3', range: [0, 100], title: { text: '% Mutations' }},
    };

    Plotly.newPlot("Ti-Tv_plot", data, layout);

}

function build_lego_plot(){

    let tri_template = {
        // each tri
        "A_A": 0,
        "A_C": 0,
        "A_G": 0,
        "A_T": 0,
        "C_A": 0,
        "C_C": 0,
        "C_G": 0,
        "C_T": 0,
        "G_A": 0,
        "G_C": 0,
        "G_G": 0,
        "G_T": 0,
        "T_A": 0,
        "T_C": 0,
        "T_G": 0,
        "T_T": 0
    }

    let tri_data = {
        // each snv
        "C>A": { ...tri_template },
        "C>G": { ...tri_template },
        "C>T": { ...tri_template },
        "T>A": { ...tri_template },
        "T>C": { ...tri_template },
        "T>G": { ...tri_template }
    }

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['tri_matrix'] != undefined) {
                for (var snv in tri_data) {
                    if (Object.prototype.hasOwnProperty.call(wes_data[i]['somatic']['tri_matrix'], snv)) {
                        // add counts from each snv to tri_data
                        tri_data[snv] = sumObjectsByKey({ ...tri_data[snv]},
                            { ...wes_data[i]['somatic']['tri_matrix'][snv]})
                    }
                }
            }
        }
    }
    // normalize tri_data in each snv
    for (var snv in tri_data) {
        let tri_sum = Object.values(tri_data[snv]).reduce((a, b) => a + b, 0)
            Object.keys(tri_data[snv])
                .forEach(function (tri) {
                    tri_data[snv][tri] = (tri_data[snv][tri] / tri_sum) * 100
                });
        }
    // subplot traces
    let data = [{
        x: Object.keys(tri_data['C>A']),
        y: Object.values(tri_data['C>A']),
        type: 'bar',
        xaxis: 'x',
        yaxis: 'y',
        name: 'C>A'
    }, {
        x: Object.keys(tri_data['C>G']),
        y: Object.values(tri_data['C>G']),
        type: 'bar',
        xaxis: 'x2',
        yaxis: 'y2',
        name: 'C>G'
    }, {
        x: Object.keys(tri_data['C>T']),
        y: Object.values(tri_data['C>T']),
        type: 'bar',
        xaxis: 'x3',
        yaxis: 'y3',
        name: 'C>T'
    }, {
        x: Object.keys(tri_data['T>A']),
        y: Object.values(tri_data['T>A']),
        type: 'bar',
        xaxis: 'x4',
        yaxis: 'y4',
        name: 'T>A'
    }, {
        x: Object.keys(tri_data['T>C']),
        y: Object.values(tri_data['T>C']),
        type: 'bar',
        xaxis: 'x5',
        yaxis: 'y5',
        name: 'T>C'
    }, {
        x: Object.keys(tri_data['T>G']),
        y: Object.values(tri_data['T>G']),
        type: 'bar',
        xaxis: 'x6',
        yaxis: 'y6',
        name: 'T>G'
    },
    ];

    let layout = {
        height: 450,
        width: 1044,
        title: 'Lego Plot',
        // subplot layout
        xaxis: {domain: [0, 0.3], anchor: 'y'},
        yaxis: {domain: [0.58, 1], anchor: 'x', title: { text: '% Mutations' }},
        yaxis2: {domain: [0.58, 1], anchor: 'x2'},
        xaxis2: {domain: [0.35, .65], anchor: 'y2'},
        xaxis3: {domain: [.7, 1], anchor: 'y3'},
        yaxis3: {domain: [.58, 1], anchor: 'x3'},
        xaxis4: {domain: [0, 0.3], anchor: 'y4'},
        yaxis4: {domain: [0, .42], anchor: 'x4', title: { text: '% Mutations' }},
        yaxis5: {domain: [0, .42], anchor: 'x5'},
        xaxis5: {domain: [0.35, .65], anchor: 'y5'},
        xaxis6: {domain: [.7, 1], anchor: 'y6'},
        yaxis6: {domain: [0, .42], anchor: 'x6'},
    };

    Plotly.newPlot("lego_plot_plot", data, layout);

}

function build_tmb_plot() {

    let tmb = [];
    let sample_ids = [];

    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            //add tmb
            tmb.push(wes_data[i]['somatic']['summary']['tmb'])
            sample_ids.push(wes_data[i]['id'])
        }
    }

    let trace1 = {
        x: sample_ids,
        y: tmb,
        type: 'bar'
    };

    let data = [trace1];

    let layout = {
        height: 450,
        width: 1044,
        xaxis: {
            title: { text: 'Sample' },
            automargin: true
        },
        yaxis: {
            title: { text: 'TMB (mut/Mb)' },
            automargin: true
        },
    };

    Plotly.newPlot("tumor_mutational_burden_plot", data, layout);

};

function build_oncoplot_plot() {

    let top_genes = [];
    // compile genes list from all samples
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['somatic']['geneList'] != undefined) {
                top_genes = top_genes.concat(wes_data[i]['somatic']['geneList'])
            }
        }
    }
    // obtain top genes
    top_genes = sortByFrequency(top_genes).slice(0, 25).reverse();
    let heatmap_z = Array.from(Array(top_genes.length), () => []);
    let sample_rank = [];
    // generate heatmap
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
                // run through top genes list
                for (let j = 0; j < top_genes.length; ++j) {
                    // evaluate if top gene is listed in sample gene list
                    if (checkAvailability(wes_data[i]['somatic']['geneList'], top_genes[j])) {
                        // hit
                        heatmap_z[j].push(1)
                        // add to sample count for later sorting
                        sample_rank.push(wes_data[i]['id'])
                    } else {
                        // non-hit
                        heatmap_z[j].push(0)
                    }
                }
            }
    }
    // sample overall abundance ranking
    sample_rank = sortByFrequency(sample_rank)
    let data = [
        {
            z: heatmap_z,
            x: current_samples,
            y: top_genes,
            type: 'heatmap',
            hoverongaps: false,
            showscale: false,
            colorscale: [
                ['0.0', 'rgb(180, 180, 180)'],
                ['1.0', 'rgb(165,0,38)']
            ],
            xgap: 1,
            ygap: 1
        }
    ];
    let layout = {
        height: 800,
        width: 800,
        title: 'Oncoplot',
        xaxis: {
            title: { text: 'Samples' },
            automargin: true,
            // sort samples by abundance ranking
            categoryorder: "array",
            categoryarray: sample_rank
        },
        yaxis: {
            title: { text: 'Genes' },
            automargin: true
        },
    };
    Plotly.newPlot('oncoplot_plot', data, layout);
}

function build_hla_oncoplot_plot() {

    let top_hla = [];
    let sample_ids = [];
    // compile hla list from all samples
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            for (let j = 0; j < sample_type.length; ++j) {
                if (typeof wes_data[i][sample_type[j]] != 'undefined') {
                    top_hla = top_hla.concat(Object.values(wes_data[i][sample_type[j]]['hla']))
                }
            }
        }
    }
    // obtain top hla counts
    top_hla = sortByFrequency(top_hla).slice(0, 25).reverse();
    let heatmap_z = Array.from(Array(top_hla.length), () => []);
    let sample_rank = [];
    // generate heatmap
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            // for normal and tumor type
            for (let j = 0; j < sample_type.length; ++j) {
                if (typeof wes_data[i][sample_type[j]] != 'undefined') {
                    // run through top hla list
                    for (let k = 0; k < top_hla.length; ++k) {
                        // evaluate if top hla is listed in sample hla lsit
                        if (checkAvailability(Object.values(wes_data[i][sample_type[j]]['hla']), top_hla[k])) {
                            // hit
                            heatmap_z[k].push(1)
                            // add to sample count for later sorting
                            sample_rank.push(wes_data[i]['id'] + sample_type_suffix[j])
                        } else {
                            // non-hit
                            heatmap_z[k].push(0)
                        }
                    }
                    sample_ids.push(wes_data[i]['id'] + sample_type_suffix[j])
                }
            }
        }
    }
    //sample abundance ranking
    sample_rank = sortByFrequency(sample_rank)
    let data = [
        {
            z: heatmap_z,
            x: sample_ids,
            y: top_hla,
            type: 'heatmap',
            hoverongaps: false,
            showscale: false,
            colorscale: [
                ['0.0', 'rgb(180, 180, 180)'],
                ['1.0', 'rgb(165,0,38)']
            ],
            xgap: 1,
            ygap: 1
        }
    ];
    let layout = {
        height: 800,
        width: 900,
        title: 'HLA Oncoplot',
        xaxis: {
            title: { text: 'Samples' },
            automargin: true,
            // sort samples by abundance ranking
            categoryorder: "array",
            categoryarray: sample_rank
        },
        yaxis: {
            title: { text: 'HLA' },
            automargin: true
        },
    };
    Plotly.newPlot('HLA_Oncoplot_plot', data, layout);
}

function build_tcellextrect_plot(){
    let sample_ids = [];
    let fraction = [];
    let hover_text = [];

    //Gets data for plot by iterating through wes_data
    for (let i = 0; i < wes_data.length; ++i) {
        if (checkAvailability(current_samples, wes_data[i]['id'])) {
            if (wes_data[i]['tcellextrect'] != undefined) {
                sample_ids.push(wes_data[i]['id']);
                fraction.push(wes_data[i]['tcellextrect']['tcell_fraction']);
                hover_text.push(
                    'Sample ID: ' + wes_data[i]['id'] + '<br>' +
                    'T-cell Fraction: ' + wes_data[i]['tcellextrect']['tcell_fraction'] + '<br>' +
                    'q-value: ' + wes_data[i]['tcellextrect']['q_value']
                );
            }
        }
    }

    // creates data and layout objects for graph
    let data = [{
        x: sample_ids,
        y: fraction,
        type: 'bar',
        hovertemplate: '%{text}<extra></extra>',
        text: hover_text

    }];
    let layout = {
        xaxis: {
        title: { text: 'Sample' },
        automargin: true,
        categoryorder: 'total descending',
        },
        yaxis: {
            title: { text: 'T-cell Fraction' },
            automargin: true
        },
        title: "T-cell Fraction by Sample"

    };

    Plotly.newPlot('tcellextrect_plot', data, layout);
}

function build_msisensor2_plot(){
  let sample_ids = [];
  let pct = [];
  let hover_text = [];

  //Gets data for plot by iterating through wes_data
  for (let i = 0; i < wes_data.length; ++i) {
      if (checkAvailability(current_samples, wes_data[i]['id'])) {
          if (wes_data[i]['msisensor2'] != undefined) {
              sample_ids.push(wes_data[i]['id']);
              pct.push(wes_data[i]['msisensor2']['percent_somatic']);
              hover_text.push(
                  'Sample ID: ' + wes_data[i]['id'] + '<br>' +
                  'Total Sites: ' + wes_data[i]['msisensor2']['total_sites'] + '<br>' +
                  'Somatic Sites: ' + wes_data[i]['msisensor2']['somatic_sites'] + '<br>' +
                  'Percent Somatic: ' + wes_data[i]['msisensor2']['percent_somatic']
              );
          }
      }
  }
  // creates data and layout objects for graph
  let data = [{
      x: sample_ids,
      y: pct,
      type: 'bar',
      hovertemplate: '%{text}<extra></extra>',
      text: hover_text

  }];

  let layout = {
      xaxis: {
          title: { text: 'Sample' },
          automargin: true,
          categoryorder: 'total descending',


      },
      yaxis: {
          title: { text: 'MSI Score' },
          automargin: true
      },
      title: "Microsatellite Score by Sample",

  };

  Plotly.newPlot('msisensor2_plot', data, layout);
}
