/* Len Taing 2020 (TGBTG) */

//var wes_resources = JSON.parse($("#wes_resources").text());

var current_samples = wes_data.map(function (x) {
    return [x][0]['id'];
});;

//Toggle filterTable show/hide btn text
$('[data-toggle="collapse"]').click(function() {
  $(this).toggleClass( "active" );
  if ($(this).hasClass("active")) {
    $(this).text("Hide\xa0Sample\xa0Selection\xa0Table");
  } else {
    $(this).text("Show\xa0Sample\xa0Selection\xa0Table");
  }
});

//Generate the filterTable-samples view first
function makeFilterTable() {
    var dset = runs_meta;
    var cols = ['id'];
    cols = cols.concat(Object.keys(dset[0]['annotations'])); //Add Annotation cols
    //var columns = cols.map(function(x){ return {'data': x}});
    var table = $('#filterTable');

    //BUILD the table up
    table.empty();
    content = "<thead><tr><th></th>"
    $.each(cols, function(i, col) { content += "<th>"+col+"</th>";});
    content += "</tr></thead>"
    $.each(dset, function(i, row) {
        content += "<tr><td></td>"
	//Add id
        content += "<td>"+row['id']+"</td>";
        for (v of Object.values(row['annotations'])) {
            content +="<td>"+v+"</td>";
         }
         content += "</tr>";
    });
    table.append(content);
    //Figure out number of searchPanes columns to use
    filter_num = $('#filterTable')[0].rows[0].cells.length-3
    if (filter_num < 6){searchpanes_col = 'columns-'+String(filter_num)} else {searchpanes_col = 'columns-6'}

    let tbl = table.DataTable({
        initComplete: function() {
            this.api().rows().select();
            $("th.select-checkbox").addClass("selected");
          },
          language: {
            select: {
                rows: ""
            }
        },
        dom: 'PSflBrtip',
        searchPanes: {
            dtOpts: {
                select: {
                    style: 'multi'
                }
            },
            layout: searchpanes_col,
            controls: false
        },
        columnDefs: [{
            orderable: false,
            className: 'select-checkbox',
            targets: 0
        }],
        select: {
            style: 'multi',
            selector: 'td:first-child'
        },
        order: [
            //[1, 'asc'] //This breaks datatables
        ],
        buttons: [
            {
                text: 'Update Sample Selection',
                action: function () {
                    current_samples = tbl.rows({ search:'applied', selected: true }).data()
                        .map(function (x) {
                            return x[1];
                        }).toArray();
                        loaded_plots = [];
                        load_section_plots(last_section);
                }
            }
        ]
    });
    tbl.on("click", "th.select-checkbox", function() {
    if ($("th.select-checkbox").hasClass("selected")) {
        tbl.rows().deselect();
        $("th.select-checkbox").removeClass("selected");
    } else {
        tbl.rows({search:'applied'}).select();
        $("th.select-checkbox").addClass("selected");
    }
}).on("select deselect", function() {
    ("Some selection or deselection going on")
    if (tbl.rows({
            selected: true
        }).count() !== tbl.rows().count()) {
        $("th.select-checkbox").removeClass("selected");
    } else {
        $("th.select-checkbox").addClass("selected");
    }

    /*
    //Version 2 -- Filter takes precedence over selections
    $('#myButton').on("click", function() {
        var searchTxt = $('.dataTables_filter input').val();
        if (searchTxt.length > 0) { //filter was applied
          //TRY to get filtered AND selected
          var rows = tbl.rows({selected:true, filter:'applied'});
          if (rows.count() == 0) {//None selected--revert to all
              rows = tbl.rows({filter:'applied'});
          }
        } else {
          //Selected only
          rows = tbl.rows({selected:true});
        }
        console.log(rows.count());
});
        */

});

}
makeFilterTable();

//Should generalize this fn- ColName, resource, handler
function handlerFactory(colName, resource, handler) {
    $(".data-coloured."+colName).on("click", function() {
	var table_id = this.closest('table').id;
	var sample_id = $(this.closest('tr')).find('th').attr('data-original-sn');
	var modal = $('#wesSubModal');
	var data = new Object();
	//#NOTE: wes_resources does not exist anymore
	data[sample_id] = wes_resources[resource][sample_id];
	handler(data, modal);
    });
    $(".data-coloured."+colName).css('cursor', 'pointer');
}
handlerFactory('Total_Mutations', 'somatic_summary', somaticSummary_submodal);
handlerFactory('TMB', 'tmb', somaticSummary_submodal);
handlerFactory('TiTv', 'ti_tv', complex_submodal);
handlerFactory('Nonsyn_Mutations', 'functional_summary', complex_submodal);

//Builds a table in the modal and displays it
function somaticSummary_submodal(data, modal) {
    //console.log(data);
    var modal_body = $(modal).find('.modal-body');
    //Clear contents
    modal_body.empty()

    //BUILD table
    var content = "<table class=\"table table-condensed mqc_table\">";
    //add header
    content += "<thead><tr>";
    //Sample names is first col
    content += "<th>&nbsp;</th>";
    var first_elm = Object.keys(data)[0];
    var hdr = Object.keys(data[first_elm]);
    $.each(hdr, function(i,val) { content += "<th>"+val+"</th>";})
    content += "</tr></thead>";
    //Add content
    const samples = Object.keys(data);
    for (s of samples) {
	content += "<tr><td>"+s+"</td>";
	var items = Object.values(data[s]);
	for (v of items) {
            content += "<td>"+v+"</td>";
	}
    }
    content += "</tr>";
    content += "</table>";
    modal_body.append(content);
    modal.modal({show:true});
}

function complex_submodal(data, modal) {
    //console.log(data);
    var modal_body = $(modal).find('.modal-body');
    //Clear contents
    modal_body.empty()

    //BUILD table
    var content = "<table class=\"table table-condensed mqc_table\">";
    //add header
    content += "<thead><tr>";
    content += "<th>&nbsp;</th>";
    var first_elm = Object.values(data)[0];
    //console.log(first_elm);
    var hdr = Object.keys(Object.values(first_elm)[0]);
    $.each(hdr, function(i,val) { content += "<th>"+val+"</th>";})
    content += "</tr></thead>";
    //Add content
    for (row of Object.keys(first_elm)) {
        content += "<tr><th>"+row+"</th>";
        var tmp = Object.values(first_elm[row]);
        for (v of tmp) { content += "<td>"+v+"</td>";}
    }
    content += "</tr>";
    content += "</table>";
    modal_body.append(content);
    modal.modal({show:true});
}

//returns table element
function buildTable(dataList) {
    var content = "<table class=\"table table-condensed\">";
    //add header
    content += "<thead><tr>";
    content += "<th></th>";
    var first_elm = dataList[0];
    var hdr = Object.keys(first_elm);
    $.each(hdr, function(i,val) { content += "<th>"+val+"</th>";})
    content += "</tr></thead>";
    //Add content
    $.each(dataList, function(i, row) {
        content += "<tr><th>"+i+"</th>";
        var tmp = Object.values(row);
        for (v of tmp) { content += "<td>"+v+"</td>";}
});
    content += "</tr>";
    content += "</table>";
    return $(content);
}

//var hlaTbl = $('#HLA_table').DataTable();
var neaontigenTbl = $('#neoantigen_table').DataTable({select: {style: 'single'}});
neoantigenTblRows = $('#neoantigen_table').find('tr');
//style the cursor
neoantigenTblRows.css({'cursor': 'pointer'});
//Add a show list btn just after the caption
var showListBtn = $('<button class=\"btn btn-primary btn-sm\">Download Full List</button>');
var neoantigenTblCap = $('#neoantigen_table').closest('.wes_table').find('.text-justify');
neoantigenTblCap.append(showListBtn);


showListBtn.on('click', function() {
    //Build a full DataTable but don't show it; invoke the save
    //table button
    var row = neaontigenTbl.rows({selected:true, filter:'applied'});
    //console.log(row);
    if (row.count() > 0) {
	var runId = row.data()[0][0];
	var neoList = wes_resources['neoantigen_table'][runId];
	//Write to file
	var tbl = buildTable(neoList);
	mb = $('#wesSubModal').find('.modal-body');
	mb.empty();
	mb.append(tbl);
	var tmp = tbl.DataTable({dom: 'Bfrtip', buttons: [{extend:'csvHtml5',title:runId+"_neoantigens"}]});
	//$('#wesSubModal').modal({show:true});
	tmp.buttons().trigger('click');
    }
});


var loaded_plots = [];

// Isolated from toggle becuase it must run when current_samples is updated
function load_section_plots(section_name){
  // builds plots only if section is not loaded
  if ($.inArray(section_name, loaded_plots) === -1){
    switch(section_name){
      case 'data_quality':
        build_data_quality();
        break;
      case 'copy_number_variation':
        build_copy_number_variation();
        break;
      case 'somatic_variants':
        build_somatic_variants();
        break;
      case 'HLA':
        build_HLA();
        break;
      case 'MISC':
        build_MISC();
        break;
    }
    loaded_plots.push(section_name);
  }
  /*else{
    console.log('already loaded');
  }*/
}

// handles all display and graphing tasks for new section when it is clicked
function sidebarSwitch(section_name){
  toggler(section_name); // hides old section, shows new section
  load_section_plots(section_name); // loads plots if needed
}

$(document).ready(function () {
    build_data_quality();
});


$(document).on('click', '.btn-export', function () {
    var plotlyDiv = $(this).parent().find('.js-plotly-plot').attr('id');
    Plotly.downloadImage(plotlyDiv, {format: 'svg', filename: plotlyDiv});
});
