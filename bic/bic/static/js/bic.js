$(document).ready(function() {
    
    $.ajaxSetup({
        data: {
        csrfmiddlewaretoken: '{{ csrf_token }}'
        },
    });

    // reload page
    $('#reset').click(function () {
        location.reload(); 
    });


    var taxonomy_level_input = "Genus"
    var cancer_type_go_input = "ACC"
    var taxonomy_level_go_input = "Genus"
    var cancer_type_coabundance_input = "ACC"
    var taxonomy_level_coabundance_input = "Genus"


    var parse_data = new FormData();
    parse_data.append('taxonomy_level_input', taxonomy_level_input);
    parse_data.append('cancer_type_go_input', cancer_type_go_input);
    parse_data.append('taxonomy_level_go_input', taxonomy_level_go_input);
    parse_data.append('cancer_type_coabundance_input', cancer_type_coabundance_input);
    parse_data.append('taxonomy_level_coabundance_input', taxonomy_level_coabundance_input);
    
    //console.log(taxonomy_level_input)
    //console.log(cancer_type_go_input)
    //console.log(taxonomy_level_go_input)

    // Autocomplete: bacteria 

    $.ajax({
        url: '/bic/analyses/bacteria_name_ajax/', 
        data: parse_data,  
        type: 'POST',
        dataType: 'json',
        processData: false,
        contentType: false,
        success: function(data)
        {
            var bacteria = data
            //console.log(bacteria)


            $("#bacteria_input").autocomplete ({
                source: function(request, response) {
                    var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                    response( $.grep( bacteria, function( item ){
                        return matcher.test( item );
                    }));
                }
            });
            
        },
        error: function()
        {
            alert('Error occurred in autocomplete during searching bacteria= =|||');
        },

    });
    
    $("#taxonomy_level_input").change(function() {
        var taxonomy_level_input = $("#taxonomy_level_input").val();
        var parse_data = new FormData();
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        
        $.ajax({
            url: '/bic/analyses/bacteria_name_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                var bacteria = data
                //console.log(bacteria)

                $("#bacteria_input").autocomplete ({
                    source: function(request, response) {
                        var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                        response( $.grep( bacteria, function( item ){
                            return matcher.test( item );
                        }));
                    }
                });

                
            },
            error: function()
            {
                alert('Error occurred in autocomplete during searching bacteria= =|||');
            },

        });
    });




    // Autocomplete: bacteria human ran network and bacteria associated gene ontology
    
    $.ajax({
        url: '/bic/analyses/bacteria_name_go_ajax/', 
        data: parse_data,  
        type: 'POST',
        dataType: 'json',
        processData: false,
        contentType: false,
        success: function(data)
        {
            var bacteria = data
            //console.log(bacteria)

            $("#bacteria_go_input").autocomplete ({
                source: function(request, response) {
                    var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex( request.term ), "i" );
                    response( $.grep( bacteria, function( item ){
                        return matcher.test( item );
                    }));
                }
            });   
        }
    });



    $("#cancer_type_go_input").change(function() {
        var cancer_type_go_input = $("#cancer_type_go_input").val();
        var taxonomy_level_go_input = $("#taxonomy_level_go_input").val();
        var parse_data = new FormData();
        parse_data.append('cancer_type_go_input', cancer_type_go_input);
        parse_data.append('taxonomy_level_go_input', taxonomy_level_go_input);

        $.ajax({
            url: '/bic/analyses/bacteria_name_go_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                var bacteria = data
                //console.log(bacteria)

                $("#bacteria_go_input").autocomplete ({
                    source: function(request, response) {
                        var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                        response( $.grep( bacteria, function( item ){
                            return matcher.test( item );
                        }));
                    }
                });

                
            },
            error: function()
            {
                alert('Error occurred in autocomplete during searching bacteria= =|||');
            },

        });
    });


    
    $("#taxonomy_level_go_input").change(function() {
        var cancer_type_go_input = $("#cancer_type_go_input").val();
        var taxonomy_level_go_input = $("#taxonomy_level_go_input").val();
        var parse_data = new FormData();
        parse_data.append('cancer_type_go_input', cancer_type_go_input);
        parse_data.append('taxonomy_level_go_input', taxonomy_level_go_input);
        
        $.ajax({
            url: '/bic/analyses/bacteria_name_go_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                var bacteria = data
                //console.log(bacteria)

                $("#bacteria_go_input").autocomplete ({
                    source: function(request, response) {
                        var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                        response( $.grep( bacteria, function( item ){
                            return matcher.test( item );
                        }));
                    }
                });                
            },
            error: function()
            {
                alert('Error occurred in autocomplete during searching bacteria= =|||');
            },

        });
    });




    // Autocomplete: bacteria coabundance
    
    $.ajax({
        url: '/bic/analyses/bacteria_name_coabundance_ajax/', 
        data: parse_data,  
        type: 'POST',
        dataType: 'json',
        processData: false,
        contentType: false,
        success: function(data)
        {
            var bacteria = data
            // console.log(bacteria)

            $("#bacteria_coabundance_input").autocomplete ({
                source: function(request, response) {
                    var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex( request.term ), "i" );
                    response( $.grep( bacteria, function( item ){
                        return matcher.test( item );
                    }));
                }
            });   
        }
    });

    $("#cancer_type_coabundance_input").change(function() {
        var cancer_type_coabundance_input = $("#cancer_type_coabundance_input").val();
        var taxonomy_level_coabundance_input = $("#taxonomy_level_coabundance_input").val();
        var parse_data = new FormData();
        parse_data.append('cancer_type_coabundance_input', cancer_type_coabundance_input);
        parse_data.append('taxonomy_level_coabundance_input', taxonomy_level_coabundance_input);

        $.ajax({
            url: '/bic/analyses/bacteria_name_coabundance_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                var bacteria = data
                // console.log(bacteria)

                $("#bacteria_coabundance_input").autocomplete ({
                    source: function(request, response) {
                        var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                        response( $.grep( bacteria, function( item ){
                            return matcher.test( item );
                        }));
                    }
                });

                
            },
            error: function()
            {
                alert('Error occurred in autocomplete during searching bacteria= =|||');
            },

        });
    });
    
    $("#taxonomy_level_coabundance_input").change(function() {
        var cancer_type_coabundance_input = $("#cancer_type_coabundance_input").val();
        var taxonomy_level_coabundance_input = $("#taxonomy_level_coabundance_input").val();
        var parse_data = new FormData();
        parse_data.append('cancer_type_coabundance_input', cancer_type_coabundance_input);
        parse_data.append('taxonomy_level_coabundance_input', taxonomy_level_coabundance_input);
        
        $.ajax({
            url: '/bic/analyses/bacteria_name_coabundance_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                var bacteria = data
                // console.log(bacteria)

                $("#bacteria_coabundance_input").autocomplete ({
                    source: function(request, response) {
                        var matcher = new RegExp( "^" + $.ui.autocomplete.escapeRegex(request.term), "i" );
                        response( $.grep( bacteria, function( item ){
                            return matcher.test( item );
                        }));
                    }
                });                
            },
            error: function()
            {
                alert('Error occurred in autocomplete during searching bacteria= =|||');
            },

        });
    });

    

    $("#progress_bar").hide();
    $('#page-result').hide();


    // analyses for bacteria abundance 

    $('#analyses_abundance_query').click(function(){

        $("#progress_bar").show();
        $('#analyses_abundance_bacteria_input_error').empty();

        var taxonomy_level_input = $("#taxonomy_level_input").val();
        var bacteria_input = $("#bacteria_input").val();
        var parse_data = new FormData();
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('bacteria_input', bacteria_input);

        $.ajax({
            url: '/bic/analyses/abundance_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {   
                data = $.parseJSON(data)
                // console.log(data)
                
                
                $("#progress_bar").hide();

                if (data["bacteria_input_error_message"] !== '')
                {
                    $('#analyses_abundance_bacteria_input_error').append(data["bacteria_input_error_message"]);
                    $("#div_analyses_abundance_plot").html('No data available in plot');
                    $('#analyses_abundance_table').DataTable().clear().draw();
                }
                
                else {

                    $('#page-result').show();

                    $("#div_analyses_abundance_plot").html(data["div"]);
                    $("head").append(data["script"]);

                    mean = $.parseJSON(data["data"])

                    $('#analyses_abundance_table').DataTable({
                        data: mean,
                        destroy: true,
                        paging: false,
                        dom: 'Bfrtip',
                        buttons: [
                            'copy', 'excel', 'pdf', 'csv'
                        ],
                        columns: [
                            {data: "cancer_type"},
                            {data: "taxonomy_level"},
                            {data: "bacteria"},
                            {data: "mean_relative_abundance"}
                        ]
                    });
                }

                
            },
            error: function()
            {
                $("#progress_bar").hide();
                alert('Something error in abundance');
            },

        });
       
    });
    

    // analyses for bacteria diversity 

    $('#analyses_diversity_query').click(function(){

        $("#progress_bar").show();

        var taxonomy_level_input = $("#taxonomy_level_input").val();
        var diversity_method_input = $("#diversity_method_input").val();
        var parse_data = new FormData();
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('diversity_method_input', diversity_method_input);

        $.ajax({
            url: '/bic/analyses/diversity_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                $('#page-result').show();
                data = $.parseJSON(data)
                //console.log(data)

                $("#div_analyses_diversity_plot").empty();
                $("#div_analyses_diversity_plot").html(data["div"]);
                $("head").append(data["script"]);

                mean = $.parseJSON(data["data"])

                $('#analyses_diversity_table').DataTable({
                    data: mean,
                    destroy: true,
                    paging: false,
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', 'excel', 'pdf', 'csv'
                    ],
                    columns: [
                        {data: "cancer_type"},
                        {data: "taxonomy_level"},
                        {data: "mean_diversity"}
                    ]
                });
                
                $("#progress_bar").hide();

            },
            error: function()
            {
                $("#progress_bar").hide();
                alert('Something error in diversity');
            },

        });
       
    });



    
    // analyses for bacteria composition 

    $('#analyses_composition_query').click(function(){

        $("#progress_bar").show();

        var cancer_type_input = $("#cancer_type_input").val();
        var taxonomy_level_input = $("#taxonomy_level_input").val();
        var top_n_expressed_taxa_input = $("#top_n_expressed_taxa_input").val();
        var parse_data = new FormData();
        parse_data.append('cancer_type_input', cancer_type_input);
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('top_n_expressed_taxa_input', top_n_expressed_taxa_input);

        $.ajax({
            url: '/bic/analyses/composition_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                $('#page-result').show();
                data = $.parseJSON(data)
                //console.log(data)

                $("#div_analyses_composition_plot").empty();
                $("#div_analyses_composition_plot").html(data["div"]);
                $("head").append(data["script"]);

                all = $.parseJSON(data["data"])

                $('#analyses_composition_table').DataTable({
                    data: all,
                    destroy: true,
                    //paging: false,
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', 'excel', 'pdf', 'csv'
                    ],
                    //responsive: true,
                    columns: [
                        {data: "cancer_type_id"},
                        {data: "taxonomy_level"},
                        {data: "name"},
                        {data: "aliquot_barcode_mirna_seq"},
                        {data: "relative_abundance_gmpr_normalized"}
                    ]
                });

                $("#progress_bar").hide();
                
            },
            error: function()
            {
                $("#progress_bar").hide();
                alert('Something error in composition');
            },

        });
       
    });




    // analyses for bacteria clinical relevance  

    $('#analyses_clinical_relevance_query').click(function(){

        $("#progress_bar").show();

        $('#analyses_clinical_relevance_bacteria_input_error').empty();
        $("#analyses_clinical_relevance_plot_bacteria_all_zero_expression").empty();

        var cancer_type_input = $("#cancer_type_input").val();
        var taxonomy_level_input = $("#taxonomy_level_input").val();
        var bacteria_input = $("#bacteria_input").val();

        var parse_data = new FormData();
        parse_data.append('cancer_type_input', cancer_type_input);
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('bacteria_input', bacteria_input);

        $.ajax({
            url: '/bic/analyses/clinical_relevance_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                data = JSON.parse(data)
                // console.log(data)

                $("#progress_bar").hide();
                $("#div_analyses_clinical_relevance_plot_survival").empty();
                $("#div_analyses_clinical_relevance_plot_TPvsTN").empty();
                $("#div_analyses_clinical_relevance_plot_stage").empty();
                $("#div_analyses_clinical_relevance_plot_race").empty();
                

                if (data["bacteria_input_error_message"] !== '')
                {
                    $('#analyses_clinical_relevance_bacteria_input_error').html(data["bacteria_input_error_message"]);
                    $("#div_analyses_clinical_relevance_plot_survival").html('No data available in plot');
                    $('#analyses_clinical_relevance_table').DataTable().clear().draw();
                }
                
                else {
                    
                    $('#page-result').show();

                    //console.log(data)

                    if (data["bacteria_expression_all_zero_message"] !== "") 
                    {
                        $("#analyses_clinical_relevance_plot_bacteria_all_zero_expression").html(data["bacteria_expression_all_zero_message"]);
                        $('#analyses_clinical_relevance_table').DataTable().clear().draw();
                    }

                    else
                    {
                        //console.log(data["survival_plot_path"])
                        //console.log(data["data"])

                        var image_not_available = '<img src="" alt="Not available">';

                        if(data["div_os"] !== "NA")
                        {
                            $("#div_analyses_clinical_relevance_plot_survival").html(data["div_os"]);
                            $("head").append(data["script_os"]);
                        }
                        else
                        {
                            $("#div_analyses_clinical_relevance_plot_survival").append(image_not_available);
                        }
                        
                        if (data["div_wilcox_TPvsTN"] !== "NA")
                        {
                            $("#div_analyses_clinical_relevance_plot_TPvsTN").html(data["div_wilcox_TPvsTN"]);
                            $("head").append(data["script_wilcox_TPvsTN"]);
                        }
                        else
                        {
                            $("#div_analyses_clinical_relevance_plot_TPvsTN").append(image_not_available);
                        }

                        if (data["div_stage_early_median_late"] !== "NA")
                        {
                            $("#div_analyses_clinical_relevance_plot_stage").html(data["div_stage_early_median_late"]);
                            $("head").append(data["script_stage_early_median_late"]);
                        }
                        else
                        {
                            $("#div_analyses_clinical_relevance_plot_stage").append(image_not_available);
                        }

                        if (data["div_race_simplified"] !== "NA")
                        {
                            $("#div_analyses_clinical_relevance_plot_race").html(data["div_race_simplified"]);
                            $("head").append(data["script_race_simplified"]);
                        }
                        else
                        {
                            $("#div_analyses_clinical_relevance_plot_race").append(image_not_available);
                        }

                        
                        

                        //all = JSON.parse(data["data"]) //=> error one = =|||
                        all = [JSON.parse(data["data"])] //=> correct one
                        //console.log(all)
            
                        $('#analyses_clinical_relevance_table').DataTable({
                            data: all,
                            destroy: true,
                            paging: false,
                            dom: 'Bfrtip',
                            buttons: [
                                'copy', 'excel', 'pdf', 'csv'
                            ],
                            //responsive: true,
                            columns: [
                                {data: "cancer_type_id"},
                                {data: "taxonomy_level"},
                                {data: "bacteria"},
                                {data: "p_value_os"},
                                {data: "p_value_wilcox_TPvsTN"},
                                {data: "p_value_stage_early_median_late"},
                                {data: "p_value_race_simplified"}
                            ]
                        });
                    }
                }
            },
            error: function()
            {
                $("#progress_bar").hide();
                alert('Something error in clinical relevance.');
            },

        });
       
    });






    // analyses for bacteria coabundance network  

    $('#analyses_coabundance_network_query').click(function(){

        $("#progress_bar").show();
        $('#analyses_coabundance_network_bacteria_input_error').empty();
        $('#analyses_coabundance_network_p_value_input_error').empty();


        var cancer_type_coabundance_input = $("#cancer_type_coabundance_input").val();
        var taxonomy_level_coabundance_input = $("#taxonomy_level_coabundance_input").val();
        var bacteria_coabundance_input = $("#bacteria_coabundance_input").val();
        var p_value_input = $("#p_value_input").val();

        var parse_data = new FormData();
        parse_data.append('cancer_type_coabundance_input', cancer_type_coabundance_input);
        parse_data.append('taxonomy_level_coabundance_input', taxonomy_level_coabundance_input);
        parse_data.append('bacteria_coabundance_input', bacteria_coabundance_input);
        parse_data.append('p_value_input', p_value_input);

        $.ajax({
            url: '/bic/analyses/coabundance_network_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                data = JSON.parse(data)
                // console.log(data)

                $("#progress_bar").hide();

                if (data["bacteria_input_error_message"] !== '' || data["p_value_cutoff_error_message"] !== '')
                {
                    $('#analyses_coabundance_network_bacteria_input_error').append(data["bacteria_input_error_message"]);
                    $('#analyses_coabundance_network_p_value_input_error').append(data["p_value_cutoff_error_message"]);
                    // clear the plot and table if no data available
                    $("#div_analyses_coabundance_network_plot").html('No data available in plot');
                    $('#analyses_coabundance_network_table').DataTable().clear().draw();
                }
                else {

                    $('#page-result').show();

                    $("#div_analyses_coabundance_network_plot").html(data["div"]);
                    $("head").append(data["script"]);


                    network_table = $.parseJSON(data["data"])
                    // console.log(network_table)

                    $('#analyses_coabundance_network_table').DataTable({
                        data: network_table,
                        destroy: true,
                        paging: false,
                        dom: 'Bfrtip',
                        buttons: [
                            'copy', 'excel', 'pdf', 'csv'
                        ],
                        //responsive: true,
                        columns: [
                            {data: "cancer_type_id"},
                            {data: "taxonomy_level_id"},
                            {data: "taxonomy_level_bacterium_1"},
                            {data: "taxonomy_level_bacterium_2"},
                            {data: "sparscc"},
                            {data: "p_value"}
                        ]
                    });
                }
            },
            error: function()
            {
                $("#progress_bar").hide();
                alert('Something error in coabundance network');
            },

        });
        
    });







    // analyses for bacteria human rna network

    $('#analyses_bacteria_human_rna_network_query').click(function(){

        $("#progress_bar").show();
        $('#analyses_bacteria_human_rna_network_bacteria_input_error').empty();
        $('#analyses_bacteria_human_rna_network_z_score_input_error').empty();
        $('#analyses_bacteria_human_rna_network_number_edge_input_error').empty();

        var cancer_type_input = $("#cancer_type_go_input").val();
        var taxonomy_level_input = $("#taxonomy_level_go_input").val();
        var bacteria_input = $("#bacteria_go_input").val();
        var z_score_input = $("#z_score_input").val();
        var number_edge_input = $("#number_edge_input").val();

        var parse_data = new FormData();
        parse_data.append('cancer_type_input', cancer_type_input);
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('bacteria_input', bacteria_input);
        parse_data.append('z_score_input', z_score_input);
        parse_data.append('number_edge_input', number_edge_input);

        $.ajax({
            url: '/bic/analyses/bacteria_human_rna_network_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                
                data = JSON.parse(data)
                // console.log(data)

                $("#progress_bar").hide();

                if (data["bacteria_input_error_message"] !== '' || data["z_score_cutoff_error_message"] !== '' || data["number_edge_input_error_message"] !== '')
                {
                    $('#analyses_bacteria_human_rna_network_bacteria_input_error').append(data["bacteria_input_error_message"]);
                    $('#analyses_bacteria_human_rna_network_z_score_input_error').append(data["z_score_cutoff_error_message"]);
                    $('#analyses_bacteria_human_rna_network_number_edge_input_error').append(data["number_edge_input_error_message"]);
                    // clear the plot and table if no data available
                    $("#div_analyses_bacteria_human_rna_network_plot").html('No data available in plot');
                    $('#analyses_bacteria_human_rna_network_table').DataTable().clear().draw();
                }
                else {
                    $('#page-result').show();
                
                    $("#div_analyses_bacteria_human_rna_network_plot").html(data["div"]);
                    $("head").append(data["script"]);


                    network_table = $.parseJSON(data["data"])
                    // console.log(network_table)

                    $('#analyses_bacteria_human_rna_network_table').DataTable({
                        data: network_table,
                        destroy: true,
                        paging: true,
                        pageLength: 15,
                        dom: 'Bfrtip',
                        buttons: [
                            'copy', 'excel', 'pdf', 'csv'
                        ],
                        //responsive: true,
                        columns: [
                            {data: "cancer_type_id"},
                            {data: "taxonomy_level_id"},
                            {data: "name"},
                            {data: "human_rna"},
                            {data: "fisher_z_transformed_scc"},
                            {data: "scc"},
                            {data: "scc_p_value"}
                        ]
                    });
                }
                
            },
            error: function()
            {
                alert('Something error in bacteria mrna network');
                $("#progress_bar").hide();
            },

        });
        
    });









    // analyses for bacteria associated gene ontology

    $('#analyses_bacteria_associated_gene_ontology_query').click(function(){
        
        $("#progress_bar").show();
        $('#analyses_bacteria_associated_gene_ontology_bacteria_input_error').empty();
        $('#analyses_bacteria_associated_gene_ontology_adjust_p_value_input_error').empty();
        $('#analyses_bacteria_associated_gene_ontology_number_term_input_error').empty();

        var cancer_type_input = $("#cancer_type_go_input").val();
        var taxonomy_level_input = $("#taxonomy_level_go_input").val();
        var bacteria_input = $("#bacteria_go_input").val();
        var adjust_p_value_input = $("#adjust_p_value_input").val();
        var number_term_input = $("#number_term_input").val();

        var parse_data = new FormData();
        parse_data.append('cancer_type_input', cancer_type_input);
        parse_data.append('taxonomy_level_input', taxonomy_level_input);
        parse_data.append('bacteria_input', bacteria_input);
        parse_data.append('adjust_p_value_input', adjust_p_value_input);
        parse_data.append('number_term_input', number_term_input);

        $.ajax({
            url: '/bic/analyses/gene_ontology_ajax/', 
            data: parse_data,  
            type: 'POST',
            dataType: 'json',
            processData: false,
            contentType: false,
            success: function(data)
            {
                data = JSON.parse(data)
                // console.log(data)
                
                $("#progress_bar").hide();

                if (data["bacteria_input_error_message"] !== '' || data["adjust_p_value_cutoff_error_message"] !== '' || data["number_term_input_error_message"] !== '')
                {
                    $('#analyses_bacteria_associated_gene_ontology_bacteria_input_error').append(data["bacteria_input_error_message"]);
                    $('#analyses_bacteria_associated_gene_ontology_adjust_p_value_input_error').append(data["adjust_p_value_cutoff_error_message"]);
                    $('#analyses_bacteria_associated_gene_ontology_number_term_input_error').append(data["number_term_input_error_message"]);
                    // clear the plot and table if no data available
                    $("#div_analyses_bacteria_associated_gene_ontology_plot").html('No data available in plot');
                    $('#analyses_bacteria_associated_gene_ontology_table').DataTable().clear().draw();
                }
                
                else {
                    
                    $('#page-result').show();

                    $("#div_analyses_bacteria_associated_gene_ontology_plot").html(data["div"]);
                    $("head").append(data["script"]);


                    data = $.parseJSON(data["data"])
                    // console.log(data)

                    $('#analyses_bacteria_associated_gene_ontology_table').DataTable({
                        data: data,
                        destroy: true,
                        paging: true,
                        pageLength: 15,
                        scrollX: true,
                        dom: 'Bfrtip',
                        buttons: [
                            'copy', 'excel', 'pdf', 'csv'
                        ],
                        //responsive: true,
                        columns: [
                            {data: "cancer_type_id"},
                            {data: "taxonomy_level_id"},
                            {data: "name"},
                            {data: "term"},
                            {data: "p_value"},
                            {data: "adjust_p_value"},
                            {data: "es"},
                            {data: "nes"},
                            {data: "n_more_extreme"},
                            {data: "size"},
                            {data: "gene_symbol", width: '100px', class: 'text-left'}
                        ]
                    });
                }  
            },
            error: function()
            {
                alert("Something error in bacteria associated gene ontology");
                $("#progress_bar").hide();
            },

        });
        
    });


    // index figures

    $('#div_index_figures').slick({
        slidesToShow: 1,
        dots: true,
        infinite: true,
        // speed: 300,
        // autoplay: true,
        // autoplaySpeed: 2000,
        adaptiveHeight: false,
        // arrows: false
        arrows: true
    });

    

    // show the statistical table in statistics html

    $('#statistics_table').DataTable({
        paging: false,
        pageLength: 40,
        dom: 'Bfrtip',
        buttons: [''],
    });

    $('#menu_documents_bacteria_abundance_processing').click(function () {
        var documents_element = document.getElementById("a_documents_bacteria_abundance_processing");

        documents_element.scrollIntoView({
            behavior: "smooth",
            block: "center",
            inline: "center"
        });
    });
    
    $('#menu_documents_analysis_bacterial_abundance').click(function () {
        var documents_element = document.getElementById("a_documents_analysis_bacterial_abundance");
        
        documents_element.scrollIntoView({
            behavior: "smooth"
            // block: "center",
            // inline: "center"
        });
    });

});

