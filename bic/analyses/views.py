from bokeh.embed.standalone import components
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, JsonResponse
from .models import *
from django.db.models import Avg, aggregates, Q
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import viridis, d3, Spectral4, Turbo256
import networkx as nx
from bokeh.models import Label, Legend, LegendItem, BoxSelectTool, Circle, EdgesAndLinkedNodes, HoverTool, MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool, CustomJSTransform, LabelSet, ColorBar, LogColorMapper, LinearColorMapper, ContinuousColorMapper, BoxZoomTool, ResetTool, SaveTool
from bokeh.plotting import from_networkx
import bokeh.layouts
from bokeh.layouts import row
from bokeh.embed import components
from bokeh.transform import transform, linear_cmap
from json import dumps
import os
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from matplotlib import pyplot as plt
from scipy.stats import ranksums, kruskal, wilcoxon
import numpy as np



# Create your views here.


## bacteria abundance 
def Analyses_bacteria_abundance(request):
    title = "analyses-bacteria-abundance"
    context = {"title": title}
    return render(request, "analyses/abundance/bacteria_abundance.html", context=context)


@csrf_exempt
def Autocomplete_search_bacteria_name(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    bacteria_name_query_list = list(Taxonomy_level_bacterium.objects.filter(taxonomy_level_id=taxonomy_level[0]).order_by("name").values_list("name", flat=True).distinct("name"))    
    return JsonResponse(bacteria_name_query_list, safe=False)


@csrf_exempt
def Autocomplete_search_bacteria_name_go(request):
    cancer_type_go = request.POST.get("cancer_type_go_input")
    taxonomy_level_go = request.POST.get("taxonomy_level_go_input")
    bacteria_name_query_list_go = list(View_bacteria_associated_gene_ontology.objects.filter(taxonomy_level_id=taxonomy_level_go[0], cancer_type_id=cancer_type_go).order_by("name").values_list("name", flat=True).distinct("name"))
    return JsonResponse(bacteria_name_query_list_go, safe=False)


@csrf_exempt
def Autocomplete_search_bacteria_name_coabundance(request):
    cancer_type_coabundance = request.POST.get("cancer_type_coabundance_input")
    taxonomy_level_coabundance = request.POST.get("taxonomy_level_coabundance_input")
    bacteria_name_query_set_coabundance = View_bacteria_coabundance.objects.filter(taxonomy_level_id=taxonomy_level_coabundance[0], cancer_type_id=cancer_type_coabundance).values()

    bacteria_name_query_set_df_coabundance = pd.DataFrame(list(bacteria_name_query_set_coabundance))
    bacteria_name_query_list_coabundance = sorted(list(set(list( bacteria_name_query_set_df_coabundance["taxonomy_level_bacterium_1"]) + list(bacteria_name_query_set_df_coabundance["taxonomy_level_bacterium_2"]) )))
    return JsonResponse(bacteria_name_query_list_coabundance, safe=False)





@csrf_exempt
def Ajax_bacteria_abundance(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    bacteria = request.POST.get("bacteria_input")

    # check if input bacteria is in our data
    bacteria_input_error_message = ""
    data = {}
    bacteria_name_available_list = list(Taxonomy_level_bacterium.objects.filter(taxonomy_level_id=taxonomy_level[0]).order_by("name").values_list("name", flat=True).distinct("name"))
    if bacteria not in bacteria_name_available_list:
        bacteria_input_error_message = "Not available"
        data["bacteria_input_error_message"] = bacteria_input_error_message
        return JsonResponse(dumps(data), safe=False)
    
    cancer_type_color_list = list(Cancer_type.objects.order_by("cancer_type_id").values_list("color_code_hex", flat=True))

    bacteria_query_set = View_bacteria_expression.objects.filter(name=bacteria, taxonomy_level_id=taxonomy_level[0]).values()
    bacteria_query_set_df = pd.DataFrame(list(bacteria_query_set))

    bacteria_cancer_type_mean_df = pd.DataFrame(bacteria_query_set_df.groupby('cancer_type_id')['relative_abundance_gmpr_normalized'].mean())
    cancer_type_list = list(bacteria_cancer_type_mean_df.index.values)

    groups = bacteria_query_set_df.loc[:,["cancer_type_id", "relative_abundance_gmpr_normalized"]].groupby('cancer_type_id')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    #print(upper)

    def outliers(cancer_type_id):
        cat = cancer_type_id.name
        #print(cat)
        return cancer_type_id[(cancer_type_id.relative_abundance_gmpr_normalized > upper.loc[cat]['relative_abundance_gmpr_normalized']) | (cancer_type_id.relative_abundance_gmpr_normalized < lower.loc[cat]['relative_abundance_gmpr_normalized'])]['relative_abundance_gmpr_normalized']
    
    out = groups.apply(outliers).dropna()

    if not out.empty:
        outx = list(out.index.get_level_values(0))
        outy = list(out.values)

    #tools="", plot_width:plot_height 16:9 or 4:3
    p = figure(plot_width=320, plot_height=180, x_range=cancer_type_list, toolbar_location="right")

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.relative_abundance_gmpr_normalized = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'relative_abundance_gmpr_normalized']),upper.relative_abundance_gmpr_normalized)]
    lower.relative_abundance_gmpr_normalized = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'relative_abundance_gmpr_normalized']),lower.relative_abundance_gmpr_normalized)]

    # stems
    p.segment(cancer_type_list, upper.relative_abundance_gmpr_normalized, cancer_type_list, q3.relative_abundance_gmpr_normalized, line_color="black")
    p.segment(cancer_type_list, lower.relative_abundance_gmpr_normalized, cancer_type_list, q1.relative_abundance_gmpr_normalized, line_color="black")

    # boxes
    p.vbar(cancer_type_list, 0.7, q2.relative_abundance_gmpr_normalized, q3.relative_abundance_gmpr_normalized, fill_color=cancer_type_color_list, line_color="black")
    p.vbar(cancer_type_list, 0.7, q1.relative_abundance_gmpr_normalized, q2.relative_abundance_gmpr_normalized, fill_color=cancer_type_color_list, line_color="black")


    # outliers
    if not out.empty:
        p.circle(outx, outy, size=0.5, color="gray", fill_alpha=0.2)

    #p.y_range.start = 0
    p.x_range.range_padding = 0.1
    p.yaxis.axis_label = "Relative abundance"
    # p.xaxis.major_label_orientation = 1
    p.xaxis.major_label_orientation = 3.414/6
    p.xgrid.grid_line_color = None
    p.sizing_mode = "scale_width"

    script, div = components(p)

    
    bacteria_cancer_type_mean_df = pd.DataFrame({'cancer_type': cancer_type_list, \
                                                 'taxonomy_level': taxonomy_level, \
                                                 'bacteria': bacteria, \
                                                 'mean_relative_abundance': bacteria_cancer_type_mean_df.iloc[:,0].values})

    bacteria_cancer_type_mean_json = bacteria_cancer_type_mean_df.to_json(orient='records')

    
    all_json = dumps(dict(script=script, div=div, data= bacteria_cancer_type_mean_json,
                          bacteria_input_error_message=bacteria_input_error_message))


    return JsonResponse(all_json, safe=False)



# #  -----------------------

## bacteria diversity 
def Analyses_bacteria_diversity(request):
    title = "analyses-bacteria-diversity"
    context = {"title": title}
    return render(request, "analyses/diversity/bacteria_diversity.html", context=context)


@csrf_exempt
def Ajax_bacteria_diversity(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    diversity_method = request.POST.get("diversity_method_input")

    cancer_type_color_list = list(Cancer_type.objects.order_by("cancer_type_id").values_list("color_code_hex", flat=True))
    # print("cancer_type_color_list: ", cancer_type_color_list)
    
    diversity_query_set = View_sample_bacteria_diversity.objects.filter(taxonomy_level_id=taxonomy_level[0], diversity_method=diversity_method).values()
    diversity_query_set_df = pd.DataFrame(list(diversity_query_set))

    diveristy_cancer_type_mean_df = pd.DataFrame(diversity_query_set_df.groupby('cancer_type_id')['diversity_value'].mean())
    cancer_type_list = list(diveristy_cancer_type_mean_df.index.values)

    groups = diversity_query_set_df.loc[:,["cancer_type_id", "diversity_value"]].groupby('cancer_type_id')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    #print(upper)

    def outliers(cancer_type_id):
        cat = cancer_type_id.name
        #print(cat)
        return cancer_type_id[(cancer_type_id.diversity_value > upper.loc[cat]['diversity_value']) | (cancer_type_id.diversity_value < lower.loc[cat]['diversity_value'])]['diversity_value']
    
    out = groups.apply(outliers).dropna()

    if not out.empty:
        outx = list(out.index.get_level_values(0))
        outy = list(out.values)

    #tools="", plot_width:plot_height 16:9 or 4:3
    p = figure(plot_width=320, plot_height=180, x_range=cancer_type_list, toolbar_location="right")

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.diversity_value = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'diversity_value']),upper.diversity_value)]
    lower.diversity_value = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'diversity_value']),lower.diversity_value)]



    # stems
    p.segment(cancer_type_list, upper.diversity_value, cancer_type_list, q3.diversity_value, line_color="black")
    p.segment(cancer_type_list, lower.diversity_value, cancer_type_list, q1.diversity_value, line_color="black")

    # boxes
    p.vbar(cancer_type_list, 0.7, q2.diversity_value, q3.diversity_value, fill_color=cancer_type_color_list, line_color="black")
    p.vbar(cancer_type_list, 0.7, q1.diversity_value, q2.diversity_value, fill_color=cancer_type_color_list, line_color="black")

    # outliers
    if not out.empty:
        p.circle(outx, outy, size=0.5, color="gray", fill_alpha=0.2)

    #p.y_range.start = 0
    p.x_range.range_padding = 0.1
    p.yaxis.axis_label = diversity_method
    # p.xaxis.major_label_orientation = 1
    p.xaxis.major_label_orientation = 3.414/6
    p.xgrid.grid_line_color = None
    p.sizing_mode = "scale_width"

    script, div = components(p)


    diveristy_cancer_type_mean_df = pd.DataFrame({'cancer_type': cancer_type_list, \
                                                  'taxonomy_level': taxonomy_level, \
                                                  'mean_diversity': diveristy_cancer_type_mean_df.iloc[:,0].values})

    diversity_cancer_type_mean_json = diveristy_cancer_type_mean_df.to_json(orient='records')

    all_json = dumps(dict(script=script, div=div, data=diversity_cancer_type_mean_json))

    return JsonResponse(all_json, safe=False)




#  -----------------------


## bacteria composition 
def Analyses_bacteria_composition(request):
    title = "analyses-bacteria-composition"
    context = {"title": title}
    return render(request, "analyses/composition/bacteria_composition.html", context=context)


@csrf_exempt
def Ajax_bacteria_composition(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    cancer_type = request.POST.get("cancer_type_input")
    top_n_expressed_taxa = int(request.POST.get("top_n_expressed_taxa_input"))
    # print(taxonomy_level)
    # print(cancer_type)
    # print(top_n_expressed_taxa)

    composition_query_set = View_bacteria_expression.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type).values()
    composition_query_set_df = pd.DataFrame(list(composition_query_set))
    #print(composition_query_set_df.columns)
    #print(composition_query_set_df.iloc[0:2, :])
    
    composition_query_subset_df = composition_query_set_df[["aliquot_barcode_mirna_seq", "name", "relative_abundance_gmpr_normalized", "cancer_type_id"]]
    composition_query_subset_df["taxonomy_level"] = taxonomy_level
    
    


    bacteria_relative_abundance_mean_df = composition_query_subset_df.groupby('name')['relative_abundance_gmpr_normalized'].mean().reset_index().sort_values(by='relative_abundance_gmpr_normalized', ascending=False)
    # print("\nbacteria_relative_abundance_mean_df.shape: ", bacteria_relative_abundance_mean_df.shape)
    # print("\bacteria_relative_abundance_mean_df.iloc[0:2, :]: ", bacteria_relative_abundance_mean_df.iloc[0:2, :])

    # bacteria_relative_abundance_mean_df = bacteria_relative_abundance_mean_df.sort_values(by='relative_abundance_gmpr_normalized', ascending=False)
    # print(bacteria_relative_abundance_mean_df)

    top_n_expressed_bacteria_name_list = bacteria_relative_abundance_mean_df.iloc[0:top_n_expressed_taxa, 0].tolist()
    # print("\ntop_n_expressed_bacteria_name_list: ", top_n_expressed_bacteria_name_list)


    others_df = composition_query_subset_df[~composition_query_subset_df.name.isin(top_n_expressed_bacteria_name_list)]
    # print("\nothers_df.shape: ", others_df.shape)

    expression_top_n_df = composition_query_subset_df[composition_query_subset_df.name.isin(top_n_expressed_bacteria_name_list)]
    # print("\nexpression_top_n_df.head(): ", expression_top_n_df.head())
    # print("\nexpression_top_n_df.shape: ", expression_top_n_df.shape)

    expression_others_df = others_df.groupby('aliquot_barcode_mirna_seq')['relative_abundance_gmpr_normalized'].sum().reset_index()
    # print("\nexpression_others_df.head(): ", expression_others_df.head())
    # print("\nexpression_others_df.shape: ", expression_others_df.shape)

    top_1_expressed_taxa_df = composition_query_subset_df[composition_query_subset_df["name"] == top_n_expressed_bacteria_name_list[0]].sort_values(by='relative_abundance_gmpr_normalized', ascending=False).reset_index()
    # print("\top_1_expressed_taxa_df.head(): ", top_1_expressed_taxa_df.head())
    # print("\top_1_expressed_taxa_df.shape: ", top_1_expressed_taxa_df.shape)

    aliquot_barcode_mirna_seq_ordered_list = top_1_expressed_taxa_df['aliquot_barcode_mirna_seq'].tolist()
    # print(aliquot_barcode_mirna_seq_ordered_list)


    top_n_expressed_taxa_dict = dict(aliquot_barcode_mirna_seq=aliquot_barcode_mirna_seq_ordered_list)
    # print(len(top_n_expressed_taxa_dict))


    for bacteria_name in top_n_expressed_bacteria_name_list:
        selected_bacteria_expr_df = composition_query_subset_df.loc[composition_query_subset_df["name"]==bacteria_name].set_index("aliquot_barcode_mirna_seq").loc[aliquot_barcode_mirna_seq_ordered_list].reset_index(inplace=False)
        top_n_expressed_taxa_dict[bacteria_name] = selected_bacteria_expr_df["relative_abundance_gmpr_normalized"].tolist()


    expression_others_df_reordered = expression_others_df.set_index("aliquot_barcode_mirna_seq").loc[aliquot_barcode_mirna_seq_ordered_list].reset_index(inplace=False)
    top_n_expressed_taxa_dict["Other"] = expression_others_df_reordered["relative_abundance_gmpr_normalized"].tolist()


    # print(len(top_n_expressed_taxa_dict))
    # print(top_n_expressed_taxa_dict)


    # for bacteria_name in top_n_expressed_bacteria_name_list:
    #     print(bacteria_name, ": ", top_n_expressed_taxa_dict[bacteria_name])

    expression_top_n_reformed_df = expression_top_n_df[["cancer_type_id", "taxonomy_level", "name", "aliquot_barcode_mirna_seq", "relative_abundance_gmpr_normalized"]]
    expression_others_reformed_df = pd.DataFrame(
        {"cancer_type_id": [cancer_type]*len(aliquot_barcode_mirna_seq_ordered_list),
        "taxonomy_level": [taxonomy_level]*len(aliquot_barcode_mirna_seq_ordered_list),
        "name": ["Other"]*len(aliquot_barcode_mirna_seq_ordered_list),
        "aliquot_barcode_mirna_seq": expression_others_df["aliquot_barcode_mirna_seq"].tolist(),
        "relative_abundance_gmpr_normalized": expression_others_df["relative_abundance_gmpr_normalized"].tolist()
        }
    )

    
    composition_output_df = pd.concat([expression_top_n_reformed_df, expression_others_reformed_df])
    # print("\ncomposition_output_df.head(): ", composition_output_df.head())
    # print("\ncomposition_output_df.shape: ", composition_output_df.shape)


    #tools="", plot_width:plot_height 16:9 or 4:3
    colors = d3['Category20'][top_n_expressed_taxa+1]
    # colors = d3['Category20'][top_n_expressed_taxa]
    # print("colors: ", colors)


    p = figure(plot_width=1040, plot_height=585, x_range=aliquot_barcode_mirna_seq_ordered_list, toolbar_location="right", 
                tools="hover", tooltips="$name's relative abundance: @$name")
    
    p.add_tools(SaveTool())
    p.add_layout(Legend(), 'right')
    p.vbar_stack(top_n_expressed_bacteria_name_list + ["Other"], x='aliquot_barcode_mirna_seq', color=colors, width=0.9, source=top_n_expressed_taxa_dict,
                    legend_label=top_n_expressed_bacteria_name_list + ["Other"])

    p.y_range.start = 0
    # p.x_range.range_padding = 0.4
    p.yaxis.axis_label = "Relative abundance"
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.outline_line_color = None
    # p.legend.location = "top_left"
    # p.legend.orientation = "horizontal"
    # p.xaxis.major_label_orientation = 1
    p.xaxis.major_label_orientation = 3.414/2.5


    script, div = components(p)

    #print("script: ", script)
    #print("div: ", div)

    composition_cancer_type_json = composition_output_df.to_json(orient='records')

    all_json = dumps(dict(script=script, div=div, data=composition_cancer_type_json))
    # print(all_json)

    return JsonResponse(all_json, safe=False)




#  -----------------------


## bacteria clinical_relevance 
def Analyses_bacteria_clinical_relevance(request):
    title = "analyses-bacteria-clinical-relevance"
    context = {"title": title}
    return render(request, "analyses/clinical_relevance/bacteria_clinical_relevance.html", context=context)


@csrf_exempt
def Ajax_bacteria_clinical_relevance(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    cancer_type = request.POST.get("cancer_type_input")
    bacteria = request.POST.get("bacteria_input")
    # print(taxonomy_level)
    # print(cancer_type)
    # print(bacteria)

    # check if input bacteria is in our data
    bacteria_input_error_message = ""
    bacteria_expression_all_zero_message = ""
    data = {}
    bacteria_name_available_list = list(Taxonomy_level_bacterium.objects.filter(taxonomy_level_id=taxonomy_level[0]).order_by("name").values_list("name", flat=True).distinct("name"))
    if bacteria not in bacteria_name_available_list:
        bacteria_input_error_message = "Not available"
        data["bacteria_input_error_message"] = bacteria_input_error_message
        return JsonResponse(dumps(data), safe=False)

    clinical_relevance_query_set = View_bacteria_clinical_relevance.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type, name=bacteria).values()
    clinical_relevance_query_set_df = pd.DataFrame(list(clinical_relevance_query_set))
    clinical_relevance_query_set_df["sample_type_id"] = clinical_relevance_query_set_df["sample_type_id"].str.replace('TN','N')
    clinical_relevance_query_set_df["sample_type_id"] = clinical_relevance_query_set_df["sample_type_id"].str.replace('TP','T')
    # print(clinical_relevance_query_set_df.columns)
    #print(clinical_relevance_query_set_df.iloc[0:2, :])
    #print(clinical_relevance_query_set_df.info())
    # print("shape of clinical_relevance_query_set_df: ", clinical_relevance_query_set_df.shape)


    p_value_os = "NA"
    p_value_wilcox_TPvsTN = "NA"
    p_value_race_simplified = "NA"
    p_value_stage_early_median_late = "NA"
    script_os = "NA"
    div_os = "NA" 
    script_wilcox_TPvsTN = "NA"
    div_wilcox_TPvsTN = "NA"
    script_race_simplified = "NA"
    div_race_simplified = "NA"
    script_stage_early_median_late = "NA"
    div_stage_early_median_late = "NA"

    # print("summazation of the taxa expression: ", sum(clinical_relevance_query_set_df["relative_abundance_gmpr_normalized"]>0))
    ## check if bacteria abundance are not all zeros
    if sum(clinical_relevance_query_set_df["relative_abundance_gmpr_normalized"]>0) > 0:
        ## overall survival (os)
        survival_df = clinical_relevance_query_set_df.copy(deep=True)
        survival_df = survival_df.loc[:, ['aliquot_barcode_mirna_seq', 'os', 'os_time', 'relative_abundance_gmpr_normalized', 'sample_type_id']]
        #print(survival_df.info)
        #print(survival_df.head())
        #print("shape of survival_df: ", survival_df.shape)

        survival_pt_df = survival_df.copy(deep=True)
        survival_pt_df = survival_pt_df.loc[survival_pt_df["sample_type_id"] == "T"]
        survival_pt_df["group"] = "low"
        #print(survival_pt_df.info)
        #print(survival_pt_df.head())
        #print("shape of survival_pt_df: ", survival_pt_df.shape)

        relative_abundance_gmpr_normalized_median = survival_pt_df["relative_abundance_gmpr_normalized"].median()

        survival_pt_df.loc[survival_pt_df["relative_abundance_gmpr_normalized"] > relative_abundance_gmpr_normalized_median, "group"] = "high"
        #print("high: ", sum(survival_pt_df["group"]=="high"))
        #print("low: ", sum(survival_pt_df["group"]=="low"))

        ## remove rows contained None in the os/os_time column
        survival_pt_omit_na_df = survival_pt_df[pd.notnull(survival_pt_df['os'])]
        survival_pt_omit_na_df = survival_pt_df[pd.notnull(survival_pt_df['os_time'])]
        #print("shape of survival_pt_omit_na_df: ", survival_pt_omit_na_df.shape)
        #print("high: ", sum(survival_pt_omit_na_df["group"]=="high"))
        #print("low: ", sum(survival_pt_omit_na_df["group"]=="low"))

        ## check if high group existed or not
        if sum(survival_pt_omit_na_df["group"]=="high") > 0:

            high_index = survival_pt_omit_na_df["group"] == "high"
            #print("high_index sum: ", sum(high_index))

            E1 = survival_pt_omit_na_df.loc[high_index]['os']
            T1 = survival_pt_omit_na_df.loc[high_index]['os_time']

            E2 = survival_pt_omit_na_df.loc[~high_index]['os']
            T2 = survival_pt_omit_na_df.loc[~high_index]['os_time']

            os_results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
            p_value_os = os_results.p_value

            if p_value_os < 0.05:
                p_value_os = "{:.3e}".format(p_value_os)
            else:
                p_value_os = "{:.3f}".format(p_value_os)

            # print(os_results.print_summary())
            

            kmf_low = KaplanMeierFitter()
            kmf_low_result = kmf_low.fit(T2, E2, label='low')
            kmf_low_survival_function = kmf_low_result.survival_function_
            kmf_low_xaxis = list(kmf_low_survival_function.index)
            kmf_low_yaxis = list(kmf_low_survival_function.iloc[:, 0])
            #print(kmf_low_survival_function.info())
            #print(kmf_low_xaxis[-10:])
            #print(kmf_low_yaxis[-10:])


            kmf_high = KaplanMeierFitter()
            kmf_high_result = kmf_high.fit(T1, E1, label='high')
            kmf_high_survival_function = kmf_high_result.survival_function_
            kmf_high_xaxis = list(kmf_high_survival_function.index)
            kmf_high_yaxis = list(kmf_high_survival_function.iloc[:, 0])
            #print(kmf_high_survival_function.info())

            p_os = figure(plot_width=400, plot_height=400)
            p_os.step(kmf_low_xaxis, kmf_low_yaxis, line_width=2, mode="after", line_color="#4472C4")
            p_os.step(kmf_high_xaxis, kmf_high_yaxis, line_width=2, mode="after", line_color="#F94472")

            p_os.y_range.start = 0
            p_os.x_range.range_padding = 0.1
            p_os.xgrid.grid_line_color = None
            p_os.axis.minor_tick_line_color = None
            p_os.outline_line_color = None
            p_os.legend.orientation = "horizontal"
            # p_os.xaxis.major_label_orientation = 0
            p_os.xaxis.major_label_orientation = 3.414/15
            p_os.xaxis.axis_label = "Days"
            p_os.yaxis.axis_label = "Survival Probability"

            p_annotaion = Label(x=int(max(kmf_high_xaxis)/2), y=0.9, \
                        text='p = '+ p_value_os )

            p_os.add_layout(p_annotaion)

            os_low_legend = LegendItem(label="low", renderers=[p_os.renderers[0]])
            os_high_legend = LegendItem(label="high", renderers=[p_os.renderers[1]])
            legend_os = Legend(items=[os_low_legend, os_high_legend], location="top_right")
            p_os.add_layout(legend_os)

            script_os, div_os = components(p_os)
        
            
            



        ## wilcox test tumor (T) vs adjacent normal (N)

        ## check if TN group existed or not
        if sum(clinical_relevance_query_set_df["sample_type_id"]=="N") > 0:

            number_TN = sum(clinical_relevance_query_set_df["sample_type_id"]=="N")
            number_TP = sum(clinical_relevance_query_set_df["sample_type_id"]=="T")

            # print("number_TN: ", number_TN)
            # print("number_TP: ", number_TP)

            sample_type_list = ["N", "T"]
            # print("sample_type_list: ", sample_type_list)
            groups_wilcox_TPvsTN = clinical_relevance_query_set_df.loc[:,["sample_type_id", "relative_abundance_gmpr_normalized"]].groupby('sample_type_id')
            # print("groups_wilcox_TPvsTN", groups_wilcox_TPvsTN)

            groups_wilcox_TPvsTN_color_list = ["#7E8BF9", "#D865AF"]
            
            TP_list = list(clinical_relevance_query_set_df.loc[clinical_relevance_query_set_df["sample_type_id"]=="T", "relative_abundance_gmpr_normalized"])
            #print("len(TP_list): ", len(TP_list))
            #print(TP_list[:5])
            TN_list = list(clinical_relevance_query_set_df.loc[clinical_relevance_query_set_df["sample_type_id"]=="N", "relative_abundance_gmpr_normalized"])
            #print("len(TN_list): ", len(TN_list))
            #print(TN_list[:5])

            stat_wilcox_TPvsTN, p_value_wilcox_TPvsTN = ranksums(TP_list, TN_list)
            # stat_wilcox_TPvsTN, p_value_wilcox_TPvsTN = wilcoxon(TP_list, TN_list, correction=True)

            # print('Statistics=%.3f, p=%.3f' % (stat_wilcox_TPvsTN, p_value_wilcox_TPvsTN))

            if p_value_wilcox_TPvsTN < 0.05:
                p_value_wilcox_TPvsTN = "{:.3e}".format(p_value_wilcox_TPvsTN)
            else:
                p_value_wilcox_TPvsTN = "{:.3f}".format(p_value_wilcox_TPvsTN)


            q1_wilcox_TPvsTN = groups_wilcox_TPvsTN.quantile(q=0.25)
            q2_wilcox_TPvsTN = groups_wilcox_TPvsTN.quantile(q=0.5)
            q3_wilcox_TPvsTN = groups_wilcox_TPvsTN.quantile(q=0.75)
            iqr_wilcox_TPvsTN = q3_wilcox_TPvsTN - q1_wilcox_TPvsTN
            upper_wilcox_TPvsTN = q3_wilcox_TPvsTN + 1.5*iqr_wilcox_TPvsTN
            lower_wilcox_TPvsTN = q1_wilcox_TPvsTN - 1.5*iqr_wilcox_TPvsTN

            # print(upper_wilcox_TPvsTN)

            def outliers_wilcox_TPvsTN(sample_type_id):
                cat = sample_type_id.name
                # print(cat)
                return sample_type_id[(sample_type_id.relative_abundance_gmpr_normalized > upper_wilcox_TPvsTN.loc[cat]['relative_abundance_gmpr_normalized']) | (sample_type_id.relative_abundance_gmpr_normalized < lower_wilcox_TPvsTN.loc[cat]['relative_abundance_gmpr_normalized'])]['relative_abundance_gmpr_normalized']
            
            out_wilcox_TPvsTN = groups_wilcox_TPvsTN.apply(outliers_wilcox_TPvsTN).dropna()

            if not out_wilcox_TPvsTN.empty:
                outx_wilcox_TPvsTN = list(out_wilcox_TPvsTN.index.get_level_values(0))
                outy_wilcox_TPvsTN = list(out_wilcox_TPvsTN.values)

            #tools="", plot_width:plot_height 16:9 or 4:3
            p_wilcox_TPvsTN = figure(plot_width=400, plot_height=400, x_range=sample_type_list, toolbar_location="right")

            # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
            qmin_wilcox_TPvsTN = groups_wilcox_TPvsTN.quantile(q=0.00)
            qmax_wilcox_TPvsTN = groups_wilcox_TPvsTN.quantile(q=1.00)
            upper_wilcox_TPvsTN.relative_abundance_gmpr_normalized = [min([x,y]) for (x,y) in zip(list(qmax_wilcox_TPvsTN.loc[:,'relative_abundance_gmpr_normalized']),upper_wilcox_TPvsTN.relative_abundance_gmpr_normalized)]
            lower_wilcox_TPvsTN.relative_abundance_gmpr_normalized = [max([x,y]) for (x,y) in zip(list(qmin_wilcox_TPvsTN.loc[:,'relative_abundance_gmpr_normalized']),lower_wilcox_TPvsTN.relative_abundance_gmpr_normalized)]

            # stems
            p_wilcox_TPvsTN.segment(sample_type_list, upper_wilcox_TPvsTN.relative_abundance_gmpr_normalized, sample_type_list, q3_wilcox_TPvsTN.relative_abundance_gmpr_normalized, line_color="black")
            p_wilcox_TPvsTN.segment(sample_type_list, lower_wilcox_TPvsTN.relative_abundance_gmpr_normalized, sample_type_list, q1_wilcox_TPvsTN.relative_abundance_gmpr_normalized, line_color="black")

            # boxes
            p_wilcox_TPvsTN.vbar(sample_type_list, 0.7, q2_wilcox_TPvsTN.relative_abundance_gmpr_normalized, q3_wilcox_TPvsTN.relative_abundance_gmpr_normalized, fill_color=groups_wilcox_TPvsTN_color_list, line_color="black")
            p_wilcox_TPvsTN.vbar(sample_type_list, 0.7, q1_wilcox_TPvsTN.relative_abundance_gmpr_normalized, q2_wilcox_TPvsTN.relative_abundance_gmpr_normalized, fill_color=groups_wilcox_TPvsTN_color_list, line_color="black")

            # outliers
            if not out_wilcox_TPvsTN.empty:
                p_wilcox_TPvsTN.circle(outx_wilcox_TPvsTN, outy_wilcox_TPvsTN, size=2, color="gray", fill_alpha=0.2)

            #p_wilcox_TPvsTN.y_range.start = 0
            p_wilcox_TPvsTN_annotaion = Label(x=1, y=max(clinical_relevance_query_set_df.loc[:, "relative_abundance_gmpr_normalized"]),
                                            text='p = ' + str(p_value_wilcox_TPvsTN) )
            p_wilcox_TPvsTN.x_range.range_padding = 0.1
            p_wilcox_TPvsTN.yaxis.axis_label = "Relative abundance"
            # p_wilcox_TPvsTN.xaxis.major_label_orientation = 1
            p_wilcox_TPvsTN.xaxis.major_label_orientation = 3.414/6
            p_wilcox_TPvsTN.xgrid.grid_line_color = None
            #p_wilcox_TPvsTN.sizing_mode = "scale_width"

            script_wilcox_TPvsTN, div_wilcox_TPvsTN = components(p_wilcox_TPvsTN)   





        ## Race Krus-Wall test tumor 
        race_pt_df = clinical_relevance_query_set_df.copy(deep=True)
        race_pt_df = race_pt_df.loc[:, ["aliquot_barcode_mirna_seq", "sample_type_id", "race_simplified", "relative_abundance_gmpr_normalized"]]
        race_pt_df = race_pt_df.loc[race_pt_df["sample_type_id"] == "T"]
        race_pt_df = race_pt_df.copy(deep=True)
        # print("dim(race_pt_df: ", race_pt_df.shape)
        
        ## there is a space in the race 'African American ' >< => re-upload the correct name (no space at the end) in bic_production version
        #race_pt_df['race_simplified'] = race_pt_df['race_simplified'].replace({'African American ' : 'African American'}, regex=False)
        # african_american_list = list(race_pt_df.loc[race_pt_df["race_simplified"]=="African American ", "relative_abundance_gmpr_normalized"])


        african_american_list = list(race_pt_df.loc[race_pt_df["race_simplified"]=="African American", "relative_abundance_gmpr_normalized"])
        asian_list = list(race_pt_df.loc[race_pt_df["race_simplified"]=="Asian", "relative_abundance_gmpr_normalized"])
        white_list = list(race_pt_df.loc[race_pt_df["race_simplified"]=="White", "relative_abundance_gmpr_normalized"])
        other_list = list(race_pt_df.loc[race_pt_df["race_simplified"]=="Other", "relative_abundance_gmpr_normalized"])
        # print("len(african_american_list_list): ", len(african_american_list))
        # print("len(asian_list): ", len(asian_list))
        # print("len(white_list): ", len(white_list))
        # print("len(other_list): ", len(other_list))

        race_simplified_list = []
        race_simplified_dict = dict()
        groups_race_simplified_color_list = []

        if (len(african_american_list)>0):
            race_simplified_list.append("African American")
            race_simplified_dict["African American"] = african_american_list
            groups_race_simplified_color_list.append("#C55A11")
        
        if (len(asian_list)>0):
            race_simplified_list.append("Asian")
            race_simplified_dict["Asian"] = asian_list
            groups_race_simplified_color_list.append("#FFC000")

        
        if (len(other_list)>0):
            race_simplified_list.append("Other")
            race_simplified_dict["Other"] = other_list
            groups_race_simplified_color_list.append("#5B9BD5")


        if (len(white_list)>0):
            race_simplified_list.append("White")
            race_simplified_dict["White"] = white_list
            groups_race_simplified_color_list.append("#44546A")


        groups_race_simplified = race_pt_df.loc[:,["race_simplified", "relative_abundance_gmpr_normalized"]].groupby("race_simplified")


        if len(race_simplified_list) > 1:

            try:
                stat_race_simplified, p_value_race_simplified = kruskal(*[race_simplified_dict[grp] for grp in race_simplified_list])
            
                # print('Statistics=%.3f, p=%.3f' % (stat_race_simplified, p_value_race_simplified))

                if p_value_race_simplified < 0.05:
                    p_value_race_simplified = "{:.3e}".format(p_value_race_simplified)
                else:
                    p_value_race_simplified = "{:.3f}".format(p_value_race_simplified)

                # print("p_value_race_simplified", p_value_race_simplified)
            
            

                q1_race_simplified = groups_race_simplified.quantile(q=0.25)
                q2_race_simplified = groups_race_simplified.quantile(q=0.5)
                q3_race_simplified = groups_race_simplified.quantile(q=0.75)
                iqr_race_simplified = q3_race_simplified - q1_race_simplified
                upper_race_simplified = q3_race_simplified + 1.5*iqr_race_simplified
                lower_race_simplified = q1_race_simplified - 1.5*iqr_race_simplified

                #print(upper_race_simplified)

                def outliers_race_simplified(race_simplified):
                    cat = race_simplified.name
                    #print(cat)
                    return race_simplified[(race_simplified.relative_abundance_gmpr_normalized > upper_race_simplified.loc[cat]['relative_abundance_gmpr_normalized']) | (race_simplified.relative_abundance_gmpr_normalized < lower_race_simplified.loc[cat]['relative_abundance_gmpr_normalized'])]['relative_abundance_gmpr_normalized']
                
                out_race_simplified = groups_race_simplified.apply(outliers_race_simplified).dropna()

                if not out_race_simplified.empty:
                    outx_race_simplified = list(out_race_simplified.index.get_level_values(0))
                    outy_race_simplified = list(out_race_simplified.values)

                #tools="", plot_width:plot_height 16:9 or 4:3
                p_race_simplified = figure(plot_width=400, plot_height=400, x_range=race_simplified_list, toolbar_location="right")

                # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
                qmin_race_simplified = groups_race_simplified.quantile(q=0.00)
                qmax_race_simplified = groups_race_simplified.quantile(q=1.00)
                upper_race_simplified.relative_abundance_gmpr_normalized = [min([x,y]) for (x,y) in zip(list(qmax_race_simplified.loc[:,'relative_abundance_gmpr_normalized']),upper_race_simplified.relative_abundance_gmpr_normalized)]
                lower_race_simplified.relative_abundance_gmpr_normalized = [max([x,y]) for (x,y) in zip(list(qmin_race_simplified.loc[:,'relative_abundance_gmpr_normalized']),lower_race_simplified.relative_abundance_gmpr_normalized)]

                # stems
                p_race_simplified.segment(race_simplified_list, upper_race_simplified.relative_abundance_gmpr_normalized, race_simplified_list, q3_race_simplified.relative_abundance_gmpr_normalized, line_color="black")
                p_race_simplified.segment(race_simplified_list, lower_race_simplified.relative_abundance_gmpr_normalized, race_simplified_list, q1_race_simplified.relative_abundance_gmpr_normalized, line_color="black")

                # boxes
                p_race_simplified.vbar(race_simplified_list, 0.7, q2_race_simplified.relative_abundance_gmpr_normalized, q3_race_simplified.relative_abundance_gmpr_normalized, fill_color=groups_race_simplified_color_list, line_color="black")
                p_race_simplified.vbar(race_simplified_list, 0.7, q1_race_simplified.relative_abundance_gmpr_normalized, q2_race_simplified.relative_abundance_gmpr_normalized, fill_color=groups_race_simplified_color_list, line_color="black")

                # outliers
                if not out_race_simplified.empty:
                    p_race_simplified.circle(outx_race_simplified, outy_race_simplified, size=2, color="gray", fill_alpha=0.2)

                #p_race_simplified.y_range.start = 0
                p_race_simplified_annotaion = Label(x=1, y=max(clinical_relevance_query_set_df.loc[:, "relative_abundance_gmpr_normalized"]),
                                                text='p = ' + str(p_value_race_simplified))
                p_race_simplified.x_range.range_padding = 0.1
                p_race_simplified.yaxis.axis_label = "Relative abundance"
                # p_race_simplified.xaxis.major_label_orientation = 1
                p_race_simplified.xaxis.major_label_orientation = 3.414/15
                p_race_simplified.xgrid.grid_line_color = None
                #p_race_simplified.sizing_mode = "scale_width"

                script_race_simplified, div_race_simplified = components(p_race_simplified)
        
            except:
                print("")







        ## stage_early_median_late Krus-Wall test tumor 
        stage_pt_df = clinical_relevance_query_set_df.copy(deep=True)
        stage_pt_df = stage_pt_df.loc[:, ["aliquot_barcode_mirna_seq", "sample_type_id", "stage_early_median_late", "relative_abundance_gmpr_normalized"]]
        stage_pt_df = stage_pt_df.loc[stage_pt_df["sample_type_id"] == "T"]
        stage_pt_df = stage_pt_df.copy(deep=True)
        #print("dim(stage_pt_df: ", stage_pt_df.shape)


        groups_stage_early_median_late = stage_pt_df.loc[:,["stage_early_median_late", "relative_abundance_gmpr_normalized"]].groupby("stage_early_median_late")
        # print("---------")
        # print(stage_pt_df.loc[:,["stage_early_median_late", "relative_abundance_gmpr_normalized"]].groupby("stage_early_median_late").sum())
        # print("---------")

        if not pd.DataFrame(groups_stage_early_median_late).empty:

            early_list = list(stage_pt_df.loc[stage_pt_df["stage_early_median_late"]=="I", "relative_abundance_gmpr_normalized"])
            median_list = list(stage_pt_df.loc[stage_pt_df["stage_early_median_late"]=="II&III", "relative_abundance_gmpr_normalized"])
            late_list = list(stage_pt_df.loc[stage_pt_df["stage_early_median_late"]=="IV", "relative_abundance_gmpr_normalized"])
            # print("len(early_list): ", len(early_list))
            # print("len(median_list): ", len(median_list))
            # print("len(late_list): ", len(late_list))

            stage_early_median_late_list = []
            stage_early_median_late_dict = {}
            groups_stage_early_median_late_color_list = []

            if (len(early_list)>0):
                stage_early_median_late_list.append("I")
                stage_early_median_late_dict["I"] = early_list
                groups_stage_early_median_late_color_list.append("#00B050")
            
            if (len(median_list)>0):
                stage_early_median_late_list.append("II&III")
                stage_early_median_late_dict["II&III"] = median_list
                groups_stage_early_median_late_color_list.append("#5B9BD5")

            if (len(late_list)>0):
                stage_early_median_late_list.append("IV")
                stage_early_median_late_dict["IV"] = late_list
                groups_stage_early_median_late_color_list.append("#7030A0")

            
            try:
                stat_stage_early_median_late, p_value_stage_early_median_late = kruskal(*[stage_early_median_late_dict[grp] for grp in stage_early_median_late_list])
                # print('Statistics=%.3f, p=%.3f' % (stat_stage_early_median_late, p_value_stage_early_median_late))

                if p_value_stage_early_median_late < 0.05:
                    p_value_stage_early_median_late = "{:.3e}".format(p_value_stage_early_median_late)
                else:
                    p_value_stage_early_median_late = "{:.3f}".format(p_value_stage_early_median_late)

                # print("p_value_stage_early_median_late", p_value_stage_early_median_late)
            

                q1_stage_early_median_late = groups_stage_early_median_late.quantile(q=0.25)
                q2_stage_early_median_late = groups_stage_early_median_late.quantile(q=0.5)
                q3_stage_early_median_late = groups_stage_early_median_late.quantile(q=0.75)
                iqr_stage_early_median_late = q3_stage_early_median_late - q1_stage_early_median_late
                upper_stage_early_median_late = q3_stage_early_median_late + 1.5*iqr_stage_early_median_late
                lower_stage_early_median_late = q1_stage_early_median_late - 1.5*iqr_stage_early_median_late

                #print(upper_stage_early_median_late)

                def outliers_stage_early_median_late(stage_early_median_late):
                    cat = stage_early_median_late.name
                    #print(cat)
                    return stage_early_median_late[(stage_early_median_late.relative_abundance_gmpr_normalized > upper_stage_early_median_late.loc[cat]['relative_abundance_gmpr_normalized']) | (stage_early_median_late.relative_abundance_gmpr_normalized < lower_stage_early_median_late.loc[cat]['relative_abundance_gmpr_normalized'])]['relative_abundance_gmpr_normalized']
                
                out_stage_early_median_late = groups_stage_early_median_late.apply(outliers_stage_early_median_late).dropna()

                if not out_stage_early_median_late.empty:
                    outx_stage_early_median_late = list(out_stage_early_median_late.index.get_level_values(0))
                    outy_stage_early_median_late = list(out_stage_early_median_late.values)

                #tools="", plot_width:plot_height 16:9 or 4:3
                p_stage_early_median_late = figure(plot_width=400, plot_height=400, x_range=stage_early_median_late_list, toolbar_location="right")

                # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
                qmin_stage_early_median_late = groups_stage_early_median_late.quantile(q=0.00)
                qmax_stage_early_median_late = groups_stage_early_median_late.quantile(q=1.00)
                upper_stage_early_median_late.relative_abundance_gmpr_normalized = [min([x,y]) for (x,y) in zip(list(qmax_stage_early_median_late.loc[:,'relative_abundance_gmpr_normalized']),upper_stage_early_median_late.relative_abundance_gmpr_normalized)]
                lower_stage_early_median_late.relative_abundance_gmpr_normalized = [max([x,y]) for (x,y) in zip(list(qmin_stage_early_median_late.loc[:,'relative_abundance_gmpr_normalized']),lower_stage_early_median_late.relative_abundance_gmpr_normalized)]

                # stems
                p_stage_early_median_late.segment(stage_early_median_late_list, upper_stage_early_median_late.relative_abundance_gmpr_normalized, stage_early_median_late_list, q3_stage_early_median_late.relative_abundance_gmpr_normalized, line_color="black")
                p_stage_early_median_late.segment(stage_early_median_late_list, lower_stage_early_median_late.relative_abundance_gmpr_normalized, stage_early_median_late_list, q1_stage_early_median_late.relative_abundance_gmpr_normalized, line_color="black")

                # boxes
                p_stage_early_median_late.vbar(stage_early_median_late_list, 0.7, q2_stage_early_median_late.relative_abundance_gmpr_normalized, q3_stage_early_median_late.relative_abundance_gmpr_normalized, fill_color=groups_stage_early_median_late_color_list, line_color="black")
                p_stage_early_median_late.vbar(stage_early_median_late_list, 0.7, q1_stage_early_median_late.relative_abundance_gmpr_normalized, q2_stage_early_median_late.relative_abundance_gmpr_normalized, fill_color=groups_stage_early_median_late_color_list, line_color="black")

                # outliers
                if not out_stage_early_median_late.empty:
                    p_stage_early_median_late.circle(outx_stage_early_median_late, outy_stage_early_median_late, size=2, color="gray", fill_alpha=0.2)

                #p_stage_early_median_late.y_range.start = 0
                p_stage_early_median_late_annotaion = Label(x=1, y=max(clinical_relevance_query_set_df.loc[:, "relative_abundance_gmpr_normalized"]),
                                                text='p = ' + str(p_value_stage_early_median_late) )
                p_stage_early_median_late.x_range.range_padding = 0.1
                p_stage_early_median_late.yaxis.axis_label = "Relative abundance"
                # p_stage_early_median_late.xaxis.major_label_orientation = 1
                p_stage_early_median_late.xaxis.major_label_orientation = 3.414/12
                p_stage_early_median_late.xgrid.grid_line_color = None
                #p_stage_early_median_late.sizing_mode = "scale_width"

                script_stage_early_median_late, div_stage_early_median_late = components(p_stage_early_median_late)\
            
            except:
                print("")
    
    else:
        bacteria_expression_all_zero_message = bacteria + " (" + taxonomy_level + ") has zero expression in all samples of " + cancer_type + "."
    
    # print("bacteria_expression_all_zero_message: ", bacteria_expression_all_zero_message)


    all_statistics = {"cancer_type_id": cancer_type, \
                      "taxonomy_level": taxonomy_level, \
                      "bacteria": bacteria, \
                      "p_value_os": p_value_os, \
                      "p_value_wilcox_TPvsTN": p_value_wilcox_TPvsTN, \
                      "p_value_stage_early_median_late": p_value_stage_early_median_late, \
                      "p_value_race_simplified": p_value_race_simplified}
    #print(all_statistics)


    all_statistics_json = dumps(all_statistics) 
    #print(all_statistics_json)

    all_json = dumps(dict(data=all_statistics_json, \
                          script_os=script_os, div_os=div_os, \
                          script_wilcox_TPvsTN=script_wilcox_TPvsTN, div_wilcox_TPvsTN=div_wilcox_TPvsTN, \
                          script_stage_early_median_late=script_stage_early_median_late, div_stage_early_median_late=div_stage_early_median_late, \
                          script_race_simplified=script_race_simplified, div_race_simplified=div_race_simplified, \
                          bacteria_input_error_message=bacteria_input_error_message, \
                          bacteria_expression_all_zero_message=bacteria_expression_all_zero_message))
    # print(all_json)

    return JsonResponse(all_json, safe=False)
    





#  -----------------------

## bacteria coabundance_network 
def Analyses_bacteria_coabundance_network(request):
    title = "analyses-bacteria-coabundance-network"
    context = {"title": title}
    return render(request, "analyses/coabundance_network/bacteria_coabundance_network.html", context=context)


@csrf_exempt
def Ajax_bacteria_coabundance_network(request):
    taxonomy_level = request.POST.get("taxonomy_level_coabundance_input")
    cancer_type = request.POST.get("cancer_type_coabundance_input")
    bacteria = request.POST.get("bacteria_coabundance_input")
    p_value_cutoff = request.POST.get("p_value_input")
    network_displaying_option = request.POST.get("network_displaying_option_input")
    # print(taxonomy_level)
    # print(cancer_type)
    # print(bacteria)
    # print(p_value_cutoff)
    # print(network_displaying_option)


    bacteria_coabundance_query_set = View_bacteria_coabundance.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type).values()
    bacteria_coabundance_query_set_df = pd.DataFrame(list(bacteria_coabundance_query_set))
    #print(bacteria_coabundance_query_set_df.columns)
    # print(bacteria_coabundance_query_set_df.shape)
    #print(bacteria_coabundance_query_set_df.iloc[0:2, :])


    # check if input bacteria is existed in our coabundance network
    bacteria_input_error_message = ""
    p_value_cutoff_error_message = ""

    data = {}
    bacteria_name_available_list = sorted(list(set(list( bacteria_coabundance_query_set_df["taxonomy_level_bacterium_1"]) + list(bacteria_coabundance_query_set_df["taxonomy_level_bacterium_2"]) )))
    if bacteria not in bacteria_name_available_list:
        bacteria_input_error_message = "Not available"
        
    # print("bacteria_input_error_message: ", bacteria_input_error_message)


    # check if p-value input correctly (0 ≤ p-value ≤ 1)
    try:
        p_value_cutoff = float(request.POST.get("p_value_input"))
        if p_value_cutoff<0 or p_value_cutoff>1:
            p_value_cutoff_error_message = "0 ≤ p-value ≤ 1"
    except:
        p_value_cutoff_error_message = "0 ≤ p-value ≤ 1"

    
    
    # print("p_value_cutoff_error_message: ", p_value_cutoff_error_message)
    
    data["bacteria_input_error_message"] = bacteria_input_error_message
    data["p_value_cutoff_error_message"] = p_value_cutoff_error_message

    if bacteria_input_error_message != "" or p_value_cutoff_error_message != "":
        return JsonResponse(dumps(data), safe=False)
    


    bacteria_coabundance_pass_df = bacteria_coabundance_query_set_df[(bacteria_coabundance_query_set_df["p_value"] < p_value_cutoff) & ((bacteria_coabundance_query_set_df["taxonomy_level_bacterium_1"] == bacteria) | (bacteria_coabundance_query_set_df["taxonomy_level_bacterium_2"] == bacteria)) ]
    
    if p_value_cutoff == 1 or p_value_cutoff == 0:
        bacteria_coabundance_pass_df = bacteria_coabundance_query_set_df[(bacteria_coabundance_query_set_df["p_value"] <= p_value_cutoff) & ((bacteria_coabundance_query_set_df["taxonomy_level_bacterium_1"] == bacteria) | (bacteria_coabundance_query_set_df["taxonomy_level_bacterium_2"] == bacteria)) ]

    if bacteria_coabundance_pass_df.empty:
        p_value_cutoff_error_message = "Too stringent!!!"

    
    ## show edges between first neighbors of the queried bacteria
    if network_displaying_option == "all":
        queried_and_first_neighbor_list = sorted(list(set(list( bacteria_coabundance_pass_df["taxonomy_level_bacterium_1"]) + list(bacteria_coabundance_pass_df["taxonomy_level_bacterium_2"]) )))
        bacteria_coabundance_pass_df = bacteria_coabundance_query_set_df[ (bacteria_coabundance_query_set_df["taxonomy_level_bacterium_1"].isin(queried_and_first_neighbor_list)) & (bacteria_coabundance_query_set_df["taxonomy_level_bacterium_2"].isin(queried_and_first_neighbor_list)) ]
        # remove rows which p-value is empty (NA)
        bacteria_coabundance_pass_df.dropna(inplace = True)
    #     p = Plot(width=1200, height=675, x_range=Range1d(-1.5,1.5), y_range=Range1d(-1.5,1.5))
    # else:
    #     p = Plot(width=800, height=450, x_range=Range1d(-1.5,1.5), y_range=Range1d(-1.5,1.5))
    
    bacteria_coabundance_pass_df["taxonomy_level_id"] = taxonomy_level
    #print(bacteria_coabundance_pass_df.columns)
    #print(bacteria_coabundance_pass_df.info())
    # print(bacteria_coabundance_pass_df.shape)
    #print(bacteria_coabundance_pass_df.iloc[0:2, :])

    if bacteria_coabundance_pass_df.empty:
        p_value_cutoff_error_message = "Too stringent!!!"


    # plot
    G = nx.from_pandas_edgelist(bacteria_coabundance_pass_df, 'taxonomy_level_bacterium_1', 'taxonomy_level_bacterium_2', ['sparscc', 'p_value'])
    # print("G nodes: ", G.nodes())
    # print("G edges: ", G.edges())

    p = Plot(width=1200, height=675, x_range=Range1d(-1.5,1.5), y_range=Range1d(-1.5,1.5))
    # p = Plot(width=800, height=450,x_range=Range1d(-1.5,1.5), y_range=Range1d(-1.5,1.5))
    p.add_tools(BoxZoomTool(), SaveTool(), ResetTool())

    if network_displaying_option == "all":
        graph_renderer = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
    else:
        graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))


    graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])

    # add node labels, code source: https://stackoverflow.com/questions/49848966/how-to-add-permanent-name-labels-not-interactive-ones-on-nodes-for-a-networkx
    source = graph_renderer.node_renderer.data_source
    source.data['names'] = [str(x) for x in source.data['index']]
    # print("source: ", source)


    # create a transform that can extract the actual x,y positions
    code = """
        var result = new Float64Array(xs.length)
        for (var i = 0; i < xs.length; i++) {
            result[i] = provider.graph_layout[xs[i]][%s]
        }
        return result
    """
    xcoord = CustomJSTransform(v_func=code % "0", args=dict(provider=graph_renderer.layout_provider))
    ycoord = CustomJSTransform(v_func=code % "1", args=dict(provider=graph_renderer.layout_provider))

    # print("xcoord: ", xcoord)
    # print("ycoord: ", ycoord)
    

    # Use the transforms to supply coords to a LabelSet 
    labels = LabelSet(x=transform('index', xcoord),
                    y=transform('index', ycoord),
                    text='names', text_font_size="12px",
                    x_offset=-35, y_offset=8,
                    source=source, render_mode='canvas')

    p.add_layout(labels)
    

    edge_color_mapper = linear_cmap(field_name='sparscc', palette=Turbo256, low=-0.2,high=0.2)
    # print("edge_color_mapper: ", edge_color_mapper)

    ## modify the line color depending on the sparcc_cor
    graph_renderer.edge_renderer.glyph = MultiLine(line_color=edge_color_mapper, line_alpha=0.8, line_width=2)
    graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=edge_color_mapper, line_width=5)
    graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=edge_color_mapper, line_width=5)

    hover_edges = HoverTool(
    tooltips=[ ('Bacteria 1','@start'), ('Bacteria 2','@end'), ('sparscc','@sparscc'), ('p-value', '@p_value') ],
    renderers=[graph_renderer.edge_renderer], line_policy="interp"
    )
    
    graph_renderer.selection_policy = NodesAndLinkedEdges()
    graph_renderer.inspection_policy = EdgesAndLinkedNodes()

    p.renderers.append(graph_renderer)
    p.add_tools(hover_edges)


    edge_colorbar_mapper = LinearColorMapper(palette=Turbo256, low=-0.2, high=0.2)
    # print("edge_colorbar_mapper: ", edge_colorbar_mapper)

    color_bar = ColorBar(title="SparCC Cor", title_text_font_size='10px', title_text_align = 'left',
                         title_standoff=20, width=15, padding=20,
                         color_mapper=edge_colorbar_mapper, label_standoff=12, location=(0,0))


    p.add_layout(color_bar, "left")



    script, div = components(p)
    #print("script: ", script)
    #print("div: ", div)


    bacteria_coabundance_pass_json = bacteria_coabundance_pass_df.to_json(orient='records')

    all_json = dumps(dict(script=script, div=div, data=bacteria_coabundance_pass_json,
                          bacteria_input_error_message=bacteria_input_error_message,
                          p_value_cutoff_error_message=p_value_cutoff_error_message))
    #print(all_json)

    return JsonResponse(all_json, safe=False)

    



## bacteria bacteria_human_rna_network 
def Analyses_bacteria_human_rna_network(request):
    title = "analyses-bacteria-human-rna-network"
    context = {"title": title}
    return render(request, "analyses/bacteria_human_rna_network/bacteria_human_rna_network.html", context=context)


@csrf_exempt
def Ajax_bacteria_human_rna_network(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    cancer_type = request.POST.get("cancer_type_input")
    bacteria = request.POST.get("bacteria_input")
    z_score = request.POST.get("z_score_input")
    number_edge = request.POST.get("number_edge_input")

    # print(taxonomy_level)
    # print(cancer_type)
    # print(bacteria)
    # print(z_score)
    # print(number_edge)


    # check if input bacteria is existed in our bacteria human rna network
    bacteria_input_error_message = ""
    z_score_cutoff_error_message = ""
    number_edge_input_error_message = ""

    data = {}
    bacteria_name_available_list = list(View_bacteria_associated_gene_ontology.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type).order_by("name").values_list("name", flat=True).distinct("name"))

    if bacteria not in bacteria_name_available_list:
        bacteria_input_error_message = "Not available"
        
    # print("bacteria_input_error_message: ", bacteria_input_error_message)


    # check if z-score input correctly (0 < z-score ≤ 20)
    try:
        z_score = float(z_score)
        if z_score<0 or z_score>20:
            z_score_cutoff_error_message = "0 ≤ z-score ≤ 20"
    except:
        z_score_cutoff_error_message = "0 ≤ z-score ≤ 20"
    
    # print("z_score_cutoff_error_message: ", z_score_cutoff_error_message)


    # check if number edge input correctly (1 < N ≤ 200)
    try:
        number_edge = int(number_edge)
        if number_edge<1 or number_edge>200:
            number_edge_input_error_message = "1 ≤ N ≤ 200"
    except:
        number_edge_input_error_message = "1 ≤ N ≤ 200"
    
    # print("number_edge_input_error_message: ", number_edge_input_error_message)
    
    
    data["bacteria_input_error_message"] = bacteria_input_error_message
    data["z_score_cutoff_error_message"] = z_score_cutoff_error_message
    data["number_edge_input_error_message"] = number_edge_input_error_message
    
    
    if bacteria_input_error_message != "" or z_score_cutoff_error_message != "" or number_edge_input_error_message != "":
        return JsonResponse(dumps(data), safe=False)
    

    bacteria_human_rna_query_set = View_bacteria_human_rna_scc.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type, name=bacteria).values()
    bacteria_human_rna_query_set_df = pd.DataFrame(list(bacteria_human_rna_query_set))
    #print(bacteria_human_rna_query_set_df.columns)
    #print(bacteria_human_rna_query_set_df.shape)
    #print(bacteria_human_rna_query_set_df.iloc[0:2, :])


    bacteria_human_rna_pass_df = bacteria_human_rna_query_set_df[ abs(bacteria_human_rna_query_set_df["fisher_z_transformed_scc"]) >= z_score ]
    bacteria_human_rna_pass_df["taxonomy_level_id"] = taxonomy_level
    #print(bacteria_human_rna_pass_df.columns)
    #print(bacteria_human_rna_pass_df.info())
    #print(bacteria_human_rna_pass_df.shape)
    #print(bacteria_human_rna_pass_df.iloc[0:2, :])


    if bacteria_human_rna_pass_df.empty:
        z_score_cutoff_error_message = "Too strigent!!!"
        data["z_score_cutoff_error_message"] = z_score_cutoff_error_message
        return JsonResponse(dumps(data), safe=False)
    
    

    else:
        bacteria_human_rna_pass_copy_df = bacteria_human_rna_pass_df.copy()
        human_rna_name = bacteria_human_rna_pass_copy_df["human_rna"].str.split("|", expand=True)
        #print("head(human_rna_name): ")
        #print(human_rna_name.head())

        bacteria_human_rna_pass_copy_df["human_rna"] = human_rna_name[0]
        #print(bacteria_human_rna_pass_copy_df.head())

        bacteria_human_rna_pass_copy_filter_df  = bacteria_human_rna_pass_copy_df[bacteria_human_rna_pass_copy_df.human_rna != 'NA']
        # print(bacteria_human_rna_pass_copy_filter_df.shape)

        min_value = min(bacteria_human_rna_pass_copy_filter_df["fisher_z_transformed_scc"])
        max_value = max(bacteria_human_rna_pass_copy_filter_df["fisher_z_transformed_scc"])
        
        edge_color_cutoff = 2
        if abs(min_value)< max_value:
            edge_color_cutoff = max_value
        else:
            edge_color_cutoff = abs(min_value)
        
        # print("min_value", min_value)
        # print("max_value", max_value)
        # print("edge_color_cutoff", edge_color_cutoff)



        bacteria_human_rna_pass_copy_filter_sort_df = bacteria_human_rna_pass_copy_filter_df.sort_values(by=["fisher_z_transformed_scc"], ascending=False)

        if bacteria_human_rna_pass_copy_filter_sort_df.shape[0] > number_edge:
            bacteria_human_rna_pass_copy_filter_sort_plot_shown_df = pd.concat([bacteria_human_rna_pass_copy_filter_sort_df.head(number_edge//2), 
                                                                                bacteria_human_rna_pass_copy_filter_sort_df.tail(number_edge//2)])
        else:
            bacteria_human_rna_pass_copy_filter_sort_plot_shown_df = bacteria_human_rna_pass_copy_filter_sort_df
        
        # print("nrow(bacteria_human_rna_pass_copy_filter_sort_plot_shown_df): ", bacteria_human_rna_pass_copy_filter_sort_plot_shown_df.shape[0])


        

        # plot
        G = nx.from_pandas_edgelist(bacteria_human_rna_pass_copy_filter_sort_plot_shown_df, 'name', 'human_rna', ['fisher_z_transformed_scc', 'scc', 'scc_p_value'])
        #print("G nodes: ", G.nodes())
        #print("G edges: ", G.edges())

        
        p = Plot(width=1200, height=675, x_range=Range1d(-1.5,1.5), y_range=Range1d(-1.5,1.5))

        p.add_tools(BoxZoomTool(), SaveTool(), ResetTool())



        graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))
        #print("graph_renderer: ", graph_renderer)

        graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
        graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])

        # add node labels, code source: https://stackoverflow.com/questions/49848966/how-to-add-permanent-name-labels-not-interactive-ones-on-nodes-for-a-networkx
        source = graph_renderer.node_renderer.data_source
        source.data['names'] = [str(x) for x in source.data['index']]
        #print("source: ", source)


        # create a transform that can extract the actual x,y positions
        code = """
            var result = new Float64Array(xs.length)
            for (var i = 0; i < xs.length; i++) {
                result[i] = provider.graph_layout[xs[i]][%s]
            }
            return result
        """
        xcoord = CustomJSTransform(v_func=code % "0", args=dict(provider=graph_renderer.layout_provider))
        ycoord = CustomJSTransform(v_func=code % "1", args=dict(provider=graph_renderer.layout_provider))

        #print("xcoord: ", xcoord)
        #print("ycoord: ", ycoord)
        

        # Use the transforms to supply coords to a LabelSet 
        labels = LabelSet(x=transform('index', xcoord),
                        y=transform('index', ycoord),
                        text='names', text_font_size="12px",
                        x_offset=-35, y_offset=8,
                        source=source, render_mode='canvas')

        p.add_layout(labels)


        
        edge_color_mapper = linear_cmap(field_name='fisher_z_transformed_scc', palette=Turbo256, low=-edge_color_cutoff,high=edge_color_cutoff)
        #print("edge_color_mapper: ", edge_color_mapper)

        ## modify the line color depending on the fisher_z_transformed_scc
        graph_renderer.edge_renderer.glyph = MultiLine(line_color=edge_color_mapper, line_alpha=0.8, line_width=2)
        graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=edge_color_mapper, line_width=5)
        graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=edge_color_mapper, line_width=5)

        hover_edges = HoverTool(
        tooltips=[ ('Bacteria','@start'), ('RNA','@end'), ('SCC z-score','@fisher_z_transformed_scc'), ('SCC p-value', '@scc_p_value') ],
        renderers=[graph_renderer.edge_renderer], line_policy="interp"
        )

        graph_renderer.selection_policy = NodesAndLinkedEdges()
        graph_renderer.inspection_policy = EdgesAndLinkedNodes()

        p.renderers.append(graph_renderer)
        p.add_tools(hover_edges)


        edge_colorbar_mapper = LinearColorMapper(palette=Turbo256, low=-edge_color_cutoff, high=edge_color_cutoff)
        #print("edge_colorbar_mapper: ", edge_colorbar_mapper)

        color_bar = ColorBar(title="Fisher’s z-transformed SCC", title_text_font_size='10px', title_text_align = 'left',
                            title_standoff=20, width=15, padding=20,
                            color_mapper=edge_colorbar_mapper, label_standoff=12, location=(0,0))

        p.add_layout(color_bar, "left")




        script, div = components(p)
        #print("script: ", script)
        #print("div: ", div)


        bacteria_human_rna_pass_json = bacteria_human_rna_pass_copy_filter_sort_df.to_json(orient='records')


        all_json = dumps(dict(script=script, div=div, data=bacteria_human_rna_pass_json, \
                            bacteria_input_error_message=bacteria_input_error_message, \
                            z_score_cutoff_error_message=z_score_cutoff_error_message, \
                            number_edge_input_error_message=number_edge_input_error_message))
        #print(all_json)

        return JsonResponse(all_json, safe=False)






## bacteria gene_ontology 
def Analyses_bacteria_associated_gene_ontology(request):
    title = "analyses-bacteria-gene-ontology"
    context = {"title": title}
    return render(request, "analyses/gene_ontology/bacteria_associated_gene_ontology.html", context=context)


@csrf_exempt
def Ajax_bacteria_associated_gene_ontology(request):
    taxonomy_level = request.POST.get("taxonomy_level_input")
    cancer_type = request.POST.get("cancer_type_input")
    bacteria = request.POST.get("bacteria_input")
    gene_set = request.POST.get("gene_set_input")
    adjust_p_value_cutoff = request.POST.get("adjust_p_value_input")
    number_term = request.POST.get("number_term_input")
    # print(taxonomy_level)
    # print(cancer_type)
    # print(bacteria)
    # print(gene_set)
    # print(adjust_p_value_cutoff)
    # print(number_term)



    # check if input bacteria is existed in our bacteria human rna network
    bacteria_input_error_message = ""
    adjust_p_value_cutoff_error_message = ""
    number_term_input_error_message = ""

    data = {}
    bacteria_name_available_list = list(View_bacteria_associated_gene_ontology.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type).order_by("name").values_list("name", flat=True).distinct("name"))

    if bacteria not in bacteria_name_available_list:
        bacteria_input_error_message = "Not available"
        
    # print("bacteria_input_error_message: ", bacteria_input_error_message)


    # check if Adjusted p-value cutoff: input correctly
    try:
        adjust_p_value_cutoff = float(adjust_p_value_cutoff)
        if adjust_p_value_cutoff<0 or adjust_p_value_cutoff>1:
            adjust_p_value_cutoff_error_message = "0 ≤ Adj p-value ≤ 1"
    except:
        adjust_p_value_cutoff_error_message = "0 ≤ Adj p-value ≤ 1"
    
    # print("adjust_p_value_cutoff_error_message: ", adjust_p_value_cutoff_error_message)


    # check if number term input correctly (1 ≤ N ≤ 500)
    try:
        number_term = int(number_term)
        if number_term<1 or number_term>500:
            number_term_input_error_message = "1 ≤ N ≤ 500"
    except:
        number_term_input_error_message = "1 ≤ N ≤ 500"
    
    # print("number_term_input_error_message: ", number_term_input_error_message)
    
    
    data["bacteria_input_error_message"] = bacteria_input_error_message
    data["adjust_p_value_cutoff_error_message"] = adjust_p_value_cutoff_error_message
    data["number_term_input_error_message"] = number_term_input_error_message
    
    
    if bacteria_input_error_message != "" or adjust_p_value_cutoff_error_message != "" or number_term_input_error_message != "":
        return JsonResponse(dumps(data), safe=False)


    bacteria_associated_gene_ontology_query_set = View_bacteria_associated_gene_ontology.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type, name=bacteria).values()
    if gene_set == "kegg_pathway":
        bacteria_associated_gene_ontology_query_set = View_bacteria_associated_kegg_pathway.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type, name=bacteria).values()
    if gene_set == "reactome_pathway":
        bacteria_associated_gene_ontology_query_set = View_bacteria_associated_reactome_pathway.objects.filter(taxonomy_level_id=taxonomy_level[0], cancer_type_id=cancer_type, name=bacteria).values()
    
    bacteria_associated_gene_ontology_query_set_df = pd.DataFrame(list(bacteria_associated_gene_ontology_query_set))
    # print(bacteria_associated_gene_ontology_query_set_df.columns)
    # print(bacteria_associated_gene_ontology_query_set_df.shape)
    # print(bacteria_associated_gene_ontology_query_set_df.dtypes)
    # print(bacteria_associated_gene_ontology_query_set_df.iloc[0:2, :])


    bacteria_associated_gene_ontology_pass_query_set_df = bacteria_associated_gene_ontology_query_set_df[ bacteria_associated_gene_ontology_query_set_df["adjust_p_value"] < adjust_p_value_cutoff ]

    if adjust_p_value_cutoff == 1 or adjust_p_value_cutoff == 0:
        bacteria_associated_gene_ontology_pass_query_set_df = bacteria_associated_gene_ontology_query_set_df[ bacteria_associated_gene_ontology_query_set_df["adjust_p_value"] <= adjust_p_value_cutoff ]
    
    # print(bacteria_associated_gene_ontology_pass_query_set_df.shape)

    bacteria_associated_gene_ontology_pass_query_set_df["taxonomy_level_id"] = taxonomy_level

    if bacteria_associated_gene_ontology_pass_query_set_df.empty:
        adjust_p_value_cutoff_error_message = "Too stringent!!!"
        data["adjust_p_value_cutoff_error_message"] = adjust_p_value_cutoff_error_message
        # print("adjust_p_value_cutoff_error_message: ", adjust_p_value_cutoff_error_message)
        return JsonResponse(dumps(data), safe=False)
    else:
        bacteria_associated_gene_ontology_pass_sort_query_set_df = bacteria_associated_gene_ontology_pass_query_set_df.sort_values(by=["nes"], ascending=False)

        go_term = bacteria_associated_gene_ontology_pass_sort_query_set_df["term"].tolist()
        nes_score = bacteria_associated_gene_ontology_pass_sort_query_set_df["nes"].tolist()

        if len(go_term) > number_term:
            go_term_plot_shown = go_term[0:number_term]
            nes_score_plot_shown = nes_score[0:number_term]
        
        else:
            go_term_plot_shown = go_term
            nes_score_plot_shown = nes_score
        

        # print("len(go_term): ", len(go_term))
        # print("len(nes_score): ", len(nes_score))
        # print("len(go_term_plot_shown): ", len(go_term_plot_shown))
        # print("len(nes_score_plot_shown): ", len(nes_score_plot_shown))

        # print("min(nes_score_plot_shown): ", min(nes_score_plot_shown))
        # print("max(nes_score_plot_shown): ", max(nes_score_plot_shown))

        width_setting = 1000
        height_setting = 15 * len(go_term_plot_shown)

        p = figure(plot_width=width_setting, plot_height=height_setting, 
                y_range=go_term_plot_shown[::-1],
                toolbar_location="right")


        p.hbar(go_term_plot_shown[::-1], right=nes_score_plot_shown[::-1], height=0.4)


        p.xaxis.axis_label = "NES"
        if nes_score_plot_shown[-1] >= 0:
            p.x_range.start = 0



        script, div = components(p)
        # print("script: ", script)
        # print("div: ", div)

        bacteria_associated_gene_ontology_pass_sort_json = bacteria_associated_gene_ontology_pass_sort_query_set_df.to_json(orient='records')

        all_json = dumps(dict(script=script, div=div, data=bacteria_associated_gene_ontology_pass_sort_json, \
                              bacteria_input_error_message=bacteria_input_error_message, \
                              adjust_p_value_cutoff_error_message=adjust_p_value_cutoff_error_message, \
                              number_term_input_error_message=number_term_input_error_message))
        # print(all_json)

        return JsonResponse(all_json, safe=False)

    









