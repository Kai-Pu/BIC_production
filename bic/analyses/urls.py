from django.urls import path
from . import views


urlpatterns = [
    path('bacteria_name_ajax/', views.Autocomplete_search_bacteria_name, name="autocomplete_search_bacteria_name"),
    path('bacteria_name_go_ajax/', views.Autocomplete_search_bacteria_name_go, name="autocomplete_search_bacteria_name_go"),
    path('bacteria_name_coabundance_ajax/', views.Autocomplete_search_bacteria_name_coabundance, name="autocomplete_search_bacteria_name_coabundance"),
    path('abundance/', views.Analyses_bacteria_abundance, name="analyses_bacteria_abundance"),
    path('abundance_ajax/', views.Ajax_bacteria_abundance, name="ajax_bacteria_abundance"),
    path('diversity/', views.Analyses_bacteria_diversity, name="analyses_bacteria_diversity"),
    path('diversity_ajax/', views.Ajax_bacteria_diversity, name="ajax_bacteria_diversity"),
    path('composition/', views.Analyses_bacteria_composition, name="analyses_bacteria_composition"),
    path('composition_ajax/', views.Ajax_bacteria_composition, name="ajax_bacteria_composition"),
    path('clinical_relevance/', views.Analyses_bacteria_clinical_relevance, name="analyses_bacteria_clinical_relevance"),
    path('clinical_relevance_ajax/', views.Ajax_bacteria_clinical_relevance, name="ajax_bacteria_clinical_relevance"),
    path('coabundance_network/', views.Analyses_bacteria_coabundance_network, name="analyses_bacteria_coabundance_network"),
    path('coabundance_network_ajax/', views.Ajax_bacteria_coabundance_network, name="ajax_bacteria_coabundance_network"),
    path('bacteria_human_rna_network/', views.Analyses_bacteria_human_rna_network, name="analyses_bacteria_human_rna_network"),
    path('bacteria_human_rna_network_ajax/', views.Ajax_bacteria_human_rna_network, name="ajax_bacteria_human_rna_network"),
    path('gene_ontology/', views.Analyses_bacteria_associated_gene_ontology, name="analyses_bacteria_associated_gene_ontology"),
    path('gene_ontology_ajax/', views.Ajax_bacteria_associated_gene_ontology, name="ajax_bacteria_associated_gene_ontology"),
]




