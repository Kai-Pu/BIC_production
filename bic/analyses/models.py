from django.db import models

# Create your models here.



class Taxonomy_level(models.Model):
    taxonomy_level_id = models.CharField(max_length=3, primary_key=True)
    level = models.CharField(max_length=10)

    class Meta:
        managed = False
        db_table = "taxonomy_level"


class Taxonomy_level_bacterium(models.Model):
    taxonomy_level_bacterium_id = models.IntegerField(primary_key=True)
    taxonomy_string_name = models.TextField()
    name = models.CharField(max_length=100)
    taxonomy_level = models.ForeignKey(Taxonomy_level, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'taxonomy_level_bacterium'


class Cancer_type(models.Model):
    cancer_type_id = models.CharField(max_length=5, primary_key=True)
    cancer_full_name = models.CharField(max_length=100)
    color_code_hex = models.CharField(max_length=10, default="")
    color_red = models.IntegerField(default=0)
    color_green = models.IntegerField(default=0)
    color_blue = models.IntegerField(default=0)

    class Meta:
        managed = False
        db_table = 'cancer_type'

class Patient(models.Model):
    patient_id = models.IntegerField(primary_key=True)
    patient_barcode = models.CharField(max_length=20)
    gender = models.CharField(max_length=10, null=True)
    age_at_initial_pathologic_diagnosis = models.IntegerField(null=True)
    race_simplified = models.CharField(max_length=30)
    stage = models.CharField(max_length=10, null=True)
    stage_early_median_late = models.CharField(max_length=20, null=True)
    stage_early_late = models.CharField(max_length=20, default="", null=True)
    os = models.IntegerField(null=True)
    os_time = models.IntegerField(null=True)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'patient'



class Sample_type(models.Model):
    sample_type_id = models.CharField(max_length=5, primary_key=True)
    sample_type = models.CharField(max_length=50)

    class Meta:
        managed = False
        db_table = 'sample_type'





class Patient_cancer_sample_type(models.Model):
    patient_cancer_sample_type_id = models.IntegerField(primary_key=True)
    sample_barcode = models.CharField(max_length=30)
    aliquot_barcode_mirna_seq = models.CharField(max_length=40)
    aliquot_barcode_rna_seq = models.CharField(max_length=40, null=True)
    sub_cluster = models.IntegerField(null=True)
    patient = models.ForeignKey(Patient, on_delete=models.CASCADE)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)
    sample_type = models.ForeignKey(Sample_type, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'patient_cancer_sample_type'



class Bacteria_expression(models.Model):
    bacteria_expression_id = models.IntegerField(primary_key=True)
    relative_abundance_gmpr_normalized = models.FloatField()
    count_gmpr_normalized = models.FloatField()
    taxonomy_level_bacterium = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    patient_cancer_sample_type = models.ForeignKey(Patient_cancer_sample_type, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'bacteria_expression'



## create view for bacteria expression
class View_bacteria_expression(models.Model):
    bacteria_expression_id = models.IntegerField(primary_key=True)
    taxonomy_level_bacterium_id = models.IntegerField()
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    taxonomy_string_name = models.TextField()
    aliquot_barcode_mirna_seq = models.CharField(max_length=40)
    cancer_type_id = models.CharField(max_length=5)
    sample_type_id = models.CharField(max_length=5)
    relative_abundance_gmpr_normalized = models.FloatField()
    count_gmpr_normalized = models.FloatField()

    class Meta:
        managed = False
        db_table = 'view_bacteria_expression'







class Sample_bacteria_diversity(models.Model):
    sample_bacteria_diversity_id = models.IntegerField(primary_key=True)
    diversity_method = models.CharField(max_length=50)
    diversity_value = models.FloatField()
    taxonomy_level = models.ForeignKey(Taxonomy_level, on_delete=models.CASCADE)
    patient_cancer_sample_type = models.ForeignKey(Patient_cancer_sample_type, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'sample_bacteria_diversity'


## create view for sample_bacteria_diversity
class View_sample_bacteria_diversity(models.Model):
    sample_bacteria_diversity_id = models.IntegerField(primary_key=True)
    taxonomy_level_id = models.CharField(max_length=3)
    aliquot_barcode_mirna_seq = models.CharField(max_length=40)
    cancer_type_id = models.CharField(max_length=5)
    sample_type_id = models.CharField(max_length=5)
    diversity_method = models.CharField(max_length=50)
    diversity_value = models.FloatField()

    class Meta:
        managed = False
        db_table = 'view_sample_bacteria_diversity'





## create view for view_bacteria_clinical_relevance
class View_bacteria_clinical_relevance(models.Model):
    bacteria_expression_id = models.IntegerField(primary_key=True)
    taxonomy_level_bacterium_id = models.IntegerField()
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    taxonomy_string_name = models.TextField()
    aliquot_barcode_mirna_seq = models.CharField(max_length=40)
    cancer_type_id = models.CharField(max_length=5)
    sample_type_id = models.CharField(max_length=5)
    relative_abundance_gmpr_normalized = models.FloatField()
    count_gmpr_normalized = models.FloatField()
    sub_cluster = models.IntegerField(null=True)
    patient_barcode = models.CharField(max_length=20)
    gender = models.CharField(max_length=10, null=True)
    age_at_initial_pathologic_diagnosis = models.IntegerField(null=True)
    race_simplified = models.CharField(max_length=30)
    stage = models.CharField(max_length=10, null=True)
    stage_early_median_late = models.CharField(max_length=20, null=True)
    stage_early_late = models.CharField(max_length=20, default="", null=True)
    os = models.IntegerField(null=True)
    os_time = models.IntegerField(null=True)

    class Meta:
        managed = False
        db_table = 'view_bacteria_clinical_relevance'




class Bacteria_coabundance(models.Model):
    bacteria_coabundance_id = models.IntegerField(primary_key=True)
    sparscc = models.FloatField(null=True)
    p_value = models.FloatField(null=True)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)
    taxonomy_level_bacterium_id_1 = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    taxonomy_level_bacterium_id_2 = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE, related_name="+")

    class Meta:
        managed = False
        db_table = 'bacteria_coabundance'



## create view for bacteria_coabundance
class View_bacteria_coabundance(models.Model):
    bacteria_coabundance_id = models.IntegerField(primary_key=True)
    cancer_type_id = models.CharField(max_length=5)
    taxonomy_level_id = models.CharField(max_length=3)
    taxonomy_level_bacterium_1 = models.CharField(max_length=100)
    taxonomy_level_bacterium_2 = models.CharField(max_length=100)
    sparscc = models.FloatField(null=True)
    p_value = models.FloatField(null=True)

    class Meta:
        managed = False
        db_table = 'view_bacteria_coabundance'





class Bacteria_human_rna_scc(models.Model):
    bacteria_human_rna_scc_id = models.IntegerField(primary_key=True)
    human_rna = models.CharField(max_length=50)
    scc = models.FloatField(null=True)
    fisher_z_transformed_scc = models.FloatField(null=True)
    scc_p_value = models.FloatField(null=True)
    taxonomy_level_bacterium = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)

    class Meta:
        managed = False
        db_table = 'bacteria_human_rna_scc'






## create view for bacteria_human_rna_scc
class View_bacteria_human_rna_scc(models.Model):
    bacteria_human_rna_scc_id = models.IntegerField(primary_key=True)
    cancer_type_id = models.CharField(max_length=5)
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    human_rna = models.CharField(max_length=50)
    scc = models.FloatField(null=True)
    fisher_z_transformed_scc = models.FloatField(null=True)
    scc_p_value = models.FloatField(null=True)

    class Meta:
        managed = False
        db_table = 'view_bacteria_human_rna_scc'




class Gene_ontology_biological_process(models.Model):
    gene_ontology_biological_process_id = models.IntegerField(primary_key=True)
    term = models.CharField(max_length=200)
    size = models.IntegerField()
    gene_symbol = models.TextField()
    
    class Meta:
        managed = False
        db_table = 'gene_ontology_biological_process'




class Bacteria_associated_gene_ontology(models.Model):
    bacteria_associated_gene_ontology_id = models.IntegerField(primary_key=True)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)
    taxonomy_level_bacterium = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    gene_ontology_biological_process = models.ForeignKey(Gene_ontology_biological_process, null=True, blank=True, on_delete=models.CASCADE)
    
    class Meta:
        managed = False
        db_table = 'bacteria_associated_gene_ontology'







## create view for bacteria_associated_gene_ontology
class View_bacteria_associated_gene_ontology(models.Model):
    bacteria_associated_gene_ontology_id = models.IntegerField(primary_key=True)
    cancer_type_id = models.CharField(max_length=5)
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    term = models.CharField(max_length=200)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    
    class Meta:
        managed = False
        db_table = 'view_bacteria_associated_gene_ontology'






class Kegg_pathway(models.Model):
    kegg_pathway_id = models.IntegerField(primary_key=True)
    term = models.CharField(max_length=200)
    size = models.IntegerField()
    gene_symbol = models.TextField()
    
    class Meta:
        managed = False
        db_table = 'kegg_pathway'




class Bacteria_associated_kegg_pathway(models.Model):
    bacteria_associated_kegg_pathway_id = models.IntegerField(primary_key=True)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)
    taxonomy_level_bacterium = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    kegg_pathway = models.ForeignKey(Kegg_pathway, null=True, blank=True, on_delete=models.CASCADE)
    
    class Meta:
        managed = False
        db_table = 'bacteria_associated_kegg_pathway'







## create view for bacteria_associated_kegg_pathway
class View_bacteria_associated_kegg_pathway(models.Model):
    bacteria_associated_kegg_pathway_id = models.IntegerField(primary_key=True)
    cancer_type_id = models.CharField(max_length=5)
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    term = models.CharField(max_length=200)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    
    class Meta:
        managed = False
        db_table = 'view_bacteria_associated_kegg_pathway'





class Reactome_pathway(models.Model):
    reactome_pathway_id = models.IntegerField(primary_key=True)
    term = models.CharField(max_length=200)
    size = models.IntegerField()
    gene_symbol = models.TextField()
    
    class Meta:
        managed = False
        db_table = 'reactome_pathway'




class Bacteria_associated_reactome_pathway(models.Model):
    bacteria_associated_reactome_pathway_id = models.IntegerField(primary_key=True)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    cancer_type = models.ForeignKey(Cancer_type, on_delete=models.CASCADE)
    taxonomy_level_bacterium = models.ForeignKey(Taxonomy_level_bacterium, on_delete=models.CASCADE)
    reactome_pathway = models.ForeignKey(Reactome_pathway, null=True, blank=True, on_delete=models.CASCADE)
    
    class Meta:
        managed = False
        db_table = 'bacteria_associated_reactome_pathway'







## create view for bacteria_associated_reactome_pathway
class View_bacteria_associated_reactome_pathway(models.Model):
    bacteria_associated_reactome_pathway_id = models.IntegerField(primary_key=True)
    cancer_type_id = models.CharField(max_length=5)
    taxonomy_level_id = models.CharField(max_length=3)
    name = models.CharField(max_length=100)
    term = models.CharField(max_length=200)
    p_value = models.FloatField(null=True)
    adjust_p_value = models.FloatField(null=True)
    es = models.FloatField(null=True)
    nes = models.FloatField(null=True)
    n_more_extreme = models.IntegerField(null=True)
    size = models.IntegerField(null=True)
    gene_symbol = models.TextField(null=True)
    
    class Meta:
        managed = False
        db_table = 'view_bacteria_associated_reactome_pathway'





class Statistics(models.Model):
    cancer_type = models.CharField(max_length=50, primary_key=True)
    number_of_tumor_sample = models.IntegerField()
    number_of_normal_sample_adjacent_to_tumor = models.IntegerField()
    survival_data_available = models.CharField(max_length=15, null=True)
    race_data_available = models.CharField(max_length=15, null=True)
    stage_data_available = models.CharField(max_length=15, null=True)

    def __unicode__(self):
            return self.name
            
    class Meta:
        managed = False
        db_table = 'statistics'


