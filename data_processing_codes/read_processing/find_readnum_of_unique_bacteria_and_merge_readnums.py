########################################################################
## usage: python3 find_readnum_of_unique_bacteria.py 
## aim: find reads info
## total_reads_unique_genus, total_reads_unique_family, total_reads_unique_order, total_reads_unique_class, total_reads_unique_phylum
## and merge into a total_reads_info.txt
## author: Kai-Pu Chen
## date: 20220908
#######################################################################

import pandas as pd
import statistics

statistical_summary_output_path = "/TCGA-BRCA/reads_info/processed_reads_statistical_summary.txt"
processed_reads_info_output_path = "/TCGA-BRCA/reads_info/processed_reads_info.txt"


# ## test of subset
# count_table_path = "/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/test3/first6column_phylum.txt"
# bacteria_count = pd.read_csv(count_table_path, sep="\t", index_col=0)
# total_counts_list = bacteria_count.sum().to_list()
# print(bacteria_count.columns)
# print("bacteria_count.iloc[0:5, 0:5]: ", bacteria_count.iloc[0:5, 0:5])
# print(bacteria_count.shape)
# print("len(total_counts_list): ", len(total_counts_list))
# print("total_counts_list[:5]: ", total_counts_list[0:5])



## turn total count into list
count_dir = "/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/"
genus_count_table_path = count_dir + "TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_Genus.txt"
family_count_table_path = count_dir + "TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_Family.txt"
order_count_table_path = count_dir + "TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_Order.txt"
class_count_table_path = count_dir + "TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_Class.txt"
phylum_count_table_path = count_dir + "TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_Phylum.txt"

total_reads_tcga_bam_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam.txt"
total_reads_tcga_bam_human_mapped_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam_human_mapped.txt"
total_reads_human_unmapped_fastq_path = "/TCGA-BRCA/reads_info/total_reads_human_unmapped_fastq.txt"
total_reads_srnanalyzer_processed_fa_path = "/TCGA-BRCA/reads_info/total_reads_srnanalyzer_processed_fa.txt"


statistical_summary = []
# statistical_summary.append(["total_reads", "mean", "sd", "median"])


genus_bacteria_count_table_df = pd.read_csv(genus_count_table_path, sep="\t", index_col=0)
print("genus_bacteria_count_table_df.iloc[0:3, 0:3]: ", genus_bacteria_count_table_df.iloc[0:3, 0:3])
print("genus_bacteria_count_table_df.shape: ", genus_bacteria_count_table_df.shape)
genus_total_counts_list = genus_bacteria_count_table_df.sum().to_list()
print("len(genus_total_counts_list): ", len(genus_total_counts_list))
print("genus_total_counts_list[:5]: ", genus_total_counts_list[0:5])
print("median: ", statistics.median(genus_total_counts_list))
print("mean: ", statistics.mean(genus_total_counts_list))
print("sd: ", statistics.stdev(genus_total_counts_list))
print("variance: ", statistics.variance(genus_total_counts_list))




family_bacteria_count_table_df = pd.read_csv(family_count_table_path, sep="\t", index_col=0)
print("family_bacteria_count_table_df.iloc[0:3, 0:3]: ", family_bacteria_count_table_df.iloc[0:3, 0:3])
print("family_bacteria_count_table_df.shape: ", family_bacteria_count_table_df.shape)
family_total_counts_list = family_bacteria_count_table_df.sum().to_list()
print("len(family_total_counts_list): ", len(family_total_counts_list))
print("family_total_counts_list[:5]: ", family_total_counts_list[0:5])
print("median: ", statistics.median(family_total_counts_list))
print("mean: ", statistics.mean(family_total_counts_list))
print("sd: ", statistics.stdev(family_total_counts_list))
print("variance: ", statistics.variance(family_total_counts_list))



order_bacteria_count_table_df = pd.read_csv(order_count_table_path, sep="\t", index_col=0)
print("order_bacteria_count_table_df.iloc[0:3, 0:3]: ", order_bacteria_count_table_df.iloc[0:3, 0:3])
print("order_bacteria_count_table_df.shape: ", order_bacteria_count_table_df.shape)
order_total_counts_list = order_bacteria_count_table_df.sum().to_list()
print("len(order_total_counts_list): ", len(order_total_counts_list))
print("order_total_counts_list[:5]: ", order_total_counts_list[0:5])
print("median: ", statistics.median(order_total_counts_list))
print("mean: ", statistics.mean(order_total_counts_list))
print("sd: ", statistics.stdev(order_total_counts_list))
print("variance: ", statistics.variance(order_total_counts_list))



class_bacteria_count_table_df = pd.read_csv(class_count_table_path, sep="\t", index_col=0)
print("class_bacteria_count_table_df.iloc[0:3, 0:3]: ", class_bacteria_count_table_df.iloc[0:3, 0:3])
print("class_bacteria_count_table_df.shape: ", class_bacteria_count_table_df.shape)
class_total_counts_list = class_bacteria_count_table_df.sum().to_list()
print("len(class_total_counts_list): ", len(class_total_counts_list))
print("class_total_counts_list[:5]: ", class_total_counts_list[0:5])
print("median: ", statistics.median(class_total_counts_list))
print("mean: ", statistics.mean(class_total_counts_list))
print("sd: ", statistics.stdev(class_total_counts_list))
print("variance: ", statistics.variance(class_total_counts_list))



phylum_bacteria_count_table_df = pd.read_csv(phylum_count_table_path, sep="\t", index_col=0)
print("phylum_bacteria_count_table_df.iloc[0:3, 0:3]: ", phylum_bacteria_count_table_df.iloc[0:3, 0:3])
print("phylum_bacteria_count_table_df.shape: ", phylum_bacteria_count_table_df.shape)
phylum_total_counts_list = phylum_bacteria_count_table_df.sum().to_list()
print("len(phylum_total_counts_list): ", len(phylum_total_counts_list))
print("phylum_total_counts_list[:5]: ", phylum_total_counts_list[0:5])
print("median: ", statistics.median(phylum_total_counts_list))
print("mean: ", statistics.mean(phylum_total_counts_list))
print("sd: ", statistics.stdev(phylum_total_counts_list))
print("variance: ", statistics.variance(phylum_total_counts_list))





## input computed total counts info of total_reads_tcga_bam, total_reads_tcga_bam_human_mapped, total_reads_human_unmapped_fastq, total_reads_srnanalyzer_processed_fa
total_reads_tcga_bam_df = pd.read_csv(total_reads_tcga_bam_path, sep="\t", header=None)
print("total_reads_tcga_bam_df.iloc[0:3]: ", total_reads_tcga_bam_df.iloc[0:3])
print("total_reads_tcga_bam_df.shape: ", total_reads_tcga_bam_df.shape)
total_reads_tcga_bam_df_list = total_reads_tcga_bam_df.iloc[:, 0].to_list()
print("len(total_reads_tcga_bam_df_list): ", len(total_reads_tcga_bam_df_list))
print("total_reads_tcga_bam_df_list[:5]: ", total_reads_tcga_bam_df_list[0:5])
# print("total_reads_tcga_bam_df_list[-5:]: ", total_reads_tcga_bam_df_list[-5:])
print("median: ", statistics.median(total_reads_tcga_bam_df_list))
print("mean: ", statistics.mean(total_reads_tcga_bam_df_list))
print("sd: ", statistics.stdev(total_reads_tcga_bam_df_list))
print("variance: ", statistics.variance(total_reads_tcga_bam_df_list))



total_reads_tcga_bam_human_mapped_df = pd.read_csv(total_reads_tcga_bam_human_mapped_path, sep="\t", header=None)
print("total_reads_tcga_bam_human_mapped_df.iloc[0:3]: ", total_reads_tcga_bam_human_mapped_df.iloc[0:3])
print("total_reads_tcga_bam_human_mapped_df.shape: ", total_reads_tcga_bam_human_mapped_df.shape)
total_reads_tcga_bam_human_mapped_df_list = total_reads_tcga_bam_human_mapped_df.iloc[:, 0].to_list()
print("len(total_reads_tcga_bam_human_mapped_df_list): ", len(total_reads_tcga_bam_human_mapped_df_list))
print("total_reads_tcga_bam_human_mapped_df_list[:5]: ", total_reads_tcga_bam_human_mapped_df_list[0:5])
# print("total_reads_tcga_bam_human_mapped_df_list[-5:]: ", total_reads_tcga_bam_human_mapped_df_list[-5:])
print("median: ", statistics.median(total_reads_tcga_bam_human_mapped_df_list))
print("mean: ", statistics.mean(total_reads_tcga_bam_human_mapped_df_list))
print("sd: ", statistics.stdev(total_reads_tcga_bam_human_mapped_df_list))
print("variance: ", statistics.variance(total_reads_tcga_bam_human_mapped_df_list))



total_reads_human_unmapped_fastq_df = pd.read_csv(total_reads_human_unmapped_fastq_path, sep="\t", header=None)
print("total_reads_human_unmapped_fastq_df.iloc[0:3]: ", total_reads_human_unmapped_fastq_df.iloc[0:3])
print("total_reads_human_unmapped_fastq_df.shape: ", total_reads_human_unmapped_fastq_df.shape)
total_reads_human_unmapped_fastq_df_list = total_reads_human_unmapped_fastq_df.iloc[:, 0].to_list()
print("len(total_reads_human_unmapped_fastq_df_list): ", len(total_reads_human_unmapped_fastq_df_list))
print("total_reads_human_unmapped_fastq_df_list[:5]: ", total_reads_human_unmapped_fastq_df_list[0:5])
# print("total_reads_human_unmapped_fastq_df_list[-5:]: ", total_reads_human_unmapped_fastq_df_list[-5:])
print("median: ", statistics.median(total_reads_human_unmapped_fastq_df_list))
print("mean: ", statistics.mean(total_reads_human_unmapped_fastq_df_list))
print("sd: ", statistics.stdev(total_reads_human_unmapped_fastq_df_list))
print("variance: ", statistics.variance(total_reads_human_unmapped_fastq_df_list))



total_reads_srnanalyzer_processed_fa_df = pd.read_csv(total_reads_srnanalyzer_processed_fa_path, sep="\t", header=None)
print("total_reads_srnanalyzer_processed_fa_df.iloc[0:3]: ", total_reads_srnanalyzer_processed_fa_df.iloc[0:3])
print("total_reads_srnanalyzer_processed_fa_df.shape: ", total_reads_srnanalyzer_processed_fa_df.shape)
total_reads_srnanalyzer_processed_fa_df_list = total_reads_srnanalyzer_processed_fa_df.iloc[:, 0].to_list()
print("len(total_reads_srnanalyzer_processed_fa_df_list): ", len(total_reads_srnanalyzer_processed_fa_df_list))
print("total_reads_srnanalyzer_processed_fa_df_list[:5]: ", total_reads_srnanalyzer_processed_fa_df_list[0:5])
# print("total_reads_srnanalyzer_processed_fa_df_list[-5:]: ", total_reads_srnanalyzer_processed_fa_df_list[-5:])
print("median: ", statistics.median(total_reads_srnanalyzer_processed_fa_df_list))
print("mean: ", statistics.mean(total_reads_srnanalyzer_processed_fa_df_list))
print("sd: ", statistics.stdev(total_reads_srnanalyzer_processed_fa_df_list))
print("variance: ", statistics.variance(total_reads_srnanalyzer_processed_fa_df_list))



number_of_sample = total_reads_srnanalyzer_processed_fa_df.shape[0]
print("\n\nnumber_of_sample: ", number_of_sample)

processed_reads_info_dict = { "patient_cancer_sample_type_id": list(range(1, number_of_sample+1, 1)),
                            "total_reads_tcga_bam": total_reads_tcga_bam_df_list,
                            "total_reads_tcga_bam_human_mapped": total_reads_tcga_bam_human_mapped_df_list,
                            "total_reads_human_unmapped_fastq": total_reads_human_unmapped_fastq_df_list,
                            "total_reads_srnanalyzer_processed_fa": total_reads_srnanalyzer_processed_fa_df_list,
                            "total_reads_phylum": phylum_total_counts_list, 
                            "total_reads_class": class_total_counts_list, 
                            "total_reads_order": order_total_counts_list, 
                            "total_reads_family": family_total_counts_list, 
                            "total_reads_genus": genus_total_counts_list }

processed_reads_info_df = pd.DataFrame(processed_reads_info_dict)
print("processed_reads_info_df.iloc[0:3]: ", processed_reads_info_df.iloc[0:3])
print("processed_reads_info_df.shape: ", processed_reads_info_df.shape)
print("processed_reads_info_df.columns: ", processed_reads_info_df.columns)

# statistical_summary: ["total_reads", "mean", "sd", "median"]
statistical_summary_dict = {"total_reads": ["mean", "sd", "median"], 
                            "tcga_bam":
                            [round(statistics.mean(total_reads_tcga_bam_df_list), 2),
                            round(statistics.stdev(total_reads_tcga_bam_df_list), 2),
                            round(statistics.median(total_reads_tcga_bam_df_list), 2)
                            ],
                            "tcga_bam_human_mapped":
                            [round(statistics.mean(total_reads_tcga_bam_human_mapped_df_list), 2),
                            round(statistics.stdev(total_reads_tcga_bam_human_mapped_df_list), 2),
                            round(statistics.median(total_reads_tcga_bam_human_mapped_df_list), 2)
                            ],
                            "human_unmapped_fastq":
                            [round(statistics.mean(total_reads_human_unmapped_fastq_df_list), 2),
                            round(statistics.stdev(total_reads_human_unmapped_fastq_df_list), 2),
                            round(statistics.median(total_reads_human_unmapped_fastq_df_list), 2)
                            ],
                            "srnanalyzer_processed_fa":
                            [round(statistics.mean(total_reads_srnanalyzer_processed_fa_df_list), 2),
                            round(statistics.stdev(total_reads_srnanalyzer_processed_fa_df_list), 2),
                            round(statistics.median(total_reads_srnanalyzer_processed_fa_df_list), 2)
                            ],
                            "Phylum":
                            [round(statistics.mean(phylum_total_counts_list), 2),
                            round(statistics.stdev(phylum_total_counts_list), 2),
                            round(statistics.median(phylum_total_counts_list), 2)
                            ],
                            "Class":
                            [round(statistics.mean(class_total_counts_list), 2),
                            round(statistics.stdev(class_total_counts_list), 2),
                            round(statistics.median(class_total_counts_list), 2)
                            ],
                            "Order":
                            [round(statistics.mean(order_total_counts_list), 2),
                            round(statistics.stdev(order_total_counts_list), 2),
                            round(statistics.median(order_total_counts_list), 2)
                            ],
                            "Family":
                            [round(statistics.mean(family_total_counts_list), 2),
                            round(statistics.stdev(family_total_counts_list), 2),
                            round(statistics.median(family_total_counts_list), 2)
                            ],
                            "Genus":
                            [round(statistics.mean(genus_total_counts_list), 2),
                            round(statistics.stdev(genus_total_counts_list), 2),
                            round(statistics.median(genus_total_counts_list), 2)
                            ]}
print("statistical_summary_dict:\n", statistical_summary_dict)

statistical_summary_df = pd.DataFrame(statistical_summary_dict)
print("statistical_summary_df.iloc[0:9]: ", statistical_summary_df.iloc[0:9])
print("statistical_summary_df.shape: ", statistical_summary_df.shape)
print("statistical_summary_df.columns: ", statistical_summary_df.columns)



processed_reads_info_df.to_csv(processed_reads_info_output_path, sep="\t", index=False)
statistical_summary_df.to_csv(statistical_summary_output_path, sep="\t", index=False)