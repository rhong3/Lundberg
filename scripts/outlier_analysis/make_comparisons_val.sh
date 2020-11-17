#!/bin/bash

gene_column_name="Hugo_Symbol"
set="TCGA_RNA"

for m in 1-yr-survival 3-yr-survival 5-yr-survival 6-month-survival
do
	location_of_py_file="compare_groups_outliers.py"
	location_of_outliers_file="/Users/HONGRUNYU/documents/Segundo_Melanoma/Results/"$set"/OLA/"$m"/ola_table.txt"

	fdr_cut_off=0.05
	blue_or_red="red"
	output_qvals="True"
	frac_filter="0.3"
	
	group1_label="alive"
	group2_label="dead"

	group1_list="/Users/HONGRUNYU/documents/Segundo_Melanoma/Results/"$set"/OLA/"$m"/G1.txt"
	group2_list="/Users/HONGRUNYU/documents/Segundo_Melanoma/Results/"$set"/OLA/"$m"/G2.txt"

	output_prefix="/Users/HONGRUNYU/documents/Segundo_Melanoma/Results/"$set"/OLA/"$m"/"$group1_label

	python3  ${location_of_py_file} \
	--outliers_table  ${location_of_outliers_file} \
	--gene_column_name ${gene_column_name} \
	--fdr_cut_off ${fdr_cut_off} \
	--output_prefix ${output_prefix} \
	--group1_label ${group1_label} \
	--group1_list ${group1_list} \
	--group2_label ${group2_label} \
	--group2_list ${group2_list} \
	--blue_or_red ${blue_or_red} \
	--output_qvals ${output_qvals} \
	--frac_filter ${frac_filter}

	output_prefix="/Users/HONGRUNYU/documents/Segundo_Melanoma/Results/"$set"/OLA/"$m"/"$group2_label

	python3  ${location_of_py_file} \
	--outliers_table  ${location_of_outliers_file} \
	--gene_column_name ${gene_column_name} \
	--fdr_cut_off ${fdr_cut_off} \
	--output_prefix ${output_prefix} \
	--group1_label ${group2_label} \
	--group1_list ${group2_list} \
	--group2_label ${group1_label} \
	--group2_list ${group1_list} \
	--blue_or_red ${blue_or_red} \
	--output_qvals ${output_qvals} \
	--frac_filter ${frac_filter}
done
