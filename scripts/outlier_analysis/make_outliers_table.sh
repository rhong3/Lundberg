#!/bin/bash


gene_column_name="gene"
set="UCEC_subtype"

updown="up"
location_of_py_file="make_outliers_table.py"
location_of_data_file="/Users/rh2740/documents/Lundberg/Results/OLA/"$set"/ola_data.csv"
iqrs_over_median=1.5 #Note 1.5 IQRs is suggested, this is just for test data.

output_prefix="/Users/rh2740/documents/Lundberg/Results/OLA/"$set"/ola_table"
sample_names_file="/Users/rh2740/documents/Lundberg/Results/OLA/"$set"/samples.txt"
aggregate=True
write_frac_table="True"

python3 ${location_of_py_file} \
--input_df ${location_of_data_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_prefix ${output_prefix} \
--sample_names_file ${sample_names_file} \
--aggregate ${aggregate} \
--up_or_down ${updown} \
--write_frac_table ${write_frac_table}
