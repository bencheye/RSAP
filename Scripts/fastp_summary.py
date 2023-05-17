#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20211229
get fastp summary report
@author: Beance
"""
import os
import sys
import pandas as pd  

def get_html_summary(sample, text):
    out_dict = {'sample' : sample}
    out_dict['length'] = text.split("mean length after filtering:</td><td class='col2'>")[1].split(',')[0]    
    fitering_value = text.split("'after_filtering_summary'>")[1]    
    out_dict['reads_num'] = fitering_value.split("total reads:</td><td class='col2'>")[1].split('</td>')[0]
    out_dict['base_num'] = fitering_value.split("total bases:</td><td class='col2'>")[1].split('</td>')[0]
    out_dict['Q20'] = fitering_value.split("Q20 bases")[1].split(' M (')[1].split(')')[0]
    out_dict['Q30'] = fitering_value.split("Q30 bases")[1].split(' M (')[1].split(')')[0]
    out_dict['GC'] = fitering_value.split("GC content:</td><td class='col2'>")[1].split('</td>')[0]
    out_dict['read_pass_rate'] = text.split("reads passed filters")[1].split('(')[1].split(')')[0]
    # out_dict['low_quality'] = text.split("reads with low quality")[1].split('(')[1].split('%')[0]
    # out_dict['N_rate'] = text.split("reads with too many N")[1].split('(')[1].split('%')[0]
    # out_dict['short_read'] = text.split("reads too short")[1].split('(')[1].split('%')[0]
    return out_dict

meta_name = sys.argv[1]
fastp_dir = sys.argv[2]
output_name = sys.argv[3]

meta_data = pd.read_csv(meta_name, sep = '\t')
samples = meta_data['Sample'].tolist()

with open(output_name,"w") as f:
    f.write("Sample\tReadNum\tBaseNum\tLength\tQ20\tQ30\tGC\tReadPassRate\n")
    for sample in samples:
        sample_name = sample + '_fastp_report.html'
        input_file = os.path.join(fastp_dir, sample_name)
        with open(input_file, 'r') as text:
            text = text.read()
        out_dict = get_html_summary(sample, text)
        out_row = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            out_dict['sample'], out_dict['reads_num'],
            out_dict['base_num'], out_dict['length'], 
            out_dict['Q20'], out_dict['Q30'],
            out_dict['GC'], out_dict['read_pass_rate'],)
        f.write(out_row)




