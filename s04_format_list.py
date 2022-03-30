# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:50:31 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

working_directory = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\Final_results"
output_directory = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\charge_match"

#working_directory = input('Enter path to working directory')
#output_directory = input('Enter path to output directory')

def get_file_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

file_query = '.txt'
all_files = (get_file_names_with_strings([file_query]))

for b in all_files:
    print(b)
    df_path = working_directory + '\\' + b
    df = pd.read_csv(df_path, sep=",",skiprows=[0], names= ['m/z','resolution','charge','intensity','MS2','Scan','status','b/y error',
                                                                               'species','theoretical ion charge','theoretical b/y ion','Sample Name'])
    df = df.drop(df[df['charge']!=df['theoretical ion charge']].index)

    b_name = b[:-4]
    file = b + '_all_tissue'
    file_name = b_name + '.txt'
    file_path = output_directory + '\\' + file_name
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            df.to_csv(filec,index=False) 

    file_name = b_name + '.csv'
    file_path = output_directory + '\\' + file_name
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            df.to_csv(filec,index=False) 
            
    file = 'all_results_combined'
    file_name = file + '.txt'
    file_path = output_directory + '\\' + file_name
    with open(file_path,'a',newline='') as filec:
            writerc = csv.writer(filec)
            df.to_csv(filec,index=False)
            
    file = 'all_results_combined'
    file_name = file + '.csv'
    file_path = output_directory + '\\' + file_name
    with open(file_path,'a',newline='') as filec:
            writerc = csv.writer(filec)
            df.to_csv(filec,index=False)