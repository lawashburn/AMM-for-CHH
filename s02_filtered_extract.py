# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:36:23 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

spectra_import = input('Enter path to .txt extract file from step 1: ')
working_directory = input('Enter path to working directory: ')
data_type = input('Enter tissue type')
trial = input('Enter trial number')
error_marg = input('Enter MS1 ppm error cutoff: ')
intensity = input('Enter intensity cutoff: ')
scan_cut = input('Enter minimum scan # cutoff: ')

#spectra_import = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\Raw_Files\Formatted_MS2\PO_3_ms2_output_list.txt"#path to spectra after RawConverter
#working_directory = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\filtered_lists"
#data_type = 'PO_'
#trial = '3_'
#error_marg = 10 #+/- ppm

trial = str(trial)

#formats spectra import values
spectra_import = pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','empty'])
#spectra_import.columns = ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','empty']

spectra_import = spectra_import.apply(pd.to_numeric)

spectra_value = pd.DataFrame()
spectra_value['m/z'] = spectra_import['m/z']
spectra_value['resolution'] = spectra_import['resolution']
spectra_value['charge'] = spectra_import['charge']
spectra_value['intensity'] = spectra_import['intensity']
spectra_value['MS2'] = spectra_import['MS2']
spectra_value['Scan #'] = spectra_import['scan_number']
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<intensity].index) #remove intensity values less than 100
spectra_value = spectra_value.drop(spectra_value[spectra_value['Scan #']<scan_cut].index) #remove intensity values less than 100
print(spectra_value)

file_name = data_type + trial + '_filtered_list.csv'
file_path = working_directory + '\\' + file_name
with open(file_path,'w',newline='') as filec:
    writerc = csv.writer(filec)
    spectra_value.to_csv(filec,index=False)
    
file_name = data_type + trial + '_filtered_list.txt'
file_path = working_directory + '\\' + file_name
with open(file_path,'w',newline='') as filec:
    writerc = csv.writer(filec)
    spectra_value.to_csv(filec,index=False)
    
print('analysis complete')