# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:23:20 2022

@author: lawashburn
"""
import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

spectra_import = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\Raw_Files\Formatted_MS2\PO_3_untarget_ms2_output_list.csv"#path to spectra after RawConverter
ion_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\ion_list.csv"
precursor_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\precursor_list.csv"
working_directory = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\num_test"
final_dir =r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\Final_results"
data_type = 'PO_'
trial = '3_'
sample_name = 'PO 3'
error_marg = 10 #+/- ppm
h_mass = 1.00784

#spectra_import = input('Enter path to formatted spectra .txt file: ')
#ion_list_import = input('Enter path to ion fragment list .csv: ')
#precursor_list_import = input('Enter path to precursor mass .csv: ')
#working_directory = input('Enter path to working directory: ')
#final_dir = input('Enter path to output directory: ')
#data_type = input('Enter tissue type: ')
#trial = input('Enter trial number: ')
#sample_name = input('Enter sample name (e.g. TG2')
#error_marg = input('Enter ppm error cutoff: ')


print('loading files', datetime.now())
#formats spectra import values
spectra_import = pd.read_csv(spectra_import, sep=",",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','empty'])
spectra_value = pd.DataFrame()
spectra_value['m/z'] = spectra_import['m/z']
spectra_value['resolution'] = spectra_import['resolution']
spectra_value['charge'] = spectra_import['charge']
spectra_value['intensity'] = spectra_import['intensity']
spectra_value['MS2'] = spectra_import['MS2']
spectra_value['Scan #'] = spectra_import['scan_number']
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<100].index) #remove intensity values less than 100
spectra_value = spectra_value.drop(spectra_value[spectra_value['Scan #']<40000].index) #remove intensity values less than 100
spectra_MS2 = spectra_value['MS2'].values.tolist()
print('sorting files', datetime.now())
precursor_list = pd.read_csv(precursor_list_import)
precursor_values = precursor_list['m/z'].values.tolist()

valid_MS2 = []
for a in precursor_values:
    for b in spectra_MS2:
        ppm_err = ((abs(a-b))/a) * 1E6
        if ppm_err<10:
            valid_MS2.append(b)

valid_MS2_short = []
for c in valid_MS2:
    if c not in valid_MS2_short:
        valid_MS2_short.append(c)

spectra_value['status'] = spectra_value['MS2'].isin(valid_MS2_short)
spectra_value = spectra_value.drop(spectra_value[spectra_value['status']==False].index)

ion_list = pd.read_csv(ion_list_import)
species_name = ion_list['Species'].values.tolist()

species_search = []
for d in species_name:
    if d not in species_search:
        species_search.append(d)
  
for e in species_search:
    print('searching species:',c,datetime.now())
    species  = []
    theo_precursor = []
    theo_ion = []
    act_precursor = []
    act_ion = []
    scan_no = []
    
    target_ions = ion_list.drop(ion_list[ion_list['Species']!=e].index)
    ion_charge1 = target_ions['1'].values.tolist()
    ion_charge2 = target_ions['2'].values.tolist()
    ion_charge3 = target_ions['3'].values.tolist()
    ion_charge4 = target_ions['4'].values.tolist()
    ion_charge5 = target_ions['5'].values.tolist()
    
    ion_charge1 = [x for x in ion_charge1 if np.isnan(x) == False]
    ion_charge2 = [x for x in ion_charge2 if np.isnan(x) == False]
    ion_charge3 = [x for x in ion_charge3 if np.isnan(x) == False]
    ion_charge4 = [x for x in ion_charge4 if np.isnan(x) == False]
    ion_charge5 = [x for x in ion_charge5 if np.isnan(x) == False]
    
    if len(ion_charge1)>=1:
        print('searching +1 ions:',c,datetime.now())
        for f in ion_charge1:
            charge = 1
            error_1 = 0.02
            spectra_c = spectra_value
            spectra_c['b/y error'] = abs(((f*charge)-(h_mass*charge)) - ((spectra_c['m/z']*charge)-(h_mass*charge)))
            spectra_d = spectra_c
            spectra_d = spectra_d.drop(spectra_d[spectra_d['b/y error']>=error_1].index) #remove charges equal to 0

            spectra_d['species'] = e
            spectra_d['theoretical ion charge'] = charge
            spectra_d['theoretical b/y ion'] = f
            spectra_d['Sample Name'] = sample_name
            if len(spectra_d)>0:
                f = str(f)
                file_name = data_type + trial + e + "_" + f + '_charge1_matches.csv'
                file_path = working_directory + '\\' + file_name
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        spectra_d.to_csv(filec,index=False)

            else:
                pass
        else:
                pass
    
    else:
            pass
        
    if len(ion_charge2)>=1:
        print('searching +2 ions:',c,datetime.now())
        for g in ion_charge2:

            charge = 2
            error_2 = 0.01
            spectra_e = spectra_value
            spectra_e['b/y error'] = abs(((g*charge)-(h_mass*charge)) - ((spectra_e['m/z']*charge)-(h_mass*charge)))
            spectra_f = spectra_e
            spectra_f = spectra_f.drop(spectra_f[spectra_f['b/y error']>error_2].index) #remove charges equal to 0
            spectra_f['species'] = c
            spectra_f['theoretical ion charge'] = charge
            spectra_f['theoretical b/y ion'] = g
            spectra_f['Sample Name'] = sample_name
            if len(spectra_f)>0:
                g = str(g)
                file_name = data_type + trial + e + "_" + g + '_charge2_matches.csv'
                file_path = working_directory + '\\' + file_name
                with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    spectra_f.to_csv(filec,index=False)
            else:
                pass           
            
    else:
        pass
    
    if len(ion_charge3)>=1:
        print('searching +3 ions:',c,datetime.now())
        for h in ion_charge3:
            charge = 3
            error_3 = 0.0067
            spectra_g = spectra_value
            spectra_g['b/y error'] = abs(((h*charge)-(h_mass*charge)) - ((spectra_g['m/z']*charge)-(h_mass*charge)))
            spectra_h = spectra_g
            spectra_h = spectra_h.drop(spectra_h[spectra_h['b/y error']>error_3].index) #remove charges equal to 0
            spectra_h['species'] = e
            spectra_h['theoretical ion charge'] = charge
            spectra_h['theoretical b/y ion'] = h
            spectra_h['Sample Name'] = sample_name
            if len(spectra_h)>0:
                h = str(h)
                file_name = data_type + trial + e + "_" + h + '_charge3_matches.csv'
                file_path = working_directory + '\\' + file_name
                with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    spectra_h.to_csv(filec,index=False)
            else:
                pass

    else:
        pass
    if len(ion_charge4)>=1:
        print('searching +4 ions:',c,datetime.now())
        for i in ion_charge4:

            charge = 4
            error_4 = 0.005
            spectra_i = spectra_value
            spectra_i['b/y error'] = abs(((i*charge)-(h_mass*charge)) - ((spectra_i['m/z']*charge)-(h_mass*charge)))
            spectra_j = spectra_i
            spectra_j = spectra_j.drop(spectra_j[spectra_j['b/y error']>error_4].index) #remove charges equal to 0
            spectra_j['species'] = e
            spectra_j['theoretical ion charge'] = charge
            spectra_j['theoretical b/y ion'] = i
            spectra_j['Sample Name'] = sample_name
            if len(spectra_j)>0:
                i = str(i)
                file_name = data_type + trial + e + "_" + i + '_charge4_matches.csv'
                file_path = working_directory + '\\' + file_name
                with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    spectra_j.to_csv(filec,index=False) 
            else:
                pass
    else:
        pass
    
    if len(ion_charge5)>=1:
        print('searching +5 ions:',c,datetime.now())
        for k in ion_charge5:

            charge = 5
            error_5 = 0.005
            spectra_k = spectra_value
            spectra_k['b/y error'] = abs(((k*charge)-(h_mass*charge)) - ((spectra_k['m/z']*charge)-(h_mass*charge)))
            spectra_l = spectra_k
            spectra_l = spectra_l.drop(spectra_l[spectra_l['b/y error']>error_5].index) #remove charges equal to 0
            spectra_l['species'] = e
            spectra_l['theoretical ion charge'] = charge
            spectra_l['theoretical b/y ion'] = k
            spectra_l['Sample Name'] = sample_name
            if len(spectra_l)>0:
                k = str(k)
                file_name = data_type + trial + e + "_" + k + '_charge5_matches.csv'
                file_path = working_directory + '\\' + file_name
                with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    spectra_l.to_csv(filec,index=False) 
            else:
                pass
    else:
        pass
    
def get_file_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

file_query = ''
all_files = (get_file_names_with_strings([file_query]))

all_files_path = []

for z in all_files:
    path = working_directory + '\\' + z
    all_files_path.append(path)

final_df = pd.concat((pd.read_csv(w) for w in all_files_path))

file_name = data_type + trial + '_results_matches.csv'
file_path = final_dir + '\\' + file_name
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        final_df.to_csv(filec,index=False) 

file_name = data_type + trial + 'results_matches.txt'
file_path = final_dir + '\\' + file_name
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        final_df.to_csv(filec,index=False) 