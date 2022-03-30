
import csv
import pandas as pd
from datetime import datetime
now = datetime.now()


working_directory = input('Enter path to working directory: ')
MS2_path = input('Enter path to RawCoverter .MS2 output: ')
tissue_type = input('Enter tissue type: ')
trial_no = input('Enter trial number: ')

#working_directory = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220318\Raw_SIM_data\MS2_format"
#MS2_path = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220318\Raw_SIM_data\RawConverter_Output\2022_0318_Brain_SIM_TR1.ms2"
#tissue_type = 'brain_'
#trial_no = 1


trial_no = str(trial_no)
print('#1: ',datetime.now())
with open(MS2_path) as input:
    lst = [line.strip() for line in input]

new_list= []
final_lst = []
final_lst.append(['m/z', 'resolution', 'charge', 'intensity', 'MS2','scan_number'])
ms2_list = []

new = lst
#print(new)
print('#2: ',datetime.now())
for i in new:
    new_list.append(i.split())
    if '@' in i:
        x = i.split()
        for y in x:
            if '@' in y:
                ms2 = y[0:y.index('@')]
                ms2_list.append(str(ms2))
print('#3: ',datetime.now())
header_list = new_list[0:25]
new_list = new_list[26:] # starts from line 26 to remove the first few header lines so that program could proceed
seperation_list = []
scan_number_list = []    
print('#4: ',datetime.now())
for i in header_list:
    if 'S' in i:
        scan_number_list.append(i[1])
print('#5: ',datetime.now())
for i in range(len(new_list)):
    #print(i, new_list[i])
    if 'RetTime' in new_list[i]:
        seperation_list.append(i-1)
    if 'PrecursorInt' in new_list[i]:
        seperation_list.append(i+2)
    if 'S' in new_list[i]:
        scan_number_list.append(new_list[i][1])
print('#6: ',datetime.now())
seperation_pairs = []
start = 0
for i in range(int(len(seperation_list)/2)):
    seperation_pairs.append((seperation_list[i+start],seperation_list[i+start+1]))
    start +=1 
print('#7: ',datetime.now())    
update_index = 0
for start,end in seperation_pairs:
    start += update_index
    end += update_index
    new_list[start:end] = '-'
    update_index -= (end-start-1)  

ms2_list_index = 0
scan_number_index = 0
for element in new_list:
    if element == '-':
        ms2_list_index+=1
        scan_number_index+=1
        continue   
    element.append(ms2_list[ms2_list_index])
    element.append(scan_number_list[scan_number_index])
    final_lst.append(element)
print('#8: ',datetime.now())
out_path = working_directory + '\\ms2_output_list.txt'
with open(out_path,'w') as output:
    for i in final_lst:
        for j in i:
            output.write(str(j + ' '))
        output.write('\n')
        
out_suf = '_ms2_output_list.txt'
out_name = working_directory + '\\' + tissue_type + trial_no + out_suf
with open(out_name,'w') as output:
    for i in final_lst:
        for j in i:
            output.write(str(j + ' '))
        output.write('\n')
csv_suf = '_untarget_ms2_output_list.csv'
csv_name =  working_directory + '\\' + tissue_type + trial_no + out_suf
read_txt = pd.read_csv(out_name)
read_txt.to_csv(csv_name, index=None)
print('#9: ',datetime.now())
read_csv_path = csv_name
SIM_result_imp = pd.read_csv(read_csv_path, sep = ' ')
SIM_result = pd.DataFrame()
SIM_result['m/z'] = SIM_result_imp['m/z']
SIM_result['resolution'] = SIM_result_imp['resolution']
SIM_result['charge'] = SIM_result_imp['charge']
SIM_result['intensity'] = SIM_result_imp['intensity']
SIM_result['MS2'] = SIM_result_imp['MS2']
SIM_result['Scan #'] = SIM_result_imp['scan_number']
print(SIM_result)

print('analysis complete')