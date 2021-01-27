#########################################################################
# peakpicking comparisator to create ven diagrams                       #
# input files are the output files from the query databases xlsx file
# files should be run through Venn Diagram maker to get some easy graph#
#########################################################################

import os
import pip
pip.main(['install','xlrd'])
import xlrd
import time

xlsx_files_path = "F:/avans/stage MM/Annotated_outpet_sherlOK/"
output = "F:/avans/stage MM/Annotated_outpet_sherlOK/venn_diagram_data.csv"
#1,2,5,7
row_number = 7

Round_down_numbers = False
deci = 2
xlsx_files = []
files_used = []

for file in os.listdir(xlsx_files_path):
    if file[-4:] == 'xlsx':
        xlsx_files.append(file)


########################################
# generating all unique MZ values found#
########################################

values = {}
mz_values_unique = []
for file in xlsx_files:
    files_used.append(file)
    #opening excel files
    file = xlsx_files_path + str(file)
    wb = xlrd.open_workbook(file)
    sheet = wb.sheet_by_index(0)
    #sheet.cell_value(0, 0)
    # reading excel files
    for i in range(sheet.nrows):
        mz_value = (str(sheet.cell(i,row_number))).replace('text:','').replace('number:','').replace("'",'')

        try:
            if Round_down_numbers == True:
                mz_value1 =str(mz_value)
                mz_value1 = float(mz_value1)
                mz_value = round(mz_value1,deci)
        except:
            print('whoops1')
        mz_value = str(mz_value)
        # reading in new found mz_values
        if mz_value not in mz_values_unique:
            mz_values_unique.append(mz_value)
print(mz_values_unique)
values = dict.fromkeys(mz_values_unique, '')
values['algorithm'] = ''
##################################
# generating the venn diagram####
#################################
counter = 0
for file in xlsx_files:
    counter += 1

    TMP_mz = []
    values['algorithm'] = str(values['algorithm']) + file + ','
    file = xlsx_files_path + str(file)
    wb = xlrd.open_workbook(file)
    sheet = wb.sheet_by_index(0)
    # sheet.cell_value(0, 0)
    # reading excel files

    for i in range(sheet.nrows):
        mz_value = (str(sheet.cell(i,row_number))).replace('text:','').replace('number:','').replace("'",'')
        if mz_value not in TMP_mz:
            TMP_mz.append(mz_value)

    for item in TMP_mz:
            # print(item)
            #print(values[item])
            try:
                if Round_down_numbers == True:
                    item = str(item)
                    item = float(item)
                    item = round(item, deci)
                    item = str(item)
                    print(item)
            except:
                print('whoops2')

            if (len(values[item])) < counter:
                values[item] = str(values[item]) + (str(1))

    for key, value in values.items():
        if len(str(value)) != counter:
            if key != 'algorithm':
                values[key] = str(value) + (str(0))

output = open(output,'w')

#removing first dic
del values[(list(values.keys())[0])]

#print(values)
for key,value in values.items():
    if key == 'algorithm':
        value = value[:-1]
        string = key + ',' + value + '\n'
        output.write(string)

for key,value in values.items():
    if key != 'algorithm':
        value = ','.join(value[i:i + 1] for i in range(0, len(value), 1))

        string = key + ',' + value + '\n'
        output.write(string)

output.close()
