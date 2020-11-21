###############################################################
# calculates the mayority vote from the venn diagram data csv##
###############################################################
import os
import pip
pip.main(['install','xlsxwriter'])
import xlsxwriter



path = 'F:/avans/stage MM/peakpicking/'
outputname = 'F:/avans/stage MM/peakpicking/Mayoiryt_vote.xlsx'



TMP_files = os.listdir(path)
files = []
for file in TMP_files:
    if 'venn_diagram_data.csv' in file:
        files.append(file)



file = outputname
wb = xlsxwriter.Workbook(file)
sheet = wb.add_worksheet()


row = -2

for file in files:
    row += 2
    column = 0
    sheet.write(row, column, file)
    file1 = path + file
    print('############################')
    print(file)
    with open(file1,'r') as f:
        #read header and create a dictionary#
        header = (f.readline()).rstrip()
        header = (header.split(','))[1:]
        Dict = {el:0 for el in header}
        print(Dict)
    with open(file1, 'r') as f:
        next(f)
        for line in f:
            TMP_line = line.rstrip()
            TMP_line = (TMP_line.split(','))[1:]
            print(TMP_line)
            TMP_value = 0
            for item in TMP_line:
                #print(item)
                TMP_value = TMP_value + int(item)

            print(TMP_value)
            # if there isnt an agreed item then the line will be dicsontinued
            if TMP_value  != 1:
                #print(TMP_line)
                counter = -1
                for number in TMP_line:
                    counter += 1
                    if int(number) != 0:
                        TMP_header = (header[counter])
                        Dict[TMP_header] += 1

    Dict = {k: v for k, v in sorted(Dict.items(), key=lambda item: item[1])}
    for x,y in Dict.items():
        column = 0
        row+= 1
        #print(x,y)
        sheet.write(row, column, x)
        column += 1
        sheet.write(row, column, y)

wb.close()
