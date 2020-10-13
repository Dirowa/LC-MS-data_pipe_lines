
import os
#########################################################
# this script reads in all peaklists in a current folder#
# with this an global data analyse can be made about how#
# well the peakpicking went in comparisation with other #
# the output data has to be runned by an R sript to     #
# generate an venn diagram                              #
#########################################################



###################
# paths ##########
#################

# folder containing all data form xcms, all data must be comma separated
# first item is feature name, second is median m/z value
folderpeakpicking = "F:/avans/stage MM/peakpicking/"

# where the preprocessed rounddown files should be
folder_preprocessed ="F:/avans/stage MM/peakpicking/pre-processed/"

#output folder of data
outputfolder = "F:/avans/stage MM/peakpicking_output/"

#folder to anouther output folder
R_input ='F:/avans/stage MM/peakpicking_output/R_input/'

#name of outout data
venn_diagram_data_output_name = 'venn_diagram_data.csv'
compressed_output = 'compressed_output.csv'
mz_overview_output = 'mz_overview_output.csv'
############
# variables#
############

# number of N after comma should be roundown
N_rounddown = 2


mz_all1 = []
files1 = []
combinations = []



######################
# pre-processing     #
######################

# this is the rounding down of the m/z values to make it more compremental when algorithms are agreeíng

#find all files
files = os.listdir(folderpeakpicking)
#print(files)


#### open all files and round down the m/z value and remove duplicates

for file in files:
    #this is to prevent from opening an output file
    if file != 'output.txt' or file != 'output.csv':
        #refining which files should be openend
        files1.append(file)
        #creating a name and pathway
        file1 = folderpeakpicking + file
        feature_tmp = []
        mz_tmp = []

        try:
            with open(file1,'r') as f:
                next(f)
                for line in f:
                    line = line.split(',')
                    #retrieving feature name
                    feature_tmp.append(line[0])
                    try:
                        #retrieving m/z value and rounding down
                        line1 = round(float(line[1]),N_rounddown)
                        mz_tmp.append(line1)
                    except:
                        print('nut a number')

            #the 2 lists are created into an directionary
            TMP_dic = dict(zip(feature_tmp, mz_tmp))
            tmp_dic1 = dict()
            temp = []

            #removing duplicates
            for key, val in TMP_dic.items():
                if val not in tmp_dic1.values():
                    tmp_dic1[key] = val

            #print(str(tmp_dic1))

            file2 = folder_preprocessed + str(file)

            #writing down the preprocessed data into the preproceseed folder
            with open(file2,'w') as f:
                for key, value in tmp_dic1.items():
                    f.write(str(key) + ',' + str(value)+"\n")
        except:
            #this is so that he only takes the files and not a folder or something and goes into an error
            print('wasnt allowed to look into ' + str(file))

#######################
# retrieve informatin##
#######################

#searching for all pre-processed data
files1 = os.listdir(folder_preprocessed)
print(files)

mz_values =dict()

#############################################################
# generate mz/ value with each algorithm containing the same#
#############################################################

for file in files1:
    if file != 'output.txt' or file != 'output.csv':
        #files1.append(file)
        file = folder_preprocessed + file

        with open(file,'r') as f:
            file = file.split('/')[5]
            for line in f:
                line = (line.split(','))[1].rstrip()

                #print(line)
                if line not in mz_values.keys():
                    mz_values[line] = str(file)
                else:
                    fileA = str(mz_values[line])
                    fileA = str(fileA)+','+str(file)
                    mz_values[line] = fileA

file = outputfolder + mz_overview_output
with open(file,'w') as f:
    for x,y in mz_values.items():
        data = str(x) + " found in these files " + str(y) + '\n'
        print(data)
        f.write(data)
####################################################
# compressing the dictionary into a a count        #
####################################################

#print(mz_values)
compressedict = dict()
for key,value in mz_values.items():
    if value not in compressedict:
        compressedict[value] = 1
    else:
        compressedict[value] +=1



file = outputfolder + compressed_output
with open(file,'w') as f:
    for x,y in compressedict.items():
        data = ((str(x) + 'amount found = ' + str(y) + '\n'))
        print(data)
        f.write(data)

print('---------------------------------')
print('will be convering to a list for R')
print('---------------------------------')
###################################################
# generating an m/z value crossover for R         #
###################################################

# this parts generates a list contaning all found m/z values
all_mz = []

for file in files1:
    if file != 'output.txt' or file != 'output.csv':
        #files1.append(file)
        file = folder_preprocessed + file

        with open(file,'r') as f:
            for line in f:

                line = (line.split(','))[1].rstrip()
                #print(line)
                if line not in all_mz:
                    all_mz.append(line)
all_mz.sort()
#print(all_mz)

##################################################
# generate dictionary for all found m.z values####
##################################################

# the list with all found m/z will be converted to a dictionary
venn = {}
for i in all_mz:
    venn[str(i)] = []
#print(venn)

# fill in the dictionary by True False if the file contains the data
files = os.listdir(folder_preprocessed)
items = 0

for file in files:
    items += 1
    if file != 'output.txt' or file != 'output.csv':
        file1 = folder_preprocessed + file

        with open(file1,'r') as f:
            for line in f:

                #print(line)
                mz_value = (line.split(','))[1].rstrip()
                #print(mz_value)


                venn[mz_value].append(str(1))
            for mz_value,value in venn.items():
                #print(mz_value,value)
                count = (len(value))
                #print(count)
                if items != count:
                    venn[mz_value].append(0)


venn_name = outputfolder + venn_diagram_data_output_name
rowname = 'Mz_value'
#rowname = ''
for file in files:
    rowname = rowname + ',' + file


with open (venn_name,'w') as f:
    f.write(rowname + ' \n')
    for x,v in venn.items():
        string = v
        string = x,v
        string = str(string)
        string = string.replace('[','').replace(']','').replace("'",'').replace("(",'').replace(")",'')
        string = string + ' \n'

        f.write(string)


#########################################
#another import thing for R ven diagrams#
#########################################


files1 = os.listdir(folder_preprocessed)

for file in files1:
    file1 = folder_preprocessed + file
    file2 = R_input + file

    mz_values = []
    with open(file1, 'r') as f:
        with open(file2,'w') as f1:
            for line in f:
                line = (line.split(',')[1]).rstrip()
                f1.write(line + ',')


