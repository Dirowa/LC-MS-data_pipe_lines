#########################################
# peaklist - database miner to xlsx file#
#########################################

# This script takes an feateure list txt file who is comma separated
# there feature names and mz values will be noted down.
# on the MZ is a range in PPM difference available in the setting to search in between
# all found databases inside the folder to databases will be used
# databases have been generated by the MetaboShiny Github package and own pubchemParser
# all data is saved in the given excel file


################
# imports ######
################
import sqlite3
import os
import pip
pip.main(['install','xlsxwriter'])
import xlsxwriter


##############################
# files & folders & Variables#
##############################

#extended_db = 'C:/Users/DVrin/OneDrive/Documenten/avans/stage_MM/databases/extended.db'
folder_to_databases = 'F:/avans/stage MM/databases/'

def MZ_query(PPM,folder_to_databases,file_to_feature_list):


    ##########################################
    # reading MZ_values and feature names in #
    ##########################################

    # opening lists
    mz_values = []

    # reading in data
    with open (file_to_feature_list,'r') as f:
        for line in f:
            line = line.split(',')
            # all feature names will be put in the possible result list to find out which fetures had no hit later on
            line = line[0:2]
            mz_values.append(line)


    mz_values = mz_values[1:]
    outputlocation_and_name = file_to_feature_list + '-PPM' + str(PPM) + '-' +'.xlsx'

    ################################################
    # check weither output file exsist if so rename#
    ################################################
    if os.path.exists(outputlocation_and_name):
        outputlocation_and_name = ((outputlocation_and_name.split('.'))[0]) +'(1).xlsx'
        TMP_counter = 2
        while os.path.exists(outputlocation_and_name):
            outputlocation_and_name = (outputlocation_and_name.split('('))[0] +'('+ str(TMP_counter) + ').xlsx'
            TMP_counter += 1

    #################
    # Database check#
    #################
    TMP_databases = os.listdir(folder_to_databases)
    databases = []
    for item in TMP_databases:
        if item[-3:] == '.db':
            if 'extended' not in item:
                databases.append(item)

    TMP_databases.clear()


    ####################################
    # querying in the extended database#
    ####################################

    #create te connection to extendend database
    extended_database = folder_to_databases + 'extended.db'




    workbook = xlsxwriter.Workbook(outputlocation_and_name)
    worksheet = workbook.add_worksheet()

    row = 0
    column = 0

    primal_header = ['Feature','inital_mz','Delta_mz','found_mz',  'fullformula', 'adduct', 'finalcharge', 'smiles', 'struct_id']
    query_items = ['compoundname', 'description', 'baseformula', 'identifier', 'charge']
    for item in primal_header:
        worksheet.write(row,column,item)
        column += 1

    for database in databases:
        #column += 1
        database = database[:-3]
        query_items =  (database + '.compound_name,' + database + '.description,' + database + '.baseformula,' + database + '.identifier,' + database + '.charge,').split(',')
        for item in query_items:
            column += 1
            worksheet.write(row, column, item)

    # database connection
    extended = sqlite3.connect(extended_database)
    c = extended.cursor()

    for item in mz_values:
        item_MZ = item[1]
        Feature = item[0]

        # the Query in the extended database
        query = 'SELECT extended.fullmz, extended.fullformula, extended.adduct, extended.finalcharge, structures.smiles, extended.struct_id FROM extended ' \
                'INNER JOIN structures ON extended.struct_id = structures.struct_id ' \
                'WHERE extended.fullmz >= '+ str(float(item[1]) - (PPM/2)) +' AND extended.fullmz <= '+ str(float(item[1]) + (PPM/2)) +';'

        # retrieve result in for loop becouse then it retrieve multiple results
        for hit in c.execute(query):
            #print(item)
            column = 3
            row += 1
            worksheet.write(row, 0, item[0])
            worksheet.write(row, 1, item[1])
            worksheet.write(row, 2, max(float(item_MZ), float(hit[0])) - (min(float(item_MZ), float(hit[0]))))
            smiles = hit[-2]

            for item_a in hit:
                worksheet.write(row,column,item_a)
                column += 1

            #column += 1
            column = 9
            counter = 0




            for database in databases:
                column += 1
                counter += 1
                database = folder_to_databases + database
                database_a = sqlite3.connect(database)
                C = database_a.cursor()
                query = ('SELECT compoundname, description, baseformula, identifier, charge FROM base '
                         ' WHERE structure = "' + smiles + '" '
                                                         ' LIMIT 10;')

                C.execute(query)
                data =  C.fetchall()

                if len(data) == 0:
                    for i in range(5):
                        worksheet.write(row, column, 'NA')
                        column +=1

                elif (len(data)) != 1:
                    compound = ''
                    description = ''
                    baseformula = ''
                    identifier = ''
                    charge = ''
                    for item_b in data:

                        compound = compound + '---' + str(item_b[0])
                        description = description + '---' + str(item_b[1])
                        baseformula = baseformula + '---' + str(item_b[2])
                        identifier = identifier + '---' + str(item_b[3])
                        charge = charge + '---' + str(item_b[4])

                    data = [compound[3:],description[3:],baseformula[3:],identifier[3:],charge[3:]]

                    for query_items in data:

                        if column / 6 == (1.5 + counter):
                         column +=1
                        worksheet.write(row, column, query_items)
                        column += 1

                else:
                    for query_items in data:
                        if column / 6 == (1.5 + counter):
                         column +=1

                        for item_c in query_items:
                            if item_c == None:
                                worksheet.write(row, column, 'NA')

                            worksheet.write(row, column, item_c)
                            column += 1

    workbook.close()

# number howmuch the MZ value is allowed to differ from the results found in other databases
PPMS = [0.0005000]

#output location file string isnt allowerd to contain (
#file containing mz values, BE SURE THE first column is feature ID and second the MZvalue! DONT REMOVE HEADER
path = 'F:/avans/stage MM/peakpicking/'
TMP_files = os.listdir(path)
files = []

for file in TMP_files:
    if '.txt' in file:
        files.append(file)

for ppm in PPMS:
    print(ppm)
    for file in files:
        print(file)
        file = path + file
        MZ_query(ppm,folder_to_databases,file)
