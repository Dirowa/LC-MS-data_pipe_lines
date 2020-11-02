################
# imports ######
import sqlite3
import os

##################
# files & folders#
##################

extended_db = 'extended database location'
folder_to_databases = 'databases locaton'

#file containing mz values, BE SURE THE first column is feature ID and second the MZvalue! DONT REMOVE HEADER
file_to_feature_list ='F:/avans/stage MM/XCMS/output/Feature_list_centwave_N5000.txt'

# nmumber
PPM = 0.0000100

#########################
# editing imported files#
#########################

################ retrieving database files ###############
#retrieving all files
databases1 = os.listdir(folder_to_databases)
#deleting extended file
databases1.remove((extended_db.split('/')[-1]))

databases = []
for file in databases1:
    if '.db' in file:
        databases.append(file)
databases1.clear()

##################### editing imported feature file ########

mz_values = []
with open (file_to_feature_list,'r') as f:
    for line in f:
        line = line.split(',')
        #print(line)
        line = line[0:2]
        mz_values.append(line)


mz_values = mz_values[1:]
print(mz_values)


###############################
# quering in extended database#
###############################


# database connection
extended = sqlite3.connect('F:/avans/stage MM/databases/extended.db')
c = extended.cursor()

# quering for items found

for item in mz_values:
    item_MZ = item[1]
    itemA = float(item_MZ) - (PPM/2)
    itemB = float(item_MZ) + (PPM/2)

    #'SELECT * FROM extended WHERE fullMZ >= 57.95000000 AND fullMZ <= 58.05 LIMIT 10;'
    query = 'SELECT * FROM extended WHERE fullMZ >= '+ str(itemA) +' AND fullMZ <= '+ str(itemB) +';'
    for row in c.execute(query):
        print(item[0], row)
    print(item[0])
