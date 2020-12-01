
########################################
# pubchem parser vor windows           #
#########################################


import os

os.system('pip install wget')
import sqlite3
import wget
import gzip
import shutil
import subprocess


#############
# Variables #
#############

Download_location = '/home/thor/Donny/databases/databases/pubchem/'
folder_to_extended_database = '/home/thor/Donny/databases/databases/extended.db'


##############################
# retrieving list to download#
##############################
items = []
identifier = 0
proc = subprocess.Popen("curl ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/", stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
out = str(out).split(' ')


def bar_custom(current, total, width=80):
    print("Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total))

for item in out:
    if 'Compound' in item:
        item = ((item.replace("\n-r--r--r--",'')))
        if item[-3:] != 'md5':
            items.append(item)

###############################################
# Downloading & unzipping & reading & deleting#
###############################################
# must be able to run parralel. this part takes arund 3 days with downloading. unzippin and exstracting information
# per run it takes maximum of 20 gb storage room and a fill will be created
# it required an already build extended db by metaboshiny


########## Download #################


print(items)
itmes_downloaded = os.listdir(Download_location)

# retrieve lasrgest identivier in extended database

extended = sqlite3.connect(folder_to_extended_database)
c = extended.cursor()
for hit in (c.execute("select MAX(struct_id) from extended  limit 1 ")):
    identifier = int(((str(hit)).split(',')[0]).replace('(',''))
print(identifier)

for item in items:

    if item not in itmes_downloaded:
        #check wether the item is already downloaded
        url = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/" + item
        location = Download_location + item
        print ("start downloading " + url + " will be send to " + Download_location)
        wget.download(url, location, bar=bar_custom)



    ########## Unzipping #############

    file = Download_location + item
    try:
        print("unzipping file " + item)

        with gzip.open(file, 'rb') as f_in:
            with open(file[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        print('finisht to unzip')


        # READING UNZIPPED FILE#
        f1 = file[:-3]
        os.remove(file)
        try:
            with open(f1, 'r') as f:

                #open extended table tsv
                name = (f1 + "Extended.txt")
                pubchem_database = open(name, 'w')

                #open structure ID base
                name1 = (f1 + "struct_ID.txt")
                pubchem_structID = open(name1, 'w')


                print('start reading ' + f1)
                # read eveyr line in
                for line in f:
                    try:
                        # parse the required items per section in
                        if "PC-Compound_charge" in line:
                            charge = (line.replace('<PC-Compound_charge>', '').replace('</PC-Compound_charge>',
                                                                                       '').rstrip())
                            # print(charge)
                            line = next(f)
                            item_names = []
                            numeric_data = []
                            while '<PC-Compound_charge>' not in line:
                                if 'PC-InfoData_value_sval' in line:
                                    line = line.replace('<PC-InfoData_value_sval>', '').replace(
                                        "</PC-InfoData_value_sval>", '')
                                    # print(line.rstrip())
                                    item_names.append(line.rstrip())


                                elif '<PC-InfoData_value_fval' in line:
                                    # line = next(f)
                                    # info_counter += 1
                                    line = line.rstrip()
                                    line = line.replace("<PC-InfoData_value_fval>", '').replace(
                                        "</PC-InfoData_value_fval>", '')
                                    # print(line)
                                    numeric_data.append(line)

                                line = next(f)

                            #####################
                            # editing found data#
                            #####################

                            UPAC_name = item_names[0]
                            Forumla = item_names[-3]
                            Structure = item_names[-1]
                            monoisotopic_mass = numeric_data[-1]
                            identifier += 1

                            try:
                                # item to write for table
                                string = (str(identifier) +  '\t'  + str(Forumla) +  '\t'  + str(charge) +  '\t'  + str(monoisotopic_mass) +  '\t'  + str("PubChem_No_Aduct") +  '\t'  + str("NA") + "\n")
                                string1 =(str(identifier) + '\t'+ str(Structure) + '\n')
                                pubchem_database.write(string)
                                pubchem_structID.write(string1)
                                # print(string)
                                string = ''
                                string1 = ''
                            except:
                                print(' unable to parse line of xml file')
                    except:
                        print(' something went wrong with a chunk of the xml file')
            pubchem_database.close()
            pubchem_structID.close()
            try:
                os.remove(f1)
                print('removed file ' + f1)
            except:
                print("unable to delete file" + f1)
        except:
            print("something went wrong with handeling the parser")
            try:
                os.remove(f1)
                print('removed file ' + f1)
            except:
                print("unable to delete file" + f1)
    except:
        print(' an error occoured')

##########################
# parsed into on tsv file#
###########################

TMP_files = os.listdir(Download_location)
types = ["Extended.txt","struct_ID.txt"]

for type in types:
    files = []
    for file in TMP_files:
        if type in file:
            files.append(file)

    TMP_files = []

    output = Download_location + type[:-4] + '.tsv'
    Output = open(output,'w')

    for file in files:
        file = Download_location + file
        with open (file, 'r') as f:
            for line in f:
                Output.write(line)
        os.remove(file)
    Output.close()

