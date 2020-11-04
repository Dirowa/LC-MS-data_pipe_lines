########################################
# pubchem parser vor windows           #
#########################################


import os

os.system('pip install wget')

import wget
import gzip
import shutil


#############
# Variables #
#############

Download_location = 'F:/avans/stage MM/databases/pubchem1/'



import subprocess
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
        item = ((item.split('\\'))[0])
        if item[-3:] != 'md5':
            items.append(item)

###############################################
# Downloading & unzipping & reading & deleting#
###############################################




########## Download #################


print(items)
for item in items:
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

                name = (f1 + ".txt")
                pubchem_database = open(name, 'w')


                header = 'identifier' + '\t' + 'compoundname' + '\t' + 'baseformula' + '\t' + 'structure' + '\t' + 'charge' + '\t' + 'description' + "\n"
                pubchem_database.write(header)


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
                            identifier_pubchem = "pubchem" + (str(identifier))
                            try:
                                # item to write for table
                                string = (str(identifier_pubchem) + '\t' + str(UPAC_name) + '\t' + str(Forumla) + '\t' + str(Structure) + '\t' + str(charge) + '\t' + str("NA") + "\n")
                                pubchem_database.write(string)
                                # print(string)

                            except:
                                print(' unable to parse line of xml file')
                    except:
                        print(' something went wrong with a chunk of the xml file')
            pubchem_database.close()

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



