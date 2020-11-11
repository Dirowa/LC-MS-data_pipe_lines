
########################################
# pubchem parser vor windows           #
#########################################


import os

os.system('pip install wget')

import wget
import gzip
import shutil
import subprocess


#############
# Variables #
#############

Download_location = '/home/thor/Donny/databases/databases/pubchem/ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/'



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



########## Download #################


print(items)
itmes_downloaded = os.listdir(Download_location)

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

                name = (f1 + ".txt")
                pubchem_database = open(name, 'w')


                header = 'identifier' + ',' + 'compoundname' +  ','  + 'baseformula' +  ','  + 'structure' +  ','  + 'charge' +  ','  + 'description' + "\n"
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
                                string = (str(identifier_pubchem) +  ','  + str(UPAC_name) +  ','  + str(Forumla) +  ','  + str(Structure) +  ','  + str(charge) +  ','  + str("NA") + "\n")
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

##########################
# parsed into on tsv file#
###########################




TMP_files = os.listdir(Download_location)
files = []
for file in TMP_files:
    if file[-3:]:
        files.append(file)

TMP_files = []

output = Download_location + output_name
Output = open(output,'w')

for file in files:
    file = Download_location + file
    with open (file, 'r') as f:
        for line in f:
            #line = line.replace('\t',',')
            #print(line)
            Output.write(line)
Output.close()

database_name = Download_location + 'pubchem.db'
conn = sqlite3.connect(database_name)

c = conn.cursor()

c.execute("CREATE TABLE Base ("
   "identifier varchar(16) NOT NULL,"
   "compoundname nvarchar(30) NOT NULL,"
   "baseformula nvarchar(20) NOT NULL,"
   "structure nvarchar(30) NOT NULL,"
   "charge int NOT NULL,"
   "description Text() NOT NULL,"
   "CONSTRAINT Base_pk PRIMARY KEY (identifier)"
");")


#data = open(output)
#data = csv.rea
#data.to_sql('Base', conn, if_exist='append', index = True)


#######################
# database creation   #
#######################

database_name = dabatase_location + 'pubchem.db'
conn = sqlite3.connect(database_name)
c = conn.cursor()
try:
    c.execute('''CREATE TABLE BASE                
                ([identifier] varchar(16) PRIMARY KEY NOT NULL, [compoundname]  varchar(30) NOT NULL, [baseformula]  varchar(20) NOT NULL, [structure] varchar(30) NOT NULL, [charge]  integer NOT NULL, [description]  Text NOT NULL)''')
except:
    print("database already excist")

with open(pubchem_tsv,'r') as file:
        next(file)
        for line in file:
            to_db = line.rstrip().replace(' ','').replace(' ','').split('\t')
            #print(to_db[0])
            Import = ("INSERT INTO BASE (identifier,compoundname,baseformula,structure,charge,description) "
                          "VALUES (" +"'"+ str(to_db[0]) +"'"','"'"+ str(to_db[1]) +"'"','"'"+ str(to_db[2]) +"'"','"'"+ str(to_db[3]) +"'"','"'"+ str(to_db[4]) +"'"','"'"+ str(to_db[5]) +"'"");")
            c.execute(Import)
