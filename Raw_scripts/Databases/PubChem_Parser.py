#############################
# pubchem parser           #
############################

# files retrieved on linux by wget $url

##############################
# variables ################
########################
import os
import xml.etree.ElementTree as ET

folder = "F:/avans/stage MM/databases/pubchem/"
file = 'F:/avans/stage MM/databases/pubchem/Compound_000000001_000500000.xml'
outputE = folder + 'output_extended.txt'
outputD = folder + 'output_database.txt'



#items to parse in folder
files = os.listdir(folder)
print(files)

identifier = 0

Extended = open(outputE,'w')
pubchem_database = open(outputD,'w')

for file in files:
    file = folder + str(file)
    print (file)
    try:
        with open(file,'r') as f:
            for line in f:
                if "PC-Compound_charge" in line:
                    charge = (line.replace('<PC-Compound_charge>','').replace('</PC-Compound_charge>','').rstrip())
                    #print(charge)
                    line = next(f)
                    item_names = []
                    numeric_data = []
                    while '<PC-Compound_charge>' not in line:
                        if 'PC-InfoData_value_sval' in line:
                                line = line.replace('<PC-InfoData_value_sval>','').replace("</PC-InfoData_value_sval>",'')
                                #print(line.rstrip())
                                item_names.append(line.rstrip())


                        elif '<PC-InfoData_value_fval' in line:
                            # line = next(f)
                           # info_counter += 1
                            line = line.rstrip()
                            line = line.replace("<PC-InfoData_value_fval>",'').replace("</PC-InfoData_value_fval>",'')
                            #print(line)
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
                    identifier_pubchem = "pubchem"+(str(identifier))


                    # item to write for table
                    string = (str(identifier_pubchem) +'\t'+ str(UPAC_name) +'\t'+ str(Forumla) +'\t'+ str(Structure) +'\t'+ str(charge) +'\t'+ str(monoisotopic_mass) +'\n')
                    pubchem_database.write(string)

                    #extend table, need to figure out how to do ADDUCT and isoprevalence
                    string = (str(identifier_pubchem)+'\t'+str(Forumla)+'\t'+str(charge)+'\t'+str(monoisotopic_mass)+'\t'+'ADDUCT\tISOPREVALENCE \n')
                    Extended.write(string)
    except:
        print ( 'could not open' + file)
