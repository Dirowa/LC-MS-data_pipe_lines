
#############################
# pubchem parser           #
############################


############## REQUIREMENTS ################

# around 350 gb of free space
# linux based program for the Wget command


# files retrieved on linux by wget $url
import os


########################
# variables to change  #
########################

location_download_pubchem = 'location where pubchem database will be downloaded'
folder = "folder who contains the XML files"


#the database tables needed
outputE = folder + 'output_extended.txt'
outputD = folder + 'output_database.txt'


###########################
#download pubchem database#
###########################

####### STILL NEEDS WORK, TO SAFE  AT CERTAIN LOCATION!
#os.system("cd" + location_download_pubchem)
# os.system( wget -m ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/)


    #items to parse in folder
files = os.listdir(folder)
print(files)

identifier = 0

Extended = open(outputE,'w')
pubchem_database = open(outputD,'w')


## Open all files
for f in files:
    #print(f)
    
    # find file extension
    ext = f[-3:]
    #print(ext)
    
    
    #we dont need md5 file so lets remove it
    if "md5" == ext:
        print(" file" + f + ' will be deleted')
        f1 = (folder + f)
        try:
            os.remove(f1)
        except:
            print(' unable to remove file ' + f1)
            
    #if it is a zipped file
    elif '.gz' == ext:

            #to find he file who needs to be unzipped
            file = folder + f
            #to find the file that will be extracted
            ext = f[:-3]
            #create adres of the extracted file
            f1 = (folder + ext)
            
            
            # all code has been put into an try loop incase an error occours at the locationg of parsing ( this happens once per file)
            
            #unzipping the file
            try:
                print("unzipping file " + f)
                os.system('gunzip -k ' + file)
                print('finisht to unzip')
                
                
                #open the created file
                try:
                    with open(f1, 'r') as f:
                                #read eveyr line in
                                for line in f:
                                    try:
                                        # parse the required items per section in
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
                                            try:
                                                # item to write for table
                                                string = (str(identifier_pubchem) +'\t'+ str(UPAC_name) +'\t'+ str(Forumla) +'\t'+ str(Structure) +'\t'+ str(charge) +'\t'+ str(monoisotopic_mass) +'\n')
                                                pubchem_database.write(string)
                                                #print(string)
                                                #extend table, need to figure out how to do ADDUCT and isoprevalence
                                                string = (str(identifier_pubchem)+'\t'+str(Forumla)+'\t'+str(charge)+'\t'+str(monoisotopic_mass)+'\t'+'ADDUCT\tISOPREVALENCE \n')
                                                Extended.write(string)
                                                #print(string)
                                            except:
                                                print(' unable to parse line of xml file')
                                    except:
                                        print(' something went wrong with a chunk of the xml file')
                    print('parsed file ' + f1 + '\n')
                    try:
                        os.remove(f1)
                        print('removed file ' + f1)
                    except:
                        print ("unable to delete file" + f1)
                except:
                    print("something went wrong with handeling the parser")
                    try:
                        os.remove(f1)
                        print('removed file ' + f1)
                    except:
                        print ("unable to delete file" + f1)
            except:
                print(' an error occoured')


Extended.close()
pubchem_database.close()

