import subprocess
import os
from os import path


current_location = os.getcwd()
current_location = current_location.replace('\\', '/')
major_settings = current_location + '/major_settings.txt'


paths = []
####### reading in global settings###
try:
    with(open(major_settings,'r')) as f:
        for line in f:
            paths.append((((line.split('='))[1]).rstrip()).replace(" ", ""))

except FileNotFoundError:
            print(' global_settings will be created')

def print_menu():  ## Your menu design here
    print(30 * "-", "MENU", 30 * "-")
    print("1.  edit Settings ")
    print("2.  Full Run")
    print("3.  Batch run ")
    print("4.  Run parts ")
    print("5.  Exit")
    print(67 * "-")
def print_menu_parts():  ## Your menu design here
    print(30 * "-", "MENU", 30 * "-")
    print("1.  Peakpicking ")
    print("2.  Batch Correction")
    print("3.  Noise filtering ")
    print("4.  Univariate testing ")
    print("5.  Multivariate testing ")
    print("6.  Feature Filtering ")
    print("7.  Database Annotation ")
    print("8.  Comparisation ")
    print("9.  Exit")
    print(67 * "-")
def print_menu_settings():  ## Your menu design here
    print(30 * "-", "EDIT SETTINGS", 30 * "-")
    print("1.  default_settings ")
    print("2.  Peakpicking ")
    print("3.  Batch Correction")
    print("4.  Noise filtering ")
    print("5.  Univariate testing ")
    print("6.  Multivariate testing ")
    print("7.  Feature Filtering ")
    print("8.  Database Annotation ")
    print("9.  Comparisation ")
    print("10.  Exit")
    print(67 * "-")
def iniatie_check(current_location):
    major_settings = current_location + '/major_settings.txt'
    if  not path.exists(major_settings):
        print('settings file is not found, will be created')
        setting = open(major_settings,'w+')

        print(30 * "-", "NOTICE", 30 * "-")
        print('You must change the folder location who contain the databases manually or in settings')
        print(30 * "-", "NOTICE", 30 * "-")


        ## list of files which should be in the settings file ##
        setting.write('Database_location = '+'PUT SOMETHING ELSE HERE'+ "\n" )
        tmp = current_location +'/scripts'
        setting.write('Scripts_location = '+  tmp + "\n")
        tmp = current_location + '/batch'
        setting.write('batch_location = '+  tmp + "\n")
        tmp = current_location + '/script_settings'
        setting.write('script_settings_location = ' + tmp + "\n")
        tmp = current_location + '/default_settings'
        setting.write('default_settings_location = ' + tmp + "\n")
        setting.write('Rscript.exe_path =  \n')

        setting.close()

    ########## generating Default settings #############
    major_settings = current_location + '/default_settings'
    if not path.exists(major_settings):
        os.mkdir(major_settings)

    ########## Peakpicking ##########
    major_settings = current_location + '/default_settings/xcms_peakpicking_default.txt'
    if not path.exists(major_settings):
        print('xcms_default_peakpicking settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('Unique_name =  Default \n')
        setting.write('path_to_data_containing folders =  <EDIT ME PLEASE> \n')
        setting.write('output_folders =  <EDIT ME PLEASE> \n')
        setting.write('sample_metadata_file_name =  <EDIT ME PLEASE> \n')
        setting.write('sample_name_column =  <EDIT ME PLEASE> \n')
        setting.write('data_file_extention =  mzML \n')
        setting.write('plot_allignment =  FALSE \n')

        setting.write('######## PEAK PICKING PARAMS #########\n')
        setting.write('ppm = 25 \n')
        setting.write('peakwidth = 20,50 \n')
        setting.write('snthresh = 10 \n')
        setting.write('mzCenterFun = WMean \n')
        setting.write('integrate = 1L \n')
        setting.write('mzdiff = -0.001 \n')
        setting.write('fitgauss   = FALSE \n')
        setting.write('noise  = 0 \n')
        setting.write('verboseColumns  = FALSE \n')
        setting.write('firstBaselineCheck  = TRUE \n')
        setting.write('extendLengthMSW = FALSE \n')

        setting.write('######## Refine found Peaks #########\n')
        setting.write('expandRt   = 2 \n')
        setting.write('expandMz   = 0 \n')
        setting.write('ppm = 10 \n')
        setting.write('minProp  = 0.75 \n')

        setting.write('######## Peak grouping #########\n')
        setting.write('bw = 30 \n')
        setting.write('minFraction = 0.5 \n')
        setting.write('minSamples  = 1 \n')
        setting.write('binSize = 0.25 \n')
        setting.write('maxFeatures = 50 \n')

        setting.write('######## Retention_time_correction #########\n')
        setting.write('binSize  = 1 \n')
        setting.write('centerSample  = 3 \n')
        setting.write('response   = 1L \n')
        setting.write('distFun  = cor_opt \n')
        setting.write('gapInit  = 0.5 \n')
        setting.write('binSize = 1 \n')
        setting.write('gapExtend = 2.5 \n')
        setting.write('factorDiag = 2 \n')
        setting.write('factorGap   = 1 \n')
        setting.write('localAlignment   = FALSE \n')
        setting.write('initPenalty  = 0 \n')

        setting.write('######## Retention_time_correction #########\n')
        setting.write('expandMz  = 0 \n')
        setting.write('expandRt    = 0 \n')
        setting.write('ppm    = 0 \n')
        setting.write('fixedMz   = 0 \n')
        setting.write('fixedRt    = 0 \n')
        setting.close()
def edit_setttings(current_location,item):
    major_settings = current_location + item
    major_settings_tmp = major_settings + '.TMP.txt'
    tmp_file = open(major_settings_tmp, 'w+')
    with open(major_settings)as f:
        for line in f:
            if '#' not in line:
                line = line.split('=')
                line1 = line[0].rstrip()
                line2 = line[1].rstrip()

                set = (input("enther path to "+ line1 +" , press enter for " + line2 + " : ") or line2)
                set = (line1 + ' = ' + set).replace('\\', '/') + "\n"
                tmp_file.write(set)
            else:
                print(line)
                tmp_file.write(line)
    tmp_file.close()
    os.remove(major_settings)
    os.renames(major_settings_tmp, major_settings)


def default_parameters_or_given(paths, item):
    loop4 = True
    while loop4:
        print("-------------Menu---------")
        print('1 = select parameters')
        print('2 = select settings file')
        print('3 = Exit')
        print("-------------Menu---------")

        choice = input("Enter your choice [1-3]: ")
        choice = str(choice)

        if choice == '1':
            print('enter unique name to keep the setting list different from the rest, will be saved in script settings')
            Unique = (input("unique name to remember by: "))

            while len(Unique) == 0:
                print('Unique name is empty, name is required!')
                Unique = (input("unique name to remember by: "))

            default_setting =  paths[4] + '/' + item
            script_settigns = paths[3] + '/' + Unique

            tmp_file = open(script_settigns, 'w+')
            with open(default_setting)as f:
                for line in f:
                    if '#' not in line:
                        line = line.split('=')
                        line1 = line[0].rstrip()
                        line2 = line[1].rstrip()

                        set = (input("enther path to " + line1 + " , press enter for " + line2 + " : ") or line2)
                        set = (line1 + ' = ' + set).replace('\\', '/') + "\n"
                        tmp_file.write(set)
                    else:
                        print(line)
                        tmp_file.write(line)
            tmp_file.close()
            return(script_settigns)



        elif choice == "2":
            file_found = False
            while not file_found:
                file = input("Type path to file: ")
                file = file.replace('\\', '/')

                if path.exists(file):
                    file_found = path.exists(file)
                    return (file)
                else:
                    print('could not find file check again')

        elif choice == "3":
            loop4 = False
        else:
            print('error in choise')
def R_script_editor_from_setting_list(path_to_settings, paths,script):
    settings = []

    with(open(path_to_settings, 'r')) as f:
        for line in f:
            if "#" not in line:
                line = ((line.split('='))[1]).rstrip().strip()
                settings.append(line)

    #print(settings)
    R = paths[1] + '/' + script
    R_tmp = R +'.TMP.R'
    counter = 0

    lines = []
    with open(R,'r') as f:
        for line in f:
            lines.append(line)


    f = open(R_tmp, 'w')
    for setting in settings:
        counter += 1
        replace_item = "!@#$%^&" + str(counter) + "&^%$#@!"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,setting)

    with open(R_tmp,'w') as f:
        for line in lines:
            f.write(line)

    return(R_tmp)

iniatie_check(current_location)
#read in major settings#
paths = []
with(open(major_settings,'r')) as f:
    for line in f:
        paths.append((((line.split('='))[1]).rstrip()).replace(" ", ""))

loop = True
while loop:  ## While loop which will keep going until loop = False
    print_menu()  ## Displays menu
    choice = input("Enter your choice [1-5]: ")
    choice = str(choice)

    if choice == '1':
        print("edit Settings has been selected")
        loop2 = True
        while loop2:
            print_menu_settings()
            choice = input("Enter your choice [1-10]: ")
            if choice == '1':
                item = '/major_settings.txt'
                edit_setttings(current_location,item)

            elif choice == '2':
                print(current_location)
                current_location1 = paths[4]
                item =  "/xcms_peakpicking_default.txt"
                edit_setttings(current_location1, item)


            elif choice == "10":
                loop2 = False

            else:
                # Any integer inputs other than values 1-5 we print an error message
                print('error in choise')

    elif choice == "2":
        print("Full Run has been selected")

    elif choice == "3":
        print("Batch run")

        ## You can add your code or functions here
    elif choice == "4":
        loop1 = True
        while loop1 == True:
            print_menu_parts()
            choice1 = input("Enter your choice [1-9]: ")
            choice1 = str(choice1)

            if choice1 == '1':
                print("Peakpicking has been selected")
                item = "xcms_peakpicking_default.txt"
                script = "xcms_peakpicking.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Peakpicking iniatiated')
                subprocess.check_call([paths[5], R_temp], shell=False)

                print('cleaning up files')
                os.remove(R_temp)

            elif choice1 == "2":
                print("Batch Correction has been selected")


            elif choice1 == "3":
                print("Noise filtering has been selected")


            elif choice1 == "4":
                print("Univariate testing  has been selected")


            elif choice1 == "5":
                print("Multivariate testing has been selected")

            elif choice1 == '6':
                print("Feature Filtering has been selected")
                ## You can add your code or functions here


            elif choice1 == "7":
                print("Database Annotation has been selected")


            elif choice1 == "8":
                print("Comparisation has been selected ")


            elif choice1 == "9":
                print("exit mode used")
                ## You can add your code or functions here
                loop1 = False  # This will make the while loop to end as not value of loop is set to False


            else:
                # Any integer inputs other than values 1-5 we print an error message
                print('error in choise')


    elif choice == "5":
        print("exit mode used")
        ## You can add your code or functions here
        loop = False  # This will make the while loop to end as not value of loop is set to False
    else:
        # Any integer inputs other than values 1-5 we print an error message
        print('error in choise')
