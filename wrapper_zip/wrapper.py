import subprocess
import os
from os import path
import platform

current_location = os.getcwd()
current_location = current_location.replace('\\', '/')
major_settings = current_location + '/major_settings.txt'

OS = platform.system()

if OS == "Windows" or "Linux":
    print('operating system ' + str(OS) + ' detected')
else:
    print('current detected OS is ' + str(OS) + ' this wrapper is not made for it \n Call the author of the script to request changes')
    exit()

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
    print("3.  Pre univariate filtering  ")
    print("4.  Univariate testing ")
    print("5.  Multivariate testing ")
    print("6.  Post univariate Filtering ")
    print("7.  Database Annotation ")
    print("8.  Comparisation ")
    print("9.  Venn diagram maker")
    print("10.  Exit")
    print(67 * "-")
def print_menu_settings():  ## Your menu design here
    print(30 * "-", "EDIT SETTINGS", 30 * "-")
    print("1.  Major settings ")
    print("2.  Peakpicking ")
    print("3.  Batch Correction")
    print("4.  Pre univariate filtering ")
    print("5.  Univariate testing ")
    print("6.  Multivariate testing ")
    print("7.  Post univariate Filtering ")
    print("8.  Database Annotation ")
    print("9.  Comparisation ")
    print("10. Venn Diagram ")
    print("11. Full run ")
    print("12. Full run _ database annotation ")
    print("13. MS2 analysis ")
    print("14. Exit")
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
        if not path.exists(tmp):
            os.mkdir(tmp)
        setting.write('batch_location = '+  tmp + "\n")
        tmp = current_location + '/script_settings'
        if not path.exists(tmp):
            os.mkdir(tmp)
        setting.write('script_settings_location = ' + tmp + "\n")
        tmp = current_location + '/default_settings'
        setting.write('default_settings_location = ' + tmp + "\n")

        if OS == "windows":
            try:
                print('searching for Rscript.exe in C://')
                Rscript_paths = []
                for r, d, f in os.walk("c:\\"):
                    for files in f:
                        if files == "Rscript.exe":
                            line = (os.path.join(r, files))
                            line = line.replace("\\", "/")
                            Rscript_paths.append(line)
                if len(Rscript_paths) == 0:
                    print('NO Rscript.exe found, withouth it script will not work, consider installing R')
            except:
                print('unable to find Rscript.exe in C: disk. if you are running this on linux dont worry, perfectly normal')
            Rscript_paths_1 = []
            for item in Rscript_paths:
                if "bin/Rscript.exe" in item:
                    Rscript_paths_1.append(item)
            sorted(Rscript_paths_1)
            Rscript_path = (Rscript_paths_1[-1])

        else:
            Rscript_path = 'NA'
        setting.write('Rscript_location = ' + str(Rscript_path) + "\n")

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
        setting.write('sample name column =  sample_name \n')
        setting.write('sample type column =  sample_type \n')
        setting.write('QC  in sample type =  QC \n')
        setting.write('blank  in sample type =  blank \n')
        setting.write('sample  in sample type =  sample \n')
        setting.write('injection order column name =  injectionOrder \n')
        setting.write('batch number column name =  batch \n')
        setting.write('sample_metadata_file_name =  sample_metadata.tsv \n')
        setting.write('data_file_extention =  mzML \n')
        setting.write('Polarity, choose between [1-2] 1-negative and 2-postive =  1 \n')

        setting.write('######## PEAK PICKING PARAMS #########\n')
        setting.write('ppm = 25 \n')
        setting.write('peakwidth = 20,50 \n')
        setting.write('snthresh = 10 \n')
        setting.write('mzCenterFun choise [1-4]. 1-wmean, 2-mean, 3-apex, 4-wmeanApex3 = 1 \n')
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
        setting.write('maxFeatures = 1000 \n')

        setting.write('######## Retention_time_correction #########\n')
        setting.write('binSize  = 1 \n')
        setting.write('centerSample  = 3 \n')
        setting.write('response   = 1L \n')
        setting.write('distFun choose [1-5] 1-cor 2-cor_opt 3-cov 4-prd 5-euc = 2 \n')
        setting.write('gapInit  = 0.5 \n')
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

    ########## batchCorrection ##########
    major_settings = current_location + '/default_settings/batch_correction_default.txt'
    if not path.exists(major_settings):
        print('batch_correction_default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('input folder =  F:/avans/stage MM/xcms_pipeline/_XCMS_default \n')
        setting.write(' DataMatrix file name=  Data_matrix_XCMS_default.tsv \n')
        setting.write('SampleMetada file name =  sample_meta_data_XCMS_default.tsv \n')
        setting.write('Variable Metadata file name =  Variable_metaData_XCMS_default.tsv \n')
        setting.write('name of output folder =  Batch_correction \n')
        setting.write('######### Best Corrected dataset ######## \n')
        setting.write('Select best parameter settings on 1-mean_var_QC_log10 2-Ratio_VAR_QC_log = 1 \n')
        setting.write('######### Graph settings ######## \n')
        setting.write('Pixelsize 1 =  20 \n')
        setting.write('Pixelsize 2 =  12 \n')
        setting.close()

    ########## pre_univariate filtering ##########
    major_settings = current_location + '/default_settings/Pre_univariate_filtering_default.txt'
    if not path.exists(major_settings):
        print('Pre_univariate_filtering default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('input folder =  F:/avans/stage MM/Sherloktest_data_2/_XCMS_Default/Batch_correction/ \n')
        setting.write('Output folder name =  pre_univriate_filterd \n')
        setting.write('DataMatrix file name=  Data_matrix_XCMS_default.tsv_batchcorrected.tsv \n')
        setting.write('SampleMetada file name =  sample_meta_data_XCMS_default.tsv_batchcorrected.tsv \n')
        setting.write('Variable Metadata file name =  Variable_metaData_XCMS_default.tsv_batchcorrected.tsv \n')
        setting.write('######## Filters #########\n')
        setting.write('Filter on isotopes (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on Relative Standard Deviation (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter Noise (TRUE or FALSE) =  TRUE \n')
        setting.write('######## Isotope filtering #########\n')
        setting.write('Deviation in isotope mass search =  0.003 \n')
        setting.write('Deviation in isotope Retention time search =  0.003 \n')
        setting.write('######## RSD filtering #########\n')
        setting.write('Relative standard deviation Treshold =  25 \n')
        setting.write('######## Noise filtering #########\n')
        setting.write('minimal amout of hits found in samples =  1 \n')
        setting.close()

    ########## Univariate ##########
    major_settings = current_location + '/default_settings/Univariate_default.txt'
    if not path.exists(major_settings):
        print('Univariate default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('input folder =  F:/avans/stage MM/pipeline_testrun/_XCMS_m/Batch_correction \n')
        setting.write('DataMatrix file name=  XCMS_default_batchcorrected_noice_reduced_matrix.tsv \n')
        setting.write('SampleMetada file name =  XCMS_default_batchcorrected_noice_reduced_sample_metadata.tsv \n')
        setting.write('Variable Metadata file name =  XCMS_default_batchcorrected_noice_reduced_variable_metadata.tsv\n')

        setting.write('######## Statistical test to perform #########\n')
        setting.write('Statistical test are 1-ttest 2-limma 3-wilcoxon 4-anova 5-kruskal 6-pearson 7-spearman 8-limma2ways 9-limma2waysInter 10-anova2ways 11-anova2waysInter =  1\n')
        setting.write('Variable(s) of interest ( comma delimited) = gender,age,bmi  \n')
        setting.write('Correct dataset according variable (NULL results in skipping) = NULL  \n')
        setting.write('Main factor of interest = gender  \n')
        setting.write('secondary factor of interest = age  \n')
        setting.write('P-value Tresh_hold = 0.05  \n')
        setting.write('max features as output? (NULL results in skipping) = NULL   \n')

        setting.write('######## Graphical #########\n')
        setting.write('graph title = univariate testing is super cool  \n')
        setting.write('prefix of reports = this makes sure that it will not overwrite other reports \n')
        setting.write('charactersize 1 = 20  \n')
        setting.write('charactersize 2 = 12  \n')

        setting.write('######## Filtering data #########\n')
        setting.write('Cutoff Hotellings P-value (lower will sample be deleted) = 0.001  \n')
        setting.write('Cutoff miss P-value (lower will sample be deleted) = 0.001  \n')
        setting.write('Cutoff deci P-value (lower will sample be deleted) = 0.001  \n')

        setting.write('######## Heathmap #########\n')
        setting.write('amount of clusters for samples = 5  \n')
        setting.write('amount of clusters for features = 5  \n')
        setting.write('heathmap statistics (1-pearson, 2-kendall, 3-spearman) choose[1:3] = 1  \n')
        setting.write('heathmap Algorithm (1-euclidean 2-maximum 3-manhattan 4-canberra 5-binary 6-minkowski 7-(1-cor) 8-(1-abs(cor)) choose[1:8] = 7  \n')

        setting.close()

    ########## Multiivariate ##########
    major_settings = current_location + '/default_settings/multivariate_default.txt'
    if not path.exists(major_settings):
        print('multivariate default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('input folder =  F:/avans/stage MM/pipeline_testrun/_XCMS_m/Batch_correction/Noise_filterd \n')
        setting.write('DataMatrix file name=  XCMS_m_batchcorrected_noice_reduced_matrix.tsv \n')
        setting.write('SampleMetada file name =  XCMS_m_batchcorrected_noice_reduced_sample_metadata.tsv \n')
        setting.write('Variable Metadata file name =  XCMS_m_batchcorrected_noice_reduced_variable_metadata.tsv \n')

        setting.write('######## Multivariate testing #########\n')
        setting.write('Variable Metadata variables_of_interest (comma delimited) =  age,gender,bmi \n')
        setting.write('Main factor of interest  =  gender \n')
        setting.write('second factor of interest  =  age \n')
        setting.write('Numerical factor of interest (or NULL) =  NULL \n')


        setting.close()

    ########## post univarite filtering ##########
    major_settings = current_location + '/default_settings/Post_univariate_filtering_default.txt'
    if not path.exists(major_settings):
        print('Post_univariate_filtering_default default default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('input folder =  F:/avans/stage MM/pipeline_testrun/_XCMS_m/Batch_correction/Noise_filterd/ttest_gender \n')
        setting.write('DataMatrix file name=  XCMS_m_batchcorrected_noice_reduced_matrix.tsv_ttest_gender.tsv \n')
        setting.write('SampleMetada file name =  XCMS_m_batchcorrected_noice_reduced_sample_metadata.tsv_ttest_gender.tsv \n')
        setting.write('Variable Metadata file name =  XCMS_m_batchcorrected_noice_reduced_variable_metadata.tsv_ttest_gender.tsv \n')
        setting.write('output folder name =  filterd \n')

        setting.write('######## Select Filter Methods #########\n')
        setting.write('Filter on significant hits (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on QC/Sample ratio features found (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on blank/Sample ratio features found (TRUE or FALSE) =  TRUE \n')

        setting.write('######## Filtering Details #########\n')
        setting.write('min qc/sample ratio =  0.1 \n')
        setting.write('max qc/sample ratio =  2.0 \n')
        setting.write('max blank/sample ratio =  0.0 \n')
        setting.write('max blank/sample ratio =  1.0 \n')
        setting.close()

    ########## database annotation ##########
    major_settings = current_location + '/default_settings/database_annotation_default.txt'
    if not path.exists(major_settings):
        print('database annotation default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('path to variable metadata =  F:/avans/stage MM/pipeline_testrun/_XCMS_m/Batch_correction/Noise_filterd/ttest_gender/filterd/noice_reduced_sample_metadata_ttest_gender_filterd_variable_metadata.tsv \n')
        setting.write('output location =  F:/avans/stage MM/Sherloktest_data_2/annotated/ \n')
        setting.write('Search range =  0.003 \n')
        setting.write('Limit output from same mass =  2 \n')
        setting.close()

    ########## annotation comparison ##########
    major_settings = current_location + '/default_settings/annotation_comparison_default.txt'
    if not path.exists(major_settings):
        print('database annotation default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Annotation_comparison #########\n')
        setting.write('path to folder containing annotated xlsx.files metadata = F:/avans/stage MM/Annotated_outpet_sherlOK/ \n')
        setting.write('output location =  F:/avans/stage MM/Annotated_outpet_sherlOK/Mjv_data \n')
        setting.write('row_number =  7 \n')
        setting.write('Round numbers down = True \n')
        setting.write('Howmany decimals in round down = 3 \n')
        setting.close()

    ########## Venn Diagram maker ##########
    major_settings = current_location + '/default_settings/Venn_diagram_maker_default.txt'
    if not path.exists(major_settings):
        print('database annotation default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Venn_diagram_maker #########\n')
        setting.write('path to file containing the anonotation comparison = F:/avans/stage MM/Annotated_outpet_sherlOK/MJV.csv \n')
        setting.write('output location=  F:/avans/stage MM/Annotated_outpet_sherlOK/ \n')
        setting.write('Pixel Size =  12 \n')
        setting.close()

    ########## Full run MS1 ##########
    major_settings = current_location + '/default_settings/Full_run_default.txt'
    if not path.exists(major_settings):
        print('Full run default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Data information #########\n')
        setting.write('path to folder containing data = F:/avans/stage MM/Sherloktest_data_2/ \n')
        setting.write('output location=  F:/avans/stage MM/full_test_run/ \n')
        setting.write('unique name = default \n')

        setting.write('######## Variable metadata information #########\n')
        setting.write('sample name header = sample_name \n')
        setting.write('sample type header = sampleType \n')
        setting.write('QC in sampletype = pool \n')
        setting.write('blank in sampletype = blank \n')
        setting.write('sample in sampletype = sample \n')
        setting.write('injection order column = injectionOrder \n')
        setting.write('batch number column = batch \n')
        setting.write('sample metadata file name = sample_metadata.tsv \n')
        setting.write('data file extenction = mzXML \n')
        setting.write('Polarity 1-negative 2-positive = 2 \n')

        setting.write('######## PEAK PICKING PARAMS #########\n')
        setting.write('ppm = 25 \n')
        setting.write('peakwidth = 20,50 \n')
        setting.write('snthresh = 10 \n')
        setting.write('mzCenterFun choise [1-4]. 1-wMean, 2-mean, 3-apex, 4-wmeanApex3 = 1 \n')
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
        setting.write('maxFeatures = 1000 \n')

        setting.write('######## Retention_time_correction #########\n')
        setting.write('binSize  = 1 \n')
        setting.write('centerSample  = 3 \n')
        setting.write('response   = 1L \n')
        setting.write('distFun choose [1-5] 1-cor 2-cor_opt 3-cov 4-prd 5-euc = 2 \n')
        setting.write('gapInit  = 0.5 \n')
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

        setting.write('######## pixelsize for graphs #########\n')
        setting.write('pixelsize1  = 20 \n')
        setting.write('pixelsize2  = 12 \n')

        setting.write('######## batch correction #########\n')
        setting.write('select best batch correction on 1-mean variation on QC 2-Ratio var total_QC  = 1 \n')

        setting.write('######## Filters #########\n')
        setting.write('Filter on isotopes (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on Relative Standard Deviation (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter Noise (TRUE or FALSE) =  TRUE \n')

        setting.write('######## Isotope filtering #########\n')
        setting.write('Deviation in isotope mass search =  0.003 \n')
        setting.write('Deviation in isotope Retention time search =  0.003 \n')

        setting.write('######## RSD filtering #########\n')
        setting.write('Relative standard deviation Treshold =  25 \n')
        setting.write('######## Noise filtering #########\n')
        setting.write('minimal amout of hits found in samples =  1 \n')

        setting.write('######## Univariate testing #########\n')

        setting.write('######## Statistical test to perform #########\n')
        setting.write('Statistical test are 1-ttest 2-limma 3-wilcoxon 4-anova 5-kruskal 6-pearson 7-spearman 8-limma2ways 9-limma2waysInter 10-anova2ways 11-anova2waysInter =  1\n')
        setting.write('Variable(s) of interest ( comma delimited) = gender,age,bmi  \n')
        setting.write('Correct dataset according variable (NULL results in skipping) = NULL  \n')
        setting.write('Main factor of interest = gender  \n')
        setting.write('secondary factor of interest = age  \n')
        setting.write('P-value Tresh_hold = 0.05  \n')
        setting.write('max features as output? (NULL results in skipping) = NULL   \n')

        setting.write('######## Graphical #########\n')
        setting.write('graph title = univariate testing is super cool  \n')
        setting.write('prefix of reports = this makes sure that it will not overwrite other reports \n')

        setting.write('######## Filtering data #########\n')
        setting.write('Cutoff Hotellings P-value (lower will sample be deleted) = 0.001  \n')
        setting.write('Cutoff miss P-value (lower will sample be deleted) = 0.001  \n')
        setting.write('Cutoff deci P-value (lower will sample be deleted) = 0.001  \n')

        setting.write('######## Heathmap #########\n')
        setting.write('amount of clusters for samples = 5  \n')
        setting.write('amount of clusters for features = 5  \n')
        setting.write('heathmap statistics (1-pearson, 2-kendall, 3-spearman) choose[1:3] = 1  \n')
        setting.write('heathmap Algorithm (1-euclidean 2-maximum 3-manhattan 4-canberra 5-binary 6-minkowski 7-(1-cor) 8-(1-abs(cor)) choose[1:8] = 7  \n')

        setting.write('######## Select Filter Methods #########\n')
        setting.write('Filter on significant hits (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on QC/Sample ratio features found (TRUE or FALSE) =  TRUE \n')
        setting.write('Filter on blank/Sample ratio features found (TRUE or FALSE) =  TRUE \n')

        setting.write('######## Filtering Details #########\n')
        setting.write('min qc/sample ratio =  0.1 \n')
        setting.write('max qc/sample ratio =  2.0 \n')
        setting.write('max blank/sample ratio =  0.0 \n')
        setting.write('max blank/sample ratio =  1.0 \n')
        setting.close()

    ########## MS2_analysis ##########
    major_settings = current_location + '/default_settings/MS2_analysis_default.txt'
    if not path.exists(major_settings):
            print('MS2_analysis default settings file is not found, will be created')
            setting = open(major_settings, 'w+')
            setting.write('######## Ms_2_analysis #########\n')
            setting.write( 'path to folder containing raw data files = F:/avans/stage MM/Sherloktest_data_2/ \n')
            setting.write( 'data file extention = mzXML \n')
            setting.write( 'path to variable metadata = F:/avans/stage MM/test1/_XCMS_default1/Variable_metaData_XCMS_default1.tsv \n')

            setting.write('path to output folder =  F:/avans/stage MM/test1 \n')
            setting.write('Unique name =  default \n')
            setting.close()

    ########## database annotation ##########
    major_settings = current_location + '/default_settings/database_annotation_default_full_run.txt'
    if not path.exists(major_settings):
        print('database annotation default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('output location =  F:/avans/stage MM/Sherloktest_data_2/annotated/ \n')
        setting.write('Search range =  0.003 \n')
        setting.write('Limit output from same mass =  2 \n')
        setting.close()

    ########## database annotation Full run ##########
    major_settings = current_location + '/default_settings/database_annotation_default_fullrun.txt'
    if not path.exists(major_settings):
        print('database annotation default settings file is not found, will be created')
        setting = open(major_settings, 'w+')
        setting.write('######## Global_settings #########\n')
        setting.write('output location =  F:/avans/stage MM/Sherloktest_data_2/annotated/ \n')
        setting.write('Search range =  0.003 \n')
        setting.write('Limit output from same mass =  2 \n')
        setting.write('######## PEAK PICKING PARAMS #########\n')
        setting.write('ppm = 25 \n')
        setting.write('peakwidth = 20,50 \n')
        setting.write('snthresh = 10 \n')
        setting.write('mzCenterFun choise [1-4]. 1-wmean, 2-mean, 3-apex, 4-wmeanApex3 = 1 \n')
        setting.write('integrate = 1L \n')
        setting.write('mzdiff = -0.001 \n')
        setting.write('fitgauss   = FALSE \n')
        setting.write('noise  = 0 \n')
        setting.write('verboseColumns  = FALSE \n')
        setting.write('firstBaselineCheck  = TRUE \n')
        setting.write('extendLengthMSW = FALSE \n')
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
            script_settigns = paths[3] + '/' + Unique +'.txt'

            tmp_file = open(script_settigns, 'w+')
            with open(default_setting)as f:
                for line in f:
                    if '#' not in line:
                        line = line.split('=')
                        line1 = line[0].rstrip()
                        line2 = line[1].rstrip()

                        set = (input("enter parameter for " + line1 + " , press enter for " + line2 + " : ") or line2)
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
                print(line)
                line = ((line.split('='))[1]).rstrip().strip()
                settings.append(line)

    print(settings)
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

    if script == "full_run.R":
        replace_item = "####scripts_location####"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[1]+'/'))


    with open(R_tmp,'w') as f:
        for line in lines:
            f.write(line)

    return(R_tmp)
def PY_script_editor_from_setting_list(path_to_settings, paths,script):
    settings = []

    with(open(path_to_settings, 'r')) as f:
        for line in f:
            if "#" not in line:
                print(line)
                line = ((line.split('='))[1]).rstrip().strip()
                settings.append(line)

    print(settings)
    #print(settings)
    R = paths[1] + '/' + script
    R_tmp = R +'.TMP.py'
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

    if script == "database_annotation.py":
        replace_item = "###database_location###"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[0]+'/'))

    if script == "database_annotation_fullrun.py":
        replace_item = "####scripts_location####"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[1]+'/'))

        replace_item = "###database_location###"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[0]+'/'))

    if script == "database_annotation_batch.py":
        replace_item = "####scripts_location####"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[1]+'/'))

        replace_item = "###database_location###"
        counterA = -1
        for line in lines:
            counterA += 1
            lines[counterA] = line.replace(replace_item,(paths[0]+'/'))


    with open(R_tmp,'w') as f:
        for line in lines:
            f.write(line)

    return(R_tmp)


iniatie_check(current_location)
#read in major settings#
paths = []
with(open(major_settings,'r')) as f:
    for line in f:
        paths.append((((line.split('='))[1]).rstrip()).strip())

Rscript = paths[5]

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

            elif choice == '3':
                print(current_location)
                current_location1 = paths[4]
                item =  "/batch_correction_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '4':
                print(current_location)
                current_location1 = paths[4]
                item =  "/Pre_univariate_filtering_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '5':
                print(current_location)
                current_location1 = paths[4]
                item =  "/Univariate_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '6':
                print(current_location)
                current_location1 = paths[4]
                item =  "/multivariate_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '7':
                print(current_location)
                current_location1 = paths[4]
                item =  "/Post_univariate_filtering_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '8':
                print(current_location)
                current_location1 = paths[4]
                item =  "/database_annotation_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '9':
                print(current_location)
                current_location1 = paths[4]
                item =  "/annotation_comparison_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '10':
                print(current_location)
                current_location1 = paths[4]
                item =  "/Venn_diagram_maker_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '11':
                print(current_location)
                current_location1 = paths[4]
                item =  "/Full_run_default.txt"
                edit_setttings(current_location1, item)

            elif choice == '12':
                print(current_location)
                current_location1 = paths[4]
                item =  "/database_annotation_default_full_run.txt"
                edit_setttings(current_location1, item)

            elif choice == '12':
                print(current_location)
                current_location1 = paths[4]
                item =  "/database_annotation_default_full_run.txt"
                edit_setttings(current_location1, item)

            elif choice == '13':
                print(current_location)
                current_location1 = paths[4]
                item =  "/MS2_analysis_default.txt"
                edit_setttings(current_location1, item)

            elif choice == "14":
                loop2 = False

            else:
                # Any integer inputs other than values 1-5 we print an error message
                print('error in choise')

    elif choice == "2":
        print("Full Run MS 1 has been selected")
        item = "Full_run_default.txt"
        script = "full_run.R"
        item1 = "database_annotation_default_full_run.txt"
        script1 = "database_annotation_fullrun.py"


        path_to_settings = default_parameters_or_given(paths, item)
        path_to_settings1 = default_parameters_or_given(paths, item1)

        R_temp = R_script_editor_from_setting_list(path_to_settings, paths, script)
        R_temp1 = PY_script_editor_from_setting_list(path_to_settings1, paths, script1)

        if OS == "Windows":
            print(R_temp)
            subprocess.check_call([Rscript, R_temp], shell=False)
            print(R_temp1)
            os.system('python ' + R_temp1)
        elif OS == "linux":
            print(R_temp)
            os.system('Rscript ' + R_temp)
            print(R_temp1)
            os.system('python3 ' + R_temp1)

        print('cleaning up files')
        try:
            os.remove(R_temp1)
        except:
            print("failed to delete " + R_temp1)

        try:
            os.remove(R_temp)
        except:
            print("failed to delete " + R_temp)

    elif choice == "3":
        print("Batch run")
        script = "full_run.R"
        script1 = "database_annotation_batch.py"


        files = (os.listdir(paths[2]))
        for file in files:
            path_to_settings = (paths[2]+'/'+file)

            R_temp = R_script_editor_from_setting_list(path_to_settings, paths, script)
            R_temp1 = PY_script_editor_from_setting_list(path_to_settings, paths, script1)

            if OS == "Windows":
                print(R_temp)
                subprocess.check_call([Rscript, R_temp], shell=False)
                print(R_temp1)
                os.system('python ' + R_temp1)
            elif OS == "linux":
                print(R_temp)
                os.system('Rscript ' + R_temp)
                print(R_temp1)
                os.system('python3 ' + R_temp1)

            print('cleaning up files')
            try:
                os.remove(R_temp1)
            except:
                print("failed to delete " + R_temp1)
            try:
                os.remove(R_temp)
            except:
                print("failed to delete " + R_temp)

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
                print('Peakpicking iniatated')
                #string = ("R-3.6.1 --slave "+ str(R_temp))
               # os.system(string)

                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)



                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)


            elif choice1 == "2":
                print("Batch Correction has been selected")
                item = "batch_correction_default.txt"
                script = "Batch_correction.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Batch_correction iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print ("failed to delete " + R_temp)




            elif choice1 == "3":
                print("Pre univariate filtering has been selected")
                item = "Pre_univariate_filtering_default.txt"
                script = "Pre_univariate_filtering.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Batch_correction iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)
            elif choice1 == "4":
                print("Univariate testing  has been selected")
                item = "Univariate_default.txt"
                script = "Univariate.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Univariate testing iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)
            elif choice1 == "5":
                print("Multivariate testing has been selected")
                item = "multivariate_default.txt"
                script = "Multivariate.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Batch_correction iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)
            elif choice1 == '6':
                print("Post Univariate Filtering has been selected")
                item = "Post_univariate_filtering_default.txt"
                script = "Post_univariate_filtering.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Batch_correction iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)
            elif choice1 == "7":
                print("Database Annotation has been selected")
                item = "database_annotation_default.txt"
                script = "database_annotation.py"

                path_to_settings = default_parameters_or_given(paths, item)
                R_temp = PY_script_editor_from_setting_list(path_to_settings, paths, script)
                print(R_temp)

                if OS == "Windows":
                    print(R_temp)
                    os.system('python ' + R_temp)
                elif OS == "linux":
                    print(R_temp)
                    os.system('python3 ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)

            elif choice1 == "8":
                print("Comparisation has been selected ")
                item = "annotation_comparison_default.txt"
                script = "annotation_comparison.py"

                path_to_settings = default_parameters_or_given(paths, item)
                R_temp = PY_script_editor_from_setting_list(path_to_settings, paths, script)
                print(R_temp)

                if OS == "Windows":
                    print(R_temp)
                    os.system('python ' + R_temp)
                elif OS == "linux":
                    print(R_temp)
                    os.system('python3 ' + R_temp)

                print('cleaning up files')

                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)

            elif choice1 == "9":
                print("Venn_diagram_maker has been selected ")
                item = "Venn_diagram_maker_default.txt"
                script = "Venn_diagram_maker.R"

                path_to_settings = default_parameters_or_given(paths,item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths,script)
                print(R_temp)
                print('Batch_correction iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)


            elif choice1 == "10":
                print("MS2 analysis has been selected ")
                item = "MS2_analysis_default.txt"
                script = "Ms2_analysis.R"

                path_to_settings = default_parameters_or_given(paths, item)
                R_temp = R_script_editor_from_setting_list(path_to_settings, paths, script)
                print(R_temp)
                print('MS2_analysis iniatated')
                if OS == "Windows":
                    print(R_temp)
                    subprocess.check_call([Rscript, R_temp], shell=False)
                elif OS == "linux":
                    print(R_temp)
                    os.system('Rscript ' + R_temp)

                print('cleaning up files')
                try:
                    os.remove(R_temp)
                except:
                    print("failed to delete " + R_temp)

            elif choice1 == "11":
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
