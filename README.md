# LC-MS2 data dynamic pipeline for metabolomic biomarker discovery

### Table of contents
* Wrapper introduction
* Wrapper function
* Full run
* Batch run
* part run
* dataset requirements
* databases
* part run scripts
  * XCMS peakpicking
  * batch correction
  * pre univariate filtering
  * univariate filtering
  * multivariate filtering
  * post univariate filtering
  * database annotation
  * annotation comparison
  * venn diagram maker
  * MS2 annotation
 * installation
 * net yet supported
## the wrapper
The LC-MS2 wrapper is not an ordinary data-pipeline.
Due the the huge amounts of variables and options is an dynamic pipeline created with an shell based interface
With the Shell based interface is it possible to analyse datasets by remote computing
The wrapper can be initialised by python and supports the Windows and Linux OS. 

#### wrapper functions

inside the wrapper are serveral inbuild functions to improve functionality
When the wrapper is excecuted for the first time it will generate a few importand setting files which worked when testing the data pipeline
These settings can be adjusted inside the Shell system also by just using copy and paste. paths will be directly put in the right format.
The settings can also be adjusted by opening the txt file and adjust the items "note keep the same layout or the wrapper will not know anymore what to do"
When running an script from the wrapper a choise menu will be opend. users can fill in their adjustments to the Default settings or press ENTER for the default setting.
After choosing the variables it will produce a txt file in the script_settings folder with the used settings to keep a track on which settings used.
These files can also be used as import as settings to the data_pipeline when running a script.

A batch run can be used. the files required for that can be found in the Batch folder with example setting file. This can be adjusted to server your needs. 

To Analyse the data a few options are made available
##### full run

with the full run the imported datasets will go through almost all available scripts.
It runs from XCMS peakpicking > batch Correction > pre univariate filtering > univariate testing > post univariate filtering > database annotation

![alt text](https://github.com/Dirowa/LC-MS-data_pipe_lines/blob/master/Full%20Run.png)

##### batch run
The batch run does the same as the full run for exception on that you will get a choice meneu to enter variables. for this the example files in the Batch folders has to be edited and possible duplicated for more runs. it is advised to test the used variables first with part run
![alt text](https://github.com/Dirowa/LC-MS-data_pipe_lines/blob/master/batch%20run.png)

##### part run 
in the part run all scripts can be manualy enterd or give a path to a file containing the settings. The Scripts are written in that way that i will accept all output formats from other scripts used for maximal flexibillity. the part runs also contains Majority Vote calculations which is needed to create venn diagrams. The input from the Majority vote script is the annotated xlsx. file generated by database annotation. used can decide on which factor it is wished to make a Majority Vote on. 

## datasets requirement
with the datasets is it important to have the following items in the Sample Meta data
It requires The sample names, Batch number, injection orders and sampletypes.
incase batch is unknown run serveral PCA plots todecide or fill in 1 with all samples.
The injection order is the order from oldest to youngest when data is analysed by the LC-MS
in the sampletypes are the following information required. Blanks, Qc (or pooled samples) and samples, Internstal standards are not yet incorperated into this data-pipeline.

## databases
the databases are generated by Using MetaboShiny.
Installing MetaboShiny can be quite troublesome, for that check "Tips_for_installing_metaboshiny_on_clean_linux.R" to install dependencies for a clean Linux system
Metaboshiny parsers small public databases and generetes from the SMILE formula the monoisotopic mass.
Later in the datapipeline compounds will be filterd for the monoisotopic mass
It is also possible to parse pubchem into metadobo Shiny. For this are 2 Python scripts created
1-  PubChem_parser_to import_to_extended.py can be run and the produces files have to be manually enterd in the extended.db from metaboshiny with SQLite3 commands.
2-  pubchem_parser_slow_window&linux.py does the same as PubChem_parser_to import_to_extended.py but the output files can be imported into Metaboshiny where it calculates the monoisotopic mass from the SMILES structure

## scripts in data pipeline.

#### XCMS peakpicking
XCMS is chosen as peakpicking algorithm due the high custombillity, Relative processing speed, RAM and CPU usage and the amount of covering significant features. 
In this script. peaks will be analysed by using the CentwaveParam, the found features will be merged together incase of splitted peaks  by the MergeNeighboringPeaksParam.
All found peaks will be corrected for their retention time based on well behaving peak (peaks found across all samples) which is done by the Obiwarpparam algorithm. After the retention time correction the FillChrompeaksparam will look back on all samples to fill in missing peak of which first was not detected. All peaks will then be grouped into features. The found Features will be analysed by the CAMERA package to find Adducts, Isotope labeling and group them into pseudospectras.
The dataframes names and values will be edited so other scripts will accepts them and not fall into an error. This reduces the amount of unnesacary filling variables.

#### Batch Correction
The generated Variable metadata, Datamatrix and sample metadata is imported from XCMS into the batch correction script.
It checks the data wether it insane or not. and produce some PCA plots before and after batch correction. All data will be Log10 transformed and missing datasets will have their NA transformed in 1/5 of the lowest found intensity across a feature. this produces a file called the new_normalised_set.csv. With that file multiple batch correction will be tested against each other. that includes : (Loess, Combat, WaveICA, EigenMS, QC_RLSC and Ancova). Based on the amount of QC samples, blanks samples, found features, injection orders and Batches an algorithm will succeed. Internal standards can also be incorperated into the script but is still on the TO DO list also becouse there isnt any good training data set available who contains internal standards. For each dataset an different batch correction is the best. To decide which one is the best, from all algorithms will the QC variance and total variance calculated and devided to each other. the user can select between lowest QC variante or lowest total/QC variance as parameter to decide which algorithm corrected the dataset the best. The best selected data matrix will be used to generate some additional plots and will be worked further in the data-pipeline. the best dataset can be recognised as data_matrix, variable_metadata, sample_metadata file with an addition of batchcorrected.tsv


#### Pre univariate filtering
Finding huge amount of features is wonderfull to have but not neccasery good. To prevent having a biased univariate testing the datasets will be filterd accordingly.
As first all isotopes will be filterd out thus that only monoisotopic mass will reside in the dataset. this is also in corrospondes with the database where the mass is calculated for the monoisotopic mass. all found features will be filterd on the Relative Standard Deviation to prevent outliers. Thereby all background noise will be filterd. For each feature in the blank sample will 2 times the intensity (normalised back to normal numbers and not Log10) used as cutoff for the features in other samples and be transformed to NA. when all intensities in a features has been reduced to NA, the feature will be deleted. if the feature isnt deleted the NA will be transformed to 1/5 of the lowest found intensity in the feature group.

#### Univariate testing
in the univariate testing the user can select a variable to test against with multiple statistical tests available. As fist the data will be inspected and filterd on only Sample intensities but also on Hotelings Pvalue, Decimal P-value and Missing P-value to get a more uniform dataset to test against. with each filtering a summary plot will be produced. The univariate testing can be done on : ("ttest", "limma", "wilcoxon", "anova", "kruskal", "pearson", "spearman","limma2ways", "limma2waysInter", "anova2ways", "anova2waysInter") and a siginificance P-value must be selected. Thereby can the user also set up multiple variables of intereset to generate more statistical plots. In this PCA and OPLS-DA modeling is implemented

#### Multivariate testing
the multivariate tool tests the dataset in multiple ways. It generates PCA plots and OPLS_DA plots.
It is also possible to predict a model with overfitting

#### Post univariate testing
in this script the user can filter the features weither it was significant.
what the QC/sample ratio is and what the Blank/sample ratio is.
This also generates boxplotes with the amount of hits in the sample, blank and QC 

#### database annotation
feature list can be annotated by an database which is done by a python script. 
In this script the users has to define the location of the databases. it is possible and advised to limit the amount of hits between the given ppm to an amount the user is comfortable. The PPM is the range the found mass will be queried in the database. where all found masses will be orderd to the difference of the initial mass to found mass. so the closests match is on top. This file produces an easy to read xlsx file.

#### Annotation comparison

when the user has multiple runs on the same datasets and want to make a comparison between various settings is it possible to import the annotated xlsx file. The users has to select on which the path containg the folder of all xlsx files. The user can dicide on wich column to make a comparison, (example choose between inital masses, closest match, SMILES structure found or compound ID. If the data is a number the user can also select to round down the number to a given amount of decimal, when the data is text the round down has to be set to False. The script will calculate the Majority votes of the given parameter which is also the input of the Venn Diagram maker to gain more insight in how much the algorithms agree with each other.

#### Venn Diagram maker

The venn diagram maker produces a venn diagramm in Euller style. this means that the size of the circles are in contrast with how big they are.

#### MS2 annotation

it is possible to annotate the dataset on MS2 level to validate the MS1 level identification. In the dataset the user has to give a feature list on MS1 level to reducome computing time. Also does it require the same CentwaveParam settings from peakpicking wich generated the feature list.
The script search with the retention time and mass of the features of intereset  on Ms2 level. when MS2 candidates are found it will be alligned. when it matches for more then 90% is is considerd valid. This script produces a MS2 feature list of the features of intereset the user has given. The Ms2 feature list can be annotated by the databases generated by MetaboShiny. So far is the option to use spectral databases not yet incorperated


## instalation

Download the Wrapper.rar and unzip in your favourite folder
Check the dependencies and install them in R
Check the tips for installing metaboshiny
launch the wrapper.py and edit the Major settings

## not yet supported
The following items are not yet supported in
* Internal standards in batch correction
* MS2 analyse with a spectrum database
* XCMS peakpicking file who will help with setting variables
