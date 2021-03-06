## Start with this staging file to set up your analysis. 
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
#setwd("")

# Set the output directory for graphs:
graph_output_dir = 'graphs'
# Set the output directory for statistics:
stats_output_dir = 'stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

####################
## File locations ##
####################
## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine. 

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

# Where is the metadata file stored on your machine?
amr_metadata_filepath = 'data/test_data/FC_meat_AMR_metadata.csv'
amr_count_matrix_filepath = 'data/test_data/strict_SNP_confirmed_AMR_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/megares_annotations_v1.03.csv'

#################################
## Microbiome - 16S or kraken? ##
#################################

# Where is the metadata file for the microbiome samples stored on your machine?
microbiome_temp_metadata_file = "data/test_data/FC_meat_metadata.csv"

# If you used the AMR++ pipeline and have the kraken2 count matrix, point to the kraken file or see below for using qiime2 results.
#kraken_temp_file = "../kraken_analytic_matrix.csv"


## First, the 16S files. 
# These are the files you'll need to export from qiime2
#qiime tools export project-taxonomy.qza --output-dir exported-biom-table-taxa
#qiime tools export project-rep-seqs.qza --output-dir exported-rep-seqs
#qiime tools export project-aligned-masked-rooted.qza --output-dir exported-tree
#qiime tools export project-dada-table-filtered.qza --output-dir exported-biom-table
#then you need to convert the biom file to "json" using qiime1
#biom convert -i feature-table.biom -o otu_table_json.biom --table-type="OTU table" --to-json

##
## If you are using qiime2 results, uncomment the four lines below and specify the location to each file
##
biom_file <- "data/test_data/exported-biom-table/otu_table_json.biom"
tre_file <- "data/test_data/exported-tree/tree.nwk"
tax_fasta <- "data/test_data/exported-rep-seqs/dna-sequences.fasta" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otus.fasta
taxa_file <- "data/test_data/exported-biom-table-taxa/taxonomy.tsv" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otu_taxonomy.txt



###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# The following is a list of analyses based on variables in 
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
# NOTE: Exploratory variables cannot be numeric. 

#AMR_exploratory_analyses = list(
#  # Analysis 1
#  # Description: 
#  list(
#    name = 'Variable1',
#    subsets = list(),
#    exploratory_var = 'Variable1'
#  ),
#  # Analysis 2
#  # Description: 
#  list(
#    name = 'Variable2_Variable1_Subset',
#    subsets = list('Variable1 == Value1'),
#    exploratory_var = 'Variable2'
#  )
#)


AMR_exploratory_analyses = list(
  # Analysis Store
  # Description: 
  list(
    name = 'Store',
    subsets = list(),
    exploratory_var = 'Blinded_Store',
    order = ''
  ),  
  # Analysis Dilution
  # Description: 
  list(
    name = 'Dilution',
    subsets = list(),
    exploratory_var = 'Dilution',
    order = ''
  ),  
  # Analysis ID
  # Description: 
  list(
    name = 'ID',
    subsets = list(),
    exploratory_var = 'ID',
    order = ''
  ),  
  # Analysis 1
  # Description: 
  list(
    name = 'City',
    subsets = list(),
    exploratory_var = 'City',
    order = ''
  ),
  # Analysis 2
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    exploratory_var = 'Treatment',
    order = ''
  ),
  # Analysis 3
  # Description:
  list(
    name = 'Packaging',
    subsets = list(),
    exploratory_var = 'Packaging',
    order = ''
  ),
  # Analysis 3
  # Description:
  list(
    name = 'sample',
    subsets = list(),
    exploratory_var = 'sample',
    order = ''
  )
)



microbiome_exploratory_analyses = list(
  # Analysis Store
  # Description: 
  list(
    name = 'Store',
    subsets = list(),
    exploratory_var = 'Blinded_Store',
    order = ''
  ), 
  # Analysis ID
  # Description: 
  list(
    name = 'ID',
    subsets = list(),
    exploratory_var = 'ID',
    order = ''
  ),  
  # Analysis 1
  # Description: 
  list(
    name = 'City',
    subsets = list(),
    exploratory_var = 'City',
    order = ''
  ),
  # Analysis 2
  # Description:
  list(
    name = 'Treatment',
    subsets = list(),
    exploratory_var = 'Treatment',
    order = ''
  ),
  # Analysis 3
  # Description:
  list(
    name = 'Packaging',
    subsets = list(),
    exploratory_var = 'Packaging',
    order = ''
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.

#AMR_statistical_analyses = list(
#  # Analysis 1
#  # Description: 
#  list(
#    name = 'Variable1',
#    subsets = list(),
#    model_matrix = '~ 0 + Variable1 + Variable2',
#    contrasts = list('Variable1Value1 - Variable1Value2'),
#    random_effect = NA
#  ),
#  # Analysis 2
#  # Description: 
#  list(
#    name = 'Variable2_Variable1_Subset',
#    subsets = list('Variable1 == Value1'),
#    model_matrix = '~ 0 + Variable2',
#    contrasts = list('Variable2Value1 - Variable2Value2',
#                     'Variable2Value1 - Variable2Value3',
#                     'Variable2Value2 - Variable2Value3'),
#    random_effect = NA
#  )
#)

AMR_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Treatment_w_store',
    subsets = list(),
    model_matrix = '~ 0 + Treatment + Blinded_Store',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 3
  # Description: 
  list(
    name = 'Treatment_w_Dilution',
    subsets = list(),
    model_matrix = '~ 0 + Treatment + Dilution',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 4
  # Description: 
  list(
    name = 'Treatment_w_store_packaging',
    subsets = list(),
    model_matrix = '~ 0 + Treatment + Packaging + Blinded_Store',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 5
  # Description: 
  list(
    name = 'Dilution',
    subsets = list(),
    model_matrix = '~ 0 + Dilution',
    contrasts = list('DilutionNone - DilutionHalf'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Treatment',
    subsets = list(),
    model_matrix = '~ 0 + Treatment ',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Treatment_w_store',
    subsets = list(),
    model_matrix = '~ 0 + Treatment + Blinded_Store',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  ),
  # Analysis 3
  # Description: 
  list(
    name = 'Treatment_w_store_packaging',
    subsets = list(),
    model_matrix = '~ 0 + Treatment + Blinded_Store + Packaging',
    contrasts = list('TreatmentCONV - TreatmentRWA'),
    random_effect = NA
  )
)


## Run the analysis
#
## Pick the correct script that handles resistome data and/or microbiome data. 
# Choose one of the following scripts
#source('scripts/metagenomeSeq_megares_kraken.R')
#source('scripts/metagenomeSeq_megares_qiime.R')
source('scripts/metagenomeSeq_qiime2.R')
source('scripts/metagenomeSeq_megaresv2.R')

# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"

# Run code to make some exploratory figures, zero inflated gaussian model, and output count matrices.
## Choose one of the following scripts:
source('scripts/print_microbiome_figures.R')
source('scripts/print_AMR_figures.R')

## For ZIG model results
source('scripts/print_microbiome_ZIG_results.R')
source('scripts/print_AMR_ZIG_results.R')


