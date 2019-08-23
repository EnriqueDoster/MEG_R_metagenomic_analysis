## Start with this staging file to set up your analysis. 
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
setwd("")

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
amr_metadata_filepath = '../metadata.csv'
amr_count_matrix_filepath = '../AMR_analytic_matrix.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/megares_annotations_v1.02.csv'

# Where is the metadata file for the microbiome samples stored on your machine?
microbiome_temp_metadata_file = "../qiime2_metadata.csv"
#kraken_temp_file = '../AMR_analytic_matrix.csv'

#################################
## Microbiome - 16S or kraken? ##
#################################

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
biom_file <- "../16S/exported-biom-table/otu_table_json.biom"
tre_file <- "../16S/exported-tree/tree.nwk"
tax_fasta <- "../16S/exported-rep-seqs/dna-sequences.fasta" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otus.fasta
taxa_file <- "../16S/exported-biom-table-taxa/taxonomy.tsv" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otu_taxonomy.txt




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
AMR_exploratory_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Farm_type',
    subsets = list(),
    exploratory_var = 'farm_type'
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Sample',
    subsets = list(),
    exploratory_var = 'sample'
  )
)

microbiome_exploratory_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Farm_type',
    subsets = list(),
    exploratory_var = 'farm_type'
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Sample',
    subsets = list(),
    exploratory_var = 'sample'
  ),
  # Analysis 3
  # Description: 
  list(
    name = 'Milk',
    subsets = list(),
    exploratory_var = 'milk'
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
AMR_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'farm_type',
    subsets = list(),
    model_matrix = '~ 0 + farm_type',
    contrasts = list('farm_typeFree_stall-farm_typeTie_stall'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'farm_type',
    subsets = list(),
    model_matrix = '~ 0 + farm_type',
    contrasts = list('farm_typeFree_stall-farm_typeTie_stall'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'milk',
    subsets = list(),
    model_matrix = '~ 0 + milk',
    contrasts = list('milkpellet-milkcream'),
    random_effect = NA
  )
)
## Run the analysis


## Run the script to convert qiime2 results into "microbiome" objects to be used with amr_plus_plus
#
## Run the script that handles resistome data and microbiome data. 
source('scripts/metagenomeSeq_megares_kraken.R')
#source('scripts/HMM_metagenomeSeq_AMR.R')

# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"

## Run code to make some exploratory figures, zero inflated gaussian model, and output count matrices.
source('scripts/print_figures.R')


all_metadata <- cbind(microbiome_metadata,metadata)
View(metadata)
