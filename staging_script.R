## Start with this staging file to set up your analysis. 

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
setwd("")

# Set the output directory for graphs:
graph_output_dir = 'MEG_analysis/graphs'
# Set the output directory for statistics:
stats_output_dir = 'MEG_analysis/stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'


# AMR analysis 
# Where is the metadata file stored on your machine?
amr_metadata_filepath = 'MEG_analysis/proj6_mapping_file.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'MEG_analysis/megares_annotations_v1.01.csv'

# Load the data, MEGARes annotations, and metadata
amr <- newMRexperiment(read.table('amr_plus_plus/AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys
amr_temp_metadata <- read.csv(amr_metadata_filepath, header=T)
amr_temp_metadata[, sample_column_id] <- make.names(amr_temp_metadata[, sample_column_id])


## First, the 16S files. 
# These are the files you'll need to export from qiime2
#qiime tools export project-taxonomy.qza --output-dir exported-biom-table-taxa
#qiime tools export project-rep-seqs.qza --output-dir exported-rep-seqs
#qiime tools export project-aligned-masked-rooted.qza --output-dir exported-tree
#qiime tools export project-dada-table-filtered.qza --output-dir exported-biom-table
#then you need to convert the biom file to "json" using qiime1
#biom convert -i feature-table.biom -o otu_table_json.biom --table-type="OTU table" --to-json

# Now, specify file location
microbiome_temp_metadata <- "16S/project_mapping_file.tsv"

biom_file <- "16S/exported-biom-table/otu_table_json.biom"
tre_file <- "16S/exported-tree/tree.nwk"
tax_fasta <- "16S/exported-rep-seqs/dna-sequences.fasta" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otus.fasta
taxa_file <- "16S/exported-biom-table-taxa/taxonomy.tsv" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otu_taxonomy.txt


## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine. 

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.


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
    name = 'Variable1',
    subsets = list(),
    exploratory_var = 'Variable1'
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Variable2_Variable1_Subset',
    subsets = list('Variable1 == Value1'),
    exploratory_var = 'Variable2'
  )
)

microbiome_exploratory_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Variable1',
    subsets = list(),
    exploratory_var = 'Variable1'
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Variable2_Variable1_Subset',
    subsets = list('Variable1 == Value1'),
    exploratory_var = 'Variable2'
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
    name = 'Variable1',
    subsets = list(),
    model_matrix = '~ 0 + Variable1 + Variable2',
    contrasts = list('Variable1Value1 - Variable1Value2'),
    random_effect = NA
  ),
  
  # Analysis 2
  # Description: 
  list(
    name = 'Variable2_Variable1_Subset',
    subsets = list('Variable1 == Value1'),
    model_matrix = '~ 0 + Variable2',
    contrasts = list('Variable2Value1 - Variable2Value2',
                     'Variable2Value1 - Variable2Value3',
                     'Variable2Value2 - Variable2Value3'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Variable1',
    subsets = list(),
    model_matrix = '~ 0 + Variable1 + Variable2',
    contrasts = list('Variable1Value1 - Variable1Value2'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'Variable2_Variable1_Subset',
    subsets = list('Variable1 == Value1'),
    model_matrix = '~ 0 + Variable2',
    contrasts = list('Variable2Value1 - Variable2Value2',
                     'Variable2Value1 - Variable2Value3',
                     'Variable2Value2 - Variable2Value3'),
    random_effect = NA
  )
)

## Run the analysis

# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')
## Run the script to convert qiime2 results into "microbiome" objects to be used with amr_plus_plus
source('scripts/qiime2_2_phyloseq.R')
## Run the script that handles resistome data and microbiome data. 
source('scripts/metagenomeSeq_analytic_template_plus_qiime.R')
# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"

## Run code to make some exploratory figures, zero inflated gaussian model, and output count matrices.
source('scripts/print_figures.R')








