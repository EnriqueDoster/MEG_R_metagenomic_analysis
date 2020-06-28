##### Set up environment
library("phyloseq")
library("dplyr")
library("ggplot2")
library("data.table")
library("tidyr")
library("forcats")
library("vegan")

# Load sample metadata
resistome_sample_metadata <-  read.delim(amr_metadata_filepath, sep = ",", stringsAsFactors=FALSE, row.names=1)
microbiome_sample_metadata <-  read.delim(microbiome_temp_metadata_file, sep = ",", stringsAsFactors=FALSE, row.names=1)
#
##
### Load qiime2 microbiome results (16S reads)
##
#

# # We import the .biom as in Lesson 1 step 4
# microbiome <- import_biom("data/Exported_16S_qiime2_results/Lesson2_16S_asv_table.biom")
# # Load the taxonomy file from qiime2
# taxa <- read.table("data/Exported_16S_qiime2_results/Lesson2_taxonomy.tsv", header=T, row.names=1, sep='\t', quote = "")
# row.names(taxa) <- paste(row.names(taxa),taxa[,1], sep= '; ')
# taxa.dt <- data.table(id=rownames(taxa)) # we'll make a column with the name "id"
# taxa.dt[, c('feature',
#             'kingdom',
#             'phylum',
#             'class',
#             'order',
#             'family',
#             'genus',
#             'species') := tstrsplit(id, '; ', type.convert = TRUE, fixed = TRUE)]
# 
# taxa.df <- as.data.frame(taxa.dt)
# taxa.df <- within(taxa.df, rm(id))
# row.names(taxa.df) <- taxa.df$feature
# taxa.df <- within(taxa.df, rm(feature))
# 
# microbiome_phylo_tree <- read_tree("./data/Exported_16S_qiime2_results/Lesson2_tree.nwk")
# 
# microbiome.ps <- merge_phyloseq(microbiome, phy_tree(microbiome_phylo_tree), tax_table(as.matrix(taxa.df)), sample_data(sample_metadata))
# 
# # Estimating richness and diversity using the easy-to-use function estimate_richness()
# microbiome_16S_diversity_values <- estimate_richness(microbiome.ps)
# 

#
##
### Loading the kraken2 microbiome results (shotgun reads)
##
#

# Load the kraken count table
kraken_microbiome_counts <- read.table(kraken_temp_file, header=T, row.names=1, sep=',')


# Convert to format that phyloseq likes with otu_table()                                      
kraken_microbiome_otu <- otu_table(kraken_microbiome_counts, taxa_are_rows = TRUE)

# Repeat similar steps to what we did with the qiime2 taxonomy
kraken_taxonomy <- data.table(id=rownames(kraken_microbiome_otu))
kraken_taxonomy[, c('domain',
                    'kingdom',
                    'phylum',
                    'class',
                    'order',
                    'family',
                    'genus',
                    'species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]

# Conver to data.frame
kraken_taxonomy <- as.data.frame(kraken_taxonomy)
# Use the id variable to rename the row.names
row.names(kraken_taxonomy) <- kraken_taxonomy$id
# Remove the "id" column
kraken_taxonomy <- within(kraken_taxonomy, rm(id))

# Create kraken phyloseq object
kraken_microbiome.ps <- merge_phyloseq(kraken_microbiome_otu, tax_table(as.matrix(kraken_taxonomy)), sample_data(microbiome_sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
microbiome_shotgun_diversity_values <- estimate_richness(kraken_microbiome.ps)


#
##
### Loading the resistome results (shotgun reads)
##
#

# Load MEGARes counts                                         
amr_counts <- read.table(amr_count_matrix_filepath, header=T, row.names=1, sep=',')

# We can convert our amr count object to the otu_table format required for phyloseq
amr_otu <- otu_table(amr_counts, taxa_are_rows = TRUE)

annotations <- read.csv(megares_annotation_filename, header=T, row.names=1)

# We can now merge these objects to make a phyloseq object
amr.ps <- merge_phyloseq(amr_otu, tax_table(as.matrix(annotations)), sample_data(resistome_sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
#amr_shotgun_diversity_values <- estimate_richness(amr.ps)


#
##
### Merge diversity values for the microbiome and resistome
##
#


# We can make a new column, named "SeqType" and give it the value of "16S"
# microbiome_16S_diversity_values <- microbiome_16S_diversity_values %>%
#   mutate(SeqType = "16S", DataType = "16S microbiome", Sample = row.names(microbiome_16S_diversity_values))

# # We can make a new column, named "SeqType" and give it the value of "shotgun"
# microbiome_shotgun_diversity_values <- microbiome_shotgun_diversity_values %>%
#   mutate(SeqType = "shotgun", DataType = "shotgun microbiome", Sample = row.names(microbiome_shotgun_diversity_values))
# 
# # We can make a new column, named "SeqType" and give it the value of "shotgun"
# amr_shotgun_diversity_values <- amr_shotgun_diversity_values %>%
#   mutate(SeqType = "shotgun", DataType = "resistome", Sample = row.names(amr_shotgun_diversity_values))
# 
# 
# # Now we can merge these tables based on identical row
# combined_diversity_values <- bind_rows(microbiome_shotgun_diversity_values, amr_shotgun_diversity_values)
# 
# # To help us better summarize the results, we can add the metadata information 
# # First, we need to add a "Sample" column in the sample_metadata like in the combined_diversity_values object
# sample_metadata$Sample <- row.names(sample_metadata)
# 
# # Now, we use left_join() to add the sample_metadata to the combined_diversity_values object
# combined_diversity_values <- left_join(combined_diversity_values, sample_metadata, by = "Sample")


# Clean up the environment
rm(amr_counts)
rm(amr_otu)
rm(annotations)
#rm(amr_shotgun_diversity_values)
rm(kraken_microbiome_counts)
rm(kraken_microbiome_otu)
rm(kraken_taxonomy)
#rm(microbiome_shotgun_diversity_values)
