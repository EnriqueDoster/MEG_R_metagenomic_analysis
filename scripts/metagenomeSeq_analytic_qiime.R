## R script template for the production of basic qualitative
## statistics on the AMR and microbiome output from AMR++

## Author: Steven Lakin
## Edited to work with 16S data processed with Qiime2
## Modified by: Enrique Doster



####################
## Automated Code ##
####################
## Modify this as necessary, though you shouldn't need to for basic use.
set.seed(154)  # Seed the RNG, necessary for reproducibility

source('scripts/qiime2_2_phyloseq.R')
##########################
## Import & Format Data ##
##########################
## These files should be standard for all analyses, as they are
## the output matrices from AMR++ nextflow.  Additionally,
## you will need to obtain the most recent megares annotations file
## from megares.meglab.org


# If subdirs for stats and exploratory variables don't exist, create them
ifelse(!dir.exists(file.path(graph_output_dir)), dir.create(file.path(graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(file.path(stats_output_dir)), dir.create(file.path(stats_output_dir), mode='777'), FALSE)

for( dtype in c('Microbiome') ) {
  ifelse(!dir.exists(file.path(graph_output_dir, dtype)),
         dir.create(file.path(graph_output_dir, dtype), mode='777'), FALSE)
  
  for( v in 1:length(microbiome_exploratory_analyses) ) {
    ifelse(!dir.exists(file.path(graph_output_dir, dtype, microbiome_exploratory_analyses[[v]]$name)),
           dir.create(file.path(graph_output_dir, dtype, microbiome_exploratory_analyses[[v]]$name), mode='777'), FALSE)
  }
  
  ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
         dir.create(file.path(stats_output_dir, dtype), mode='777'), FALSE)
  
  for( a in 1:length(microbiome_statistical_analyses) ) {
    ifelse(!dir.exists(file.path(stats_output_dir, dtype, microbiome_statistical_analyses[[a]]$name)),
           dir.create(file.path(stats_output_dir, dtype, microbiome_statistical_analyses[[a]]$name), mode='777'), FALSE)
  }
}

ifelse(!dir.exists(file.path('microbiome_matrices')), dir.create(file.path('microbiome_matrices'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('microbiome_matrices/sparse_normalized')), dir.create(file.path('microbiome_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('microbiome_matrices/normalized')), dir.create(file.path('microbiome_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('microbiome_matrices/raw')), dir.create(file.path('microbiome_matrices/raw'), mode='777'), FALSE)


##### optional edits ###
## add step to remove samples with 1 or 0 features
#microbiome_sparseFeatures = which(colSums(MRcounts(microbiome) > 0) < 2) ## just counts how many rows have less than 2 hits
# microbiome = microbiome[, -microbiome_sparseFeatures]

# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(microbiome)

# Extract the normalized counts into data tables for aggregation
microbiome_norm <- data.table(MRcounts(microbiome, norm=T))
microbiome_raw <- data.table(MRcounts(microbiome, norm=F))

## Make objects for aggregating microbiome counts to different taxa levels
microbiome_norm[, id :=(rownames(microbiome)), ]
setkey(microbiome_norm, id)
microbiome_norm <- microbiome_taxonomy[microbiome_norm]  # left outer join

microbiome_raw[, id :=(rownames(microbiome)), ]
setkey(microbiome_raw, id)
microbiome_raw <- microbiome_taxonomy[microbiome_raw]  # left outer join

# Group the microbiome data by level for analysis, removing NA entries
microbiome_domain <- microbiome_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
microbiome_domain_analytic <- newMRexperiment(counts=microbiome_domain[, .SD, .SDcols=!'Domain'])
rownames(microbiome_domain_analytic) <- microbiome_domain$Domain

microbiome_domain_raw <- microbiome_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
microbiome_domain_raw_analytic <- newMRexperiment(counts=microbiome_domain_raw[, .SD, .SDcols=!'Domain'])
rownames(microbiome_domain_raw_analytic) <- microbiome_domain_raw$Domain

microbiome_phylum <- microbiome_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
microbiome_phylum_analytic <- newMRexperiment(counts=microbiome_phylum[, .SD, .SDcols=!'Phylum'])
rownames(microbiome_phylum_analytic) <- microbiome_phylum$Phylum

microbiome_phylum_raw <- microbiome_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
microbiome_phylum_raw_analytic <- newMRexperiment(counts=microbiome_phylum_raw[, .SD, .SDcols=!'Phylum'])
rownames(microbiome_phylum_raw_analytic) <- microbiome_phylum_raw$Phylum

microbiome_class <- microbiome_norm[!is.na(Class) & Class != 'NA' & Class != '', lapply(.SD, sum), by='Class', .SDcols=!1:8]
microbiome_class_analytic <- newMRexperiment(counts=microbiome_class[, .SD, .SDcols=!'Class'])
rownames(microbiome_class_analytic) <- microbiome_class$Class

microbiome_class_raw <- microbiome_raw[!is.na(Class) & Class != 'NA' & Class != '', lapply(.SD, sum), by='Class', .SDcols=!1:8]
microbiome_class_raw_analytic <- newMRexperiment(counts=microbiome_class_raw[, .SD, .SDcols=!'Class'])
rownames(microbiome_class_raw_analytic) <- microbiome_class_raw$Class

microbiome_order <- microbiome_norm[!is.na(Order) & Order != 'NA' & Order != '', lapply(.SD, sum), by='Order', .SDcols=!1:8]
microbiome_order_analytic <- newMRexperiment(counts=microbiome_order[, .SD, .SDcols=!'Order'])
rownames(microbiome_order_analytic) <- microbiome_order$Order

microbiome_order_raw <- microbiome_raw[!is.na(Order) & Order != 'NA' & Order != '', lapply(.SD, sum), by='Order', .SDcols=!1:8]
microbiome_order_raw_analytic <- newMRexperiment(counts=microbiome_order_raw[, .SD, .SDcols=!'Order'])
rownames(microbiome_order_raw_analytic) <- microbiome_order_raw$Order

microbiome_family <- microbiome_norm[!is.na(Family) & Family != 'NA' & Family != '', lapply(.SD, sum), by='Family', .SDcols=!1:8]
microbiome_family_analytic <- newMRexperiment(counts=microbiome_family[, .SD, .SDcols=!'Family'])
rownames(microbiome_family_analytic) <- microbiome_family$Family

microbiome_family_raw <- microbiome_raw[!is.na(Family) & Family != 'NA' & Family != '', lapply(.SD, sum), by='Family', .SDcols=!1:8]
microbiome_family_raw_analytic <- newMRexperiment(counts=microbiome_family_raw[, .SD, .SDcols=!'Family'])
rownames(microbiome_family_raw_analytic) <- microbiome_family_raw$Family

microbiome_genus <- microbiome_norm[!is.na(Genus) & Genus != 'NA' & Genus != '', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
microbiome_genus_analytic <- newMRexperiment(counts=microbiome_genus[, .SD, .SDcols=!'Genus'])
rownames(microbiome_genus_analytic) <- microbiome_genus$Genus

microbiome_genus_raw <- microbiome_raw[!is.na(Genus) & Genus != 'NA' & Genus != '', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
microbiome_genus_raw_analytic <- newMRexperiment(counts=microbiome_genus_raw[, .SD, .SDcols=!'Genus'])
rownames(microbiome_genus_raw_analytic) <- microbiome_genus_raw$Genus

microbiome_species <- microbiome_norm[!is.na(Species) & Species != 'NA' & Species != '', lapply(.SD, sum), by='Species', .SDcols=!1:8]
microbiome_species_analytic <- newMRexperiment(counts=microbiome_species[, .SD, .SDcols=!'Species'])
rownames(microbiome_species_analytic) <- microbiome_species$Species

microbiome_species_raw <- microbiome_raw[!is.na(Species) & Species != 'NA' & Species != '', lapply(.SD, sum), by='Species', .SDcols=!1:8]
microbiome_species_raw_analytic <- newMRexperiment(counts=microbiome_species_raw[, .SD, .SDcols=!'Species'])
rownames(microbiome_species_raw_analytic) <- microbiome_species_raw$Species


# Make long data frame for plotting with ggplot2
microbiome_melted_analytic <- rbind(melt_dt(MRcounts(microbiome_domain_analytic), 'Domain'),
                                melt_dt(MRcounts(microbiome_phylum_analytic), 'Phylum'),
                                melt_dt(MRcounts(microbiome_class_analytic), 'Class'),
                                melt_dt(MRcounts(microbiome_order_analytic), 'Order'),
                                melt_dt(MRcounts(microbiome_family_analytic), 'Family'),
                                melt_dt(MRcounts(microbiome_genus_analytic), 'Genus'),
                                melt_dt(MRcounts(microbiome_species_analytic), 'Species'))
microbiome_melted_raw_analytic <- rbind(melt_dt(MRcounts(microbiome_domain_raw_analytic), 'Domain'),
                                    melt_dt(MRcounts(microbiome_phylum_raw_analytic), 'Phylum'),
                                    melt_dt(MRcounts(microbiome_class_raw_analytic), 'Class'),
                                    melt_dt(MRcounts(microbiome_order_raw_analytic), 'Order'),
                                    melt_dt(MRcounts(microbiome_family_raw_analytic), 'Family'),
                                    melt_dt(MRcounts(microbiome_genus_raw_analytic), 'Genus'),
                                    melt_dt(MRcounts(microbiome_species_raw_analytic), 'Species'))

# Ensure that the metadata entries match the factor order of the MRexperiments
microbiome_metadata <- data.table(microbiome_temp_metadata[match(colnames(MRcounts(microbiome_phylum_analytic)), microbiome_temp_metadata[, sample_column_id]), ])
setkeyv(microbiome_metadata, sample_column_id)

microbiome_analytic_data <- c(microbiome_domain_analytic,
                          microbiome_phylum_analytic,
                          microbiome_class_analytic,
                          microbiome_order_analytic,
                          microbiome_family_analytic,
                          microbiome_genus_analytic,
                          microbiome_species_analytic)
microbiome_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
microbiome_raw_analytic_data <- c(microbiome_domain_raw_analytic,
                              microbiome_phylum_raw_analytic,
                              microbiome_class_raw_analytic,
                              microbiome_order_raw_analytic,
                              microbiome_family_raw_analytic,
                              microbiome_genus_raw_analytic,
                              microbiome_species_raw_analytic)
microbiome_raw_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')


## Add diversity values to metadata object

## Microbiome diversity
microbiome_metadata$microbiome_domain_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[1]])))
microbiome_metadata$microbiome_domain_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[1]])))

microbiome_metadata$microbiome_phylum_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[2]])))
microbiome_metadata$microbiome_phylum_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[2]])))

microbiome_metadata$microbiome_class_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[3]])))
microbiome_metadata$microbiome_class_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[3]])))

microbiome_metadata$microbiome_order_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[4]])))
microbiome_metadata$microbiome_order_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[4]])))

microbiome_metadata$microbiome_family_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[5]])))
microbiome_metadata$microbiome_family_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[5]])))

microbiome_metadata$microbiome_genus_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[6]])))
microbiome_metadata$microbiome_genus_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[6]])))

microbiome_metadata$microbiome_species_Richness = specnumber(t(MRcounts(microbiome_analytic_data[[7]])))
microbiome_metadata$microbiome_species_Shannon = diversity(t(MRcounts(microbiome_analytic_data[[7]])))


for( l in 1:length(microbiome_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(microbiome_analytic_data[[l]])), microbiome_metadata[[sample_column_id]])
    pData(microbiome_analytic_data[[l]]) <- data.frame(
        microbiome_metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(microbiome_analytic_data[[l]])) <- microbiome_metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(microbiome_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(microbiome_analytic_data[[l]])))
    rownames(fData(microbiome_analytic_data[[l]])) <- rownames(MRcounts(microbiome_analytic_data[[l]]))
}

for( l in 1:length(microbiome_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(microbiome_raw_analytic_data[[l]])), microbiome_metadata[[sample_column_id]])
    pData(microbiome_raw_analytic_data[[l]]) <- data.frame(
        microbiome_metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(microbiome_raw_analytic_data[[l]])) <- microbiome_metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(microbiome_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(microbiome_raw_analytic_data[[l]])))
    rownames(fData(microbiome_raw_analytic_data[[l]])) <- rownames(MRcounts(microbiome_raw_analytic_data[[l]]))
}


