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


# We usually filter out genes with wild-type potential.  If you want to include these
# in your analysis, comment this vector out
snp_regex = c('ACRR',
              'CATB',
              'CLS',
              'DFRC',
              'DHFR',
              'DHFRIII',
              'DHFRIX',
              'EMBA',
              'embB',
              'EMBB',
              'EMBC',
              'EMBR',
              'ETHA',
              'FOLP',
              'GIDB',
              'GYRA',
              'gyrB',
              'GYRB',
              'INHA',
              'INIA',
              'INIC',
              'KASA',
              'LIAFSR',
              'LMRA',
              'MARR',
              'MEXR',
              'MEXZ',
              'mprF',
              'MPRF',
              'NDH',
              'omp36',
              'OMP36',
              'OMPF',
              'OPRD',
              'PARC',
              'parE',
              'PARE',
              'PGSA',
              'phoP',
              'PHOP',
              'PNCA',
              'POR',
              'PORB',
              'RAMR',
              'rpoB',
              'RPOB',
              'RPOC',
              'RPSL',
              'SOXS',
              'tetR',
              'TETR',
              'TLYA',
              'TUFAB')


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

for( dtype in c('AMR') ) {
    ifelse(!dir.exists(file.path(graph_output_dir, dtype)),
           dir.create(file.path(graph_output_dir, dtype), mode='777'), FALSE)
    
    for( v in 1:length(AMR_exploratory_analyses) ) {
        ifelse(!dir.exists(file.path(graph_output_dir, dtype, AMR_exploratory_analyses[[v]]$name)),
               dir.create(file.path(graph_output_dir, dtype, AMR_exploratory_analyses[[v]]$name), mode='777'), FALSE)
    }
    
    ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
           dir.create(file.path(stats_output_dir, dtype), mode='777'), FALSE)
    
    for( a in 1:length(AMR_statistical_analyses) ) {
        ifelse(!dir.exists(file.path(stats_output_dir, dtype, AMR_statistical_analyses[[a]]$name)),
               dir.create(file.path(stats_output_dir, dtype, AMR_statistical_analyses[[a]]$name), mode='777'), FALSE)
    }
}

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


ifelse(!dir.exists(file.path('amr_matrices')), dir.create(file.path('amr_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('microbiome_matrices')), dir.create(file.path('microbiome_matrices'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('amr_matrices/sparse_normalized')), dir.create(file.path('amr_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/normalized')), dir.create(file.path('amr_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/raw')), dir.create(file.path('amr_matrices/raw'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('microbiome_matrices/sparse_normalized')), dir.create(file.path('microbiome_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('microbiome_matrices/normalized')), dir.create(file.path('microbiome_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('microbiome_matrices/raw')), dir.create(file.path('microbiome_matrices/raw'), mode='777'), FALSE)



####################
##       AMR      ##
####################
# Load the data, MEGARes annotations, and metadata
amr <- newMRexperiment(read.table(amr_count_matrix_filepath, header=T, row.names=1, sep=','))
amr <- newMRexperiment(round(MRcounts(amr),0))
amr_temp_metadata <- read.csv(amr_metadata_filepath, header=T)
amr_temp_metadata[, sample_column_id] <- make.names(amr_temp_metadata[, sample_column_id])
# Annotation for regular AMR++ analysis
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys



##### optional edits ###
## add step to remove samples with 1 or 0 features
#amr_sparseFeatures = which(colSums(MRcounts(amr) > 0) < 2) ## just counts how many rows have less than 2 hits
#amr = amr[, -amr_sparseFeatures]
# microbiome_sparseFeatures = which(colSums(MRcounts(microbiome) > 0) < 2) ## just counts how many rows have less than 2 hits
# microbiome = microbiome[, -microbiome_sparseFeatures]

# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(amr)

# Extract the normalized counts into data tables for aggregation
amr_norm <- data.table(MRcounts(amr, norm=T))
amr_raw <- data.table(MRcounts(amr, norm=F))

# Aggregate the normalized counts for AMR using the annotations data table, SQL
# outer join, and aggregation with vectorized lapply
amr_norm[, header :=( rownames(amr) ), ]
setkey(amr_norm, header)
amr_norm <- annotations[amr_norm]  # left outer join

amr_raw[, header :=( rownames(amr) ), ]
setkey(amr_raw, header)
amr_raw <- annotations[amr_raw]  # left outer join

# subset groups that correspond to potentially wild-type genes
amr_snp_confirm <- amr_raw[group %in% snp_regex, ]

# Remove groups that correspond to potentially wild-type genes
amr_raw <- amr_raw[!(group %in% snp_regex), ]
amr_norm<- amr_norm[!(group %in% snp_regex), ]

# Group the AMR data by level for analysis
amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
rownames(amr_class_analytic) <- amr_class$class

amr_class_raw <- amr_raw[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_raw_analytic <- newMRexperiment(counts=amr_class_raw[, .SD, .SDcols=!'class'])
rownames(amr_class_raw_analytic) <- amr_class_raw$class

amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_analytic) <- amr_mech$mechanism

amr_mech_raw <- amr_raw[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_raw_analytic <- newMRexperiment(counts=amr_mech_raw[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_raw_analytic) <- amr_mech_raw$mechanism

amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
rownames(amr_group_analytic) <- amr_group$group

amr_group_raw <- amr_raw[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_raw_analytic <- newMRexperiment(counts=amr_group_raw[, .SD, .SDcols=!'group'])
rownames(amr_group_raw_analytic) <- amr_group_raw$group

amr_gene_analytic <- newMRexperiment(
    counts=amr_norm[!(group %in% snp_regex),
                    .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
amr_gene_raw_analytic <- newMRexperiment(
    counts=amr_raw[!(group %in% snp_regex),
                   .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])

rownames(amr_gene_analytic) <- amr_norm$header
rownames(amr_gene_raw_analytic) <- amr_raw$header


# Make long data frame for plotting with ggplot2
amr_melted_analytic <- rbind(melt_dt(MRcounts(amr_class_analytic), 'Class'),
                             melt_dt(MRcounts(amr_mech_analytic), 'Mechanism'),
                             melt_dt(MRcounts(amr_group_analytic), 'Group'),
                             melt_dt(MRcounts(amr_gene_analytic), 'Gene'))
amr_melted_raw_analytic <- rbind(melt_dt(MRcounts(amr_class_raw_analytic), 'Class'),
                                 melt_dt(MRcounts(amr_mech_raw_analytic), 'Mechanism'),
                                 melt_dt(MRcounts(amr_group_raw_analytic), 'Group'),
                                 melt_dt(MRcounts(amr_gene_raw_analytic), 'Gene'))



####################
##   Microbiome   ##
####################
# Aggregate the kraken data using the rownames:
# this set of commands splits the rownames into their taxonomic levels and
# fills empty values with NA.  We then join that taxonomy data table with
# the actual data and aggregate using lapply as before.
temp_microbiome <- read.table(kraken_temp_file, header=T, row.names=1, sep=',')
microbiome <- newMRexperiment(temp_microbiome[rowSums(temp_microbiome) > 0, ])
# Load microbiome metadata
microbiome_temp_metadata = read.delim(microbiome_temp_metadata_file, sep = ",", stringsAsFactors=FALSE, row.names=NULL)

microbiome_taxonomy <- data.table(id=rownames(microbiome))
setDT(microbiome_taxonomy)[, c('Domain',
                               'Kingdom',
                               'Phylum',
                               'Class',
                               'Order',
                               'Family',
                               'Genus',
                               'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
setkey(microbiome_taxonomy, id)


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
microbiome_domain <- microbiome_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:9]
microbiome_domain_analytic <- newMRexperiment(counts=microbiome_domain[, .SD, .SDcols=!'Domain'])
rownames(microbiome_domain_analytic) <- microbiome_domain$Domain

microbiome_domain_raw <- microbiome_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:9]
microbiome_domain_raw_analytic <- newMRexperiment(counts=microbiome_domain_raw[, .SD, .SDcols=!'Domain'])
rownames(microbiome_domain_raw_analytic) <- microbiome_domain_raw$Domain

microbiome_phylum <- microbiome_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:9]
microbiome_phylum_analytic <- newMRexperiment(counts=microbiome_phylum[, .SD, .SDcols=!'Phylum'])
rownames(microbiome_phylum_analytic) <- microbiome_phylum$Phylum

microbiome_phylum_raw <- microbiome_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:9]
microbiome_phylum_raw_analytic <- newMRexperiment(counts=microbiome_phylum_raw[, .SD, .SDcols=!'Phylum'])
rownames(microbiome_phylum_raw_analytic) <- microbiome_phylum_raw$Phylum

microbiome_class <- microbiome_norm[!is.na(Class) & Class != 'NA' & Class != '', lapply(.SD, sum), by='Class', .SDcols=!1:9]
microbiome_class_analytic <- newMRexperiment(counts=microbiome_class[, .SD, .SDcols=!'Class'])
rownames(microbiome_class_analytic) <- microbiome_class$Class

microbiome_class_raw <- microbiome_raw[!is.na(Class) & Class != 'NA' & Class != '', lapply(.SD, sum), by='Class', .SDcols=!1:9]
microbiome_class_raw_analytic <- newMRexperiment(counts=microbiome_class_raw[, .SD, .SDcols=!'Class'])
rownames(microbiome_class_raw_analytic) <- microbiome_class_raw$Class

microbiome_order <- microbiome_norm[!is.na(Order) & Order != 'NA' & Order != '', lapply(.SD, sum), by='Order', .SDcols=!1:9]
microbiome_order_analytic <- newMRexperiment(counts=microbiome_order[, .SD, .SDcols=!'Order'])
rownames(microbiome_order_analytic) <- microbiome_order$Order

microbiome_order_raw <- microbiome_raw[!is.na(Order) & Order != 'NA' & Order != '', lapply(.SD, sum), by='Order', .SDcols=!1:9]
microbiome_order_raw_analytic <- newMRexperiment(counts=microbiome_order_raw[, .SD, .SDcols=!'Order'])
rownames(microbiome_order_raw_analytic) <- microbiome_order_raw$Order

microbiome_family <- microbiome_norm[!is.na(Family) & Family != 'NA' & Family != '', lapply(.SD, sum), by='Family', .SDcols=!1:9]
microbiome_family_analytic <- newMRexperiment(counts=microbiome_family[, .SD, .SDcols=!'Family'])
rownames(microbiome_family_analytic) <- microbiome_family$Family

microbiome_family_raw <- microbiome_raw[!is.na(Family) & Family != 'NA' & Family != '', lapply(.SD, sum), by='Family', .SDcols=!1:9]
microbiome_family_raw_analytic <- newMRexperiment(counts=microbiome_family_raw[, .SD, .SDcols=!'Family'])
rownames(microbiome_family_raw_analytic) <- microbiome_family_raw$Family

microbiome_genus <- microbiome_norm[!is.na(Genus) & Genus != 'NA' & Genus != '', lapply(.SD, sum), by='Genus', .SDcols=!1:9]
microbiome_genus_analytic <- newMRexperiment(counts=microbiome_genus[, .SD, .SDcols=!'Genus'])
rownames(microbiome_genus_analytic) <- microbiome_genus$Genus

microbiome_genus_raw <- microbiome_raw[!is.na(Genus) & Genus != 'NA' & Genus != '', lapply(.SD, sum), by='Genus', .SDcols=!1:9]
microbiome_genus_raw_analytic <- newMRexperiment(counts=microbiome_genus_raw[, .SD, .SDcols=!'Genus'])
rownames(microbiome_genus_raw_analytic) <- microbiome_genus_raw$Genus

microbiome_species <- microbiome_norm[!is.na(Species) & Species != 'NA' & Species != '', lapply(.SD, sum), by='Species', .SDcols=!1:9]
microbiome_species_analytic <- newMRexperiment(counts=microbiome_species[, .SD, .SDcols=!'Species'])
rownames(microbiome_species_analytic) <- microbiome_species$Species

microbiome_species_raw <- microbiome_raw[!is.na(Species) & Species != 'NA' & Species != '', lapply(.SD, sum), by='Species', .SDcols=!1:9]
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
metadata <- data.table(amr_temp_metadata[match(colnames(MRcounts(amr_class_analytic)), amr_temp_metadata[, sample_column_id]), ])
setkeyv(metadata, sample_column_id)

microbiome_metadata <- data.table(microbiome_temp_metadata[match(colnames(MRcounts(microbiome_phylum_analytic)), microbiome_temp_metadata[, sample_column_id]), ])
setkeyv(microbiome_metadata, sample_column_id)

# Vector of objects for iteration and their names
AMR_analytic_data <- c(amr_class_analytic,
                       amr_mech_analytic,
                       amr_group_analytic,
                       amr_gene_analytic)
AMR_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
AMR_raw_analytic_data <- c(amr_class_raw_analytic,
                           amr_mech_raw_analytic,
                           amr_group_raw_analytic,
                           amr_gene_raw_analytic)
AMR_raw_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
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
## Resistome diversity
metadata$AMR_class_Richness = specnumber(t(MRcounts(AMR_raw_analytic_data[[1]])))
metadata$AMR_class_Shannon = diversity(t(MRcounts(AMR_raw_analytic_data[[1]])))

metadata$AMR_mech_Richness = specnumber(t(MRcounts(AMR_raw_analytic_data[[2]])))
metadata$AMR_mech_Shannon = diversity(t(MRcounts(AMR_raw_analytic_data[[2]])))

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

for( l in 1:length(AMR_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_analytic_data[[l]])))
    rownames(fData(AMR_analytic_data[[l]])) <- rownames(MRcounts(AMR_analytic_data[[l]]))
}

for( l in 1:length(AMR_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_raw_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_raw_analytic_data[[l]])))
    rownames(fData(AMR_raw_analytic_data[[l]])) <- rownames(MRcounts(AMR_raw_analytic_data[[l]]))
}


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


