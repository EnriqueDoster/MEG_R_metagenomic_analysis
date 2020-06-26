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

#source('scripts/qiime2_2_phyloseq.R')

# We usually filter out genes with wild-type potential.  If you want to include these
# in your analysis, comment this vector out
snp_regex = c("A16S",
"ACRR",
"AXYZ",
"BASRS",
"CAP16S",
"CATB",
"CATP",
"CDSA",
"CLS",
"DFRC",
"DHFR",
"DHFRIII",
"EATAV",
"EMBA",
"EMBB",
"EMBC",
"EMBR",
"ETHA",
"FABG",
"FOLC",
"FOLP",
"FUSA",
"FUSE",
"GIDB",
"GLPT",
"GSHF",
"GYRA",
"GYRB",
"GYRBA",
"GYRC",
"ILES",
"INHA",
"INIA",
"INIB",
"INIC",
"KASA",
"KATG",
"LIAFSR",
"LMRA",
"LPXA",
"LPXC",
"MDR23S",
"MENA",
"MEXS",
"MEXT",
"MEXZ",
"MLS23S",
"MPRFD",
"MTRR",
"MURA",
"MURG",
"NALC",
"NALD",
"NDH",
"NFXB",
"O23S",
"OMP36",
"OMPF",
"OMPFB",
"OPRD",
"P16S",
"P23S",
"PARC",
"PARE",
"PAREF",
"PBP1A",
"PBP2",
"PBP2B",
"PBP2X",
"PBP3",
"PGSA",
"PH23S",
"PHOB",
"PHOP",
"PHOQ",
"PNCA",
"POR",
"PTSL",
"RAMR",
"RIBD",
"RPOB",
"RPOBL",
"RPOC",
"RPOCL",
"RPSA",
"RPSJ",
"RPSL",
"RRSA",
"RRSC",
"RRSH",
"SOXR",
"SOXS",
"TET16S",
"TETR",
"TETRM",
"THYA",
"TLYA",
"TUFAB",
"UHPT",
"VAN",
"WALK",
"YKKCL",
"YYBT")

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

ifelse(!dir.exists(file.path('amr_matrices')), dir.create(file.path('amr_matrices'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('amr_matrices/sparse_normalized')), dir.create(file.path('amr_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/normalized')), dir.create(file.path('amr_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/raw')), dir.create(file.path('amr_matrices/raw'), mode='777'), FALSE)


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

amr_sparseFeatures = which(colSums(MRcounts(amr) > 0) < 2) ## just counts how many rows have less than 2 hits
if (length(amr_sparseFeatures) == 0) {
  print("No AMR features removed")
} else {
  amr = amr[, -amr_sparseFeatures]
  print(paste0(length(amr_sparseFeatures)," features removed. Listed below:"))
  print(names(amr_sparseFeatures))
}

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


# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(amr_temp_metadata[match(colnames(MRcounts(amr_class_analytic)), amr_temp_metadata[, sample_column_id]), ])
setkeyv(metadata, sample_column_id)

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


## Add diversity values to metadata object
## Resistome counts

new_AMR_metadata <- as.data.table(pData(AMR_raw_analytic_data[[1]]), keep.rownames = TRUE)
names(new_AMR_metadata)[names(new_AMR_metadata) == "rn"] <- "ID"

new_AMR_metadata$resistome_raw_mapped_reads = colSums(MRcounts(AMR_raw_analytic_data[[1]]))
new_AMR_metadata$resistome_CSS_counts = colSums(MRcounts(AMR_analytic_data[[1]]))

## Resistome diversity
new_AMR_metadata$AMR_class_Richness = specnumber(t(MRcounts(AMR_raw_analytic_data[[1]])))
new_AMR_metadata$AMR_class_Shannon = diversity(t(MRcounts(AMR_raw_analytic_data[[1]])))
#metadata$AMR_class_InvSimpson = diversity(t(MRcounts(AMR_raw_analytic_data[[1]])),"invsimp")

new_AMR_metadata$AMR_mech_Richness = specnumber(t(MRcounts(AMR_raw_analytic_data[[2]])))
new_AMR_metadata$AMR_mech_Shannon = diversity(t(MRcounts(AMR_raw_analytic_data[[2]])))

metadata <- new_AMR_metadata



