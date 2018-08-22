## Qiime2 export to phyloseq and metagenomeSeq
## Author: Enrique Doster
## Date: Spring 2018
## This code takes exported files from Qiime2 and converts them to a phyloseq object. Additionally, R objects are created to work with Dr. Zaid Abdo's analytic code. 

# Make OTU table piece
OTU_table <- import_biom(biom_file,tre_file,tax_fasta) # this has count table and phylogenetic tree

## Taxonomic information. The taxonomy.tsv file must be edited. Erase the "confidence" column". data.table to split up the taxonomic classification by ";" and rename each column by it's taxonomic level
taxa <- read.delim(taxa_file, sep = "\t", stringsAsFactors=FALSE, row.names=1)
row.names(taxa) <- paste(row.names(taxa),taxa[,1], sep= '; ')
taxa.dt <- as.data.table(taxa)
fixed_taxonomy <- data.table(id=rownames(taxa))
setDT(fixed_taxonomy)[, c('feature',
                          'kingdom',
                          'phylum',
                          'class',
                          'order',
                          'family',
                          'genus',
                          'species') := tstrsplit(id, '; ', type.convert = TRUE, fixed = TRUE)]
fixed_taxonomy[, id:= NULL]

taxa.df <- as.data.frame(fixed_taxonomy,row.names = fixed_taxonomy$feature)
row.names(taxa.df) <- fixed_taxonomy$feature
taxa.df <- taxa.df[,-1]
taxa <- as.matrix(taxa.df)
TAXA <- tax_table(taxa)

# Now, it's sample metadata
sampleDF = read.delim(microbiome_temp_metadata, sep = "\t", stringsAsFactors=FALSE, row.names=1)
all(rownames(sampleDF) %in% sample_names(OTU_table)) # Check that the rownames match the sample names
METADATA = sample_data(sampleDF)# Convert to "sample_data" class in phyloseq

# Now merge into final phyloseq object
data.ps.original = merge_phyloseq(OTU_table,TAXA, METADATA,refseqFunction = parse_taxonomy_greengenes) #refseqFunction = parse_taxonomy_greengenes
taxa.df <- as.data.frame(as(tax_table(data.ps.original),"matrix"))
data.df <- t(as.data.frame(as(otu_table(data.ps.original),"matrix")))
meta.df <- as.data.frame(as(sample_data(data.ps.original),"matrix"))

data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))

## At this point we have all the files we need:
# phyloseq object: data.ps
# 3 dataframe objects
# taxa.df 
# data.df 
# meta.df 

## Now, need to make the taxonomy object set up for use with amr_plus_plus. 
data.df <- t(data.df)
microbiome <- newMRexperiment(data.df)

## Domain is incorrect, just used to work with amrplusplus
microbiome_taxonomy <- data.table(id= rownames(microbiome))

microbiome_taxonomy$Domain <- gsub("[a-z]__$","unclassified", taxa.df$kingdom)
microbiome_taxonomy$Phylum <- gsub("[a-z]__$","unclassified", taxa.df$phylum)
microbiome_taxonomy$Class <- gsub("[a-z]__$","unclassified", taxa.df$class)
microbiome_taxonomy$Order <- gsub("[a-z]__$","unclassified", taxa.df$order)
microbiome_taxonomy$Family <- gsub("[a-z]__$","unclassified", taxa.df$family)
microbiome_taxonomy$Genus <- gsub("[a-z]__$","unclassified", taxa.df$genus)
microbiome_taxonomy$Species <- gsub("[a-z]__$","unclassified", taxa.df$species)

microbiome_taxonomy$Domain <- gsub("[a-z]__","", taxa.df$kingdom)
microbiome_taxonomy$Phylum <- gsub("[a-z]__","", taxa.df$phylum)
microbiome_taxonomy$Class <- gsub("[a-z]__","", taxa.df$class)
microbiome_taxonomy$Order <- gsub("[a-z]__","", taxa.df$order)
microbiome_taxonomy$Family <- gsub("[a-z]__","", taxa.df$family)
microbiome_taxonomy$Genus <- gsub("[a-z]__","", taxa.df$genus)
microbiome_taxonomy$Species <- gsub("[a-z]__","", taxa.df$species)

for( v in 1:nrow(microbiome_taxonomy) ) {
  if (!is.na(microbiome_taxonomy[v]$Species) && microbiome_taxonomy[v]$Species !="" ) {
    microbiome_taxonomy[v]$Species <- paste(microbiome_taxonomy[v]$Genus,microbiome_taxonomy[v]$Species, sep=" ")
  }
}

setkey(microbiome_taxonomy, id)


