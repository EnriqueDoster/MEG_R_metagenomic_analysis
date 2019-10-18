## Phyloseq for unifrac distance comparison
microbiome.ps <- data.ps
#taxa_names(phylum_microbiome) <- fData(AMR_analytic_data[[1]])
otu_table(microbiome.ps) <- otu_table(MRcounts(microbiome, norm=TRUE), taxa_are_rows = TRUE)
phy_tree(microbiome.ps) <- phy_tree(OTU_table)

filtered_microbiome.ps <- tax_glom(microbiome.ps, "phylum")
#filtered_microbiome.ps <- filter_taxa(filtered_microbiome.ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
filtered_microbiome.ps <- prune_samples(sample_sums(filtered_microbiome.ps)>=20, filtered_microbiome.ps)


##################### Only for testing for both #######################

### Testing various distance metrics
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods = dist_methods[-which(dist_methods=="ANY")]

# # These require tree
# dist_methods[(1:3)]
# # Remove them from the vector
# dist_methods <- dist_methods[-(1:3)]
# # This is the user-defined method:
# dist_methods["designdist"]
# # Remove the user-defined distance
# dist_methods = dist_methods[-which(dist_methods=="ANY")]

# Loop through each distance method, save each plot to a list, called plist.
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(filtered_microbiome.ps, method=i)
  # Calculate ordination
  iMDS  <- ordinate(filtered_microbiome.ps, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(filtered_microbiome.ps, iMDS, color="Group")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

# Combine results and shade according to a metadata variable
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color="Group", shape="filtered_microbiome.ps"))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p