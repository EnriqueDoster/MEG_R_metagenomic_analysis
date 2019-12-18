## Adonis script
##
###
####
#####
## AMR with all samples
#####
####
###
##


#AMR_all_phylum <- AMR_raw_analytic_data[[1]]
##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
#p = cumNormStatFast(AMR_all_phylum)
### We adjust our container to contain the normalization pecentile using cumNorm()
#AMR_all_phylum.ms = cumNorm(AMR_all_phylum,p)
### We find the normalization factors per sample using normFactors()
#nf = normFactors(AMR_all_phylum.ms)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(AMR_analytic_data[[1]])
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)


meta.df = as.data.table(pData( AMR_raw_analytic_data[[1]])[,c("Project","Total_ADD","Feed_ADD","Parenteral_ADD","DOF","total_tetracycline_ADD","total_MLS_ADD",
                                                              "feed_MLS_ADD","feed_tetracycline_ADD","parenteral_tetracycline_ADD","parenteral_MLS_ADD",
                                                              "parenteral_phenicol_ADD","parenteral_sulfonamide_ADD")])

proj3_metadata <- meta.df[Project == "Project_3"]
proj4_metadata <- meta.df[Project == "Project_4"]

proj3_metadata[,("Project") := NULL]
proj4_metadata[,("Project") := NULL]

scaled_proj3_metadata <- as.data.frame(scale(proj3_metadata))
scaled_proj4_metadata <- as.data.frame(scale(proj4_metadata))

# error during scaling, parenteral_betalactams_ADD, parenteral_fluoroquinolones_ADD
hist(AMR_raw_analytic_data[[1]]$DDD,breaks = 20)

# check that we get mean of 0 and sd of 1
colMeans(scaled_proj4_metadata)
apply(scaled_proj3_metadata, 2, sd)


## Add time and feedlot categories
meta.df$Time <-  pData(AMR_raw_analytic_data[[1]])$Time
meta.df$Feedlot_categories <-  pData(AMR_raw_analytic_data[[1]])$Feedlot_categories
#meta.df$Group <-  pData(AMR_raw_analytic_data[[1]])$Group
## 
combined_scaled_metadata <- rbind(scaled_proj3_metadata, scaled_proj4_metadata)
combined_scaled_metadata$Time <- meta.df$Time
combined_scaled_metadata$Feedlot_categories <- meta.df$Feedlot_categories
#combined_scaled_metadata$Group <- meta.df$Group
##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ ., data =  combined_scaled_metadata)
anova(mod1, by = "term", perm = 1000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data =  combined_scaled_metadata)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm", direction="both")
anova(mod, by = "term", perm = 1000)
anova(mod,perm = 1000)
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <-ordiplot(mod,type="n")
colvec <- c("black","grey51","blue","red","green") 
#pchvec <-c(10,2)
points(mod, "sites", pch=21, col=as.factor(combined_scaled_metadata$Time), bg=as.factor(combined_scaled_metadata$Time))
groupz <- sort(unique(combined_scaled_metadata$Time))
for(i in seq(groupz)) {ordispider(mod,  combined_scaled_metadata$Time,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legv <- sort(unique(combined_scaled_metadata$Time))
#legend("right", legend = levels(combined_scaled_metadata$Group),col=legv, title = "Sample groups",bty = "n", pt.bg = legv)
#legend("topright", legend = levels(TreatmentArrival),col=colvec[legv], pch=pchvec[legv], title = "Arrival: Mechanism",bty = "n", pt.bg = colvec[legv])
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red","green"), cex=2, pos = 4)




# Apriori model
mod_apriori <- rda(hell.norm.mt ~ DOF +  Feed_ADD  + Feedlot_categories , data =  combined_scaled_metadata)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <-ordiplot(mod_apriori,type="n")
colvec <- c("black","grey51","blue","red","green") 
#pchvec <-c(10,2)
points(mod, "sites", pch=21, col=as.factor(combined_scaled_metadata$Time), bg=as.factor(combined_scaled_metadata$Time))
groupz <- sort(unique(combined_scaled_metadata$Time))
for(i in seq(groupz)) {ordispider(mod,  combined_scaled_metadata$Time,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legv <- sort(unique(combined_scaled_metadata$Time))
#legend("right", legend = levels(combined_scaled_metadata$Group),col=legv, title = "Sample groups",bty = "n", pt.bg = legv)
#legend("topright", legend = levels(TreatmentArrival),col=colvec[legv], pch=pchvec[legv], title = "Arrival: Mechanism",bty = "n", pt.bg = colvec[legv])
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red"), cex=2, pos = 4)



#text(fig, "species", col="red", cex=0.9)



##
###
####
#####
## AMR only rehandling samples
#####
####
###
##

rehandling_samp_num <- which(pData(AMR_analytic_data[[1]])$Time != "Arrival")
rehandling_class_ms <- AMR_analytic_data[[1]][,rehandling_samp_num]

#AMR_all_phylum <- AMR_raw_analytic_data[[1]]
##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
#p = cumNormStatFast(AMR_all_phylum)
### We adjust our container to contain the normalization pecentile using cumNorm()
#AMR_all_phylum.ms = cumNorm(AMR_all_phylum,p)
### We find the normalization factors per sample using normFactors()
#nf = normFactors(AMR_all_phylum.ms)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(rehandling_class_ms)
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)


meta.df = as.data.table(pData( rehandling_class_ms)[,c("Project","Total_ADD","Feed_ADD","Parenteral_ADD","DOF","total_tetracycline_ADD","total_MLS_ADD",
                                                "feed_MLS_ADD","feed_tetracycline_ADD","parenteral_tetracycline_ADD","parenteral_MLS_ADD",
                                                "parenteral_phenicol_ADD","parenteral_sulfonamide_ADD")])

proj3_metadata <- meta.df[Project == "Project_3"]
proj4_metadata <- meta.df[Project == "Project_4"]

proj3_metadata[,("Project") := NULL]
proj4_metadata[,("Project") := NULL]

scaled_proj3_metadata <- as.data.frame(scale(proj3_metadata))
scaled_proj4_metadata <- as.data.frame(scale(proj4_metadata))

# error during scaling, parenteral_betalactams_ADD, parenteral_fluoroquinolones_ADD
#hist(AMR_raw_analytic_data[[1]]$DDD,breaks = 20)

# check that we get mean of 0 and sd of 1
colMeans(scaled_proj4_metadata)
apply(scaled_proj3_metadata, 2, sd)


## Add time and feedlot categories
meta.df$Time <-  pData(rehandling_class_ms)$Time
meta.df$Feedlot_categories <-  pData(rehandling_class_ms)$Feedlot_categories
meta.df$Sample_type <-  pData(rehandling_class_ms)$Sample_type
meta.df$Group <-  pData(rehandling_class_ms)$Group
## 
combined_scaled_metadata <- rbind(scaled_proj3_metadata, scaled_proj4_metadata)
combined_scaled_metadata$Time <- meta.df$Time
combined_scaled_metadata$Feedlot_categories <- meta.df$Feedlot_categories
combined_scaled_metadata$Sample_type <- meta.df$Sample_type
combined_scaled_metadata$Group <- meta.df$Group
##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ ., data =  combined_scaled_metadata)
anova(mod1, by = "term", perm = 1000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data =  combined_scaled_metadata)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm", direction="both")
anova(mod, by = "term", perm = 1000)
anova(mod,perm = 1000)
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <- ordiplot(mod,type="n")
points(mod, "sites", pch=21, col=as.factor(combined_scaled_metadata$Group), bg=as.factor(combined_scaled_metadata$Group))
text(fig, "species", col="red", cex=0.9)
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red"), cex=2, pos = 4)




# Apriori model
mod_apriori <- rda(hell.norm.mt ~ Time + Feedlot_categories , data =  combined_scaled_metadata)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <- ordiplot(mod_apriori,type="n")
points(mod_apriori, "sites", pch=21, col=as.factor(combined_scaled_metadata$Feedlot_categories), bg=as.factor(combined_scaled_metadata$Feedlot_categories))
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red"), cex=2.5, labels = c("Arrival","Rehandling"),adj = c(-.3,-3))


#text(fig, "species", col="red", cex=0.9)










##
###
####
#####
## AMR split by Time of collection
#####
####
###
##

AMR_arrival_samples = which(pData(AMR_analytic_data[[1]])$Time == "Arrival")
AMR_exit_samples = which(pData(AMR_analytic_data[[1]])$Time == "Exit")

pData(AMR_analytic_data[[1]])$Days_since_tx <- as.integer(pData(AMR_analytic_data[[1]])$Days_since_tx)

AMR_class_arrival <- AMR_analytic_data[[1]][, AMR_arrival_samples]
AMR_mech_arrival <- AMR_analytic_data[[2]][, AMR_arrival_samples]
AMR_class_exit <- AMR_analytic_data[[1]][, AMR_exit_samples]
AMR_mech_exit <- AMR_analytic_data[[2]][, AMR_exit_samples]


##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
p = cumNormStatFast(AMR_class_exit)
### We adjust our container to contain the normalization pecentile using cumNorm()
AMR_class_exit = cumNorm(AMR_class_exit,p)
### We find the normalization factors per sample using normFactors()
nf = normFactors(AMR_class_exit)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(AMR_class_exit)
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)

meta.df = pData(AMR_class_exit)[,c("Total_ADD","Feed_ADD","Parenteral_ADD","DOF","num_tx","Days_since_tx","total_tetracycline_ADD","total_MLS_ADD",
                                   "feed_MLS_ADD","feed_tetracycline_ADD","parenteral_tetracycline_ADD","parenteral_MLS_ADD",
                                   "parenteral_phenicol_ADD","parenteral_sulfonamide_ADD")]
str(meta.df)
meta.df$Days_since_tx <- as.integer(meta.df$Days_since_tx)
scaled.meta.df <- as.data.frame(scale(meta.df))
# error during scaling, parenteral_betalactams_ADD, parenteral_fluoroquinolones_ADD


# check that we get mean of 0 and sd of 1
colMeans(scaled.meta.df)  # faster version of apply(scaled.dat, 2, mean)
apply(scaled.meta.df, 2, sd)

scaled.meta.df$Feedlot_categories <-  pData(AMR_class_exit)$Feedlot_categories
#scaled.meta.df$PEN.ID <-  pData(AMR_class_exit)$PEN.ID

meta.df$Feedlot_categories <-  pData(AMR_class_exit)$Feedlot_categories
#meta.df$PEN.ID <-  pData(AMR_class_exit)$PEN.ID


##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ ., data = scaled.meta.df)
anova(mod1, by = "term", perm = 10000)
anova(mod1, perm = 10000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data = scaled.meta.df)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm")
anova(mod, by = "term", perm = 999)
anova(mod, perm = 1000)
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <- ordiplot(mod,type="n")
points(mod, "sites", pch=21, col="red", bg="black")
text(fig, "biplot", col="red", cex=1)


### A priori test model ###
setDT(meta.df)
meta.df[,.(MLS_feed_total = sum(feed_MLS_ADD), Tetracycline_feed_total = sum(feed_tetracycline_ADD))]

mod_apriori <- rda(hell.norm.mt ~ DOF  + feed_MLS_ADD + feed_tetracycline_ADD  , data =  scaled.meta.df)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <- ordiplot(mod_apriori,type="n")
points(mod_apriori, "sites", pch=21, col=as.factor(meta.df$Time), bg="black")
text(fig, "biplot", col="red", cex=0.9)
#text(fig, "species", col="red", cex=0.9)


##
###
####
#####
## microbiome All samples
#####
####
###
##


pData(microbiome_analytic_data[[2]])$Days_since_tx <- as.integer(pData(microbiome_analytic_data[[1]])$Days_since_tx)


##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
p = cumNormStatFast(microbiome_analytic_data[[2]])
### We adjust our container to contain the normalization pecentile using cumNorm()
microbiome_analytic_data[[2]] = cumNorm(microbiome_analytic_data[[2]],p)
### We find the normalization factors per sample using normFactors()
nf = normFactors(microbiome_analytic_data[[2]])
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(microbiome_analytic_data[[2]])
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)

meta.df = pData(microbiome_analytic_data[[2]])[,c("Total_ADD","Feed_ADD","Parenteral_ADD","DOF","num_tx","total_tetracycline_ADD","total_MLS_ADD",
                                                  "feed_MLS_ADD","feed_tetracycline_ADD","parenteral_tetracycline_ADD","parenteral_MLS_ADD",
                                                  "parenteral_phenicol_ADD","parenteral_sulfonamide_ADD")]
setDT(meta.df)
#meta.df <- meta.df[Days_since_tx !="None"]


str(meta.df)
scaled.meta.df <- as.data.frame(scale(meta.df))
# error during scaling, parenteral_betalactams_ADD, parenteral_fluoroquinolones_ADD
meta.df$Time <- pData(microbiome_analytic_data[[2]])$Time

# check that we get mean of 0 and sd of 1
colMeans(scaled.meta.df)  # faster version of apply(scaled.dat, 2, mean)
apply(scaled.meta.df, 2, sd)

#scaled.meta.df$Feedlot_categories <-  pData(microbiome_phylum_exit)$Feedlot_categories
scaled.meta.df$Time <-  pData(microbiome_analytic_data[[2]])$Time

#meta.df$Feedlot_categories <-  pData(microbiome_phylum_exit)$Feedlot_categories
#meta.df$PEN.ID <-  pData(microbiome_class_exit)$PEN.ID


##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ ., data = scaled.meta.df)
anova(mod1, by = "term", perm = 10000)
anova(mod1, perm = 10000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data = scaled.meta.df)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm")
anova(mod, by = "term", perm = 999)
anova(mod, perm = 1000)
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <- ordiplot(mod,type="n")
points(mod, "sites", pch=21, col=as.factor(scaled.meta.df$Time), bg=as.factor(scaled.meta.df$Time))
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red"), cex=2.5, labels = c("Arrival","Exit"),adj = c(-.3,-3))




### A priori test model ###
setDT(meta.df)
meta.df[,.(MLS_feed_total = sum(feed_MLS_ADD), Tetracycline_feed_total = sum(feed_tetracycline_ADD))]

mod_apriori <- rda(hell.norm.mt ~ DOF  +Time , Total_ADD  , data =  scaled.meta.df)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <- ordiplot(mod_apriori,type="n")
points(mod_apriori, "sites", pch=21, col=as.factor(meta.df$Time), bg="black")
text(fig, "biplot", col="red", cex=0.9)
#text(fig, "species", col="red", cex=0.9)






##
###
####
#####
## microbiome split by Time of collection
#####
####
###
##

microbiome_arrival_samples = which(pData(microbiome_analytic_data[[2]])$Time == "Arrival")
microbiome_exit_samples = which(pData(microbiome_analytic_data[[2]])$Time == "Exit")

pData(microbiome_analytic_data[[1]])$Days_since_tx <- as.integer(pData(microbiome_analytic_data[[1]])$Days_since_tx)

microbiome_phylum_arrival <- microbiome_analytic_data[[2]][, microbiome_arrival_samples]
microbiome_class_arrival <- microbiome_analytic_data[[3]][, microbiome_arrival_samples]
microbiome_phylum_exit <- microbiome_analytic_data[[2]][, microbiome_exit_samples]
microbiome_exit_samples_days <- which(pData(microbiome_phylum_exit)$Days_since_tx != "None")
microbiome_phylum_exit <- microbiome_phylum_exit[,microbiome_exit_samples_days]

microbiome_class_exit <- microbiome_analytic_data[[3]][, microbiome_exit_samples]


##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
p = cumNormStatFast(microbiome_phylum_exit)
### We adjust our container to contain the normalization pecentile using cumNorm()
microbiome_phylum_exit = cumNorm(microbiome_phylum_exit,p)
### We find the normalization factors per sample using normFactors()
nf = normFactors(microbiome_phylum_exit)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(microbiome_phylum_exit)
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)

meta.df = pData(microbiome_phylum_exit)[,c("Total_ADD","Feed_ADD","Parenteral_ADD","DOF","num_tx","Days_since_tx","total_tetracycline_ADD","total_MLS_ADD",
                                           "feed_MLS_ADD","feed_tetracycline_ADD","parenteral_tetracycline_ADD","parenteral_MLS_ADD",
                                           "parenteral_phenicol_ADD","parenteral_sulfonamide_ADD")]
setDT(meta.df)
#meta.df <- meta.df[Days_since_tx !="None"]

meta.df$Days_since_tx <- as.integer(meta.df$Days_since_tx)
str(meta.df)
scaled.meta.df <- as.data.frame(scale(meta.df))
# error during scaling, parenteral_betalactams_ADD, parenteral_fluoroquinolones_ADD


# check that we get mean of 0 and sd of 1
colMeans(scaled.meta.df)  # faster version of apply(scaled.dat, 2, mean)
apply(scaled.meta.df, 2, sd)

#scaled.meta.df$Feedlot_categories <-  pData(microbiome_phylum_exit)$Feedlot_categories
#scaled.meta.df$PEN.ID <-  pData(microbiome_class_exit)$PEN.ID

#meta.df$Feedlot_categories <-  pData(microbiome_phylum_exit)$Feedlot_categories
#meta.df$PEN.ID <-  pData(microbiome_class_exit)$PEN.ID


##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ ., data = scaled.meta.df)
anova(mod1, by = "term", perm = 10000)
anova(mod1, perm = 10000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data = scaled.meta.df)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm")
anova(mod, by = "term", perm = 999)
anova(mod, perm = 1000)
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <- ordiplot(mod,type="n")
points(mod, "sites", pch=21, col="red", bg="black")
text(fig, "biplot", col="red", cex=1)


### A priori test model ###
setDT(meta.df)
meta.df[,.(MLS_feed_total = sum(feed_MLS_ADD), Tetracycline_feed_total = sum(feed_tetracycline_ADD))]

mod_apriori <- rda(hell.norm.mt ~ DOF  + feed_MLS_ADD + feed_tetracycline_ADD  , data =  scaled.meta.df)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <- ordiplot(mod_apriori,type="n")
points(mod_apriori, "sites", pch=21, col=as.factor(meta.df$Time), bg="black")
text(fig, "biplot", col="red", cex=0.9)
#text(fig, "species", col="red", cex=0.9)



