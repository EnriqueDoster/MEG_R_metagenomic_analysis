## proj 4 bayesian regression
## need to match sample names, remove 1 from AMR
#knitr::opts_chunk$set(echo = TRUE)

#https://github.com/lakinsm/meg-bayesian-regression
#https://github.com/duncanwadsworth/dmbvs

library(data.table)
library(RCurl)
library(metagenomeSeq)
require(MASS)
require(dirmult)

## Publicly hosted scripts with helper functions
eval(parse(text=getURL("https://raw.githubusercontent.com/duncanwadsworth/dmbvs/master/code/helper_functions.R", ssl.verifypeer=FALSE)))
eval(parse(text=getURL("https://raw.githubusercontent.com/duncanwadsworth/dmbvs/master/code/wrapper.R", ssl.verifypeer=FALSE)))

## Local paths
# setwd("~/Dropbox/Reference_resources/Useful_code/R/dmbvs/")
# system("cd code; make")
executable_location = "~/Dropbox/Reference_resources/Useful_code/R/dmbvs/code/dmbvs.x"
save_prefix = "resistome"

## Data run from AmrPlusPlus analytic template
mech <- amr_mech

kraken_genus <- microbiome_genus
kraken_genus_analytic <- microbiome_genus_analytic
rownames(kraken_genus_analytic) <- kraken_genus$Genus
kraken_genus_analytic <- MRcounts(kraken_genus_analytic)

kraken_phylum <- microbiome_phylum
kraken_phylum_analytic <- microbiome_phylum_analytic
rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum
kraken_phylum_analytic <- MRcounts(kraken_phylum_analytic)



XX <- scale(t(mech[, .SD, .SDcols=!"mechanism"]), center=TRUE, scale=TRUE)
colnames(XX) <- mech$mechanism
YY <- t(kraken_phylum_analytic)
cat("Dimensions of Taxa Matrix: ", dim(YY), "\n")
cat("Dimensions of Covariate Matrix: ", dim(XX), "\n")

# MCMC and hyperparameters
# these values are reasonable for the data simulated here but should be changed
# depending on the characteristics of other datasets
#GG = 301L; thin = 2L; burn = 101L; # fast, for testing
GG = 11001L; thin = 10L; burn = 1001L; # good defaults, in this case
# reasonable default parameters, see further discussion in the manuscript
bb_alpha = 0.02; bb_beta = 2 - bb_alpha
proposal_alpha = 0.5; proposal_beta = 0.5
slab_variance = 10; intercept_variance = 10

# description
cat("Beta-Binomial mean:", bb_alpha/(bb_alpha + bb_beta), "\n")
cat("Number of kept iterations:", (GG - burn)/thin, "\n")

# Run MCMC, comment out after running 
results = dmbvs(XX = XX, YY = YY, intercept_variance = intercept_variance,
                slab_variance = slab_variance, bb_alpha = bb_alpha,
                bb_beta = bb_beta, GG = GG, thin = thin, burn = burn,
                init_beta = "warmstart", init_alpha = "warmstart",
                proposal_alpha = proposal_alpha, proposal_beta = proposal_beta,
                exec = executable_location, selection_type = "ss",
                output_location = "/home/enrique/Dropbox/Projects/NIFA_project345/analysis/Proj4_pulldown_analysis/amr_plus_plus/meg-bayesian-regression/")


params = data.frame(GG, burn, thin, intercept_variance,
                    slab_variance, bb_alpha, bb_beta,
                    proposal_alpha, proposal_beta)
save(results, params, XX, YY,
     file = paste0("results-", save_prefix, "-", Sys.Date(), ".RData"))

## Check results
mppi = colMeans((results$beta != 0) + 0)
(blfdrate = bfdr(mppi, threshold = 0.1)$threshold)
MPPI = data.frame(expand.grid(covariates = colnames(results$hyperparameters$inputdata$XX),
                              taxa = colnames(results$hyperparameters$inputdata$YY)),
                  mppi = mppi,
                  beta = colMeans(results$beta))
plot(mppi, type = "h", ylab = "MPPI",
     xlab = "beta index", main = "Manhattan plot")

# active variable traceplot
plot.ts(rowSums((results$beta != 0) + 0), main = "Active variables traceplot",
        ylab = "number of betas in the model", xlab = "iteration")

# some of the selected beta traceplots
selected = which(mppi > 0.5)
fortraces = selected[sample(length(selected), 4)] # use to be 10
plot.ts(results$beta[,fortraces], main = "Some selected beta traceplots",
        xlab = "iteration", ylab = "")

# visualize the associations
png(file="association_plot.png", width=1200, height=1200, units="px")
association_plot(MPPI, graph_layout = "bipartite", main = "Sample Results")
dev.off()

mm = subset(MPPI, mppi > 0.5)
mppi2 = mm[order(mm$covariates, mm$taxa),]
write.csv(mppi2, file="associations.csv", row.names=FALSE)

# visualize subset associations
png(file="association_plot.png", width=1200, height=1200, units="px")
association_plot(mppi2, graph_layout = "bipartite", main = "Sample Results")
dev.off()

