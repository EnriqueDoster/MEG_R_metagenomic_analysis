## Analysis for ZIG model 
library(statmod)

## Filter Data
AMR_raw_analytic_data <- meg_filter_data(AMR_raw_analytic_data,filter_min_threshold = 0.15)
microbiome_analytic_data <- meg_filter_data(microbiome_analytic_data,filter_min_threshold = 0.15)


## Metadata variables for analysis
Time = pData(AMR_analytic_data[[1]])$Time
Treatment = pData(AMR_analytic_data[[1]])$Treatment
Animal = pData(AMR_analytic_data[[1]])$Animal
Batch = pData(AMR_analytic_data[[1]])$Batch

## By time
AMR_arrival_samples = which(pData(AMR_analytic_data[[1]])$Time == "Arrival")
AMR_day11_samples = which(pData(AMR_analytic_data[[1]])$Time == "Day11")

AMR_class_arrival <- AMR_raw_analytic_data[[1]][, AMR_arrival_samples]
AMR_mech_arrival <- AMR_raw_analytic_data[[2]][, AMR_arrival_samples]
AMR_class_day11 <- AMR_raw_analytic_data[[1]][, AMR_day11_samples]
AMR_mech_day11 <- AMR_raw_analytic_data[[2]][, AMR_day11_samples]

## By treatment
AMR_treated_samples = which(pData(AMR_analytic_data[[1]])$Treatment == "Treated")
AMR_untreated_samples = which(pData(AMR_analytic_data[[1]])$Treatment == "Untreated")

AMR_class_treated <- AMR_analytic_data[[1]][, AMR_treated_samples]
AMR_mech_treated <- AMR_analytic_data[[2]][, AMR_treated_samples]
AMR_class_untreated <- AMR_analytic_data[[1]][, AMR_untreated_samples]
AMR_mech_untreated <- AMR_analytic_data[[2]][, AMR_untreated_samples]



## AMR group subsets
AMR_arrival_treated_samples = which(pData(AMR_analytic_data[[1]])$Group == "Arrival-Treated")
AMR_arrival_untreated_samples = which(pData(AMR_analytic_data[[1]])$Group == "Arrival-Untreated")
AMR_day11_treated_samples = which(pData(AMR_analytic_data[[1]])$Group == "Day11-Treated")
AMR_day11_untreated_samples = which(pData(AMR_analytic_data[[1]])$Group == "Day11-Untreated")
#Class
AMR_class_arrival_treated <- AMR_raw_analytic_data[[1]][, AMR_arrival_treated_samples]
AMR_class_arrival_untreated <- AMR_raw_analytic_data[[1]][, AMR_arrival_untreated_samples]
AMR_class_day11_treated <- AMR_raw_analytic_data[[1]][, AMR_day11_treated_samples]
AMR_class_day11_untreated <- AMR_raw_analytic_data[[1]][, AMR_day11_untreated_samples]
#Mechanism
AMR_mech_arrival_treated <- AMR_raw_analytic_data[[2]][, AMR_arrival_treated_samples]
AMR_mech_arrival_untreated <- AMR_raw_analytic_data[[2]][, AMR_arrival_untreated_samples]
AMR_mech_day11_treated <- AMR_raw_analytic_data[[2]][, AMR_day11_treated_samples]
AMR_mech_day11_untreated <- AMR_raw_analytic_data[[2]][, AMR_day11_untreated_samples]



settings = zigControl(maxit=20, verbose=TRUE)





## By time
zero_mod <- model.matrix(~0+log(libSize(amr)))
designTimeTreatment = model.matrix(~0 + Time)

dup_resTimeTreatment <- duplicateCorrelation(MRcounts(AMR_analytic_data[[1]]), block=Animal, design= designTimeTreatment)

resAll_Class_TimeTreatment = fitZig(obj= AMR_analytic_data[[1]][2,], mod = designTimeTreatment, control = settings,zeroMod=zero_mod,useCSSoffset=FALSE, useMixedModel=TRUE, block=Animal)
zigFit_Class_TimeTreatment = resAll_Class_TimeTreatment$fit
finalMod_Class_TimeTreatment = resAll_Class_TimeTreatment$fit$design
contrast_Class_TimeTreatment= makeContrasts(TimeDay11-TimeDay1, levels=finalMod_Class_TimeTreatment)
res2_Class_TimeTreatment = contrasts.fit(zigFit_Class_TimeTreatment, contrast_Class_TimeTreatment)
resEB_Class_TimeTreatment = eBayes(res2_Class_TimeTreatment )


fz_Mechtrt_bh <- topTable(resEB_Class_TimeTreatment, coef=1, adjust.method="BH",number = 1000)
fz_Mechtrt <- topTable(resEB_Class_TimeTreatment, coef=1,number = 1000)
fz_Mechtrt_MDR_bh <- topTable(resEB_Class_TimeTreatment, coef=1, adjust.method="BH",number = 1000)
fz_Mechtrt_MDR <- topTable(resEB_Class_TimeTreatment, coef=1,number = 1000)
fz_Mechtrt_MLS_bh <- topTable(resEB_Class_TimeTreatment, coef=1, adjust.method="BH",number = 1000)
fz_Mechtrt_MLS <- topTable(resEB_Class_TimeTreatment, coef=1,number = 1000)


View(MRcounts(AMR_analytic_data[[1]][2,]))



## Treatment groups over time
AMR_class_treated


resAll_Class_TimeTreatment = fitZig(obj= AMR_class_treated, mod = designTimeTreatment, control = settings,zeroMod=zero_mod,useCSSoffset=FALSE, useMixedModel=dup_resTimeTreatment$consensus)
zigFit_Class_TimeTreatment = resAll_Class_TimeTreatment$fit
finalMod_Class_TimeTreatment = resAll_Class_TimeTreatment$fit$design
contrast_Class_TimeTreatment= makeContrasts(TimeDay11-TimeArrival, levels=finalMod_Class_TimeTreatment)
res2_Class_TimeTreatment = contrasts.fit(zigFit_Class_TimeTreatment, contrast_Class_TimeTreatment)
resEB_Class_TimeTreatment = eBayes(res2_Class_TimeTreatment )


zero_mod <-model.matrix(~1+log(libSize(MicroTreated)))

designAllMicroTreated = model.matrix(~0+TimeTreated + BatchTreated)
dup_resAllMicroTreated <- duplicateCorrelation(MRcounts(MicroTreated),block=Animal_Treated, design= designAllMicroTreated)
resAllMicroTreated = fitZig(obj= MicroTreated, mod = designAllMicroTreated, control = settings, zeroMod=zero_mod,useCSSoffset=TRUE, useMixedModel=dup_resAllMicroTreated$consensus)
zigFitAllMicroTreated = resAllMicroTreated$fit
finalModAllMicroTreated = resAllMicroTreated$fit$design
contrastAllMicroTreated = makeContrasts(TimeTreatedDay11-TimeTreatedArrival, levels=finalModAllMicroTreated)
resAllMicro2Treated = contrasts.fit(zigFitAllMicroTreated, contrastAllMicroTreated)
resAllMicro2EBTreated = eBayes(resAllMicro2Treated)


zero_mod <-model.matrix(~1+log(libSize(MicroUntreated)))

designAllMicroUntreated = model.matrix(~0+TimeUntreated + BatchUntreated)
dup_resAllMicroUntreated <- duplicateCorrelation(MRcounts(MicroUntreated),block=Animal_Untreated, design= designAllMicroUntreated)
resAllMicroUntreated = fitZig(obj = MicroUntreated, mod = designAllMicroUntreated, control = settings, zeroMod=zero_mod,useCSSoffset=TRUE, useMixedModel=dup_resAllMicroUntreated$consensus)
zigFitAllMicroUntreated = resAllMicroUntreated$fit
finalModAllMicroUntreated = resAllMicroUntreated$fit$design
contrastAllMicroUntreated = makeContrasts(TimeUntreatedDay11-TimeUntreatedArrival, levels=finalModAllMicroUntreated)
resAllMicro2Untreated = contrasts.fit(zigFitAllMicroUntreated, contrastAllMicroUntreated)
resAllMicro2EBUntreated = eBayes(resAllMicro2Untreated)

zero_mod <-model.matrix(~1+log(libSize(MicroArrival)))

designAllMicroArrival = model.matrix(~0+TreatmentArrival)
resAllMicroArrival = fitZig(obj = MicroArrival, mod = designAllMicroArrival, control = settings, useCSSoffset=TRUE)
zigFitAllMicroArrival = resAllMicroArrival$fit
finalModAllMicroArrival = resAllMicroArrival$fit$design
contrastAllMicroArrival = makeContrasts(TreatmentArrivalUntreated-TreatmentArrivalTreated, levels=finalModAllMicroArrival)
resAllMicro2Arrival = contrasts.fit(zigFitAllMicroArrival, contrastAllMicroArrival)
resAllMicro2EBArrival = eBayes(resAllMicro2Arrival)

zero_mod <-model.matrix(~1+log(libSize(MicroDay11)))

designAllMicroDay11 = model.matrix(~0+TreatmentDay11)
resAllMicroDay11 = fitZig(obj = MicroDay11, mod = designAllMicroDay11, control = settings, useCSSoffset=TRUE)
zigFitAllMicroDay11 = resAllMicroDay11$fit
finalModAllMicroDay11 = resAllMicroDay11$fit$design
contrastAllMicroDay11 = makeContrasts(TreatmentDay11Untreated-TreatmentDay11Treated, levels=finalModAllMicroDay11)
resAllMicro2Day11 = contrasts.fit(zigFitAllMicroDay11, contrastAllMicroDay11)
resAllMicro2EBDay11 = eBayes(resAllMicro2Day11)



fz_AllMicro_Time <- topTable(resAllMicro2EBTime, coef=1, adjust.method="BH",number = 1000)
fz_AllMicro_Treatment <- topTable(resAllMicro2EBTx, coef=1, adjust.method="BH",number = 1000)
fz_AllMicro_Treated <- topTable(resAllMicro2EBTreated, coef=1, adjust.method="BH",number = 1000)
fz_AllMicro_Untreated <- topTable(resAllMicro2EBUntreated, coef=1, adjust.method="BH",number = 1000)
fz_AllMicro_Arrival <- topTable(resAllMicro2EBArrival, coef=1, adjust.method="BH",number = 1000)
fz_AllMicro_Day11 <-topTable(resAllMicro2EBDay11, coef=1, adjust.method="BH",number = 1000)



##### Microbiome level
## Microbiome subsets
microbiome_arrival_treated_samples = which(pData(microbiome_analytic_data[[1]])$Group == "Arrival-Treated")
microbiome_arrival_untreated_samples = which(pData(microbiome_analytic_data[[1]])$Group == "Arrival-Untreated")
microbiome_day11_treated_samples = which(pData(microbiome_analytic_data[[1]])$Group == "Day11-Treated")
microbiome_day11_untreated_samples = which(pData(microbiome_analytic_data[[1]])$Group == "Day11-Untreated")
#phylum
microbiome_phylum_arrival_treated <- microbiome_analytic_data[[2]][, microbiome_arrival_treated_samples]
microbiome_phylum_arrival_untreated <- microbiome_analytic_data[[2]][, microbiome_arrival_untreated_samples]
microbiome_phylum_day11_treated <- microbiome_analytic_data[[2]][, microbiome_day11_treated_samples]
microbiome_phylum_day11_untreated <- microbiome_analytic_data[[2]][, microbiome_day11_untreated_samples]
#class
microbiome_class_arrival_treated <- microbiome_analytic_data[[3]][, microbiome_arrival_treated_samples]
microbiome_class_arrival_untreated <- microbiome_analytic_data[[3]][, microbiome_arrival_untreated_samples]
microbiome_class_day11_treated <- microbiome_analytic_data[[3]][, microbiome_day11_treated_samples]
microbiome_class_day11_untreated <- microbiome_analytic_data[[3]][, microbiome_day11_untreated_samples]
#order
microbiome_order_arrival_treated <- microbiome_analytic_data[[4]][, microbiome_arrival_treated_samples]
microbiome_order_arrival_untreated <- microbiome_analytic_data[[4]][, microbiome_arrival_untreated_samples]
microbiome_order_day11_treated <- microbiome_analytic_data[[4]][, microbiome_day11_treated_samples]
microbiome_order_day11_untreated <- microbiome_analytic_data[[4]][, microbiome_day11_untreated_samples]
#family
microbiome_family_arrival_treated <- microbiome_analytic_data[[5]][, microbiome_arrival_treated_samples]
microbiome_family_arrival_untreated <- microbiome_analytic_data[[5]][, microbiome_arrival_untreated_samples]
microbiome_family_day11_treated <- microbiome_analytic_data[[5]][, microbiome_day11_treated_samples]
microbiome_family_day11_untreated <- microbiome_analytic_data[[5]][, microbiome_day11_untreated_samples]
#genus
microbiome_genus_arrival_treated <- microbiome_analytic_data[[6]][, microbiome_arrival_treated_samples]
microbiome_genus_arrival_untreated <- microbiome_analytic_data[[6]][, microbiome_arrival_untreated_samples]
microbiome_genus_day11_treated <- microbiome_analytic_data[[6]][, microbiome_day11_treated_samples]
microbiome_genus_day11_untreated <- microbiome_analytic_data[[6]][, microbiome_day11_untreated_samples]


