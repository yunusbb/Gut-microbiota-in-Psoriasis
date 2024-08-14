# ACTUAL: Quasi-pair Wilcoxon test########
# with CLR transformation

################################################################
# NEAREST NEIGHBOR MATCHING - WITHOUT computing Propensity Score

###########################################################################
# *****  PSEUDO-PAIR CONSTRUCTION based on GENUS-level PCs
# 

install.packages('MatchIt') # R package
install.packages('optmatch')
install.packages("PairedData")


# load

library(MatchIt)
library(optmatch)
library(PairedData)


library("phyloseq"); packageVersion("phyloseq")
#[1] ‘1.48.0’

library("ggplot2"); packageVersion("ggplot2")
#[1] ‘3.5.1’

library("plyr"); packageVersion("plyr")
#1] ‘1.8.9’


########################################################################################
#
# *********** Converting to relative abundance - compositions
#
#************ Aggregating at different taxonomic levels
#

# ("Kingdom", "Phylum", "Class","Order", "Family","Genus","Species")

######################################################
# convert to compositional, i.e. relative abundance


library(microbiome)
library(dplyr)

id_diagnosis<-sample_data(psq3)[,c('ID','Diagnosis')]

psor_id <- c(id_diagnosis[id_diagnosis$Diagnosis=='psoriasis',1])$ID


### at phylum level: #####


psq_phylum <- aggregate_taxa(psq3, level = "Phylum")


# STEP 1:
# replace 0 values with the count zero multiplicative method
# this function expects the samples to be in rows

otu_phylum_CZM <- t(cmultRepl(t(otu_table(psq_phylum)), method="CZM", output="p-counts"))


# CONVERT to PROPORTIONS by sample (columns) using the apply function

otu_phylum_CZM_prop <- apply(otu_phylum_CZM, 2, function(x){x/sum(x)})

# STEP 2: CLR Transformation:

clr_phylum <- t(apply(otu_phylum_CZM_prop, 2, function(x){log(x) - mean(log(x))}))



# NOTE: sample_id order will be identical at aggregated phyloseq objects

# 1) take the transformed abundance table, 
# 2) extract row-names, i.e. sample ids:
# 3) find matching ids for psoriasis patients, assign 'green4' if match OR 'red' if not - controls

mycols=data.frame(c(ifelse(rownames(clr_phylum) %in% psor_id,'red','green4')))

colnames(mycols) <- "colors"

# Singular value decompositon: PCA 

pcx_phylum <- prcomp(clr_phylum)

# Sum the total variance

d.mvar <- sum(pcx_phylum$sdev^2)

# modify donor IDS - rownames - make shorter to look nice on the plot:

row.names(pcx_phylum$x) <- gsub("_profile","", row.names(pcx_phylum$x))


# make a covariance biplot of the data with compositions function

# BIG PLOT

coloredBiplot(pcx_phylum, col="grey66", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              arrow.len=0.05,
              xlab=paste("PC1: ", round(sum(pcx_phylum$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_phylum$sdev[2]^2)/d.mvar, 3)),
              expand=0.8,var.axes=T, scale=1, main="Biplot")


#dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Phylum_7_7_inch.pdf",width=7,height=7)



# REMOVE Species names

par(mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,cex.sub=0.7,cex=0.7)

coloredBiplot(pcx_phylum, col="#ffffff00", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              var.axes=F,
              xlab=paste("PC1: ", round(sum(pcx_phylum$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_phylum$sdev[2]^2)/d.mvar, 3)),
              expand=0.8, scale=1, main="Biplot")


#dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Phylum_TRANSPARENT_3.5_3.5_inch.pdf",width=3.5,height=3.5)

# extract eigenvalues

round(sum(pcx_phylum$sdev[1]^2)/d.mvar, 3)
[1] 0.423

round(sum(pcx_phylum$sdev[2]^2)/d.mvar, 3)
[1] 0.31

round(sum(pcx_phylum$sdev[3]^2)/d.mvar, 3)
[1] 0.224

# At phylum level, top PCs explain:
#  (0.232+0.179+0.15)*100 = 56%

### at class level: ####

psq_class <- aggregate_taxa(psq3, level = "Class")

# ******* at CLASS level:

# replace 0 values with the count zero multiplicative method

otu_class_CZM <- t(cmultRepl(t(otu_table(psq_class)), method="CZM", output="p-counts"))

# CONVERT to PROPORTIONS by sample (columns) using the apply function

otu_class_CZM_prop <- apply(otu_class_CZM, 2, function(x){x/sum(x)})

# STEP 2: CLR Transformation:

clr_class <- t(apply(otu_class_CZM_prop, 2, function(x){log(x) - mean(log(x))}))


# Singlular value decompositon: PCA 

pcx_class <- prcomp(clr_class)

# Sum the total variance

d.mvar <- sum(pcx_class$sdev^2)

# trim donor IDS - rownames - make shorter to look nice on the plot:

row.names(pcx_class$x) <- gsub("_profile","", row.names(pcx_class$x))


# make a covariance biplot of the data with compositions function

# BIG PLOT

coloredBiplot(pcx_class, col="grey66", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              arrow.len=0.05,
              xlab=paste("PC1: ", round(sum(pcx_class$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_class$sdev[2]^2)/d.mvar, 3)),
              expand=0.8,var.axes=T, scale=1, main="Biplot")


#dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Class_7_7_inch.pdf",width=7,height=7)


# REMOVE Species names

par(mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,cex.sub=0.7,cex=0.7)

coloredBiplot(pcx_class, col="#ffffff00", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              var.axes=F,
              xlab=paste("PC1: ", round(sum(pcx_class$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_class$sdev[2]^2)/d.mvar, 3)),
              expand=0.8, scale=1, main="Biplot")


#dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Class_TRANSPARENT_3.5_3.5_inch.pdf",width=3.5,height=3.5)

# extract eigenvalues

round(sum(pcx_class$sdev[1]^2)/d.mvar, 3)
[1] 0.234

round(sum(pcx_class$sdev[2]^2)/d.mvar, 3)
[1] 0.17
> 
  round(sum(pcx_class$sdev[3]^2)/d.mvar, 3)
[1] 0.117


# At class level, PCs explain:
# PC1 and PC2 -  (0.10+0.08+0.072)*100 = 25%


### at order level: ####

psq_order <- aggregate_taxa(psq3, level = "Order")

# at ORDER level:

# replace 0 values with the count zero multiplicative method

otu_order_CZM <- t(cmultRepl(t(otu_table(psq_order)), method="CZM", output="p-counts"))


# CONVERT to PROPORTIONS by sample (columns) using the apply function

otu_order_CZM_prop <- apply(otu_order_CZM, 2, function(x){x/sum(x)})

# STEP 2: CLR Transformation:

clr_order <- t(apply(otu_order_CZM_prop, 2, function(x){log(x) - mean(log(x))}))


# Singular value decompositon: PCA 

pcx_order <- prcomp(clr_order)

# Sum the total variance

d.mvar <- sum(pcx_order$sdev^2)

# trim donor IDS - rownames - make shorter to look nice on the plot:

row.names(pcx_order$x) <- gsub("_profile","", row.names(pcx_order$x))


# make a covariance biplot of the data with compositions function


# BIG PLOT

coloredBiplot(pcx_order, col="grey66", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              arrow.len=0.05,
              xlab=paste("PC1: ", round(sum(pcx_order$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_order$sdev[2]^2)/d.mvar, 3)),
              expand=0.8,var.axes=T, scale=1, main="Biplot")


dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Order_7_7_inch.pdf",width=7,height=7)



# REMOVE Species names

par(mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,cex.sub=0.7,cex=0.7)

coloredBiplot(pcx_order, col="#ffffff00", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              var.axes=F,
              xlab=paste("PC1: ", round(sum(pcx_order$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_order$sdev[2]^2)/d.mvar, 3)),
              expand=0.8, scale=1, main="Biplot")


dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Order_TRANSPARENT_3.5_3.5_inch.pdf",width=3.5,height=3.5)



# extract eigenvalues

round(sum(pcx_order$sdev[1]^2)/d.mvar, 3)
[1] 0.179
round(sum(pcx_order$sdev[2]^2)/d.mvar, 3)
[1] 0.139
round(sum(pcx_order$sdev[3]^2)/d.mvar, 3)
[1] 0.108


# At class level, PCs explain:
# PC1 and PC2 -  0.097+0.083 = 18%


### at family level: ####

psq_family <- aggregate_taxa(psq3, level = "Family")

# at FAMILY level:

# replace 0 values with the count zero multiplicative method

otu_family_CZM <- t(cmultRepl(t(otu_table(psq_family)), method="CZM", output="p-counts"))


# CONVERT to PROPORTIONS by sample (columns) using the apply function

otu_family_CZM_prop <- apply(otu_family_CZM, 2, function(x){x/sum(x)})

# STEP 2: CLR Transformation:

clr_family <- t(apply(otu_family_CZM_prop, 2, function(x){log(x) - mean(log(x))}))


# Singular value decompositon: PCA 

pcx_family <- prcomp(clr_family)

# Sum the total variance

d.mvar <- sum(pcx_family$sdev^2)

# trim donor IDS - rownames - make shorter to look nice on the plot:

row.names(pcx_family$x) <- gsub("_profile","", row.names(pcx_family$x))


# make a covariance biplot of the data with compositions function


# BIG PLOT

coloredBiplot(pcx_family, col="grey66", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              arrow.len=0.05,
              xlab=paste("PC1: ", round(sum(pcx_family$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_family$sdev[2]^2)/d.mvar, 3)),
              expand=0.8,var.axes=T, scale=1, main="Biplot")

#dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Family_7_7_inch.pdf",width=7,height=7)

# REMOVE Species names

par(mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,cex.sub=0.7,cex=0.7)

coloredBiplot(pcx_family, col="#ffffff00", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              var.axes=F,
              xlab=paste("PC1: ", round(sum(pcx_family$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_family$sdev[2]^2)/d.mvar, 3)),
              expand=0.8, scale=1, main="Biplot")

dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Family_TRANSPARENT_3.5_3.5_inch.pdf",width=3.5,height=3.5)

# extract eigenvalues

round(sum(pcx_family$sdev[1]^2)/d.mvar, 3)
[1] 0.136
round(sum(pcx_family$sdev[2]^2)/d.mvar, 3)
[1] 0.119
round(sum(pcx_family$sdev[3]^2)/d.mvar, 3)
[1] 0.077



### at genus level: ####

psq_genus <- aggregate_taxa(psq3, level = "Genus")


# prepare vector of colors to paint the samples:


############################################
# at GENUS level:

# replace 0 values with the count zero multiplicative method

otu_genus_CZM <- t(cmultRepl(t(otu_table(psq_genus)), method="CZM", output="p-counts"))


# CONVERT to PROPORTIONS by sample (columns) using the apply function

otu_genus_CZM_prop <- apply(otu_genus_CZM, 2, function(x){x/sum(x)})

###########################
# STEP 2: CLR Transformation:

clr_genus <- t(apply(otu_genus_CZM_prop, 2, function(x){log(x) - mean(log(x))}))

# Singular value decompositon: PCA 

pcx_genus <- prcomp(clr_genus)

# Sum the total variance

d.mvar <- sum(pcx_genus$sdev^2)

# trim donor IDS - rownames - make shorter to look nice on the plot:

row.names(pcx_genus$x) <- gsub("_profile","", row.names(pcx_genus$x))


# make a covariance biplot of the data with compositions function


# BIG PLOT

coloredBiplot(pcx_genus, col="grey66", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              arrow.len=0.05,
              xlab=paste("PC1: ", round(sum(pcx_genus$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_genus$sdev[2]^2)/d.mvar, 3)),
              expand=0.8,var.axes=T, scale=1, main="Biplot")

dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Genus_7_7_inch.pdf",width=7,height=7)


# REMOVE Species names

par(mgp=c(2,0.5,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,cex.sub=0.7,cex=0.7)

coloredBiplot(pcx_genus, col="#ffffff00", cex=c(0.6, 0.5), xlabs.col=mycols$colors,
              var.axes=F,
              xlab=paste("PC1: ", round(sum(pcx_genus$sdev[1]^2)/d.mvar, 3)),
              ylab=paste("PC2: ", round(sum(pcx_genus$sdev[2]^2)/d.mvar, 3)),
              expand=0.8, scale=1, main="Biplot")


dev.copy2pdf(file="PC_biplots/Tripleton_filtered/Biplot_Genus_TRANSPARENT_3.5_3.5_inch.pdf",width=3.5,height=3.5)





# extract eigenvalues


round(sum(pcx_genus$sdev[1]^2)/d.mvar, 3)
[1] 0.136
round(sum(pcx_genus$sdev[2]^2)/d.mvar, 3)
[1] 0.083
round(sum(pcx_genus$sdev[3]^2)/d.mvar, 3)
[1] 0.061


##############
###########################################################################
# *****  PSEUDO-PAIR CONSTRUCTION based on GENUS-level PCs
# 

install.packages('MatchIt') # R package

install.packages('optmatch')

install.packages("PairedData")


# load

library(MatchIt)

library(optmatch)


library(PairedData)


library("phyloseq"); packageVersion("phyloseq")

library("ggplot2"); packageVersion("ggplot2")

library("plyr"); packageVersion("plyr")


# ("Kingdom", "Phylum", "Class","Order", "Family","Genus","Species")


###########################################################################
# 
# ***********  MATCHED-PAIR CONSTRUCTION based on FAMILY-level PCs
#                   (a.k.a Quasi-pair construction)
#

# packages:

library(MatchIt)

library(optmatch)

library("phyloseq"); packageVersion("phyloseq")

library("ggplot2"); packageVersion("ggplot2")

library("plyr"); packageVersion("plyr")


eigs_family <- ord_family$values

eigs_family$Eigenvalues[1]/sum(eigs_family$Eigenvalues)*100
69.37967

eigs_family$Eigenvalues[2]/sum(eigs_family$Eigenvalues)*100
12.29686

# 69.37967+12.29686 = 81.67653


id_diagnosis <- sample_data(comp_psq)[,c('ID','Diagnosis')]

class(id_diagnosis)
[1] "sample_data"
attr(,"package")
[1] "phyloseq"

class(id_diagnosis) <- c("data.frame")

PC_coords <- ord_family$vectors

head(id_diagnosis)
ID Diagnosis
D9_profile   D9_profile   control
D99_profile D99_profile   control


ids_pcs <- data.frame(IDs=c(row.names(PC_coords)),PC_coords[,1:3])

ids_pcs[1:3,]
IDs     Axis.1      Axis.2       Axis.3
D9_profile   D9_profile -0.1525980 -0.07106407  0.000130203
D99_profile D99_profile -0.2198335  0.15713327  0.013817484


PCs_status <- merge(x=id_diagnosis, y=ids_pcs, by.x='ID', by.y='IDs',sort=F)


PCs_status[1:3,1:5]
ID Diagnosis     Axis.1      Axis.2       Axis.3
1  D9_profile   control -0.1525980 -0.07106407  0.000130203
2 D99_profile   control -0.2198335  0.15713327  0.013817484


PCs_status$Diagnosis <- as.factor(PCs_status$Diagnosis)

#########################################
# 1) Check Imbalance

psom.out0 <- matchit(Diagnosis ~ Axis.1 + Axis.2, data = PCs_status, method = NULL, distance = "glm")

# Checking balance prior to matching

summary(psom.out0)
Call:
  matchit(formula = Diagnosis ~ Axis.1 + Axis.2, data = PCs_status, 
          method = NULL, distance = "glm")

Summary of Balance for All Data:
  Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
distance        0.5309        0.5290          0.0859     1.0417    0.0556   0.1690
Axis.1         -0.0060        0.0067         -0.0408     0.9931    0.0306   0.1224
Axis.2          0.0047       -0.0053          0.0774     0.9690    0.0377   0.0943

Sample Sizes:
  Control Treated
All            47      53
Matched        47      53
Unmatched       0       0
Discarded       0       0


# Quite well balanced

####
#################################################################
# NEAREST NEIGHBOR MATCHING - WITHOUT computing Propensity Score
#

# RAW Euclidean distances are used to find matches:

psomNN_euclidean <- matchit(Diagnosis ~ Axis.1 + Axis.2, data = PCs_status, method = "nearest", distance = "euclidean", replace=T)

nrow(psomNN_euclidean$match.matrix)
[1] 53


NODUPS_psomNN_euclidean <- matchit(Diagnosis ~ Axis.1 + Axis.2, data = PCs_status, method = "nearest", distance = "euclidean", replace=F)

uniq_matches <- NODUPS_psomNN_euclidean$match.matrix[c(!is.na(NODUPS_psomNN_euclidean$match.matrix)),]



################################################################
# extract pairs - NN on Euclidean distance based on PC1 and PC2
# psomNN_euclidean$match.matrix


library(MatchIt)

library(optmatch)


library(PairedData)


library(PairedData)

# RAW Euclidean distances are used to find matches:

Euclid_NN_pair_index <- data.frame(treated=row.names(psomNN_euclidean$match.matrix), control=psomNN_euclidean$match.matrix)

Uniq_Euclid_NN_pair_index <- data.frame(treated=names(uniq_matches), control=uniq_matches)

nrow(Euclid_NN_pair_index)
[1] 53


Euclid_NN_pair_index[1:3,]
treated control
21      21      34
22      22      77
24      24      34


nrow(Uniq_Euclid_NN_pair_index)
[1] 47


Uniq_Euclid_NN_pair_index[1:3,]
treated control
21      21      34
22      22      77
24      24      12


# replace with sample IDs

Eucl_NN_pairIDs=cbind(PCs_status$ID[as.numeric(Euclid_NN_pair_index$treated)],PCs_status$ID[as.numeric(Euclid_NN_pair_index$control)])

Eucl_NN_pairIDs[1:5,]
[,1]          [,2]         
[1,] "D81_profile" "D6_profile" 
[2,] "D80_profile" "D30_profile"
[3,] "D79_profile" "D6_profile" 
...

# All pairs are maximally close to each other!!!
# Good solution


Uniq_Eucl_NN_pairIDs<- cbind(PCs_status$ID[as.numeric(Uniq_Euclid_NN_pair_index$treated)],PCs_status$ID[as.numeric(Uniq_Euclid_NN_pair_index$control)])

Uniq_Eucl_NN_pairIDs[1:5,]
[,1]          [,2]         
[1,] "D81_profile" "D6_profile" 
[2,] "D80_profile" "D30_profile"
[3,] "D79_profile" "D8_profile" 
...


####
############################################################################
#
# PREPARE TWO versions of the input dataset for wilcoxon matched-paired rank sum test
#
# Version 1: separate phyloseq elements: metadata and otu_table, where samples are grouped by pair

# ->> Input for microbiotaPair function 
# Version 2: phyloseq with metadata, where samples are grouped by case-control
# ->> Input for R base wilcox.test function

# Version 1: phyloseq with metadata, where samples are grouped by pair 
# ->> Input for microbiotaPair function 

# IMPORTANT: make sure that sample IDs are stored in 'SampleID' column
# this is what is being used the R wrapper function:

library("varhandle")

######################################################
# TRY with unique pairs:


uniq_samples_by_pairs=data.frame(SampleID=c(Uniq_Eucl_NN_pairIDs[,1],Uniq_Eucl_NN_pairIDs[,2]), ppairs=c(paste("pair",1:nrow(Uniq_Eucl_NN_pairIDs),sep="_"),paste("pair",1:nrow(Uniq_Eucl_NN_pairIDs),sep="_")), group=c(rep("psoriasis",length(Uniq_Eucl_NN_pairIDs[,1])),rep("control",length(Uniq_Eucl_NN_pairIDs[,1]))))


head(uniq_samples_by_pairs,6)
SampleID ppairs     group
1 D81_profile pair_1 psoriasis
2 D80_profile pair_2 psoriasis
3 D79_profile pair_3 psoriasis
4 D78_profile pair_4 psoriasis
5 D77_profile pair_5 psoriasis
6 D76_profile pair_6 psoriasis



uniq_samples_by_pairs <-uniq_samples_by_pairs[order(uniq_samples_by_pairs$ppairs),]

head(uniq_samples_by_pairs,6)

head(uniq_samples_by_pairs,6)
SampleID  ppairs     group
1  D81_profile  pair_1 psoriasis
48  D6_profile  pair_1   control
10 D72_profile pair_10 psoriasis
57 D29_profile pair_10   control
11 D71_profile pair_11 psoriasis
58 D17_profile pair_11   control


comp_psq3 <- microbiome::transform(psq3, "compositional")
comp_psq3

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 249 taxa and 100 samples ]
sample_data() Sample Data:       [ 100 samples by 57 sample variables ]
tax_table()   Taxonomy Table:    [ 249 taxa by 7 taxonomic ranks ]

# copy/paste 'Wilcoxon_Sign_Rank_Test.R' function

W_test_res_occ_0.01 <- Wilcoxon_Sign_Rank_Test(physeq=comp_psq3, phen=uniq_samples_by_pairs,PID="ppairs",GROUP="group", grp1="control",grp2="psoriasis",paired=T,occ=0.01)

#Было 50 или более предупреждений (введите warnings() чтобы просмотреть первые 50)

# export select columns
write.table(W_test_res_occ_0.01[which(W_test_res_occ_0.01$Pvalue<=0.10),c(1,5,8)],file="quasi-pair-wilcox-p_0.1_occur_0.01_sp.txt",quote=F,row.names=F,sep="\t")

# export all
write.table(W_test_res_occ_0.01[which(W_test_res_occ_0.01$Pvalue<=0.10),],file="quasi-pair-wilcox-p_0.1_occur_0.01_ALL_RESULTS.txt",quote=F,row.names=F,sep="\t")
save.image(file="quasi-pair_analysis.RData")

write.table(uniq_samples_by_pairs, file="Matched_Pairs_IDs.txt", quote=F)




# SUCCESSFUL
# Note: Not all species survived internal filtering to run Wilcoxon test


dim(W_test_res_occ_0.1)
[1] 172  16


W_test_res_occ_0.1[which(W_test_res_occ_0.1$Pvalue<=0.05),c(1,5,7:10)]
type     Pvalue   Enrich_median Enrich_rank control_occurence psoriasis_occurence
1    Eubacterium_sp_CAG_180 0.00602998       psoriasis   psoriasis            0.2128              0.5106
2      Megasphaera_elsdenii 0.00971774 No significance   psoriasis            0.1277              0.4043
3 Ruminococcus_obeum_CAG_39 0.01078740 No significance     control            0.1915              0.1064
4   Roseburia_inulinivorans 0.01603560         control     control                 1              0.8936
5       Eubacterium_eligens 0.01882620         control     control            0.8936              0.6809
6  Haemophilus_sp_HMSC71H05 0.01893960         control     control            0.5106              0.2979
7        Clostridium_leptum 0.02042690 No significance   psoriasis            0.2979              0.4468
8       Rothia_mucilaginosa 0.02337640 No significance   psoriasis            0.0851              0.1702
9      Bacteroides_vulgatus 0.02793330         control     control            0.9787              0.9574




########################

###################################################################
# 
# PLOT for top species (14) matched-pair Wilcoxon test
#

load(file="quasi-pair_analysis.RData")

library(microbiome)
library(gcookbook)
library(tidyverse)



# Take a subset of species increased in psoriasis

top7_comp_psq3 <- subset_taxa(comp_psq3, Species == "Megasphaera_elsdenii" | Species == "Eubacterium_sp_CAG_180" | Species == "Rothia_mucilaginosa" | Species == "Bacteroides_xylanisolvens" | Species == "Clostridium_leptum" | Species == "Catenibacterium_mitsuokai" | Species == "Prevotella_sp_AM42_24")


# Keep only samples with Calprotectin data - remove NA's

top7_comp_psq3_calp <- subset_samples(top7_comp_psq3, Calprotectin != "NA")

top7_comp_psq3_calp
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 7 taxa and 95 samples ]

##################################
tax_table(top7_comp_psq3_calp)[,2]
Taxonomy Table:     [7 taxa by 1 taxonomic ranks]:
  Phylum          
Rothia_mucilaginosa       "Actinobacteria"
Bacteroides_xylanisolvens "Bacteroidetes" 
Prevotella_sp_AM42_24     "Bacteroidetes" 
Eubacterium_sp_CAG_180    "Firmicutes"    
Clostridium_leptum        "Firmicutes"    
Catenibacterium_mitsuokai "Firmicutes"    
Megasphaera_elsdenii      "Firmicutes"    



# Set your custom colors

rownames(tax_table(top7_comp_psq3_calp))
[1] "Rothia_mucilaginosa"       "Bacteroides_xylanisolvens" "Prevotella_sp_AM42_24"     "Eubacterium_sp_CAG_180"   
[5] "Clostridium_leptum"        "Catenibacterium_mitsuokai" "Megasphaera_elsdenii" 

colors_bad_7 <- c("black","grey55","grey77","cornflowerblue","violet","lightseagreen","mediumblue")

names( colors_bad_7) <- c(rownames(tax_table(top7_comp_psq3_calp)))


fp <- top7_comp_psq3_calp %>% plot_composition(sample.sort = "Calprotectin", otu.sort = "abundance", x.label = "shortID", group_by="Diagnosis") 

fp + theme_classic() + scale_fill_manual("Species",values = colors_bad_7) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(face = "italic", size=6)) + 
  theme(legend.key.size = unit(0.1, "inch"))

dev.copy2pdf(file="Matched-Pair-Wilcox-Top-7-bad-candidates-Calprotectin-sorted.pdf",width=7,height=3.0)


# print your colors
colors_bad_7
Rothia_mucilaginosa Bacteroides_xylanisolvens     Prevotella_sp_AM42_24    Eubacterium_sp_CAG_180 
"black"                  "grey55"                  "grey77"          "cornflowerblue" 
Clostridium_leptum Catenibacterium_mitsuokai      Megasphaera_elsdenii 
"violet"           "lightseagreen"              "mediumblue"


