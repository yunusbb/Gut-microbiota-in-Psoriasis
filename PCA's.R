############################################################################
#
# *********** Converting to relative abundance - compositions
#

# convert to compositional, i.e. relative abundance


library(microbiome)
library(dplyr)

# convert to relative abundance.

comp_psq <- microbiome::transform(psq1, "compositional")


# Aggregate at different taxonomic levels:
#  ("Kingdom", "Phylum", "Class","Order", "Family","Genus","Species")

# at phylum level:
comp_psq_phylum <- aggregate_taxa(comp_psq, level = "Phylum")

# at class level:
comp_psq_class <- aggregate_taxa(comp_psq, level = "Class")

# at order level:
comp_psq_order <- aggregate_taxa(comp_psq, level = "Order")

# at family level:
comp_psq_family <- aggregate_taxa(comp_psq, level = "Family")

# at genus level:
comp_psq_genus <- aggregate_taxa(comp_psq, level = "Genus")


#########################################################################
# *************** PCA at different taxonomic levels:
#
##########################################
# ******* at PHYLUM level:
ord_phylum <- ordinate(comp_psq_phylum, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq_phylum, ord_phylum, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom")

pcaplot + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_PHYLUM_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)
#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_PHYLUM_wh_5_5_inch.pdf",width=5,height=5)


# extract eigenvalues

eigs_phylum <- ord_phylum$values

eigs_phylum$Eigenvalues[1]/sum(eigs_phylum$Eigenvalues)*100

eigs_phylum$Eigenvalues[2]/sum(eigs_phylum$Eigenvalues)*100

# At phylum level, PCs explain:
# PC1 and PC2 -  84.75929+11.796 = 96.6%
#
# For comparison at Species level:
# first two PCs explain 34%+8% = 42% of total variation


##########################################
# ******* at CLASS level:


ord_class <- ordinate(comp_psq_class, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq_class, ord_class, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom") + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_CLASS_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_CLASS_wh_5_5_inch.pdf",width=5,height=5)

eigs_class <- ord_class$values

eigs_class$Eigenvalues[1]/sum(eigs_class$Eigenvalues)*100

eigs_class$Eigenvalues[2]/sum(eigs_class$Eigenvalues)*100

# At Class level, PCs explain:
# PC1 and PC2 -  75.90892+13.27988 = 89.2%
#
# For comparison at Species level:
# first two PCs explain 34%+8% = 42% of total variation


##########################################################
# at ORDER level:
ord_order <- ordinate(comp_psq_order, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq_order, ord_order, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom") + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_ORDER_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)
#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_ORDER_wh_5_5_inch.pdf",width=5,height=5)

eigs_order <- ord_order$values

eigs_order$Eigenvalues[1]/sum(eigs_order$Eigenvalues)*100
eigs_order$Eigenvalues[2]/sum(eigs_order$Eigenvalues)*100

# At Class level, PCs explain:
# PC1 and PC2 -  74.97379+13.30784 = 88.3%
#
# For comparison at Species level:
# first two PCs explain 34%+8% = 42% of total variation

##########################################################
# at FAMILY level:

ord_family <- ordinate(comp_psq_family, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq_family, ord_family, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom") + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_FAMILY_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)
#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_FAMILY_wh_5_5_inch.pdf",width=5,height=5)


#######################################
# TRY ENCODING Calprotectin with colors

min(sample_data(comp_psq_genus)$Calprotectin,na.rm=T)
max(sample_data(comp_psq_genus)$Calprotectin,na.rm=T)

round_any(exp(seq(log(1), log(409.35), length = 6)), 5)

my_breaks = c(0, 5, 10, 35, 125, 410)

pcaplot <- plot_ordination(comp_psq_family, ord_family, color = "Calprotectin", label="shortID",shape="Diagnosis") + geom_point(size = 3) + theme_classic()

pcaplot + scale_colour_gradient2()

# no apparent correlation with PC1 or PC2 axes!

eigs_family <- ord_family$values

eigs_family$Eigenvalues[1]/sum(eigs_family$Eigenvalues)*100

eigs_family$Eigenvalues[2]/sum(eigs_family$Eigenvalues)*100

# 69.37967+12.29686 = 81.67653
#
# For comparison at Species level:
# first two PCs explain 34%+8% = 42% of total variation

############################################
# at GENUS level:

ord_genus <- ordinate(comp_psq_genus, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq_genus, ord_genus, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom") + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_GENUS_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)
#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_GENUS_wh_5_5_inch.pdf",width=5,height=5)

eigs_genus <- ord_genus$values

eigs_genus$Eigenvalues[1]/sum(eigs_genus$Eigenvalues)*100
eigs_genus$Eigenvalues[2]/sum(eigs_genus$Eigenvalues)*100

#  66.89841+11.19056 = 78.08897

ord_species <- ordinate(comp_psq, "PCoA", "bray")

pcaplot <- plot_ordination(comp_psq, ord_species, color = "Diagnosis", label="shortID",shape="Diagnosis") + geom_point(size=0.5) + theme_classic() + theme( legend.position = "bottom") + scale_shape_manual(values=c(3,4)) + scale_colour_manual(values = c("control"="green4","psoriasis" ="red"))

# one copy there:

#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_307_species_wh_3.5_3.5_inch-PAPER.pdf",width=3.5,height=3.5)
#dev.copy2pdf(file="PCA_rel_abund_Bray-Curtis_distance_based/PCoA_bray-curtis_307_species_wh_5_5_inch-PAPER.pdf",width=5,height=5)


# extract eigenvalues - various summary about eigenvalues (principal components)

eigsummary <- ord_species$values

eigsummary$Eigenvalues[1]/sum(eigsummary$Eigenvalues)*100
eigsummary$Eigenvalues[2]/sum(eigsummary$Eigenvalues)*100

# first PC explains 34% of total variation
# second PC explains 8%

# # ("Kingdom", "Phylum", "Class","Order", "Family","Genus","Species")