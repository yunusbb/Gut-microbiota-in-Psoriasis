########################################
# LOAD PACKAGES:


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

BiocManager::install("ALDEx2")

install.packages('devtools')
devtools::install_github('ggloor/CoDaSeq/CoDaSeq')

devtools::install_github("tpq/propr")


install.packages("./propr_4.2.6.tar.gz", repos = NULL, type = "source")

install.packages('fastcluster')


writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

pkgbuild::check_build_tools(debug = TRUE)

BiocManager::install("PERFect")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PERFect")



#

library(phyloseq)


library(dplyr)

library(tidyr)

library(stringr)


# Compositional

library(ade4)

library(knitr)

library(ALDEx2)

library(CoDaSeq)

library(zCompositions)

library(compositions)

###

###

library(igraph)

library(car)



library(vegan)

####

library(PERFect)

library(propr)