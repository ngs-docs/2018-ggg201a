install.packages("tidyverse")
install.packages("vegan")
install.packages("pheatmap")
install.packages("sf")
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("GenomicAlignments", "VariantAnnotation", 
			  "Biostrings", "GenomicFeatures"))
