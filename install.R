if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocInstaller::install(c("GenomicAlignments", "VariantAnnotation", "Biostrings", "GenomicFeatures"))

install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("vegan")
install.packages("pheatmap")
