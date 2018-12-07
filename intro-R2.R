##############
# 12/07/18
# GGG 201a
# Intro to R 2
##############

# this lesson material is available at:
# https://github.com/ngs-docs/2018-ggg201a

# Installing packages in R ------------------------------------------------

# If the package lives on CRAN
install.packages("devtools")
install.packages(c("RColorBrewer", "ggthemes"))

# If the package lives on Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicFeatures", version = "3.8")
#BiocManager::install(c("ggbio", "GenomicAlignments", "VariantAnnotation", "Biostrings"), version = "3.8")

source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("ggbio")

# If the package lives on github
library(devtools)
install_github('arendsee/rhmmer')

# Once a package is installed, use library() to load its functionality
library()

# Organizing data in R ----------------------------------------------------
library(dplyr)
library(tidyr)

ecoli <- read.csv("data2/Ecoli_metadata_composite.csv")
head(ecoli)
dim(ecoli)

# filter out sequences for which mutator is not known
ecoli_mut <- ecoli %>%
              filter(!is.na(mutator))
dim(ecoli_mut)

# select specific columns
ecoli_sub <- ecoli %>%
              select(strain, generation, clade, mutator, cit)
dim(ecoli_sub)

# string these functions together
ecoli_keep <- ecoli %>%
                select(strain, generation, clade, mutator, cit) %>%
                filter(!is.na(mutator))

# create a new column for the number of basepairs in the dataset
ecoli <- ecoli %>%
          mutate(total_bp = read_length * sequencing_depth)

# separate author from year
ecoli <- separate(data = ecoli, col = reference, into = c("author", "year"), sep = "\\.")

# put the period back at the end of "et al"
ecoli$author <- gsub("et al", "et al.", ecoli$author)

# how many samples came from each study?
ecoli %>%
  group_by(author) %>%
  tally()

# How many different generations exist in the data?
ecoli %>%
  group_by(generation) %>%
  tally() %>%
  nrow()

# How many citrate plus mutants occurred per generation?
ecoli %>%
  filter(cit == "plus") %>%
  group_by(cit, generation) %>%
  tally()

# Working with sequencing data in R ---------------------------------------
library(Biostrings)

# download two E. coli genomes
# rel 606
url <- "https://github.com/ngs-docs/2018-ggg201a/blob/master/data2/GCF_000017985.1_ASM1798v1_protein.faa.gz?raw=true"
download.file(url = url, destfile = "data2/GCF_000017985.1_ASM1798v1_protein.faa.gz")

# k-12
url <- "https://github.com/ngs-docs/2018-ggg201a/blob/master/data2/GCA_000005845.2_ASM584v2_protein.faa.gz?raw=true"
download.file(url = url, destfile = "data2/GCA_000005845.2_ASM584v2_protein.faa.gz")

# read in the files
rel606 <- readAAStringSet(filepath = "data2/GCF_000017985.1_ASM1798v1_protein.faa.gz")
k12 <- readAAStringSet(filepath = "data2/GCA_000005845.2_ASM584v2_protein.faa.gz")

length(rel606)
head(names(rel606))
rel606[10]

rel606 %in% k12
rel606 %in% k12
length(rel606[rel606 %in% k12])

rel606[duplicated(rel606)]

alphabetFrequency(rel606)

vcountPattern(pattern = "VIL", subject = rel606)

names(rel606[10])
grep("hydrogenase 1 maturation protease", x = names(k12))
aln <- pairwiseAlignment(rel606[10], k12[887])
pid(aln)

# Visualizing data with ggplot2 -------------------------------------------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(sf)
library(ggbio)
#library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(VariantAnnotation)

## ggplot2, ggthemes, and RColorBrewer
ggplot(ecoli, aes(x = cit, y = generation)) + 
  geom_point()

ggplot(ecoli %>%
         filter(cit != ""), 
       aes(x = cit, y = generation)) + 
  geom_jitter()

ggplot(ecoli %>%
         filter(cit != ""), 
       aes(x = cit, y = generation)) + 
  geom_jitter(size = 5)

ggplot(ecoli %>%
         filter(cit != ""),
       aes(x = cit, y = generation)) + 
  geom_jitter(aes(size = read_length))

ggplot(ecoli %>%
         filter(cit != ""), 
       aes(x = cit, y = generation, color = clade)) + 
  geom_jitter(aes(size = read_length))

ggplot(ecoli%>%
         filter(cit != "") %>%
         filter(clade != ""), 
       aes(x = cit, y = generation, color = clade)) + 
  geom_jitter(aes(size = read_length))

ggplot(ecoli%>%
         filter(cit != "") %>%
         filter(clade != ""), 
       aes(x = cit, y = generation, color = clade)) + 
  geom_jitter(aes(size = read_length)) +  
  theme_wsj()

ggplot(ecoli%>%
         filter(cit != "") %>%
         filter(clade != ""), 
       aes(x = cit, y = generation, color = clade)) + 
  geom_jitter(aes(size = read_length)) +  
  theme_wsj() + 
  scale_color_brewer(palette = "Paired")

## plotting a map


nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

ggplot(nc) +
  geom_sf(aes(fill = AREA))

ggplot(nc) +
  geom_sf(aes(fill = AREA)) +
  geom_point(data = data.frame(lon = -79, lat = 35.5), 
             aes(x = lon, y = lat))

ggplot(nc) +
  geom_sf(aes(fill = AREA)) +
  geom_point(data = data.frame(lon = -79, lat = 35.5), 
             aes(x = lon, y = lat), color = "white", shape = 2)


## ggbio
# download some more data
url = "http://biocluster.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_12_16_2013/Rgraphics/data.zip"
download.file(url = url, destfile = "data2/data.zip")
unzip("data2/data.zip", exdir = "data2")
list.files("data2")

#library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(VariantAnnotation)

ga <- readGAlignments(file = "data2/data/SRR064167.fastq.bam", 
                      use.names=TRUE, 
                      param=ScanBamParam(which=GRanges("Chr5", IRanges(4000, 8000))))


p1 <- autoplot(ga)
p1

p2 <- autoplot(ga, geom = "line", stat = "coverage")
p2

vcf <- readVcf(file="data2/data/varianttools_gnsap.vcf", 
               genome="ATH1")
vcf

p3 <- autoplot(vcf[seqnames(vcf)=="Chr5"], type = "fixed") + 
        xlim(4000, 8000) + 
        theme(legend.position = "none", 
              axis.text.y = element_blank(), 
              axis.ticks.y=element_blank())

txdb <- makeTxDbFromGFF(file="./data/TAIR10_GFF3_trunc.gff", 
                        format="gff3")

p4 <- autoplot(txdb, 
               which = GRanges("Chr5", IRanges(4000, 8000)), 
               names.expr = "gene_id")
p4

tracks(Reads=p1, 
       Coverage=p2, 
       Variant=p3, 
       Transcripts=p4, 
       heights = c(0.3, 0.2, 0.1, 0.35)) + 
  ylab("")


# Runnng R on the command line --------------------------------------------




# About the datasets ------------------------------------------------------

# ecoli:
# - The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment).
# 
# - The experiment was designed to assess adaptation in *E. coli*. A population (designated **Ara-3**) were propagated for more than 40,000 generations in a glucose-limited minimal medium (in most conditions glucose is the best carbon source for *E. coli*, providing faster growth than other sugars). This medium was supplemented with citrate which *E. coli* cannot metabolize in the aerobic conditions of the experiment. Sequencing of the populations at regular time points reveals that spontaneous citrate-using variant (**Cit+**) appeared between 31,000 and 31,500 generations causing an increase in population size and diversity. In addition, this experiment showed hypermutability in certain regions. Hypermutability is important and can help accelerate adaptation to novel environments, but also can be selected against in well-adapted populations.
# 
# - To see a timeline of the experiment to date, check out this [figure](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment#/media/File:LTEE_Timeline_as_of_May_28,_2016.png), and this paper [Blount et al. 2008: Historical contingency and the evolution of a key innovation in an experimental population of *Escherichia coli*](http://www.pnas.org/content/105/23/7899).


# Arabadopsis:
# Illumina sequencing of Arabidopsis thaliana Translatome at flower stage 4, replicate 2

# Yeast:
# “How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?” Schurch et al., 2016.

# Where to get help -------------------------------------------------------
# Davis R Users Group
# https://d-rug.github.io/
# join list serve: https://groups.google.com/d/forum/davis-rug
# Fall 2018: Meets Wed, 10-12pm, 360 Sheilds Library (DSI classroom)

# Meet and Analyze Data
# join list serve: training+subscribe@dib-lab.groups.io
# Wed 3-5pm, usually in a room in Valley Hall of CCAH

# Other resources ---------------------------------------------------------

# ANGUS materials
# https://angus.readthedocs.io/en/2018/

# to find examples/explanation about code, run:
browseVignettes("ggplot2")

