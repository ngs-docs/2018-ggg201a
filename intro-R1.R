##############
# 11/30/18
# GGG 201a
# Intro to R 1
##############

# This lesson is adapted from:
# http://evomics.org/learning/programming/introduction-to-r/
# this lesson material is available at:
# https://github.com/ngs-docs/2018-ggg201a

library(tidyverse)
library(vegan)
library(pheatmap)

# Download data files

url <- "https://raw.githubusercontent.com/ngs-docs/2018-ggg201a/master/data/myoviridae.csv"
download.file(url = url, destfile = "myoviridae.csv")

# Read data into R
myov <- read.csv("myoviridae.csv", row.names = 1)

# look at the top of each file
head(myov)
summary(myov)

# Download metadata file
url <- "https://raw.githubusercontent.com/ngs-docs/2018-ggg201a/master/data/metadata.csv"
download.file(url = url, destfile = "metadata.csv")

# Read data into R
metadata <- read.csv("metadata.csv")

# look at the top of the file
head(metadata)
summary(metadata)

# look at the first column of the metadata file
metadata[ , 1]

# plot the age of healthy and sick patients
ggplot(data = metadata, aes(x = Diagnosis, y = Age)) +
  geom_boxplot() + 
  geom_jitter() +
  ggtitle("Comparison of Age Between Groups")

# transform the count data
myov_total <- decostand(myov, method="total")
?decostand()

# move the rownames into a column
myov_total_labeled$X = rownames(myov_total)

# change the format of the count data
myov_long <- gather(data = myov_total_labeled, key = species, value = count, -X)
head(myov_long)
# join the count data with metadata
all_data <- left_join(myov_long, metadata, by = "X")


# pull out key species in myov data
key_species <- c("Tevenvirinae", "Punalikevirus", "Clostridium_phase_c.st", "PhiCD119likevirus")

key_species_data <- all_data %>%
                      filter(species %in% key_species)

# plot the abundances of key species in the myov object
ggplot(data = key_species_data, aes(x = Diagnosis, y = count)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ species)

# plot a heatmap of all species
pheatmap(myov_total, cellwidth=8, cellheight=8, main="Healthy vs. Sick")

# export the transformed data that went into the heatmap
write.csv(myov_total, "myoviridae_total.csv", quote = F)
