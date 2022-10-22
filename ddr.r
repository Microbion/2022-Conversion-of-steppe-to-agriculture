library(tidyverse)
library(phyloseq)
library(vegan)
source("scripts/dist_3col.r")
# use phyloseq to align sample names in profile table and meta table
ps <- phyloseq(read.delim("data/kegg_KO_percent.xls", row.names = 1) %>%
               t %>%
               decostand("hellinger") %>%
               otu_table(taxa_are_rows = FALSE),
read.delim("data/env.txt", row.names = 1) %>% sample_data)

# split data in to agr(agricultual) and nat(natural)
is_natural <- sample_names(ps) %>% str_starts("N")

otu_agr <- subset_samples(ps, !is_natural) %>% otu_table %>% data.frame
otu_nat <- subset_samples(ps, is_natural) %>% otu_table %>% data.frame
env_agr <- subset_samples(ps, !is_natural) %>% sample_data %>% data.frame
env_nat <- subset_samples(ps, is_natural) %>% sample_data %>% data.frame

# geo distance
geo_dist_nat <- geoXY(env_nat$lat, env_nat$lon, unit = 1000) %>%
vegdist("euclidean") %>%
dist_3col
geo_dist_agr <- geoXY(env_agr$lat, env_agr$lon, unit = 1000) %>%
vegdist("euclidean") %>%
dist_3col

# env distance
env_nat %>% scale %>% vegdist("euclidean") %>% dist_3col
env_agr %>% scale %>% vegdist("euclidean") %>% dist_3col

# bio distance, default is dissimilarity. For, bray-curtis distance, Simmilarity = 1 - dissimilarity.
otu_dist_nat <- otu_nat %>% vegdist("bray") %>% dist_3col
otu_dist_agr <- otu_agr %>% vegdist("bray") %>% dist_3col