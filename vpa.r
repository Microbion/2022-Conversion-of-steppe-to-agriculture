library(tidyverse)
library(phyloseq)
library(vegan)
source("scripts/pcnm.r")
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

varpart.nat <- pcnm_vpa(otu_nat, env_nat)
varpart.agr <- pcnm_vpa(otu_agr, env_agr)

# plot
var_agr <- varpart(otu_agr,varpart.agr[[1]],varpart.agr[[2]],varpart.agr[[3]])
var_agr$part$indfract$Adj.R.square * 100
plot(var_agr)
var_nat <- varpart(otu_nat,varpart.nat[[1]],varpart.nat[[2]],varpart.nat[[3]])
var_nat$part$indfract$Adj.R.square * 100
plot(var_nat)
