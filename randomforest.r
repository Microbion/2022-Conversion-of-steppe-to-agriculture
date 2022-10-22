library(tidyverse)
library(phyloseq)
library(vegan)
# use phyloseq to align sample names in profile table and meta table
ps <- phyloseq(read.delim("data/kegg_KO_profile.xls", row.names = 1) %>%
               select(-Total, -Description, -Hyperlink) %>%
               otu_table(taxa_are_rows = TRUE),
read.delim("data/env.txt", row.names = 1) %>% sample_data)
ps <- rarefy_even_depth(ps, rngseed = TRUE, verbose = FALSE)
ko_alpha <- estimate_richness(ps, measures=c("Observed", "Shannon"))
ps <- merge_phyloseq(ps, ko_alpha)

# generate RF models
env_nat <- subset_samples(ps, is_natural) %>% sample_data %>% data.frame
env_agr <- subset_samples(ps, !is_natural) %>% sample_data %>% data.frame
train_set_nat <- sample(nrow(env_nat), nrow(env_nat) * 0.7)
train_set_agr <- sample(nrow(env_agr), nrow(env_agr) * 0.7)
train.nat <- rfPermute(env_nat[train_set_nat,4]~., data=env_nat[train_set_nat, c("pH", "AP", "TC")], importance=TRUE)
train.agr <- rfPermute(env_agr[train_set_agr,4]~., data=env_agr[train_set_agr, c("pH", "TC.TN", "TC")], importance=TRUE)
pred.nat <- predict(train.nat, env_nat[-train_set_nat,])
pred.agr <- predict(train.agr, env_agr[-train_set_agr,])
pred_all.nat <- predict(train.nat, env_nat)
pred_all.agr <- predict(train.agr, env_agr)

cor(env_nat[,4], pred_all.nat)^2
cor(env_nat[-train_set_nat,4], pred.nat)^2
cor(env_agr[,4], pred_all.agr)^2
cor(env_agr[-train_set_agr,4], pred.agr)^2

importance(train.nat)
importance(train.agr)