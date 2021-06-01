library(disruptr)
library(dplyr)
library(tidyr)
library(magrittr)
cache = "G:/My Drive/MIR_Combo_Targeting/code/miRNA_Targeting/data_files/"

g <- crosstalkr::prep_biogrid(cache = cache)
#isolate a single

load(paste0(cache, "stjude_counts.Rda"))
jude_df <- jude_df %>% filter(sj_diseases == "EWS")
jude_df <- pivot_wider(jude_df, id_cols = Geneid, names_from = sample_name,
                       values_from = log_count)

jude_dfnp <- compute_np(cache = cache, experiment_name = "sj", ppi = "biogrid",
                        exp_mat= jude_df, ncores = 4)
