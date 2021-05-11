library(disruptr)
df <- read.csv(system.file("test_data/rld_Counts.csv", package = "disruptr"))
df_small <- dplyr::slice_sample(df, prop = 0.02)

test1 <- disruptr::compute_np(cache = "G:/My Drive/data/mir_paper/", exp_mat = df_small,
           mir_paper = TRUE, ncores = 4, experiment_name = "test")

df_short <- df[,1:2] #this grabs exactly one sample
test2 <- disruptr::compute_np(cache = "G:/My Drive/data/mir_paper/", exp_mat = df_short,
                              mir_paper = TRUE, ncores = 4, experiment_name = "test")

Position(function(x) x == "AADACP1", x = test2$gene_name)
#Trying to figure out what gives with all the NAs
df_short <- tidy_expression(df_short)
exp <- df_short$expression
names(exp) <- df_short$gene_name

g <- crosstalkr::prep_biogrid(cache = "G:/My Drive/data/mir_paper/")



test_that("compute_")
