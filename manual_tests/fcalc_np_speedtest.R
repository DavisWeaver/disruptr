#grab ppi to use for the scaffold
g <- crosstalkr::prep_stringdb(cache = cache, min_score = 600)

#clean gdsc
df_list <- clean_gdsc(cache = cache)

#grab expression data
df_exp <- df_list[[3]]

#get list of cell lines to iterate over
cell_lines <- unique(df_exp$cell_line)
cell_line_test <- cell_lines[1]
df <- dplyr::filter(df_exp, cell_line == cell_line_test)

exp <- df$log_expression
names(exp) <- df$gene_symbols

microbenchmark::microbenchmark(calc_np_all(exp = exp, g = g), calc_np_all2(exp = exp, g = g), times = 1)
#calc np for that expression vector

np <- calc_np_all(exp, g)
np2 <- calc_np_all2(exp, g)
