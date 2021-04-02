load(system.file("test_data/toy_graph.Rda", package = "disruptr"))
#np_all
seeds <- c("OLR1", "APP", "VAV2", "ITGAV", "JAG1", "APOH")
toy_exp <- c(4.9, 9.9, 1.0, 0.2, 7.5, 8.4)
names(toy_exp) <- seeds

expected_np_old <- c(15.7, -4.8, Inf, -0.94, -2.23, -1.38)
names(expected_np_old) <- seeds
v_check <- seeds[5]
g <- igraph::induced_subgraph(g, seeds)

expected_np_diff <- c(0.0, 6.21, NaN, 0.082, 0.0)
names(expected_np_diff) <- seeds[seeds != v_check]


test_that("node_repression produces expected results", {
  expect_equal(node_repression(g, v_rm = v_check,
                  state_function = calc_np_all, exp = toy_exp),
               expected_np_diff, tolerance = 0.1)

})

