test_that("consensus GGM and plotting return expected objects", {
  data("toy_er")

  ens <- ensemble_ggm(
    dat = toy_er$dat[, 1:10, drop = FALSE],
    num_iteration = 2,
    sample_class = NULL,
    estimate_CI = TRUE,
    boot = TRUE,
    method = "D-S_NW_SL",
    alpha = 0.05,
    global = TRUE,
    n_cores = 1
  )

  cons <- consensus_net_ggm(
    ggm_networks = ens,
    filter = "fdr",
    threshold = 0.2
  )

  testthat::expect_type(cons, "list")
  testthat::expect_true("consensus_network" %in% names(cons))
  testthat::expect_s3_class(cons$consensus_network, "igraph")

  vis <- plot_cn(
    ig = cons$consensus_network,
    vertex_symbol = "name",
    edge_label = "pcor",
    CI_show = TRUE
  )
  testthat::expect_type(vis, "list")
  testthat::expect_true(all(c("p", "data") %in% names(vis)))
})

test_that("consensus CGBN works with both supported methods", {
  g1 <- igraph::graph_from_literal(A -+ B, B -+ C, C -+ D)
  g2 <- igraph::graph_from_literal(A -+ B, A -+ C, C -+ D)
  g3 <- igraph::graph_from_literal(B -+ C, A -+ D, C -+ D)
  networks <- list(g1, g2, g3)

  avg_res <- consensus_net_cgbn(networks = networks, method = "average", cut = 0.34)
  all_res <- consensus_net_cgbn(networks = networks, reference_network = g1, method = "all", cut = 0.34)

  testthat::expect_type(avg_res, "list")
  testthat::expect_s3_class(avg_res$consensus_network, "igraph")
  testthat::expect_equal(avg_res$method, "average")

  testthat::expect_type(all_res, "list")
  testthat::expect_s3_class(all_res$consensus_network, "igraph")
  testthat::expect_equal(all_res$method, "all")
})
