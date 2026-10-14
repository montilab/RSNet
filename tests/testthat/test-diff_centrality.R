test_that("bootstrap differential centrality returns three columns per statistic", {
  node_names <- paste0("node_", seq_len(60))
  g1 <- igraph::make_ring(60)
  g2 <- igraph::add_edges(
    g1, as.vector(rbind(seq_len(60), (seq_len(60) + 4L) %% 60L + 1L))
  )
  igraph::V(g1)$name <- node_names
  igraph::V(g2)$name <- node_names
  obs <- list(one = g1, two = g2)
  nulls <- list(
    null_networks = rep(list(obs), 10),
    tag = "null_ggm"
  )
  statistics <- c("degree", "strength", "eigen", "betweenness", "closeness", "pagerank")
  expected_names <- as.vector(rbind(
    paste0(statistics, "_diff"),
    paste0(statistics, "_pval"),
    paste0(statistics, "_95%_CI")
  ))
  expected_diff <- as.matrix(centrality(g2)) - as.matrix(centrality(g1))

  for (alternative in c("two.sided", "greater", "less")) {
    out <- diff_centrality(
      obs, dat = NULL, null_networks = nulls,
      shuffle_method = "bootstrap", shuffle_iter = 10,
      alternative = alternative
    )

    expect_s3_class(out, "data.frame")
    expect_identical(dim(out), c(60L, 18L))
    expect_named(out, expected_names)
    expect_identical(rownames(out), node_names)
    for (j in seq_along(statistics)) {
      expect_equal(out[[3L * j - 2L]], unname(expected_diff[, j]))
      expect_type(out[[3L * j - 1L]], "double")
      expect_type(out[[3L * j]], "character")
      expect_false(anyNA(out[[3L * j]]))
    }
  }
})

test_that("permutation differential centrality retains two columns per statistic", {
  g1 <- igraph::make_ring(60)
  g2 <- igraph::make_full_graph(60)
  obs <- list(one = g1, two = g2)
  nulls <- list(null_networks = rep(list(obs), 10), tag = "null_ggm")

  out <- diff_centrality(
    obs, dat = NULL, null_networks = nulls,
    shuffle_method = "permutation", shuffle_iter = 10
  )

  statistics <- colnames(centrality(g1))
  expect_identical(dim(out), c(60L, 12L))
  expect_named(out, as.vector(rbind(
    paste0(statistics, "_diff"), paste0(statistics, "_pval")
  )))
})
