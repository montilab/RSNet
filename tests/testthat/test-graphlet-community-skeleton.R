test_that("gdvm/gcm functions return expected structures", {
  g <- model_er(p = 18, e_prob = 0.2, DAG = FALSE)
  igraph::V(g)$name <- paste0("V", seq_len(igraph::vcount(g)))

  gdvm <- gdvm_gcm(g, level = "4", include_gcm = FALSE)
  both <- gdvm_gcm(g, level = "4", include_gcm = TRUE)

  testthat::expect_true(is.matrix(gdvm))
  testthat::expect_equal(nrow(gdvm), igraph::vcount(g))
  testthat::expect_type(both, "list")
  testthat::expect_named(both, c("gdvm", "gcm"))
  testthat::expect_equal(gcm_distance(both$gcm, both$gcm), 0)

  intra <- intra_gdv_distance(gdvm)
  testthat::expect_equal(dim(intra), c(nrow(gdvm), nrow(gdvm)))

  pair <- paired_gdv_distance(gdvm, gdvm)
  testthat::expect_true(is.numeric(pair))
  testthat::expect_equal(length(pair), nrow(gdvm))
})

test_that("signed graphlet and community tools work on small graphs", {
  g <- igraph::make_ring(6)
  igraph::E(g)$sign <- rep(c(1, -1), length.out = igraph::ecount(g))
  igraph::V(g)$name <- paste0("N", seq_len(igraph::vcount(g)))

  sgdv <- signed_gdvm_gcm(g, n_cores = 1, include_gcm = FALSE)
  testthat::expect_true(is.matrix(sgdv))
  testthat::expect_equal(nrow(sgdv), igraph::vcount(g))

  g_iso <- g + igraph::vertices("isolated")
  igraph::V(g_iso)$name <- c(paste0("N", 1:6), "isolated")
  out_iso <- remove_isolated(g_iso, mode = "total", vertex_symbol = "name")
  testthat::expect_type(out_iso, "list")
  testthat::expect_true("ig" %in% names(out_iso))

  cd <- community_detection(g, method = "walk_trap", steps = 4)
  testthat::expect_type(cd, "list")
  testthat::expect_true(all(c("ig", "community_partition") %in% names(cd)))

  g_col <- colorize_community(cd$ig, community_attr_name = "community")
  testthat::expect_true("color" %in% igraph::vertex_attr_names(g_col))
})
