test_that("model generators and simulation utilities return expected shapes", {
  g_er <- model_er(p = 20, e_prob = 0.1, DAG = FALSE)
  g_sw <- model_sw(p = 20, dim = 1, nei = 2, rw_prob = 0.05)
  g_pa <- model_pa(p = 20, power = 1, z_appeal = 1, DAG = FALSE)
  g_mpa <- model_mpa(p = 60, m = 4, q_hub = 0.8, m_links = 3, power = 1.2, z_appeal = 1)
  g_sbm <- model_sbm(p = 20, m = 4, sizes = c(5, 5, 5, 5), p_in = 0.7, p_out = 0.1)

  testthat::expect_s3_class(g_er, "igraph")
  testthat::expect_s3_class(g_sw, "igraph")
  testthat::expect_s3_class(g_pa, "igraph")
  testthat::expect_s3_class(g_mpa, "igraph")
  testthat::expect_s3_class(g_sbm, "igraph")
  testthat::expect_equal(igraph::vcount(g_mpa), 60)

  sim <- model_sim(n = 30, ig = g_er, seed = 1)
  testthat::expect_s3_class(sim, "data.frame")
  testthat::expect_equal(nrow(sim), 30)
  testthat::expect_equal(ncol(sim), igraph::vcount(g_er))
})

test_that("resampling and matrix utility functions behave consistently", {
  sample_data <- data.frame(
    ID = paste0("S", 1:10),
    class = rep(c("A", "B"), each = 5),
    cluster = rep(1:5, each = 2),
    stringsAsFactors = FALSE
  )

  b <- resample(sample_data = sample_data[, c("ID", "class")], sample_class = sample_data$class, boot = TRUE)
  s <- resample(sample_data = sample_data[, c("ID", "class")], sample_class = sample_data$class, boot = FALSE, sub_ratio = 0.5)
  c_res <- resample_cluster(sample_data = sample_data[, c("ID", "cluster")], id_col = "ID", cluster_col = "cluster", boot = TRUE)

  testthat::expect_type(b, "character")
  testthat::expect_true(length(s) > 0)
  testthat::expect_type(c_res, "list")
  testthat::expect_named(c_res, c("bagging", "bagging_cluster"))

  p_mat <- matrix(runif(16), nrow = 4)
  p_mat <- (p_mat + t(p_mat)) / 2
  diag(p_mat) <- 1
  q_mat <- matrix_p_adjust(p_mat)
  testthat::expect_true(is.matrix(q_mat))
  testthat::expect_equal(q_mat, t(q_mat))

  upper <- c(0.1, 0.2, 0.3)
  rec <- upper_tri_to_matrix(upper, variable_names = c("A", "B", "C"), diagl = 1)
  testthat::expect_equal(dim(rec), c(3, 3))
  testthat::expect_equal(diag(rec), rep(1, 3))
})

test_that("graph metrics and misc helpers return valid outputs", {
  d <- gdv_distance(c(0, 1, 2), c(1, 1, 3))
  s <- gdv_distance(c(0, 1, 2), c(1, 1, 3), similarity = TRUE)
  testthat::expect_true(is.numeric(d) && length(d) == 1)
  testthat::expect_true(is.numeric(s) && length(s) == 1)

  j1 <- jaccard_similarity(c(1, 2, 3), c(2, 3, 4))
  j2 <- jaccard_similarity(diag(3), diag(3))
  testthat::expect_true(j1 >= 0 && j1 <= 1)
  testthat::expect_equal(j2, 1)

  cols <- distinct_colors(25)
  vals <- normalize_range(c(1, 2, 3), a = -1, b = 1)
  colorized <- colorize(c(0.1, 0.5, 0.9))
  captured <- capture_all({ cat("hidden output"); 42 })

  testthat::expect_length(cols, 25)
  testthat::expect_equal(range(vals), c(-1, 1))
  testthat::expect_length(colorized, 3)
  testthat::expect_equal(captured, 42)
})
