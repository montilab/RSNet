test_that("null_ggm and differential wrappers return tabular statistics", {
  data("toy_er")

  dat <- as.data.frame(toy_er$dat[, 1:8, drop = FALSE])
  dat$phenotype <- rep(c("A", "B"), each = nrow(dat) / 2)

  nn <- null_ggm(
    dat = dat,
    group_col = "phenotype",
    inference_method = "D-S_NW_SL",
    shuffle_method = "permutation",
    shuffle_iter = 2,
    balanced = TRUE,
    filter = "none",
    n_cores = 1,
    seed = 1
  )

  testthat::expect_type(nn, "list")
  testthat::expect_named(nn, c("null_networks", "tag"))
  testthat::expect_equal(nn$tag, "null_ggm")
  testthat::expect_equal(length(nn$null_networks), 2)

  obs <- nn$null_networks[[1]]

  dc <- diff_centrality(
    obs_networks = obs,
    dat = dat,
    group_col = "phenotype",
    null_networks = nn,
    shuffle_method = "permutation"
  )
  dg <- diff_gdv(
    obs_networks = obs,
    dat = dat,
    group_col = "phenotype",
    null_networks = nn,
    sign = FALSE
  )

  testthat::expect_s3_class(dc, "data.frame")
  testthat::expect_s3_class(dg, "data.frame")
  testthat::expect_true(any(grepl("_pval$", colnames(dc))))
  testthat::expect_true(any(grepl("_pval$", colnames(dg))))
})
