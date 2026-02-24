test_that("ensemble_cgbn is available when RHugin is installed", {
  testthat::skip_if_not_installed("RHugin")
  data("toy_cgbn")

  res <- ensemble_cgbn(
    dat = toy_cgbn,
    discrete_variable = paste0("D", 1:5),
    num_iteration = 1,
    boot = FALSE,
    sub_ratio = 0.8,
    hugin = FALSE,
    n_cores = 1
  )

  testthat::expect_type(res, "list")
  testthat::expect_true(all(c("ig_networks", "ig_skeletons", "baggings") %in% names(res)))
  testthat::expect_equal(length(res$ig_networks), 1)
})
