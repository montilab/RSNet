test_that("ensemble_ggm returns expected output",{

  data('toy_er')

  ## Learning via bootstrap
  ensemble_networks <- ensemble_ggm(dat = toy_er$dat,
                                    num_iteration = 10,
                                    sample_class =  NULL,
                                    estimate_CI = TRUE,
                                    boot = TRUE,
                                    sub_ratio = NULL,
                                    correlated = FALSE,
                                    method = "D-S_NW_SL",
                                    alpha = 0.05,
                                    global = TRUE,
                                    n_cores = 2)



  testthat::expect_length(ensemble_networks, 8L)
  testthat::expect_type(ensemble_networks, "list")
  testthat::expect_named(ensemble_networks,
                         c("method","estimate_CI","n","p","partialCor_mat",
                           "avg_z_score_partialCor","vids","baggings"))


  testthat::expect_type(ensemble_networks$partialCor_mat,"double")
  testthat::expect_type(ensemble_networks$avg_z_score_partialCor,"double")
  testthat::expect_type(ensemble_networks$bagging,"list")
  testthat::expect_length(ensemble_networks$bagging, 10L)

  ## Consensus Network
  cons_network <- consensus_net_ggm(ggm_networks = ensemble_networks,
                                    filter = "fdr",
                                    threshold = 0.05)

  testthat::expect_length(cons_network, 4L)
  testthat::expect_type(cons_network, "list")
  testthat::expect_named(cons_network,c("consensus_network","pcor_CI","CI","method"))

})
