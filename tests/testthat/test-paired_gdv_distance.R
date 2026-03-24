test_that("paired_gdv_distance accepts data frames and preserves shared node names", {
  gdvm1 <- data.frame(
    orbit_1 = c(1, 2, 3),
    orbit_2 = c(4, 5, 6),
    row.names = c("node_a", "node_b", "node_c")
  )

  gdvm2 <- data.frame(
    orbit_1 = c(2, 4, 8),
    orbit_2 = c(1, 3, 5),
    row.names = c("node_b", "node_c", "node_d")
  )

  out <- paired_gdv_distance(gdvm1, gdvm2)

  testthat::expect_type(out, "double")
  testthat::expect_named(out, c("node_b", "node_c"))
  testthat::expect_length(out, 2L)
})
