# tests/testthat.R
if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check("RSNet")
} else {
  message("Skipping tests because 'testthat' is not installed.")
}

