if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat package is required to run tests")
}

testthat::test_dir("tests/testthat")
