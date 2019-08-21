test_that("get_variable_regions() working properly", {

  mat <- matrix(sample.int(15, 9*100, TRUE), 9, 100) %>% as.data.frame()
  mat <- mat %>%
  tibble::rowid_to_column(var = "start") %>%
  dplyr::mutate(end = start + 1000) %>%
  dplyr::mutate(chr = "chr1") %>%
  dplyr::select(chr, start, end, dplyr::everything())

  res <- get_variable_regions(enrichments_df = mat, num_regions = 50)

  expect_equal(nrow(res), 50)

})
