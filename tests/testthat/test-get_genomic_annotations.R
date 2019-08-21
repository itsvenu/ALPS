test_that("get_genomic_annotations gives correct error", {

  sample_data_table <- data.frame(bw_path = c("bw1.bw", "bw2.bw"),
                                  sample_id = c("bw1", "bw2"))

  expect_error(get_genomic_annotations(data_table = sample_data_table))

})
