test_that("plot_browser_tracks.R gives correct error", {

  sample_data_table <- data.frame(bw_path = c("bw1.bw", "bw2.bw"),
                                  sample_id = c("bw1", "bw2"))

  expect_error(plot_browser_tracks(data_table = sample_data_table))

})
