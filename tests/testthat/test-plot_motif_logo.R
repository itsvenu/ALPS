test_that("plot_motif_logo() works fine", {

  ex_transfac <- system.file("extdata/motifs", "MA0147.2.transfac", package = "ALPS", mustWork = TRUE)

  plot_motif_logo(motif_path = ex_transfac, database = "transfac", plot_type = "logo")

})
