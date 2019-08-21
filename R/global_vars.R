# defining global variables and functions to appease R CMD Check

utils::globalVariables(
  names = c(
    ".",
    "PO",
    "V1",
    "V2",
    "V3",
    "avg",
    "bw_path",
    "chr",
    "end",
    "group",
    "plt_val",
    "pos",
    "region",
    "sample_id",
    "score",
    "seqnames",
    "start",
    "value",
    "variable",
    "color_code",
    ":=",
    "sample_status",
    "bed_path",
    "Frequency",
    "Feature"),
  package = "ALPS",
  add = FALSE
)
