getChrStartEnd <- function(x){

  x_chr <- gsub(":.*", "", x) %>% as.character()
  x_start <- gsub(".*:", "", x) %>% gsub("-.*", "", .) %>% as.numeric()
  x_end <- gsub(".*:", "", x) %>% gsub(".*-", "", .) %>% as.numeric()

  x_res <- list("chr" = x_chr, "from" = x_start, "to" = x_end)

  return(x_res)
}

#' UCSC Genome browser like plots
#' @description Function to plot genome browser like plots given a very minimal information such as a set of bigwig files and a genomics region.
#' Tracks are arranged as they are in given input data.frame. Function uses highly customizable Gviz R/bioconductor package to plot browser like plots.
#'
#' @param data_table a dataframe that contains \code{bw_path}, \code{sample_id} and \code{color_code}. Tracks are colored according to \code{color_code}. Default \code{NULL}
#' @param gene_range genomic region for which browser-like plots needed in format \code{chr:start-end}. Default \code{NULL}
#' @param ref_gen reference genome, to get gene annotations, currently supports  \code{hg19} and \code{hg38}. Default, \code{hg38}
#' @param cex.axis axis label size, default \code{0.5}
#' @param cex.title axis title size, default \code{0.8}
#' @param ... additional arguments to change the appearence of a plot. All arguments that can be passed to base R graphics are supported
#'
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom dplyr filter pull
#' @importFrom Gviz GeneRegionTrack DataTrack plotTracks
#'
#' @return plot of genome browser tracks
#'
#' @export
#'
#' @examples
#'
#' ## load example data
#'
#' chr21_data_table <- system.file("extdata/bw", "ALPS_example_datatable.txt", package = "ALPS", mustWork = TRUE)
#'
#' ## attach path to bw_path and bed_path
#' d_path <- dirname(chr21_data_table)
#'
#' chr21_data_table <- read.delim(chr21_data_table, header = TRUE)
#' chr21_data_table$bw_path <- paste0(d_path, "/", chr21_data_table$bw_path)
#' chr21_data_table$bed_path <- paste0(d_path, "/", chr21_data_table$bed_path)
#'
#' gene_range = "chr21:45643725-45942454"
#'
#' plot_browser_tracks(data_table = chr21_data_table,
#' gene_range = gene_range, ref_gen = "hg38")

plot_browser_tracks <- function(data_table=NULL, gene_range = NULL,
                                ref_gen = "hg38", cex.axis = 0.5,
                                cex.title = 0.8, ...){

  assertthat::assert_that(!is.null(data_table), msg = "Please provide `data_table`")
  assertthat::assert_that(assertthat::has_name(data_table, "bw_path"))
  assertthat::assert_that(assertthat::has_name(data_table, "sample_id"))
  assertthat::assert_that(assertthat::has_name(data_table, "color_code"))

  assertthat::assert_that(!is.null(gene_range), msg = "Please provide `gene_range` in format; chr:start-end")

  gene_range_split <- getChrStartEnd(x = gene_range)

  ## build genetrack
  if(ref_gen == "hg38"){

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  } else {

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }

  ##
  grtrack <- Gviz::GeneRegionTrack(txdb, genome = ref_gen,
                             chromosome = gene_range_split$chr,
                             window = -1, name = NULL,
                             background.title = "white",
                             col.frame = "white", col.axis = "black",
                             col = "black", col.title = "black",
                             geneSymbol=TRUE, showId=TRUE,
                             from = gene_range_split$from, to = gene_range_split$to)

  ## build enrichment tracks
  dt_list <- list()
  all_sample_ids <- data_table$sample_id %>% as.character()
  for(i in 1:length(all_sample_ids)){

    x_id <- all_sample_ids[i]
    x_dat <- data_table %>%
      dplyr::filter(sample_id == get("x_id"))

    x_bw <- x_dat %>% dplyr::pull(bw_path) %>% as.character()
    x_col <- x_dat %>% dplyr::pull(color_code) %>% as.character()

    x_dt <- Gviz::DataTrack(range = x_bw, genome = ref_gen, chromosome = gene_range_split$chr,
                      name = x_id, type = "hist", window=-1,
                      fill.histogram = x_col, col.histogram = "NA",
                      background.title = "white",
                      col.frame = "white", col.axis = "black",
                      col = "black", col.title = "black")

    dt_list[[i]] <- x_dt

  }

  gtrack_pos <- length(dt_list) + 1
  dt_list[[gtrack_pos]] <- grtrack

  Gviz::plotTracks(dt_list, from = gene_range_split$from, to = gene_range_split$to, cex.axis = cex.axis, cex.title = cex.title, ...)

}
