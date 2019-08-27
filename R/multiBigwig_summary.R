#' Enrichments at genomics regions
#'
#' @description \code{multiBigwig_summary} is a function to calculate
#' enrichments from a set of given bigwig files
#' and a set of genomics regions.
#' This function is similar to
#' \code{deeptools multiBigwigSummary} python package.
#'
#' @param data_table a dataframe that contains \code{bw_path}, \code{sample_id}.
#' \code{sample_id} ids will be used in the final result matrix
#' @param summary_type whether to calculate mean, median, min or max
#' for each genomic region in within consensus peak-set
#' from the bigwig files. Default is \code{mean}
#' @param parallel logical. Whether to parallelize the calculation process,
#' default is \code{TRUE}
#'
#' @importFrom dplyr mutate select rename rename left_join
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer BigWigFileList summary
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr separate
#'
#' @return data.frame of enrichments within given genomic regions
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
#' enrichments <- multiBigwig_summary(data_table = chr21_data_table,
#'                                    summary_type = "mean",
#'                                    parallel = FALSE)

multiBigwig_summary <- function(data_table = NULL,
                                summary_type = "mean",
                                parallel = TRUE){

  assertthat::assert_that(!is.null(data_table), msg = "Please provide `data_table`")
  assertthat::assert_that(assertthat::has_name(data_table, "bw_path"))
  assertthat::assert_that(assertthat::has_name(data_table, "sample_id"))
  assertthat::assert_that(assertthat::has_name(data_table, "bed_path"))

  ## create a consensus peak-set from all files
  all_beds <- data_table$bed_path %>% as.character()
  peaks_gr <- merge_GR(x = all_beds)

  ## bw list
  all_bw_files <- data_table$bw_path %>% as.character()
  names(all_bw_files) <- data_table$sample_id %>% as.character()

  bwL <- rtracklayer::BigWigFileList(all_bw_files)

  if(parallel){

    .par_fun <- function(x){

      x_bw <- bwL[[x]]

      x_res <- rtracklayer::summary(x_bw, peaks_gr, type = summary_type) %>%
        as.data.frame() %>%
        dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
        tibble::column_to_rownames(var = "region") %>%
        dplyr::select(score) %>%
        dplyr::rename(!!x := score)

      return(x_res)
    }

    all_pid <- bwL %>% names

    all_pid_reslst <- parallel::mclapply(all_pid, .par_fun)
    count_mat <- all_pid_reslst %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "region") %>%
      tidyr::separate(region, into = c("chr", "start", "end"))

  } else {

    count_mat <- peaks_gr %>% as.data.frame() %>%
      dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
      dplyr::select(region)

    for(i in 1:length(bwL)){

      sample_id <- bwL[i] %>% names
      sample_path <- bwL[[i]]

      sample_res <- rtracklayer::summary(sample_path, peaks_gr, type = summary_type) %>%
        as.data.frame() %>%
        dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
        dplyr::select(region, score) %>%
        dplyr::rename(!!sample_id := score)

      suppressMessages(count_mat <- dplyr::left_join(count_mat, sample_res, by = "region"))
    }
    count_mat <- count_mat %>% tidyr::separate(region, into = c("chr", "start", "end"))
  }
  return(count_mat)
}


