#' Merge genomic regions
#'
#' @description merge overlaping genomic regions from multiple peak/bed files
#' @param x a character vector of bed file paths
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#'
#' @return GRanges object
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
#' x <- as.character(chr21_data_table$bed_path)
#'
#' merge_GR(x = x)

merge_GR <- function(x){

  bed_df <- NULL

  for(i in x){

    x_df <- data.table::fread(i) %>% unique()
    bed_df <- rbind(bed_df, x_df)

  }

  bed_df <- unique(bed_df)

  bed_gr <- GenomicRanges::makeGRangesFromDataFrame(df = bed_df,
                                                    seqnames.field = "V1",
                                                    start.field = "V2",
                                                    end.field = "V3",
                                                    keep.extra.columns = FALSE)
  bed_merged <- GenomicRanges::reduce(bed_gr)

  return(bed_merged)
}


#' Annotate genomic regions
#'
#' @description annotate genomic regions by simultaneosuly merging overlaping regions and
#' preparing consensus set of genomic regions from multiple files
#'
#' @param data_table a data.frame of sample table, as is for \code{multiBigwig_summary} input table, default NULL
#' @param ref_gen reference genome, either \code{hg38 or hg19}, default \code{hg38}
#' @param tss_region bp ± TSS to define promoter regions
#' @param merge_level either \code{all, group_level} or \code{none}.
#' \code{all} prepares a consensus set of peaks from all files present in \code{data_table}.
#' \code{group_level} prepares consensus set of peaks from all samples present each group separately.
#' \code{none} does not prepare any consensus peak set, annotates each peak file separately. Default \code{all}
#'
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom dplyr filter pull
#' @importFrom ChIPseeker annotatePeak
#' @importFrom data.table dcast
#' @import org.Hs.eg.db
#'
#' @return a data.frame of annotations and percentages
#'
#' @export
#'
#' @examples
#' chr21_data_table <- system.file("extdata/bw", "ALPS_example_datatable.txt", package = "ALPS", mustWork = TRUE)
#'
#' ## attach path to bw_path and bed_path
#' d_path <- dirname(chr21_data_table)
#'
#' chr21_data_table <- read.delim(chr21_data_table, header = TRUE)
#' chr21_data_table$bw_path <- paste0(d_path, "/", chr21_data_table$bw_path)
#' chr21_data_table$bed_path <- paste0(d_path, "/", chr21_data_table$bed_path)
#'
#' get_genomic_annotations(data_table = chr21_data_table,
#' merge_level = "group_level")

get_genomic_annotations <- function(data_table = NULL,
                                    ref_gen = "hg38",
                                    tss_region = c(-1000, 1000),
                                    merge_level = "all"){

  ## check ref_gen
  if(ref_gen == "hg38"){

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else {

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }

  ## check input_data
  assertthat::assert_that(!is.null(data_table), msg = "Please provide `data_table` as a data.frame bed paths! See examples!")
  assertthat::assert_that(assertthat::has_name(data_table, "bed_path"))

  ## prepare annotatePeaks input here depending on option from `merge_level`

  if(merge_level == "all"){

    all_bed <- data_table$bed_path %>% as.character()
    all_bed_gr <- merge_GR(x = all_bed)

  } else if(merge_level == "group_level"){

    assertthat::assert_that(assertthat::has_name(data_table, "group"))

    all_groups <- data_table$group %>% as.character() %>% unique

    all_bed_gr <- list()

    for(i in all_groups){

      i_beds <- data_table %>%
        dplyr::filter(group == get("i")) %>%
        dplyr::pull(bed_path) %>% as.character()

      i_beds_gr <- merge_GR(x = i_beds)

      all_bed_gr[[i]] <- i_beds_gr
    }

  } else {

    assertthat::assert_that(assertthat::has_name(data_table, "sample_id"))

    all_pid <- data_table$sample_id %>% as.character()

    all_bed_gr <- data_table$bed_path %>% as.character()
    names(all_bed_gr) <- all_pid
  }

  ## if not list obj, run without lapply
  if(!is.list(all_bed_gr) && !is.character(all_bed_gr)){

    suppressMessages(all_bed_ap <- ChIPseeker::annotatePeak(peak = all_bed_gr,
                                                            tssRegion = tss_region,
                                                            TxDb = txdb, annoDb = "org.Hs.eg.db",
                                                            verbose = FALSE))

    return_dat <- all_bed_ap@annoStat %>% as.data.frame() %>%
      tibble::remove_rownames()

  } else {

    all_bed_ap <- lapply(all_bed_gr, ChIPseeker::annotatePeak,
                         TxDb = txdb, tssRegion = tss_region,
                         verbose=FALSE)

    return_dat <- lapply(all_bed_ap, function(x) x@annoStat %>%
                           as.data.frame() %>%
                           tibble::remove_rownames()) %>%
      plyr::ldply() %>%
      data.table::dcast(Feature ~ .id, value.var = "Frequency")
  }

  return(return_dat)
}
