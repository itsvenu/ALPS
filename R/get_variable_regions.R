#' Get variable regions
#' @description given a data.frame of genomic regions enrichments, the function returns the number of
#' variable regions across all samples. The resulting matrix can be directly used with \code{PCAtools} or \code{ComplexHeatmap}
#' for further downstream explorative analysis e.g. unsupervised clustering
#'
#' @param enrichments_df a data.frame o enrichments at genomic regions. Output of
#' \code{multiBigwig_summary} or a similar format is compatible
#' @param log_transform logical, whether to log2 transform the counts, default \code{TRUE}
#' @param scale logical, whether to \code{scale} the variable regions before returning the results, default \code{TRUE}
#' @param num_regions number of variable regions to return, default \code{500}
#'
#' @importFrom genefilter rowVars
#' @importFrom dplyr mutate select
#' @importFrom tibble column_to_rownames
#'
#' @return a data.frame of scaled variable regions
#'
#' @export
#'
#' @examples
#'
#' mat <- matrix(sample.int(15, 9*100, TRUE), 9, 100) %>% as.data.frame()
#' mat <- mat %>%
#' tibble::rowid_to_column(var = 'start') %>%
#' dplyr::mutate(end = start + 1000) %>%
#' dplyr::mutate(chr = 'chr1') %>%
#' dplyr::select(chr, start, end, dplyr::everything())
#'
#' get_variable_regions(enrichments_df = mat, num_regions = 50)

get_variable_regions <- function(enrichments_df,
                                 log_transform = TRUE,
                                 scale = TRUE,
                                 num_regions = 500) {


    enrichments_fmt <- enrichments_df %>%
      dplyr::mutate(region = paste(chr, start, end, sep = "_")) %>%
      dplyr::select(-c(chr:end)) %>%
      tibble::column_to_rownames(var = "region")

    ## log transform
    if (log_transform) {

        enrichments_lg <- log2(enrichments_fmt + 1)

    } else {

        enrichments_lg <- enrichments_fmt

    }

    ## select var regions
    rv <- genefilter::rowVars(enrichments_lg)
    idx <- order(-rv)[1:num_regions]
    enrichments_var <- enrichments_lg[idx,]

    ## scaling
    if (scale) {

        enrichments_sc = t(apply(enrichments_var, 1, scale))
        colnames(enrichments_sc) <- colnames(enrichments_var)

    } else {

        enrichments_sc <- enrichments_var

    }

    return(enrichments_sc)

}
