#' Correlations among replicates/groups
#'
#' @description Function to calculate correlations of ChIP/ATAC-seq
#' enrichment among replicates/samples or groups.
#' The function is compatible with the output of
#' \code{\link{multiBigwig_summary}} or any custom data.frame
#' with the similar format.
#'
#' @param enrichments_df dataframe of enrichments, usually in the form of the output from function \code{\link{multiBigwig_summary}}, default \code{NULL}
#' @param log_transform logical, whether to log2 transform \code{enrichments_df}
#' @param method method to calculate correlation coefficient, one of \code{'pearson'} (default), \code{'spearman'} or \code{'kendall'}
#' @param plot_type whether to plot \code{replicate_level} correlations or \code{group_level} correlations. \code{replicate_level} plot represents the pairwise correlation values among all the samples/columns, where as the \code{group_level} plot is a paired plot of genomic regions after averaging of all samples in a group. Default \code{replicate_level}
#' @param select_var logical, whether to plot only a subset of columns from \code{enrichments_df}, if TRUE then colvar and rowvar must also be specified. Can also be used to change order of plotted variables.
#' @param colvar When \code{select_var} = TRUE, define the columns from \code{enrichments_df} that should be used as column variables in the correlation plot. Enter a single number, a range (e.g. 1:3), or selected columns (e.g. c(3,1,5))
#' @param rowvar When \code{select_var} = TRUE, define the columns from \code{enrichments_df} that should be used as row variables in the correlation plot. Enter a single number, a range (e.g. 1:3), or selected columns (e.g. c(3,1,5))
#' @param sample_metadata a data.frame required if \code{plot_type = 'group_level'}. The data.frame must contain columns \code{sample_id} and \code{group}
#' @param ... additional arguments either to \code{corrplot::corrplot} or \code{GGally::ggpairs} depending on arg \code{plot_type}
#'
#' @import corrplot
#' @importFrom stats cor
#' @importFrom dplyr mutate select filter pull rename
#' @importFrom tibble column_to_rownames
#' @importFrom GGally ggpairs
#'
#' @return corrplot or ggplot2 object
#'
#' @export
#'
#' @examples
#'
#' ## load example data
#' ## load example data
#'
#' chr21_data_table <- system.file('extdata/bw', 'ALPS_example_datatable.txt', package = 'ALPS', mustWork = TRUE)
#'
#' ## attach path to bw_path and bed_path
#' d_path <- dirname(chr21_data_table)
#'
#' chr21_data_table <- read.delim(chr21_data_table, header = TRUE)
#' chr21_data_table$bw_path <- paste0(d_path, '/', chr21_data_table$bw_path)
#' chr21_data_table$bed_path <- paste0(d_path, '/', chr21_data_table$bed_path)
#'
#' enrichments <- multiBigwig_summary(data_table = chr21_data_table,
#'                                    summary_type = 'mean',
#'                                    parallel = TRUE)
#'
#' ## replicate_level correlation plot
#' plot_correlation(enrichments_df = enrichments,
#' log_transform = TRUE, plot_type = 'replicate_level',
#' sample_metadata = chr21_data_table)
#'
#' ## group_level correlation plot
#' plot_correlation(enrichments_df = enrichments,
#' log_transform = TRUE, plot_type = 'group_level',
#' sample_metadata = chr21_data_table)

plot_correlation <- function(enrichments_df,
                             log_transform = TRUE,
                             method = "pearson",
                             plot_type = "replicate_level",
                             select_var=FALSE,
                             colvar,
                             rowvar,
                             sample_metadata, ...) {
  
  assertthat::assert_that(length(colnames(enrichments_df)[-1:-3]) == length(sample_metadata$sample_id),
                          msg = "colnames in `enrichments_df` & sample_id in `sample_metadata` are not matching!")
  
  assertthat::assert_that(assertthat::has_name(enrichments_df, "chr"))
  assertthat::assert_that(assertthat::has_name(enrichments_df, "start"))
  assertthat::assert_that(assertthat::has_name(enrichments_df, "end"))
  
  enrichments_df <- enrichments_df %>%
    dplyr::mutate(region = paste(chr, start, end, sep = "_")) %>%
    dplyr::select(-c(chr, start, end)) %>%
    tibble::column_to_rownames(var = "region")
  
  ## sanity check
  length(colnames(enrichments_df)) == length(sample_metadata$sample_id)
  
  if (plot_type == "replicate_level") {
    
    if (log_transform) {
      
      enrichments_mat <- log2(enrichments_df + 1)
      
    } else {
      
      enrichments_mat <- enrichments_df
      
    }
    
    if (select_var) {
      
      enrichments_cor <- cor(enrichments_mat[rowvar], enrichments_mat[colvar], method = method)
      
    } else {enrichments_cor <- cor(enrichments_mat, enrichments_mat, method = method)
    
    }
    
    corrplot::corrplot(enrichments_cor, tl.col="black", tl.srt=45,
                       cl.pos="b", addrect=2, ...)
    
  } else {
    
    assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))
    assertthat::assert_that(assertthat::has_name(sample_metadata, "group"))
    
    all_groups <- sample_metadata$group %>%
      as.character %>% unique()
    
    .mean_fun <- function(x) {
      
      x_sampleids <- sample_metadata %>%
        dplyr::filter(group == get("x")) %>%
        dplyr::pull(sample_id) %>%
        as.character()
      
      x_res <- enrichments_df %>%
        dplyr::select(dplyr::one_of(x_sampleids)) %>%
        rowMeans() %>% as.data.frame() %>%
        dplyr::rename(`:=`(!!x, "."))
      
      return(x_res)
    }
    
    group_df <- lapply(all_groups, .mean_fun) %>%
      as.data.frame()
    
    if (log_transform) {
      
      group_mat <- log2(group_df + 1)
    } else {
      
      group_mat <- group_df
    }
    
    ## paired plot
    GGally::ggpairs(group_mat, upper=list(continuous=GGally::wrap("cor", size=7)), ...) +
      theme_bw(base_size = 12) +
      theme(strip.background=element_rect(colour="white", fill="white"),
            strip.text = element_text(color="black", size=14),
            panel.grid.minor = element_blank())
  }
  
}
