#' Enrichment plots
#'
#' @description Function to plot enrichments from ChIP-seq/ATAC-seq at
#' genomics regions either as an individual groups or as
#' paired condition e.g untreated-treated
#'
#' @param enrichments_df enrichments at genomics regions from all samples, as in the format of output from \code{\link{multiBigwig_summary}}
#' @param log_transform logical. Should the data be \code{log2} transformed? Default is \code{TRUE}
#' @param plot_type either \code{separate} or \code{overlap}
#' @param sample_metadata metadata associated with the columns present in \code{enrichments_df}, information in this table will be used depending on the option in \code{plot_type}
#' @param box_alpha alpha/transparency to use for box plots, default 0.8
#' @param violin_alpha alpha/transparecny to use for violin plots, default 0.8
#' @param x_order ordering of levels on x-axis in resulting plot, default \code{NULL}
#' @param overlap_order ordering of overlaying if \code{plot_type = 'overlap'}. E.g. overlaying treatment data on top of untreated data, default \code{NULL}
#'
#' @importFrom dplyr mutate select rename left_join group_by summarise
#' @importFrom reshape2 melt
#' @importFrom gghalves geom_half_violin
#'
#' @return \code{ggplot2} object
#'
#' @export
#'
#' @examples
#'
#' ## load example data
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
#' ## plot_type == 'separate'
#' plot_enrichments(enrichments_df = enrichments, log_transform = TRUE,
#' plot_type = 'separate', sample_metadata = chr21_data_table)
#'
#' ## plot_type == 'overlap'
#' enrichemnts_4_overlapviolins <- system.file('extdata/overlap_violins', 'enrichemnts_4_overlapviolins.txt', package = 'ALPS', mustWork = TRUE)
#' enrichemnts_4_overlapviolins <- read.delim(enrichemnts_4_overlapviolins, header = TRUE)
#'
#' ## metadata associated with above enrichments
#' data_table_4_overlapviolins <- system.file('extdata/overlap_violins', 'data_table_4_overlapviolins.txt', package = 'ALPS', mustWork = TRUE)
#' data_table_4_overlapviolins <- read.delim(data_table_4_overlapviolins, header = TRUE)
#'
#' plot_enrichments(enrichments_df = enrichemnts_4_overlapviolins, log_transform = FALSE,
#' plot_type = 'overlap', sample_metadata = data_table_4_overlapviolins)

plot_enrichments <- function(enrichments_df = NULL,
                             log_transform = TRUE, plot_type = "separate",
                             sample_metadata, box_alpha = 0.8, violin_alpha = 0.8,
                             x_order = NULL, overlap_order = NULL) {

    assertthat::assert_that(assertthat::has_name(enrichments_df, "chr"))
    assertthat::assert_that(assertthat::has_name(enrichments_df, "start"))
    assertthat::assert_that(assertthat::has_name(enrichments_df, "end"))

    assertthat::assert_that(assertthat::has_name(sample_metadata, "group"))
    assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))

    ## melt, log_transform and join with
    ## metadata before plotting
    suppressMessages(enrichments_mlt <- enrichments_df %>%
                       dplyr::mutate(region = paste(chr, start, end, sep = "_")) %>%
                       dplyr::select(-c(chr, start, end)) %>%
                       reshape2::melt())

    ## log transform
    if (log_transform) {

        enrichments_lg <- enrichments_mlt %>%
            dplyr::mutate(plt_val = log2(value + 1))

    } else {

        enrichments_lg <- enrichments_mlt %>%
          dplyr::mutate(plt_val = value)

    }

    ## join with metadata
    suppressWarnings(enrichments_jn <- enrichments_lg %>%
                       dplyr::rename(sample_id = variable) %>%
                       dplyr::left_join(sample_metadata,
                                        by = "sample_id"))

    ## check color pals within `plot_type`
    ## condition

    if (plot_type == "separate") {

        # average of each region from all samples
        # of a group
        enrichments_mlt_jn_avg <- enrichments_jn %>%
          dplyr::select(region, group, plt_val) %>%
          dplyr::group_by(.dots = c("region", "group")) %>%
          dplyr::summarise(avg = mean(plt_val))

        ## orderings
        if (!is.null(x_order)) {

            enrichments_final <- enrichments_mlt_jn_avg$group <- factor(enrichments_mlt_jn_avg$group, levels = x_order)

        } else {

            enrichments_final <- enrichments_mlt_jn_avg

        }

        ## check col argument
        if ("color_code" %in% colnames(sample_metadata)) {

            col_pal <- sample_metadata$color_code %>% as.character() %>% unique()
            names(col_pal) <- sample_metadata$group %>% as.character() %>% unique()

            plt <- ggplot(enrichments_final, aes(x = group, y = avg, fill = group)) +
              gghalves::geom_half_violin(side = "r",
                                         position = position_nudge(x = 0.25, y = 0),
                                         adjust = 2, trim = FALSE,
                                         alpha = violin_alpha) +
              geom_boxplot(width = 0.2,
                           outlier.shape = NA,
                           alpha = box_alpha) +
              scale_fill_manual(values = col_pal)

        } else {

            plt <- ggplot(enrichments_final, aes(x = group, y = avg, fill = group)) +
              gghalves::geom_half_violin(side = "r",
                                         position = position_nudge(x = 0.25, y = 0),
                                         adjust = 2, trim = FALSE, alpha = violin_alpha) +
              geom_boxplot(width = 0.2, outlier.shape = NA, alpha = box_alpha)

        }

        plt <- plt + theme_classic(base_size = 14) +
          ylab("normalized read density") +
          xlab("") + theme(axis.text = element_text(color = "black"),
                           legend.position = "none",
                           axis.ticks.length = unit(0.2,"cm"),
                           axis.ticks = element_line(color = "black"))

    } else {

        assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_status"))
        assertthat::assert_that(assertthat::has_name(sample_metadata, "group"))

        ## summarise based on 3 features in this
        ## case
        enrichments_mlt_jn_avg <- enrichments_jn %>%
          dplyr::group_by(.dots = c("region", "sample_status", "group")) %>%
          dplyr::summarise(avg = mean(plt_val))

        ## orderings X
        if (!is.null(x_order)) {

            enrichments_final <- enrichments_mlt_jn_avg
            enrichments_final$group <- factor(enrichments_final$group, levels = x_order)

        } else {

            enrichments_final <- enrichments_mlt_jn_avg

        }

        ## orderings of violins
        if (!is.null(overlap_order)) {

            enrichments_finalx <- enrichments_final
            enrichments_finalx$sample_status <- factor(enrichments_final$sample_status, levels = overlap_order)

        } else {

            enrichments_finalx <- enrichments_final

        }

        ## check col_pal
        if ("color_code" %in% colnames(sample_metadata)) {

            col_pal <- sample_metadata$color_code %>% as.character() %>% unique()
            names(col_pal) <- sample_metadata$sample_status %>% as.character() %>% unique()

            plt <- ggplot(enrichments_finalx, aes(x = group, y = avg, fill = sample_status)) +
              gghalves::geom_half_violin(side = "r",
                                         position = position_nudge(x = 0.25, y = 0),
                                         adjust = 2, trim = FALSE, alpha = violin_alpha) +
              geom_boxplot(width = 0.2, outlier.shape = NA, alpha = box_alpha) +
              scale_fill_manual(values = col_pal)

        } else {

            plt <- ggplot(enrichments_finalx, aes(x = group, y = avg, fill = sample_status)) +
              gghalves::geom_half_violin(side = "r",
                                         position = position_nudge(x = 0.25, y = 0),
                                         adjust = 2, trim = FALSE, alpha = violin_alpha) +
              geom_boxplot(width = 0.2, outlier.shape = NA, alpha = box_alpha)

        }

        plt <- plt + theme_classic(base_size = 14) +
          ylab("normalized read density") +
          xlab("") +
          theme(axis.text = element_text(color = "black"),
                axis.ticks.length = unit(0.2,"cm"),
                axis.ticks = element_line(color = "black"))

    }

    return(plt)
}

