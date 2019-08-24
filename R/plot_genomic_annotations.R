#' Plot genomic annotations
#'
#' @description plot the annotations of genomic regions
#' either as \code{stacked bar} or \code{heatmap}.
#' The function takes the output of \code{\link{get_genomic_annotations}}
#' directly or it is also compatible with a similar data.frame.
#'
#' @param annotations_df a data.frame of genomic annotations.
#' It can either be of one sample or of multiple samples as a data.frame
#' @param plot_type either \code{bar} plot or \code{heatmap},
#' default \code{bar}
#' @param col vector of colors for each feature in \code{annotations_df}.
#' If provided these colors are used, if not a custom set of distinct
#' colors will be used. Default \code{NULL}
#'
#' @importFrom reshape2 melt
#' @importFrom utils head
#'
#' @return \code{ggplot2} plot
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
#' g_annotations <- get_genomic_annotations(data_table = chr21_data_table,
#' merge_level = "group_level")
#'
#' plot_genomic_annotations(annotations_df = g_annotations, plot_type = "heatmap")

plot_genomic_annotations <- function(annotations_df = NULL,
                                     plot_type = "bar",
                                     col = NULL){

  ## colors
  if(!is.null(col)){

    col_pal <- col
  } else {

    dist_cols <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8",
                   "#f58231", "#911eb4", "#42d4f4", "#f032e6",
                   "#bfef45", "#fabebe", "#469990", "#e6beff",
                   "#9A6324", "#fffac8", "#800000", "#aaffc3",
                   "#808000", "#ffd8b1", "#000075", "#a9a9a9",
                   "#ffffff", "#000000")

    feature_names <- annotations_df$Feature %>%
      as.character() %>% unique

    num_features <- feature_names %>% length

    col_pal <- dist_cols %>% utils::head(n = num_features)
    names(col_pal) <- feature_names
  }

  ## single barplot
  if(ncol(annotations_df) == 2){

    message("  `annotation_df` has only 2 columns, plotting a barplot")

    annotations_stack <- annotations_df %>%
      dplyr::mutate(group = "group")

    plt <- ggplot(annotations_stack, aes(x = group, y = Frequency, fill = Feature))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values = col_pal)+
      theme_classic(base_size = 14)+
      scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
      theme(axis.title.x = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank())+
      ylab("% of peaks in feature")

  } else {

    ## if > 1 samples in annotations_df
    suppressMessages(annotations_mlt <- annotations_df %>% reshape2::melt())

    if(plot_type == "bar"){

      plt <- ggplot(annotations_mlt, aes(x = variable, y = value, fill = Feature))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = col_pal)+
        theme_classic(base_size = 14)+
        scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
        theme(axis.title.x = element_blank(),
              axis.text = element_text(color = "black"),
              axis.ticks.length = unit(.25, "cm"),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank())+
        ylab("% of peaks in feature")


    } else{
      ## heatmap
      plt <- ggplot(annotations_mlt, aes(x = variable, y = Feature, fill = value))+
        geom_tile(size = 0.65, color = "white")+
        theme_minimal(base_size = 14)+
        theme(axis.ticks = element_blank())+
        scale_fill_viridis_c(name = "% of peaks")+
        theme(axis.text = element_text(color = "black"),
              panel.grid = element_blank())+
        xlab("")
    }
  }

  return(plt)
}
