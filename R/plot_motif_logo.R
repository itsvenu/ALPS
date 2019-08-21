#' Process meme format
#' @param x path to meme format file
#'
#' @importFrom data.table fread
#' @return data.frame
#'
#' @export
#'
#' @examples
#' myc_meme <- system.file("extdata/motifs", "MA0147.2.meme", package = "ALPS", mustWork = TRUE)
#' myc_df <- process_meme(x = myc_meme)

process_meme <- function(x){

  suppressWarnings(x_rd <- data.table::fread(x, skip = 11))
  colnames(x_rd) <- c("A", "C", "G", "T")

  return(x_rd)
}

#' Process jaspar format
#' @param x path to jaspar format file
#'
#' @importFrom stringr str_squish str_split
#' @return data.frame
#'
#' @export
#'
#' @examples
#' myc_jaspar <- system.file("extdata/motifs", "MA0147.2.jaspar", package = "ALPS", mustWork = TRUE)
#' myc_df <- process_jaspar(x = myc_jaspar)


process_jaspar <- function(x){

  x_rd <- readLines(x) %>% .[-1] %>%
    gsub("*.\\[ ", "", .) %>%
    gsub("\\ ].*", "", .) %>%
    stringr::str_squish()

  x_rd2 <- x_rd %>% stringr::str_split(pattern = " ")

  x_df <- matrix(nrow = length(x_rd2[[1]]) - 1, ncol =4)

  for(i in 1:length(x_rd2)){

    x_df[, i] <- x_rd2[[i]] %>% .[-1] %>% as.numeric()

  }

  x_df <- x_df %>% as.data.frame() %>%
    dplyr::rename(A = "V1", C = "V2", G = "V3", T = "V4")

  x_df <- x_df/rowSums(x_df)

  return(x_df)
}

#' Process transfac format
#' @param x path to transfac format file
#'
#' @importFrom data.table fread
#' @return data.frame
#'
#' @export
#'
#' @examples
#' myc_transfac <- system.file("extdata/motifs", "MA0147.2.transfac", package = "ALPS", mustWork = TRUE)
#' my_df <- process_transfac(x = myc_transfac)

process_transfac <- function(x){

  suppressWarnings(x_rd <- data.table::fread(x, skip = 5))
  x_rd <- x_rd %>% dplyr::select(-c(PO))

  x_rd <- x_rd/rowSums(x_rd)

  return(x_rd)
}

#' Process homer format
#' @param x path to homer format file
#'
#' @importFrom data.table fread
#' @return data.frame
#'
#' @export
#'
#' @examples
#' myc_homer <- system.file("extdata/motifs", "cmyc.homer", package = "ALPS", mustWork = TRUE)
#' myc_df <- process_homer(x = myc_homer)

process_homer <- function(x){

  x_rd <- data.table::fread(x, skip = 1) %>%
    dplyr::rename(A = "V1", C = "V2", G = "V3", T = "V4")

  return(x_rd)
}

#' Process PFM format
#' @param x path to PFM format file
#'
#' @importFrom data.table fread
#' @return data.frame
#'
#' @export
#'
#' @examples
#' myc_pfm <- system.file("extdata/motifs", "MA0147.2.pfm", package = "ALPS", mustWork = TRUE)
#' myc_df <- process_pfm(x = myc_pfm)


process_pfm <- function(x){

  x_rd <- data.table::fread(x, skip = 1)

  x_rd <- x_rd %>% t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    dplyr::rename(A = "V1", C = "V2", G = "V3", T = "V4")

  x_rd <- x_rd/rowSums(x_rd)

  return(x_rd)
}

#' Plot sequence motifs
#'
#' @description Function to plot transscription factor motif logos in two different styles, \code{bar} plot or \code{logo} plot.
#' It supports motif formats from different databases e.g. \code{JASPAR}, \code{MEME}, \code{TRANSFAC}, \code{HOMER} and \code{PFM}
#'
#' @param motif_path path to motif file, default \code{NULL}
#' @param database database name from which motif has taken, default \code{NULL}
#' @param plot_type either \code{bar} or \code{logo}, default \code{bar}
#'
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @importFrom ggseqlogo ggseqlogo
#'
#' @return \code{ggplot2} object
#'
#' @export
#'
#' @examples
#'
#' ## examplr motif file paths
#' myc_meme <- system.file("extdata/motifs", "MA0147.2.meme", package = "ALPS", mustWork = TRUE)
#' myc_jaspar <- system.file("extdata/motifs", "MA0147.2.jaspar", package = "ALPS", mustWork = TRUE)
#' myc_transfac <- system.file("extdata/motifs", "MA0147.2.transfac", package = "ALPS", mustWork = TRUE)
#' myc_homer <- system.file("extdata/motifs", "cmyc.homer", package = "ALPS", mustWork = TRUE)
#' myc_pfm <- system.file("extdata/motifs", "MA0147.2.pfm", package = "ALPS", mustWork = TRUE)
#'
#' ## plot motifs
#' plot_motif_logo(motif_path = myc_homer, database = "homer", plot_type = "logo")

plot_motif_logo <- function(motif_path = NULL,
                            database = NULL,
                            plot_type = "bar"){

  assertthat::assert_that(file.exists(motif_path), msg = "File doesn't exist. Please check!")

  ## nt colors
  nt_cols <- c("#109648", "#255c99", "#F7B32B", "#D62839")
  names(nt_cols) <- c("A", "C", "G", "T")

  if(database == "jaspar"){

    motif_df <- process_jaspar(x = motif_path)

  } else if(database == "meme"){

    motif_df <- process_meme(x = motif_path)

  } else if(database == "transfac"){

    motif_df <- process_transfac(x = motif_path)

  } else if(database == "homer"){

    motif_df <- process_homer(x = motif_path)

  } else {

    motif_df <- process_pfm(x = motif_path)
  }

  ## plot
  if(plot_type == "bar"){

    motif_probs <- motif_df
    motif_probs$pos = as.numeric(as.character(rownames(motif_probs)))

    motif_probs$height <- apply(motif_probs[,c('A', 'C','G','T')], MARGIN=1,
                                FUN=function(x){2-sum(log(x^x,base=2))})

    motif_logo_df <- data.frame(A = motif_probs$A*motif_probs$height, C = motif_probs$C*motif_probs$height,
                                G = motif_probs$G*motif_probs$height, T = motif_probs$T*motif_probs$height,
                                pos = motif_probs$pos)


    motif_logo_mlt <- suppressMessages(motif_logo_df %>%
                                         dplyr::mutate(pos = as.factor(pos)) %>%
                                         reshape2::melt())

    motif_logo_mlt$variable <- factor(motif_logo_mlt$variable, levels = c("A", "C", "G", "T"))

    plt <- ggplot(data = motif_logo_mlt, aes(x = pos, y = value))  +
      geom_bar(aes(fill = variable), position='stack', stat='identity') +
      scale_fill_manual(values = nt_cols)+
      theme_classic(base_size = 18)+
      theme(axis.text = element_text(color = "black"),
            legend.title = element_blank(),
            axis.ticks.length = unit(.2, "cm"),
            axis.ticks = element_line(color = "black"))+
      xlab("Position") + ylab("Information")

  } else {

    ## logo
    plt <- ggseqlogo::ggseqlogo(t(motif_df))+
      theme_classic(base_size = 18)+
      theme(axis.text = element_text(color = "black"),
            legend.title = element_blank(),
            axis.ticks.length = unit(.2, "cm"),
            axis.ticks = element_line(color = "black"))+
      xlab("Position") + ylab("Bits")

  }

  return(plt)

}

