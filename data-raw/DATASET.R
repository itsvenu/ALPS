## datasets for ALPS package

## code to prepare `DATASET` dataset goes here

library(dplyr)

usethis::use_directory("inst/extdata/bw")

## copy chr21 bigwig files into extdata
# all_files <- list.files(path = "~/Desktop/work/TCGA_ATAC/bw/chr21_bw5/", pattern = "*bw", full.names = TRUE)
#
# for(i in 1:length(all_files)){
#
#   file.copy(from = all_files[i], to = "inst/extdata/bw")
# }

#chr21_tar <- "~/Desktop/work/TCGA_ATAC/bw/chr21_bw5/chr21.bw.tar.gz"

#file.copy(from = chr21_tar, to = "inst/extdata/bw/chr21.bw.tar.gz")

## data_table

bw_files <- list.files(path = "~/Desktop/work/TCGA_ATAC/bw/chr21_bw5", pattern = "bw$")

data_table <- bw_files %>%
  as.data.frame() %>%
  dplyr::rename(bw_path = ".")

group_info <- basename(bw_files) %>% gsub("_.*", "", .)

data_table$group <- group_info

group_info <- paste0(group_info, "_", rep(1:5, 12))

data_table$sample_id <- group_info

## add colors
data_table <- data_table %>%
  dplyr::mutate(color_code = dplyr::case_when(grepl("ACCx", group) ~ "#8DD3C7",
                                              grepl("BLCA", group) ~ "#FFFFB3",
                                              grepl("BRCA", group) ~ "#FFED6F",
                                              grepl("CESC", group) ~ "#BEBADA",
                                              grepl("CHOL", group) ~ "#FB8072",
                                              grepl("COAD", group) ~ "#80B1D3",
                                              grepl("ESCA", group) ~ "#FDB462",
                                              grepl("GBMx", group) ~ "#B3DE69",
                                              grepl("HNSC", group) ~ "#FCCDE5",
                                              grepl("KIRC", group) ~ "#D9D9D9",
                                              grepl("KIRP", group) ~ "#BC80BD",
                                              grepl("LGGx", group) ~ "#CCEBC5"))

## modify this data-table to include bed files
## bed files
bed_data = "/Users/venu/Desktop/work/TCGA_ATAC/bw/chr21_bed5"

bed_data_list <- list.files(bed_data)

data_table2 <- bed_data_list %>%
  as.data.frame() %>%
  dplyr::rename(bed_path = ".") %>%
  dplyr::mutate(sample_id = basename(as.character(bed_path))) %>%
  dplyr::mutate(sample_id = gsub("\\.bed", "", sample_id)) %>%
  dplyr::left_join(data_table, by = "sample_id") %>%
  dplyr::select(sample_id, group, color_code, bw_path, bed_path)

write.table(data_table2, file = "inst/extdata/bw/ALPS_example_datatable.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## sample set for genomic regions
file.copy(from = "../VISITOR/ACCx_specific_chr21.bed", to = "inst/extdata/bw/ALPS_ACCx_example_GenomicRegions.bed")

## subset above data_table to 20 samples to use as an example for `plot_browser_tracks`
set.seed(123)
data_table_bt <- data_table %>%
  dplyr::sample_n(size = 20) %>%
  dplyr::arrange(group)

write.table(data_table_bt, file = "inst/extdata/bw/ALPS_example_datatable_browsertracks.txt", sep = "\t", row.names = F, quote = F)

## motif files
usethis::use_directory("inst/extdata/motifs")

motif_files <- list.files(path = "..", full.names = TRUE)
motif_files <- motif_files[grepl("MA0147|cmyc", motif_files)]

for(i in 1:length(motif_files)){

  file.copy(from = motif_files[i], to = "inst/extdata/motifs")

}

## overlap violins example enrichments and metadata
enrichemnts_4_overlapviolins <- data.table::fread("../VISITOR/enrichments_overlapviolins_example.txt", header = TRUE)
data_table_4_overlapviolins <- read.delim("../VISITOR/enrichments_overlapviolins_example_metadata.txt", header = TRUE)

data_table_4_overlapviolins$bw_path <- "path.bw"
data_table_4_overlapviolins$color_code <- c(rep("gray", 4), rep("red", 4))

data_table_4_overlapviolins <- data_table_4_overlapviolins %>%
  dplyr::select(bw_path, sample_id, sample_status, group, color_code)

usethis::use_directory("inst/extdata/overlap_violins")

## write these to inst/extdata/overlap-violins
write.table(enrichemnts_4_overlapviolins, file = "inst/extdata/overlap_violins/enrichemnts_4_overlapviolins.txt", sep = "\t", row.names = F, quote = F)
write.table(data_table_4_overlapviolins, file = "inst/extdata/overlap_violins/data_table_4_overlapviolins.txt", sep = "\t", row.names = F, quote = F)

#usethis::use_data(enrichemnts_4_overlapviolins, overwrite = TRUE)
#usethis::use_data(data_table_4_overlapviolins, overwrite = TRUE)

## in future add, deeptools|ngs.plot|homer metagene example files



