## datasets for ALPS package

## code to prepare `DATASET` dataset goes here

library(dplyr)

usethis::use_directory("inst/extdata/bw")

bw_files <- list.files(path = "~/Desktop/work/TCGA_ATAC/bw/ALPS_example/example_bw", pattern = "bw$")

data_table <- bw_files %>%
  as.data.frame() %>%
  dplyr::rename(bw_path = ".")

group_info <- basename(bw_files) %>% gsub("_.*", "", .)

data_table$group <- group_info

group_info <- paste0(group_info, "_", rep(1:2, 4))

data_table$sample_id <- group_info

## add colors
data_table <- data_table %>%
  dplyr::mutate(color_code = dplyr::case_when(grepl("ACCx", group) ~ "#8DD3C7",
                                              grepl("BRCA", group) ~ "#FFED6F",
                                              grepl("GBMx", group) ~ "#B3DE69",
                                              grepl("LGGx", group) ~ "#CCEBC5"))

## modify this data-table to include bed files
## bed files
bed_data = "~/Desktop/work/TCGA_ATAC/bw/ALPS_example/example_bed"

bed_data_list <- list.files(bed_data)

data_table2 <- bed_data_list %>%
  as.data.frame() %>%
  dplyr::rename(bed_path = ".") %>%
  dplyr::mutate(sample_id = basename(as.character(bed_path))) %>%
  dplyr::mutate(sample_id = gsub("\\.bed", "", sample_id)) %>%
  dplyr::left_join(data_table, by = "sample_id") %>%
  dplyr::select(sample_id, group, color_code, bw_path, bed_path)

## copy bw and beds to extdata/bw
for(i in bw_files){

  file.copy(from = paste0("~/Desktop/work/TCGA_ATAC/bw/ALPS_example/example_bw/", i), to = "inst/extdata/bw")
}

##
for(i in bed_data_list){

  file.copy(from = paste0("~/Desktop/work/TCGA_ATAC/bw/ALPS_example/example_bed/", i), to = "inst/extdata/bw")
}

## use this table to prepare table with paths depending on the OS
write.table(data_table2, file = "inst/extdata/bw/ALPS_example_datatable.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
write.table(enrichemnts_4_overlapviolins, file = "inst/extdata/overlap_violins/enrichemnts_4_overlapviolins.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data_table_4_overlapviolins, file = "inst/extdata/overlap_violins/data_table_4_overlapviolins.txt", sep = "\t", row.names = FALSE, quote = FALSE)




