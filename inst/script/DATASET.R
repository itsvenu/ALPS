## datasets for ALPS package

## Following is a documentation of how the data
## for ALPS package is prepared


# Step - 1: Download data ----------------------------------------------------------------

## 1a) download data from TCGA GDC
## 1b) https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
## 1c) Download bigwig files related to entities ACC, BRCA, GBM, LGG & unzip all into a directory 'all_bw'
## 1d) Download entity specific peak sets: `All cancer type-specific peak sets``

# Step - 2: create example bed files ------------------------------------------------------

## 2a) creata a random region on chr21
## 2b) echo -e "chr21\t45639550\t46668243" > chr21_roi.bed
## 2c) create example entity specific peaksets (for ACC, BRCA, HGG, LGG)
##     using data from `Step 1d` and chr21_roi.bed from `Step 2b` using
##    `bedtools intersect` utility
##
## bash script for above task
## ls ../../../{ACCx,BRCA,GBMx,LGGx}_log2norm.txt | while read -r FILE; do OF=$(basename "$FILE" / | sed -e 's/_log2norm.txt//'); cat $FILE | cut -f1-3 | grep -v 'seqnames' | intersectBed -a - -b chr21_roi.bed > "$OF"_1.bed; cp "$OF"_1.bed "$OF"_2.bed; done

# Step - 3: create example bigwig files ---------------------------------------------------

## 3a) download the script to subset bigwig files: https://gist.github.com/itsvenu/70bbdaa9410aa87f63a223618d7e5e51
## 3b) subset 2 bigwig files from entities ACC, BRCA, GBM, LGG with chr21_roi.bed to create example bigwig files
##
## bash script for the whole task
## wget https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes
##
## ls all_bw/*{025FE5F8,000CFD9F,01112370,09C0DCE7,116B070A}* | grep -v '000CFD9F.*T1' | grep -v '01112370.*T2'  | while read -r BW;
## do
##
## OF1=$(basename "$BW" / | cut -d "_" -f1,2);
## OF2=$(basename "$BW" / | rev | cut -d "_" -f2 | rev);
## python subsetBW.py $BW "$OF1"_"$OF2".chr21.bw chr21 --headerToo;
## bigWigToBedGraph "$OF1"_"$OF2".chr21.bw "$OF1"_"$OF2".chr21.bg;
## cat "$OF1"_"$OF2".chr21.bg | intersectBed -a - -b chr21_roi.bed -wa > "$OF1"_"$OF2".roi.bed;
## bedGraphToBigWig "$OF1"_"$OF2".roi.bed hg38.chrom.sizes "$OF1"_"$OF2".bw;
## rm "$OF1"_"$OF2".chr21.bw "$OF1"_"$OF2".chr21.bg;
##
## done

# Step - 4: prepare sample_table for bigwig files ------------------------------------------

library(dplyr)

usethis::use_directory("inst/extdata/bw")

## copy above generated example bigwig and bed files to `inst/extdata/bw` directory

bw_files <- list.files(path = "inst/extdata/bw", pattern = "bw$")

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
bed_data = "inst/extdata/bw"

bed_data_list <- list.files(bed_data, pattern = "bed$")

data_table2 <- bed_data_list %>%
  as.data.frame() %>%
  dplyr::rename(bed_path = ".") %>%
  dplyr::mutate(sample_id = basename(as.character(bed_path))) %>%
  dplyr::mutate(sample_id = gsub("\\.bed", "", sample_id)) %>%
  dplyr::left_join(data_table, by = "sample_id") %>%
  dplyr::select(sample_id, group, color_code, bw_path, bed_path)

## use this table to prepare table with paths depending on the OS
write.table(data_table2, file = "inst/extdata/bw/ALPS_example_datatable.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Step - 5: get raw motif files -------------------------------------------

usethis::use_directory("inst/extdata/motifs")

## motif files are taken from respective databases & copied to `inst/extdata/motifs`
## homer: http://homer.ucsd.edu/homer/
## MEME, TRANSFAC, JASPAR, PFM: http://hocomoco11.autosome.ru/downloads_v11


# Step - 6: data for overlap violin plots ---------------------------------

usethis::use_directory("inst/extdata/overlap_violins")

xx_meta2 <- data.frame(sample_id = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"),
                       sample_status = c(rep("untreated", 4), c(rep("treated", 4))),
                       group = c("tf1", "tf2", "tf1", "tf2", "tf1", "tf2", "tf1", "tf2"))

xx_meta2$bw_path <- "path.bw"
xx_meta2$color_code <- c(rep("gray", 4), rep("red", 4))

xx_meta2 <- xx_meta2 %>%
  dplyr::select(bw_path, sample_id, sample_status, group, color_code)

## matrix = enrichments
xx2 <- matrix(runif(800), ncol=8) %>% as.data.frame()
colnames(xx2) <- xx_meta2$sample_id %>% as.character()
xx2 <- cbind(xx %>% dplyr::select(chr:end), xx2)

xx2$s5 <- xx2$s5 - 0.2
xx2$s7 <- xx2$s7 - 0.22

xx2$s6 <- xx2$s6 + 0.15
xx2$s8 <- xx2$s8 + 0.19

## write these to inst/extdata/overlap_violins
write.table(xx2, file = "inst/extdata/overlap_violins/enrichemnts_4_overlapviolins.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(xx_meta2, file = "inst/extdata/overlap_violins/data_table_4_overlapviolins.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## END



