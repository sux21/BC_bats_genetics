# record allele sizes

source("scripts/helper_functions.R")

# load FSA data. All data will store in a list once the scripts are working.

osiris_out <- read.table("output/osiris_out/microsat_data/fragment_analysis_mobix_2024-02-02_plate1_1.tab"
                         , sep = "\t"
                         , header = TRUE
                         , na.strings = "")

# make a subset for one sample
osiris_subset <- subset(osiris_out, File.Name == "E05_0725-05--M2_020_5759")

df <- split(osiris_subset, seq_len(nrow(osiris_subset)))

allele_4channel <- lapply(seq_along(df), function (i) {chr_to_num(df[[i]]$Allele)})

intensity_4channel <- lapply(seq_along(df), function(i) {chr_to_num(df[[i]]$RFU)})


strong_signal_index <- lapply(seq_along(intensity_4channel), function(i) {which(intensity_4channel[[i]] >= 27000)}) 
# return integer index which intensity >= 27000

strong_signal <- lapply(seq_along(allele_4channel), function(i) {allele_4channel[[i]][strong_signal_index[[i]]]})
# return peaks which intensity >= 27000

index <- vector(mode = "list", length = 4)

# continue to work on this. We have the index. Need to find out how to put in the a list
for (i in seq_along(strong_signal)) {
  for (j in seq_along(strong_signal[[i]])) {
    for (k in seq_along(strong_signal)[-i]) {
      pos <- which(allele_4channel[[k]] >= (strong_signal[[i]][j] - 0.5) & allele_4channel[[k]] <= (strong_signal[[i]][j] + 0.5))
      index[[k]] <- append(NULL, pos)
      print(paste0(i,",",j,",",k))
      print(pos)
    }
  }
}
# return integer index of pull-up peaks in each channel

which(allele_4channel[[3]] >= (strong_signal[[4]][3] - 0.5) & allele_4channel[[3]] <= (strong_signal[[4]][3] + 0.5))
# function 1: remove pull-up signals

# function 2: average split peaks





