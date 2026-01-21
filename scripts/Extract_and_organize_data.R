# Data extraction and formatting

# In the WCS_data folder,  
# file: 2017-2025 capture data for reference nov 19 2025.xlsx:
# metadata of bats from Vancouver and surrounding areas from WCS Canada. Hasn't been formatted. 

# file: Lillooet_Creston_metadata.csv
# metadata of bats from Lillooet and Creston areas from WCS Canada. Hasn't been formatted. 

# BC_MYLU_MYYU_genetics_SampleList_XS_2025-11-21.csv: 
# file with a list of sample names waiting for joining the metadata. 
# Capture date is included in case same batID has multiple dates.

library(tidyverse)
library(readxl)

sample_list <- read.csv("E:/BC_bats_genetics/data/BC_MYLU_MYYU_genetics_SampleList_XS_2025-11-21.csv")

South_coast_WCS_data <- readxl::read_xlsx("E:/BC_bats_genetics/data/WCS_data/2017-2025 capture data for reference nov 19 2025.xlsx",
                                        sheet = "2017-2025 CLEAN Nov 13 2025")

Lillooet_Creston_WCS_data <- read.csv("E:/BC_bats_genetics/data/WCS_data/Lillooet_Creston_metadata.csv")


# format variable names
South_coast_WCS_data <- janitor::clean_names(South_coast_WCS_data)

Lillooet_Creston_WCS_data <- janitor::clean_names(Lillooet_Creston_WCS_data)

# format variable types
South_coast_WCS_data$capture_time <- as.character(South_coast_WCS_data$capture_time)

South_coast_WCS_data$capture_date <- as.character(South_coast_WCS_data$capture_date)

South_coast_WCS_data$start_time <- as.character(South_coast_WCS_data$start_time)

South_coast_WCS_data$end_time <- as.character(South_coast_WCS_data$end_time)

South_coast_WCS_data$moon_phase <- as.character(South_coast_WCS_data$moon_phase)

South_coast_WCS_data$mass_g <- as.numeric(South_coast_WCS_data$mass_g)

South_coast_WCS_data$forearm_length_mm <- as.numeric(South_coast_WCS_data$forearm_length_mm)

South_coast_WCS_data$number_of_biopsies <- as.character(South_coast_WCS_data$number_of_biopsies)

# variables which South coast metadata has but Lillooet/Creston metadata does not
colnames(South_coast_WCS_data)[!colnames(South_coast_WCS_data) %in% colnames(Lillooet_Creston_WCS_data)]

# variables which Lillooet/Creston metadata has but in South coast metadata does not
colnames(Lillooet_Creston_WCS_data)[!colnames(Lillooet_Creston_WCS_data) %in% colnames(South_coast_WCS_data)]


# Correct batIDs in the sample list. The following batIDs have been corrected in the South_coast_WCS_data:
# "22AL0022"   "22AL0371"  "22AL0453"   "22AL0909"   "22BV00XX"  
# "22SL0130"   "22SL0836"  "22SL0843"   "22SL0887" 

correct_batID_df <- data.frame(old = c("22AL0022", "22AL0371", "22AL0453", "22AL0909", "22BV00XX",
                                      "22SL0130", "22SL0836", "22SL0843", "22SL0887"),
                              new = c("22AL0021", "22AL1059",  "22AL0003", "22AL0102", "22DEA19-069",
                                      "22SL0369", "22SL0846", "22SL0850", "22SL0901"))

# 22AL0371 may be 22AL1059 or 22AL0951

sample_list$batID[sample_list$batID %in% correct_batID_df$old] <- correct_batID_df$new #correct the wrong batIDs



# format and extract our sample data from the South_coast_WCS_data
## 1. append last two numbers of capture_year to bat_id
South_coast_WCS_data$bat_id <- with(South_coast_WCS_data, paste0(substr(capture_year, 3, 4), bat_id))

## 2. rename duplicated bat_id such that each bat_id is unique
South_coast_WCS_data$bat_id <- make.unique(South_coast_WCS_data$bat_id, sep = ".")

## 3. extract our sample data
South_coast_metadata <- South_coast_WCS_data[South_coast_WCS_data$bat_id %in% sample_list$batID,]


# For Lillooet Creston data, use area name as site name
Lillooet_Creston_WCS_data$site_name <- Lillooet_Creston_WCS_data$area

## 4. Lillooet_Creston_WCS_data only has our sample data so no need for data extraction
##    combine South_coast_metadata and the Lillooet_Creston_WCS_data

metadata_BC_MYLU_MYYU_genetics <- dplyr::bind_rows(South_coast_metadata, Lillooet_Creston_WCS_data)





# some investigations on the metadata
sample_size <- with(metadata_BC_MYLU_MYYU_genetics, table(field_species_id, metadata_BC_MYLU_MYYU_genetics$site_name))

# sum sample sizes per size
margin.table(sample_size, margin = 2)

# sum of sample sizes per species
margin.table(sample_size, margin = 1)

# total sample sizes
margin.table(sample_size)

# with(metadata_BC_MYLU_MYYU_genetics, table(field_species_id, capture_month, capture_year))








