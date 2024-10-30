## to run this script, data from zenodo is required

library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(cowplot)

# load all txt files in the Data folder into a big data frame, except log.txt
files <- list.files("Data/n2/", pattern = "*.txt", full.names = TRUE)
# remove Log.txt
files <- files[!basename(files) %in% "Log.txt"]

# load all into big data frame
allData <- do.call(rbind, lapply(files, function(x) {
  temp <- read.table(x, header = TRUE, fill = TRUE);
  temp$pos <- 1:nrow(temp);
  temp}))

# remove rows with NA
allData <- allData[complete.cases(allData), ]

# normalise intensity
allData$Norm <- (allData$Intensity - allData$BG) / (allData$MX - allData$BG)
allData$Signal <- allData$Intensity - allData$BG

# extract T and Norm to new data frame
allResults <- allData[, c("T", "Norm", "Signal", "pos")]

# load log file
logFile <- readLines("Data/n2/Log.txt")
lookup <- data.frame(file = logFile,
                     "T" = seq(1, length(logFile)))
# remove everything up to and including "LD446/" in the file column
lookup$file <- sub(".*LD446/", "", lookup$file)
# if file column contains "noTG", then set cond to ctrl; if it contains "withTG", then set cond to TG
lookup$cond <- ifelse(grepl("noTG", lookup$file), "ctrl", "TG")

# merge the lookup data frame with allResults
allResults <- merge(allResults, lookup, by = "T")

# save the dataframes we used for plotting
write.csv(allResults, "Output/Data/allResults_n2.csv", row.names = FALSE)
