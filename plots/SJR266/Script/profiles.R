# the goal of this script is to load in line profile data from ImageJ and
# the fraction of the line that has signal for GFP.
# The data is the output from lineProfiler.ijm with an mCherry-FRB signal in
# ch1 and FKBP-GFP signal in ch2. The data is in tab-separated format in many
# files in the directory Data/. It is hierarchical, with the top level being
# the organelle, then the protein, then the condition. Files are generated from
# one cell, multiple rois and 3 channels for each.

library(tidyverse)
library(ggforce)
library(cowplot)

# Functions ----

# this function helps to assemble all the profile data
# all the data is in many tab-separated format files in the directory Data/
# load all into a data frame, append the filename as a column when doing so
read_profile <- function(file) {
  # if file is ends 3_data.txt, skip it
  if (grepl("3_data.txt$", file)) {
    return(NULL)
  }
  # read in data
  df <- read.delim(file)
  if (nrow(df) < 3) {
    return(NULL)
  }
  # deduce the parent folder name
  df$cond <- basename(dirname(file))
  # deduce the folder above that
  df$prot <- basename(dirname(dirname(file)))
  # deduce the folder above that
  df$organelle <- basename(dirname(dirname(dirname(file))))
  # add file name to data frame
  df$file <- basename(file)
  # filename is of the form foo_1_1_data.txt extract foo into column called cell
  df$cell <- sub("_([[:digit:]]+)_([[:digit:]]+)_data.txt$", "", df$file)
  # extract the 1st and 2nd regmatch from the filename - the roi and channel
  df$roi <- as.numeric(sub(".*_([[:digit:]]+)_([[:digit:]]+)_data.txt$", "\\1", df$file))
  df$ch <- as.numeric(sub(".*_([[:digit:]]+)_([[:digit:]]+)_data.txt$", "\\2", df$file))
  # add row number as a column
  df$row <- seq_len(nrow(df))
  # subtract BG from Intensity
  df$Intensity <- df$Intensity - df$BG
  # normalise using MX - BG
  df$Intensity <- df$Intensity / (df$MX - df$BG)
  return(df)
}

# Variables ----

# Nikon SoRa 60X, 2.8X objective, 0.065131666298805 um per pixel
pxsize <- 0.065131666298805
threshold <- 0.4

# Script ----

# list all csv files in subfolders of "Data" folder
files <- list.files("Data", pattern = "*data.txt",
                    recursive = TRUE, full.names = TRUE)

df_all <- data.frame()
for (i in seq_len(length(files))) {
  cat(paste0("\nFile: ", files[i], "\n"))
  # fetch the calculated data
  temp_df <- read_profile(files[i])
  if (!is.null(temp_df)) {
    df_all <- rbind(df_all, temp_df)
  }
}

# pivot_wider so that each channel is a column, drop the C and file column first
df_all <- df_all %>%
  select(-C, -file, -BG, -MX) %>%
  pivot_wider(names_from = ch, names_prefix = "ch", values_from = Intensity)

# find the fraction of Norm values for each file that are greater than 1
# fractionAboveThresh <- aggregate(ch2 ~ cell + cond + prot + organelle, data = df_all, FUN = function(x) sum(x > threshold) / length(x))
fractionAboveThresh <- aggregate(ch2 ~ roi + cell + cond + prot + organelle, data = df_all, FUN = function(x) sum(x > threshold) / length(x))

# take the rapa condition
sub_df <- fractionAboveThresh %>%
  filter(cond == "Rapa")

summary_df <- sub_df %>%
  group_by(prot, organelle, cond) %>%
  summarise(mean = mean(ch2), sd = sd(ch2))

organelle_list = c("LDs", "Endosomes", "Lysosomes", "Mitochondria")

for(i in organelle_list) {
  sub_df %>%
    filter(organelle == i) %>%
    ggplot(aes(x = prot, y = ch2)) +
    geom_sina(colour = "#00a651", size = 0.5) +
    geom_errorbar(data = summary_df %>% filter(organelle == i), aes(x = prot, y = mean, ymin = mean - sd, ymax = mean + sd), colour = "black", linewidth = 0.2, width = 0) +   #errorbar
    geom_crossbar(data = summary_df %>% filter(organelle == i), aes(x = prot, y = mean, ymin = mean, ymax = mean), colour = "black", linewidth = 0.4, size = 1 ) +   #mean line
    labs(x = "", y = "MCS Fraction") +
    lims(y = c(-0.15, 1.15)) +
    theme_cowplot(8)
    ggsave(paste0("Output/Plots/Fractions_", i, ".pdf"), width = 50, height = 40, units = "mm")
}

# ttest to compare prots at each Organelle
t.test(ch2 ~ prot, data = filter(sub_df, organelle == "LDs"))
t.test(ch2 ~ prot, data = filter(sub_df, organelle == "Endosomes"))
t.test(ch2 ~ prot, data = filter(sub_df, organelle == "Lysosomes"))
t.test(ch2 ~ prot, data = filter(sub_df, organelle == "Mitochondria"))

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}

F3D <- fractionAboveThresh %>%
  filter(organelle == "Mitochondria") %>%
  select(organelle, prot, ch2)
write.csv(F3D, "Output/Data/F3D.csv", row.names = FALSE)
F5B <- fractionAboveThresh %>%
  filter(organelle != "Mitochondria") %>%
  select(organelle, prot, ch2)
write.csv(F5B, "Output/Data/F5B.csv", row.names = FALSE)
