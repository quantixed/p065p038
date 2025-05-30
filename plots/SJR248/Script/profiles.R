# the goal of this script is to load in line profile data from ImageJ and plot.
# The data is the output from lineProfiler.ijm with an mCherry-FRB signal in
# ch1 and FKBP-GFP signal in ch2. The data is in tab-separated format in many
# files in the directory Data/. It is hierarchical, with the top level being
# the organelle, then the protein, then the condition. Files are generated from
# one cell, multiple rois and 3 channels for each.

library(tidyverse)

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
  # filename is of the form foo_1_1_data.txt
  df$cell <- sub("_([0-9])_([0-9])_data.txt", "", df$file)
  # use last 12 characters of the filename as suffix
  suffix <- substr(df$file, nchar(df$file) - 15, nchar(df$file) - 4)
  # extract the first number from underscore separated string
  df$roi <- as.numeric(strsplit(suffix, "_")[[1]][2])
  # extract the second number from underscore separated string
  df$ch <- as.numeric(strsplit(suffix, "_")[[1]][3])
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

# Script ----

# list all csv files in subfolders of "Data" folder
files <- list.files("Data", pattern = "*.txt",
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

# now we need to work on sub-data frames for each roi
# we need to find the max for each channel and make a column for ch1 and ch2
# normalised to their respective max
df_all <- df_all %>%
  group_by(organelle, prot, cond, cell) %>%
  mutate(max_ch1 = max(ch1),
         max_ch2 = max(ch2)) %>%
  ungroup() %>%
  mutate(ch1 = ch1 / max_ch1,
         ch2 = ch2 / max_ch2)
# and also offset the row column so that the maximum for ch2 is at the midpoint
#of row, the row should be wrapped around so that values greater than the
# maximum row move to the beginning and values less than the minimum row
# move to the end
df_all <- df_all %>%
  group_by(organelle, prot, cond, cell, roi) %>%
  mutate(max_row = max(row)) %>%
  mutate(row = row - which.max(ch2) + median(row)) %>%
  mutate(rerow = ifelse(row > max_row, row - max_row, row)) %>%
  mutate(rerow = ifelse(rerow < 0, rerow + max_row, rerow)) %>%
  mutate(real = rerow * pxsize) %>%
  ungroup()

ggplot(df_all) +
  geom_line(aes(x = real, y = ch1), color = "#da70d6") +
  geom_line(aes(x = real, y = ch2), color = "#00A651") +
  theme_bw(8) +
  labs(y = "Intensity", x = "Length (µm)") +
  lims(y = c(0, 1)) +
  facet_grid(organelle + prot + cond ~ roi, scales = "free")
# ggsave("Output/Plots/profiles.png", width = 10, height = 20, dpi = 300)

ggplot(df_all) +
  geom_histogram(aes(x = log2(ch2 / ch1)), fill = "grey", bins = 60) +
  theme_bw(8) +
  labs(x = "Ch1 / Ch2 (Log2)") +
  facet_grid(prot ~ organelle + cond, scales = "free_y")
# ggsave("Output/Plots/histo.png", width = 10, height = 10, dpi = 300)

# for each organelle, make a plot of ch1 and ch2 vs real
for (orgnl in unique(df_all$organelle)) {
  df_all %>%
    filter(organelle == orgnl) %>%
    ggplot() +
    geom_line(aes(x = real, y = ch1), color = "#da70d6") +
    geom_line(aes(x = real, y = ch2), color = "#00A651") +
    theme_bw(8) +
    labs(y = "Intensity", x = "Length (µm)") +
    lims(y = c(0, 1)) +
    facet_grid(prot + cond ~ roi, scales = "free")
  # ggsave(paste0("Output/Plots/", orgnl, "_profiles.png"),
  #        width = 10, height = 10, dpi = 300)
}

# Figure ----
# for Figure panels we will use 1 roi for all each prot and cond
selected_rois <- data.frame(prot =
                              rep(rep(c("LBR", "Sec61"), each = 2),
                                  times = 4),
                            cond =
                              rep(c("Control", "Rapa"), times = 8),
                            organelle =
                              rep(c("Endosomes", "LDs", "Lysosomes", "Mitochondria"),
                                  each = 4),
                            roi = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 0, 2, 2, 2, 2, 2))
# filter df_all based on the selected_rois
df_selected <- df_all %>%
  inner_join(selected_rois, by = c("prot", "cond", "organelle", "roi"))

# for each organelle, make a plot of ch1 and ch2 vs real
for (orgnl in unique(df_selected$organelle)) {
  df_selected %>%
    filter(organelle == orgnl) %>%
    ggplot() +
    geom_line(aes(x = real, y = ch1), color = "#da70d6") +
    geom_line(aes(x = real, y = ch2), color = "#00A651") +
    theme_bw(8) +
    labs(y = "Intensity", x = "Length (µm)") +
    lims(y = c(0, 1)) +
    facet_wrap(~prot + cond, scales = "free")
  ggsave(paste0("Output/Plots/", orgnl, "_selprofiles.pdf"),
         width = 80, height = 80, units = "mm", dpi = 300)
}

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}

F3C <- df_selected %>%
  filter(organelle == "Mitochondria") %>%
  select(organelle, prot, cond, real, ch1, ch2)
write.csv(F3C, "Output/Data/F3C.csv", row.names = FALSE)
S9 <- df_selected %>%
  filter(organelle != "Mitochondria") %>%
  select(organelle, prot, cond, real, ch1, ch2)
write.csv(S9, "Output/Data/S9.csv", row.names = FALSE)
