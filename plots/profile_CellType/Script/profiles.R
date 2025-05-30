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
  # if file is ends 2_data.txt, skip it
  if (grepl("2_data.txt$", file)) {
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
  df$organelle <- basename(dirname(dirname(file)))
  # deduce the folder above that
  df$cellline <- basename(dirname(dirname(dirname(file))))
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
# if there's a column called ch3, rename it to ch2
if ("ch3" %in% colnames(df_all)) {
  df_all <- df_all %>%
    rename(ch2 = ch3)
}
# remove the ch1 and ch2 columns that are NA
df_all <- df_all %>%
  filter(!is.na(ch1) & !is.na(ch2))


# now we need to work on sub-data frames for each roi
# we need to find the max for each channel and make a column for ch1 and ch2
# normalised to their respective max
df_all <- df_all %>%
  group_by(cellline, organelle, cond, cell) %>%
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
  group_by(cellline, organelle, cond, cell, roi) %>%
  mutate(max_row = max(row)) %>%
  mutate(row = row - which.max(ch2) + median(row)) %>%
  mutate(rerow = ifelse(row > max_row, row - max_row, row)) %>%
  mutate(rerow = ifelse(rerow < 0, rerow + max_row, rerow)) %>%
  mutate(real = rerow * pxsize) %>%
  ungroup()

# loop through all organelles and make a plot for each
celllines <- unique(df_all$cellline)
organelles <- unique(df_all$organelle)
for (cl in celllines) {
  for (orgnl in organelles) {
    df_all %>%
      filter(organelle == orgnl) %>% 
      filter(cellline == cl) %>%
      ggplot() +
      geom_line(aes(x = real, y = ch2), color = "#da70d6") +
      geom_line(aes(x = real, y = ch1), color = "#00A651") +
      theme_bw(8) +
      labs(y = "Intensity", x = "Length (µm)") +
      lims(y = c(0, 1)) +
      facet_grid(cond + cell ~ roi, scales = "free")
#    ggsave(paste0("Output/Plots/profiles_",cl, "_" , orgnl,".png"), width = 10, height = 5, dpi = 300)
  }
}


# Figure ----
# for Figure panels we will use 1 roi for all each prot and cond
selected_rois <- data.frame(cellline =
                              rep(c("Cos7", "HeLa", "RPE1"), each = 6),
                            cond =
                              rep(c("Cntl", "Rapa"), times = 9),
                            organelle =
                              rep(rep(c("Golgi", "Mitochondria", "PM"),
                                  each = 2), times = 3),
                            roi = c(0, 2, 3, 2, 2, 2,
                                  0, 3, 1, 2, 1, 2,
                                  2, 1, 6, 3, 2, 1),
                            cellno = c(1, 2, 2, 2, 2, 2,
                                    1, 1, 1, 1, 1, 1,
                                    2, 1, 1, 2, 1, 3))
# in the data frame above, the cell column refers to the cell that is referred to
# extract this value and add as a column to the selected_rois data frame
# to do this find the unique cell values for each cellline, organelle, cond
master_list <- df_all %>%
  select(cellline, organelle, cond, cell) %>%
  distinct() %>%
  arrange(cellline, organelle, cond)
# now we need a column that shows the instance of each cell value for combination of cellline, organelle and cond
# but keep the cell column as it is
labelled_list <- master_list %>%
  group_by(cellline, organelle, cond) %>%
  mutate(cellno = row_number()) %>%
  ungroup()
# join the master)list labelled_list
master_list <- master_list %>%
  left_join(labelled_list, by = c("cellline", "organelle", "cond", "cell"))

# now we can join the master_list to the selected_rois data frame
selected_rois <- merge(selected_rois, master_list, by = c("cellline", "organelle", "cond", "cellno"), all.x = TRUE)


# filter df_all based on the selected_rois
df_selected <- df_all %>%
  inner_join(selected_rois, by = c("cellline", "organelle", "cond", "roi", "cell"))
# filter df_selected for cond == Rapa
df_selected <- df_selected %>%
  filter(cond == "Rapa")

# for each organelle, make a plot of ch1 and ch2 vs real
for (cl in unique(df_selected$cellline)) {
    df_selected %>%
      filter(cellline == cl) %>%
      ggplot() +
      geom_line(aes(x = real, y = ch2), color = "#da70d6") +
      geom_line(aes(x = real, y = ch1), color = "#00A651") +
      theme_bw(8) +
      labs(y = "Intensity", x = "Length (µm)") +
      lims(x = c(0, 4), y = c(0, 1)) +
      facet_wrap(~organelle, scales = "free")
    ggsave(paste0("Output/Plots/", cl, "_selprofiles.pdf"),
           width = 90, height = 50, units = "mm", dpi = 300)
}

# offset real by -2 but only in Cos7 Golgi or PM
df_selected <- df_selected %>%
  mutate(real = ifelse(cellline == "Cos7" & organelle %in% c("Golgi", "PM"), real - 2, real))
# order of plots should be RPE1, Cos7, HeLa
df_selected$cellline <- factor(df_selected$cellline,
                                 levels = c("RPE1", "Cos7", "HeLa"))
# order of organelles should be PM, Mitochondria, Golgi
df_selected$organelle <- factor(df_selected$organelle,
                                  levels = c("PM", "Mitochondria", "Golgi"))


df_selected %>%
  ggplot() +
  geom_line(aes(x = real, y = ch2), color = "#da70d6") +
  geom_line(aes(x = real, y = ch1), color = "#00A651") +
  theme_bw(9) +
  labs(y = "Intensity", x = "Length (µm)") +
  lims(x = c(0, 4), y = c(0, 1)) +
  facet_wrap(cellline~organelle, scales = "free")
ggsave(paste0("Output/Plots/all_selprofiles.pdf"),
       width = 120, height = 100, units = "mm", dpi = 300)

