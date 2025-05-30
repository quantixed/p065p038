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
  df$prot <- basename(dirname(file))
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

# add a column called cond which is Control if the cell contains neither "cntlnorapa" or "-RAP"
# and it is Rapa otherwise
df_all <- df_all %>%
  mutate(cond = ifelse(grepl("cntlnorapa", cell), "Control", "Rapa")) %>%
  mutate(cond = ifelse(grepl("-RAP", cell), "Control", cond))

# now we need to work on sub-data frames for each roi
# we need to find the max for each channel and make a column for ch1 and ch2
# normalised to their respective max
df_all <- df_all %>%
  group_by(prot) %>%
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
  group_by(prot, cell) %>%
  mutate(max_row = max(row)) %>%
  mutate(row = row - which.max(ch2) + median(row)) %>%
  mutate(rerow = ifelse(row > max_row, row - max_row, row)) %>%
  mutate(rerow = ifelse(rerow < 0, rerow + max_row, rerow)) %>%
  mutate(real = rerow * pxsize) %>%
  ungroup()

# make a graph
df_all %>%
  ggplot() +
  geom_line(aes(x = real, y = ch1), color = "#da70d6") +
  geom_line(aes(x = real, y = ch2), color = "#00A651") +
  theme_bw(8) +
  labs(y = "Intensity", x = "Length (µm)") +
  lims(x = c(18,24), y = c(0, 1)) +
  facet_grid(cond ~ prot, scales = "free")
# ggsave(paste0("Output/Plots/region_profiles.png"), width = 10, height = 5, dpi = 300)

# reorder prots
df_all$prot <- factor(df_all$prot,
                       levels = c("LBRWTCterm_mitosis",
                                  "LBRN547D_mitosis",
                                  "LBRR583Q_mitosis",
                                  "LBRWTNterm_mitosis",
                                  "LBR1_288NS_mitosis",
                                  "LBR1_245S_mitosis"))

my_labels <-  as_labeller(
  c("LBRWTCterm_mitosis" = "LBR-FKBP-GFP",
  "LBRN547D_mitosis" = "LBR(N547D)-FKBP-GFP",
  "LBRR583Q_mitosis" = "LBR(R583Q)-FKBP-GFP",
  "LBRWTNterm_mitosis" = "GFP-FKBP-LBR",
  "LBR1_288NS_mitosis" = "GFP-FKBP-LBR(1-288)",
  "LBR1_245S_mitosis" = "FKBP-LBR(1-245)-GFP",
  "Control" = "Control",
  "Rapa" = "Rapamycin"))

# make a series of graphs for each prot
df_all %>%
  mutate(real = real - 19) %>% 
  ggplot() +
  geom_line(aes(x = real, y = ch1), color = "#da70d6") +
  geom_line(aes(x = real, y = ch2), color = "#00A651") +
  theme_bw(8) +
  labs(y = "Intensity", x = "Length (µm)") +
  lims(x = c(0, 5), y = c(0, 1)) +
  facet_grid(cond ~ prot, labeller = my_labels)
ggsave(paste0("Output/Plots/region_profiles.pdf"), width = 170, height = 60, units = "mm")



