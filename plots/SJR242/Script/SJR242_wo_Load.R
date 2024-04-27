library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrastr)
library(scales)

# functions ----

# this function helps to assemble all the spot data
read_and_distinguish <- function(file) {
  # read in data
  df <- read.csv(file)
  if(nrow(df) < 2) {
    return(NULL)
  }
  # add file name to data frame
  df$file <- basename(file)
  
  return(df)
}

# wrapper to collate all the spot data
read_and_collate <- function(filelist) {
  df <- data.frame()
  for(i in 1:length(filelist)) {
    # fetch the spot data
    temp_df <- read_and_distinguish(filelist[i])
    if(!is.null(temp_df)) {
      df <- rbind(df, temp_df)
    }
  }
  return(df)
}

id_mito <- function(df_spots_cell, df_mito, i) {
  cx <- floor(df_spots_cell$CX..pix.[i])
  cy <- floor(df_spots_cell$CY..pix.[i])
  cz <- floor(df_spots_cell$CZ..pix.[i])
  # clip a 11x11x1 image around the centroid
  clip <- df_mito[df_mito$X > cx - 9 & df_mito$X <= cx + 9 & df_mito$Y > cy - 9 & df_mito$Y <= cy + 9 & df_mito$Z > cz - 1 & df_mito$Z <= cz + 2,]
  id <- NA
  if(nrow(clip) == 0) {
    id <- NA
  } else {
    id <- unique(clip$Label)
    if(length(id) > 1) {
      id <- find_nearest_mito(cx, cy, cz, clip)
    }
  }
  return(id)
}

# disambiguate the mitochondria ID for each contact
find_nearest_mito <- function(xx, yy, zz, df) {
  # calculate the distance between each mitochondria and the contact
  df$dist <- sqrt((df$X - xx)^2 + (df$Y - yy)^2 + (df$Z - zz)^2)
  # find the row with the minimum distance and return the corresponding mitochondria ID
  return(df$Label[which.min(df$dist)])
}

# Script ----

# load the data - all contacts
df_spots <- read.csv("Output/Data/Spots.csv")

# rename cells 21 81 and 101 as Control and the rest as Rapamycin
df_spots$cell_label <- ifelse(df_spots$cell == "cell21" | df_spots$cell == "cell81" | df_spots$cell == "cell101", "Control", "Rapamycin")

# get unique list of cell ids
cell_ids <- unique(df_spots$cell)

# for each mitochondrion, calculate the number of contacts and the total volume of contacts
# so this is contacts summarised by mitochondrion but there is no mitochondrial 3D info here yet
mito_summary <- df_spots %>%
  group_by(cond, cell, mito_id) %>%
  summarise(n = n(),
            total_vol = sum(Vol..unit.),
            total_sa = sum(Surf..unit. / 2))
# rename cells 21, 81 and 101 as Control and the others as Rapamycin
mito_summary$cell_label <- ifelse(mito_summary$cell == "cell21" | mito_summary$cell == "cell81" | mito_summary$cell == "cell101", "Control", "Rapamycin")

# discard rows where n > 100
mito_summary <- mito_summary[mito_summary$n <= 100,]
mito_summary <- mito_summary[!is.na(mito_summary$mito_id),]

# plot the total surface area of all contacts per mitochondrion
p5 <- ggplot(mito_summary, aes(x = cond, y = total_sa)) + 
  geom_quasirandom(colour = "#00a651", size = 0.5, alpha = 0.5) + #individual points spaced with density
  scale_y_log10(limits = c(1000,NA)) +
  facet_wrap(cell_label ~ cell) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Contact Area (nm"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_sa_mito_facet.png", p5, width = 8, height = 6, dpi = 300, bg = "white")

# load in the mitochondrial summary data from 3D Image Suite
mito_files <- c()
for(cell_id in cell_ids) {
  # find the corresponding mitochondria file
  mito_file <- list.files("Data", pattern = paste0("M_mito_stats_", cell_id, ".csv"), full.names = TRUE)
  if(length(mito_file) == 0) {
    next
  }
  # add mito_file to the list
  mito_files <- c(mito_files, mito_file)
}
# get all the data
df_mitos <- read_and_collate(mito_files)

df_mitos$file <- sapply(strsplit(df_mitos$file, "\\."), "[", 1)
df_mitos$cell <- sapply(strsplit(df_mitos$file, "_"), "[", 4)

# merge the data using cell and mito_id in mito_summary and Label in df_mito
# df_mito has the contact info mapped to mitos and the corresponding mito 3D info
df_mito <- merge(mito_summary, df_mitos, by.x = c("cell", "mito_id"), by.y = c("cell","Label"))

# plot the total surface area of all contacts per mitochondrion vs the mitochondrion surface area
p6 <- ggplot(df_mito, aes(x = Surf..unit., y = total_sa, colour = cell_label)) + 
  geom_point_rast(shape = ".", alpha = 0.5) + 
  scale_x_log10(limits = c(1e3,1.1e8)) +
  scale_y_log10(limits = c(1e3,5e6)) +
  scale_colour_manual(values = c("#7f7f7f","#00a651")) +
  facet_wrap(cell_label ~ cond) +
  labs(x = expression(paste("Mitochondria surface (nm"^"2", ")")), y = expression(paste("Contact surface area (nm"^"2", ")"))) +
  theme_cowplot(9) +
  coord_fixed() +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_vs_mito_sa.png", p6, width = 8, height = 6, dpi = 300, bg = "white")

# ggplot of surface area for each cond faceted by cell
p7 <- ggplot(df_mito, aes(x = cond, y = total_sa)) + 
  geom_boxplot(aes(group = cond)) +
  geom_quasirandom(colour = "#00a651", size = 0.2, alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ cell_label) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Total surface area (nm"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_per_mito_facet.png", p7, width = 8, height = 6, dpi = 300, bg = "white")


# ggplot of contact area per mitochondrion sa for each cond faceted
p8 <- ggplot(df_mito, aes(x = cond, y = total_sa / Surf..unit.)) + 
  geom_boxplot(aes(group = cond)) +
  geom_quasirandom(colour = "#00a651", size = 0.2, alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ cell_label) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Total contact area / Mitochondrion surface"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_per_mito_facet_norm.png", p8, width = 8, height = 6, dpi = 300, bg = "white")

# now redo the original plots after excluding the contacts that we dismissed above
# we will keep rows in df_spots where the cell and mito_id are in mito_summary
# concatenate the cell and mito_id columns in both dataframes
df_spots$cell_mito <- paste(df_spots$cell, df_spots$mito_id, sep = "_")
mito_summary$cell_mito <- paste(mito_summary$cell, mito_summary$mito_id, sep = "_")
# keep only the rows in df_spots where cell_mito is in mito_summary
unique_cell_mito <- unique(mito_summary$cell_mito)
df_spots <- df_spots[df_spots$cell_mito %in% unique_cell_mito,]

# ggplot of sa for each cond faceted by cell_label
p1 <- ggplot(df_spots, aes(x = cond, y = Surf..unit. / 2)) + 
  geom_boxplot(aes(group = cond)) +
  geom_quasirandom(colour = "#00a651", size = 0.2, alpha = 0.2) + #individual points spaced with density
  scale_y_log10(limits = c(1000,NA)) +
  facet_wrap(cell_label ~ cell) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Contact surface area (nm"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_sa_facet.png", p1, width = 8, height = 6, dpi = 300, bg = "white")


# regenerate the summary plots
df_summary <- df_spots %>%
  group_by(cond, cell_label) %>%
  summarise(n = n(),
            mean_vol = mean(Vol..unit.),
            median_vol = median(Vol..unit.),
            mean_sa = mean(Surf..unit. / 2),
            median_sa = median(Surf..unit. / 2))

p2 <- ggplot(df_summary, aes(x = cell_label, y = n, fill = as.factor(cond))) + 
  geom_col(position = "dodge") +
  labs(x = "", y = "Number of contacts") +
  theme_cowplot(8) +
  theme(legend.position = "none")

p3 <- ggplot(df_summary, aes(x = cell_label, y = median_vol, fill = as.factor(cond))) + 
  geom_col(position = "dodge") +
  labs(x = "", y = expression(paste("Median volume (nm"^"3", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

p4 <- ggplot(df_summary, aes(x = cell_label, y = median_sa, fill = as.factor(cond))) + 
  geom_col(position = "dodge") +
  labs(x = "", y = expression(paste("Median surface (nm"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

q1 <- p1 / (p2 | p3 | p4)
ggsave("Output/Plots/contacts_summary.png", q1, width = 190, height = 160, dpi = 300, units = "mm")

# this plot shows the "evolution" of contacts per mitochondrion as the search area is expanded - pretty pointless!
p5 <- ggplot(mito_summary, aes(x = cond, y = total_sa, group = mito_id, colour = as.factor(mito_id))) + 
  geom_path(alpha = 0.2) + 
  scale_y_log10(limits = c(1000,NA)) +
  facet_wrap(~ cell_label) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Surface area (nm"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_sa_mito_distance.png", p5, width = 8, height = 6, dpi = 300, bg = "white")


# calculate summaries per cell
cell_summary <- df_spots %>%
  group_by(cell_label, cell, cond,) %>%
  summarise(n = n(),
            mean_vol = mean(Vol..unit.),
            median_vol = median(Vol..unit.),
            mean_sa = mean(Surf..unit. / 2),
            median_sa = median(Surf..unit. / 2))
cell_summary_overlay <- cell_summary %>% 
  group_by(cell_label, cond) %>%
  summarise(mean_ol = mean(mean_sa),
            sd_ol = sd(mean_sa))

# ggplot showing the mean_sa with cell_label  this is showing the average contact size per cell (not total per mito per cell)
dodge <- position_dodge(width = 8)
p9 <- ggplot(cell_summary, aes(x = cond, y = mean_sa / 1e6, group = cell_label, colour = cell_label, shape = cell_label)) + 
  geom_quasirandom(dodge.width = 8, size = 1) +
  scale_colour_manual(values = c("#7f7f7f","#f59331")) +
  scale_shape_manual(values = c(1, 16)) +
  geom_errorbar(data = cell_summary_overlay, aes(y = mean_ol / 1e6, ymin = (mean_ol - sd_ol) / 1e6, ymax = (mean_ol + sd_ol) / 1e6), colour = "black", linewidth = 0.2, width = 0, position = dodge) +
  geom_crossbar(data = cell_summary_overlay, aes(y = mean_ol / 1e6, ymin = mean_ol / 1e6, ymax = mean_ol / 1e6), colour = "black", linewidth = 0.4, width = 6, position = dodge) +
  scale_y_continuous(limits = c(0, 0.035)) +
  labs(x = "Mitochondria-ER distance (nm)", y = expression(paste("Contact s.a. (",mu,"m"^"2", ")"))) +
  theme_cowplot(9) +
  theme(legend.position = c(0.01, 0.9), legend.title = element_blank(), legend.text = element_text(size = 9))
ggsave("Output/Plots/contacts_sa_cell_distance.png", p9, width = 8, height = 6, dpi = 300, bg = "white")

# redo p6 from above
p6 <- ggplot(df_mito, aes(x = Surf..unit. / 1e6, y = total_sa / 1e6, colour = cell_label)) + 
  geom_point_rast(shape = ".", alpha = 0.5) + 
  scale_x_log10(limits = c(1e-3,1.1e2), labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(limits = c(1e-3,5e0), labels = label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = c("#7f7f7f","#f59331")) +
  facet_wrap(cell_label ~ cond) +
  labs(x = expression(paste("Mitochondrion surface area (",mu,"m"^"2", ")")), y = expression(paste("Total contact surface area per mitochondrion (",mu,"m"^"2", ")"))) +
  theme_cowplot(9) +
  coord_fixed() +
  theme(legend.position = "none")
ggsave("Output/Plots/contacts_vs_mito_sa.png", p6, width = 8, height = 6, dpi = 300, bg = "white")


# for the figure
ggsave("Output/Plots/contacts_vs_mito_sa.pdf", p6, width = 8, height = 4.2, dpi = 300, bg = "white")
ggsave("Output/Plots/contacts_sa_cell_distance.pdf", p9, width = 69, height = 47, units = "mm", bg = "white")

## stats on mean contact s.a.
for(i in unique(cell_summary$cond)){
  temp <- cell_summary[cell_summary$cond == i, c("cell_label", "mean_sa")]
  outcome <- temp$mean_sa
  treatment <- temp$cell_label
  test1 <- t.test(outcome~treatment)
  print(paste("cond:", i, "p-value:", test1$p.value))
}
