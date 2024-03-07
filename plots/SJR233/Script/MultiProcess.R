library(ggplot2)
library(ggdist)
library(alphashape3d)
library(ggbeeswarm)
library(cowplot)
library(dplyr)
library(patchwork)

# functions ----

# function to calculate area of a triangle in 3D space
heron <- function(a,b,c){
  s <- (a + b + c) / 2
  area <- sqrt(s*(s-a) * (s-b)*(s-c))
  return(area)
}

distance3d <- function(x1,y1,z1,x2,y2,z2){
  a <- (x1-x2)^2+(y1-y2)^2 + (z1-z2)^2
  d <- sqrt(a)
  return(d)
}

areatriangle3d <- function(x1,y1,z1,x2,y2,z2,x3,y3,z3){
  a <- distance3d(x1,y1,z1,x2,y2,z2)
  b <- distance3d(x2,y2,z2,x3,y3,z3)
  c <- distance3d(x3,y3,z3,x1,y1,z1)
  A <- heron(a,b,c)
  # print(paste("area of triangle is",A))
  return(A)
}

# this function calculates the area of each triangle
calculate_area <- function(x, triangles) {
  # x is a matrix of 3D coordinates
  # triangles is a data frame of triangle indices
  # calculate the area of each triangle
  
  a <- areatriangle3d(x[triangles$tr1, 1], x[triangles$tr1, 2], x[triangles$tr1, 3],
                      x[triangles$tr2, 1], x[triangles$tr2, 2], x[triangles$tr2, 3],
                      x[triangles$tr3, 1], x[triangles$tr3, 2], x[triangles$tr3, 3])
  return(a)
}

# this functions does all the calculations, prints the results and returns data frame
calculations_3d <- function(file, cond, iter) {
  # read in data
  df <- read.csv(file)
  if(nrow(df) < 3) {
    return(NULL)
  }
  # convert to matrix
  mat <- as.matrix(df[,8:10])
  alphamat <- ashape3d(mat, 5)
  plot(alphamat, col = c(2,4,1))
  bg3d(color = "white")
  par3d(windowRect = c(20, 30, 800, 800))
  
  # take a picture of the xquartz window and save it
  rgl.snapshot(paste0("Output/Plots/",cond,"_",iter,"_","3dplot.png"), fmt = "png")
  # volume of alpha shape
  cat(paste0("Volume: ",volume_ashape3d(alphamat),"\n"))
  # surface area of alpha shape
  selection <- alphamat$triang[, 9] >= 2
  triangles <- alphamat$triang[selection, c("tr1", "tr2", "tr3")]
  triangles <- as.data.frame(triangles)
  triangles$area <- calculate_area(alphamat$x, triangles)
  # sum the area of all triangles
  cat(paste0("Surface area: ",sum(triangles$area),"\n"))
  # number of spots that went in
  cat(paste0("Number of spots: ",nrow(mat),"\n"))
  # density of spots is the number of spots divided by the surface area
  cat(paste0("Density of spots: ",nrow(mat) / sum(triangles$area),"\n"))
  # how many spots are in the alpha shape? the number of unique triangle indices
  cat(paste0("Number of spots in alpha shape: ",length(unique(c(triangles$tr1, triangles$tr2, triangles$tr3))),"\n"))
  # density of the spots in alpha shape
  cat(paste0("Density of spots in alpha shape: ",length(unique(c(triangles$tr1, triangles$tr2, triangles$tr3))) / sum(triangles$area),"\n"))
  # package each of these values into a data frame
  df <- data.frame(Vol = volume_ashape3d(alphamat),
                   Surf = sum(triangles$area),
                   n = nrow(mat),
                   n_in = length(unique(c(triangles$tr1, triangles$tr2, triangles$tr3))))
  # add file name to data frame
  df$file <- basename(file)
  # add condition to data frame
  df$cond <- cond
  
  return(df)
}

# this function helps to assemble all the spot data
read_spots <- function(file, cond) {
  # read in data
  df <- read.csv(file)
  if(nrow(df) < 3) {
    return(NULL)
  }
  # add file name to data frame
  df$file <- basename(file)
  # add condition to data frame
  df$cond <- cond
  
  return(df)
}

# Script ----

# list all csv files in subfolders of "Data" folder
files <- list.files("Data", pattern = "*.csv", recursive = TRUE, full.names = TRUE)

df_all <- data.frame()
df_spots <- data.frame()
for(i in 1:length(files)) {
  # get condition from file name
  cond <- strsplit(files[i], "/")[[1]][2]  
  if(basename(files[i]) == cond) {
    next
  }
  cat(paste0("\nFile: ",files[i],"\n"))
  # fetch the calculated data
  temp_df <- calculations_3d(files[i], cond, i)
  if(!is.null(temp_df)) {
    df_all <- rbind(df_all, temp_df)
  }
  # fetch the spot data
  temp2_df <- read_spots(files[i], cond)
  if(!is.null(temp_df)) {
    df_spots <- rbind(df_spots, temp2_df)
  }
}

df_all$density <- df_all$n / df_all$Surf
df_all$density_in <- df_all$n_in / df_all$Surf


# plots ----
# ggplot of density
p1 <- ggplot(df_all, aes(x = cond)) + 
  geom_quasirandom(aes(y = density), colour = "#00a651", size = 0.5) + #individual points spaced with density
  geom_errorbar(data = df_all %>%
                  group_by(cond) %>%
                  summarise(mean = mean(density), sd = sd(density)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black", linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = df_all %>%
                  group_by(cond) %>%
                  summarise(mean = mean(density), sd = sd(density)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black", linewidth = 0.4, size = 1 ) +   #mean line
  scale_color_manual(values = c("#00a651", "#00a651")) +
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "", y = expression(paste("Density (", mu, "m"^"-2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

p2 <- ggplot(df_all, aes(x = cond)) + 
  geom_quasirandom(aes(y = density_in), colour = "#00a651", size = 0.5) + #individual points spaced with density
  geom_errorbar(data = df_all %>%
                  group_by(cond) %>%
                  summarise(mean = mean(density_in), sd = sd(density_in)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black", linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = df_all %>%
                  group_by(cond) %>%
                  summarise(mean = mean(density_in), sd = sd(density_in)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black", linewidth = 0.4, size = 1 ) +   #mean line
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "", y = expression(paste("Density (", mu, "m"^"-2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

r1 <- p1 | p2
ggsave("Output/Plots/contacts_density.png", r1, width = 85, height = 50, dpi = 300, units = "mm")

# ggplot of volume
p3 <- ggplot(data = df_spots, aes(x = cond, y = Vol..unit., fill = cond)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.5, slab_linewidth = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "", y = expression(paste("Volume (", mu, "m"^"3", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none") 

ggsave("Output/Plots/contacts_vol.png", p3, width = 3, height = 4, dpi = 300, bg = "white")

# ggplot of surface area

p4 <- ggplot(data = df_spots, aes(x = cond, y = Surf..unit./2, fill = cond)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.5, slab_linewidth = NA) +
  scale_fill_manual(values = rep("#00a651",2)) +
  scale_y_log10() +
  labs(x = "", y = expression(paste("Contact area (", mu, "m"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none") 

ggsave("Output/Plots/contacts_area.png", p4, width = 3, height = 6, dpi = 300, bg = "white")

# ggplot of volume faceted by file
p5 <- ggplot(data = df_spots, aes(x = 0, y = Vol..unit., fill = cond)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_linewidth = NA) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ file) +
  labs(x = "", y = expression(paste("Volume (", mu, "m"^"3", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none") 

ggsave("Output/Plots/contacts_vol_facet.png", p5, width = 8, height = 8, dpi = 300, bg = "white")

# ggplot of surface area faceted by file
p6 <- ggplot(data = df_spots, aes(x = 0, y = Surf..unit./2, fill = cond)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_linewidth = NA) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ file) +
  labs(x = "", y = expression(paste("Contact area (", mu, "m"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

ggsave("Output/Plots/contacts_area_facet.png", p6, width = 8, height = 8, dpi = 300, bg = "white")

# summary data for all files
summary_df <- df_spots %>%
  group_by(cond, file) %>%
  summarise(medianVol = median(Vol..unit.), medianArea = median(Surf..unit./2))
p7 <- ggplot(data = summary_df, aes(x = cond, y = medianVol)) + 
  geom_quasirandom(colour = "#00a651", size = 0.5) + #individual points spaced with density
  geom_errorbar(data = summary_df %>%
                  group_by(cond) %>%
                  summarise(mean = mean(medianVol), sd = sd(medianVol)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black",linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = summary_df %>%
                  group_by(cond) %>%
                  summarise(mean = mean(medianVol), sd = sd(medianVol)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black",linewidth = 0.4, size = 1 ) +   #mean line
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "", y = expression(paste("Median volume (", mu, "m"^"3", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

p8 <- ggplot(data = summary_df, aes(x = cond, y = medianArea)) + 
  geom_quasirandom(colour = "#00a651", size = 0.5) + #individual points spaced with density
  geom_errorbar(data = summary_df %>%
                  group_by(cond) %>%
                  summarise(mean = mean(medianArea), sd = sd(medianArea)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black",linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = summary_df %>%
                  group_by(cond) %>%
                  summarise(mean = mean(medianArea), sd = sd(medianArea)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black",linewidth = 0.4, width = 1) +   #errorbar
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "", y = expression(paste("Median contact area (", mu, "m"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

r2 <- p7 | p8
ggsave("Output/Plots/contact_vol_area_per_cell.png", r2, width = 85, height = 50, dpi = 300, units = "mm")

# we will use the density measure (not density_in) this is the total number of spots divided by the surface area = p1
# we will use the median contact area per cell = p8
# p4 shows all contact data (surface area)

r3 <- (p4 | p8) / (plot_spacer() | p1)
ggsave("Output/Plots/all_plots.pdf", r3, width = 85, height = 90, dpi = 300, units = "mm")

## stats ----

# t test of medianArea with cond as factor from summary_df
t.test(medianArea ~ cond, data = summary_df)
# t test of density with cond as factor from summary_df
t.test(density ~ cond, data = df_all)
# print number of spots in each condition from df_spots
df_spots %>%
  group_by(cond) %>%
  summarise(n = n())
# print number of files in each condition from df_all
df_all %>%
  group_by(cond) %>%
  summarise(n = n())

## Utility ----

# We want to inspect individual files to see if the alpha shape is working correctly
# This function plots the alpha shape and the points that went into it
plot_ashape <- function(file) {
  df <- read.csv(file, header = TRUE)
  # convert to matrix
  mat <- as.matrix(df[,8:10])
  alphamat <- ashape3d(mat, 2)
  plot(alphamat, col = c(2,4,1))
  par3d(windowRect = c(20, 30, 800, 800))
  bg3d(color = "white")
}

inspect_ashape <- function(filelist, number) {
  mat <- read.csv(filelist[number], header = TRUE)
  # get condition from file name
  cond <- strsplit(files[number], "/")[[1]][2]  
  if(basename(files[number]) == cond) {
    next
  }
  cat(paste0("\nFile: ",filelist[number],"\n"))
  plot_ashape(filelist[number])
}

# inspect_ashape(files, 1)

# if a disc is 0.5 um^2 then its diameter is:
sqrt(0.5/pi)*2
