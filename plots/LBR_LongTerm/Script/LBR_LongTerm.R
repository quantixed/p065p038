library(ggplot2)
library(ggdist)
library(alphashape3d)
library(ggbeeswarm)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyr)


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
  # add file name to data frame
  df$file <- basename(file)
  # add condition to data frame
  df$cond <- cond
  
  # convert to matrix
  mat <- as.matrix(df[,8:10])
  alphamat <- ashape3d(mat, 5)
  # plot(alphamat, col = c(2,4,1))
  # bg3d(color = "white")
  # par3d(windowRect = c(20, 30, 800, 800))
  # 
  # # take a picture of the xquartz window and save it
  # rgl.snapshot(paste0("Output/Plots/",cond,"_",iter,"_","3dplot.png"), fmt = "png")
  
  # process alpha shape
  selection <- alphamat$triang[, 9] >= 2
  triangles <- alphamat$triang[selection, c("tr1", "tr2", "tr3")]
  triangles <- as.data.frame(triangles)
  surfpoints <- unique(c(triangles$tr1, triangles$tr2, triangles$tr3))
  # filter df for the clusters that are on the surface of the alpha shape
  df <- df[df$Label %in% surfpoints, ]
  
  # assemble first data frame
  df1 <- data.frame(n = nrow(mat),
                    n_in = length(surfpoints))
  alphamat <- ashape3d(mat, 15)
  selection <- alphamat$triang[, 9] >= 2
  triangles <- alphamat$triang[selection, c("tr1", "tr2", "tr3")]
  triangles <- as.data.frame(triangles)
  triangles$area <- calculate_area(alphamat$x, triangles)
  # package each of these values into a data frame
  df2 <- data.frame(Vol = volume_ashape3d(alphamat),
                    Surf = sum(triangles$area))
  # add file name to data frame
  df2$file <- basename(file)
  # add condition to data frame
  df2$cond <- cond
  summary <- cbind(df1, df2)
  
  both <-  list(summary,df)
  
  return(both)
}


# Script ----

# list all csv files in subfolders of "Data" folder
files <- list.files("Data", pattern = "*.csv", recursive = TRUE, full.names = TRUE)

df_all <- data.frame()
df_spots <- data.frame()
for(i in 1:length(files)) {
  # get condition from file name
  cond <- strsplit(files[i], "_")[[1]][6]  
  if(basename(files[i]) == cond) {
    next
  }
  
  cat(paste0("\nFile: ",files[i],"\n"))
  # fetch the calculated data
  temp_obj <- calculations_3d(files[i], cond, i)
  temp_df <- temp_obj[[1]]
  if(!is.null(temp_df)) {
    df_all <- rbind(df_all, temp_df)
  }
  # fetch the spot data
  temp2_df <- temp_obj[[2]]
  if(!is.null(temp_df)) {
    df_spots <- rbind(df_spots, temp2_df)
  }
}

df_all$density <- df_all$n / df_all$Surf
df_all$density_in <- df_all$n_in / df_all$Surf

# from the file names, extract the cell name which is the 7th element
df_all$cell <- sapply(strsplit(df_all$file, "_"), function(x) x[7])

# plots ----
# ggplot of density
p1 <- ggplot(df_all, aes(x = cond, y = density_in)) + 
  geom_line(aes(group = cell), color = "#00a6513f") +
  geom_point(colour = "#00a651", size = 0.5) + #individual points spaced with density
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
  scale_y_continuous(limits = c(0,0.4)) +
  labs(x = "", y = expression(paste("Density (", mu, "m"^"-2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

# ggplot of surface area

p4 <- ggplot(data = df_spots, aes(x = cond, y = Surf..unit./2, fill = cond)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.5, slab_linewidth = NA) +
  scale_fill_manual(values = rep("#00a651",8)) +
  scale_y_log10() +
  labs(x = "", y = expression(paste("Cluster area (", mu, "m"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none") 

# summary data for all files
summary_df <- df_spots %>%
  group_by(cond, file) %>%
  summarise(medianVol = median(Vol..unit.), medianArea = median(Surf..unit./2))
# add the cell name column by parsing the file for the 7th element
summary_df$cell <- sapply(strsplit(summary_df$file, "_"), function(x) x[7])

p8 <- ggplot(data = summary_df, aes(x = cond, y = medianArea)) + 
  geom_line(aes(group = cell), color = "#00a6513f") +
  geom_point(colour = "#00a651", size = 0.5) + #individual points spaced with density
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
  scale_y_continuous(limits = c(0,0.8)) +
  labs(x = "", y = expression(paste("Median cluster area (", mu, "m"^"2", ")"))) +
  theme_cowplot(8) +
  theme(legend.position = "none")

r1 <- (p4 | p8 | p1)
ggsave("Output/Plots/all_plots_longterm.pdf", r1, width = 176, height = 58, units = "mm")

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}
# extract Protein and IntDen columns to a new dataframe and save as csv
S8B <- df_spots %>%
  select(cond, Surf..unit.) %>% 
  mutate(area = Surf..unit. / 2) %>%
  select(cond, area)
write.csv(S8B, "Output/Data/S8B.csv", row.names = FALSE)
S8C <- summary_df %>%
  select(cond, cell, medianArea)
write.csv(S8C, "Output/Data/S8C.csv", row.names = FALSE)
S8D <- df_all %>%
  select(cond, cell, density_in)
write.csv(S8D, "Output/Data/S8D.csv", row.names = FALSE)
