library(ggplot2)
library(ggdist)
library(ggbeeswarm)
library(cowplot)
library(dplyr)
library(patchwork)
library(foreach)
library(doParallel)

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

# list all ER contact csv files in "Data" folder
er_files <- list.files("Data", pattern = "^M_stats", recursive = TRUE, full.names = TRUE)
# read into a big data frame
df_spots <- read_and_collate(er_files)

# filenames are things like "M_stats_10_cell21.csv"
# we want to extract the condition from the filename
# we do this by splitting the filename by "_"
# and then taking the second element of the resulting vector
# filename without extension
df_spots$file <- sapply(strsplit(df_spots$file, "\\."), "[", 1)
df_spots$cond <- sapply(strsplit(df_spots$file, "_"), "[", 3)
df_spots$cell <- sapply(strsplit(df_spots$file, "_"), "[", 4)

# we will just look at a subset of the cells
df_spots <- df_spots[df_spots$cell %in% c("cell7", "cell21", "cell31", "cell41", "cell81", "cell91", "cell101"),]

# next we will figure out proximity of each contact to each mitochondrial ID
# do this by identifying the cell and then loading the corresponding mitochondria data
# then we can calculate the distance between each contact and each mitochondria
# and then find the minimum distance for each contact
# for each cell id in df_spots, find the corresponding mitochondria file. Mitochondria files are csvs with the name V_mito_cellXX.csv
# get unique list of cell ids
cell_ids <- unique(df_spots$cell)

#setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# for each cell id, load the corresponding mitochondria data
df_spots_all <- data.frame()

for(cell_id in cell_ids) {
  # find the corresponding mitochondria file
  mito_file <- list.files("Data/", pattern = paste0("V_mito_vx_", cell_id, ".csv"), full.names = TRUE)
  if(length(mito_file) == 0) {
    next
  }
  # read in the mitochondria data
  df_mito <- read.csv(mito_file)
  # create new df which is a subset of df_spots for the current cell
  df_spots_cell <- df_spots[df_spots$cell == cell_id,]
  spots <- nrow(df_spots_cell)
  df_spots_cell$mito_id <- foreach(i=1:spots, .combine = rbind) %dopar% {
    id <- id_mito(df_spots_cell, df_mito, i)
  }
  df_spots_all <- rbind(df_spots_all, df_spots_cell)
}
#stop cluster
stopCluster(cl)

# free some memory
rm(df_mito)
# rename the final column of df_spots_all as mito_id
colnames(df_spots_all)[ncol(df_spots_all)] <- "mito_id"
# replace df_spots with df_spots all
df_spots <- df_spots_all
rm(df_spots_all)
# save the data
write.csv(df_spots, "Output/Data/Spots.csv", row.names = FALSE)


