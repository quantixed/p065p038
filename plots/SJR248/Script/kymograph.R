library(EBImage)
library(abind)

# In profiles.R we selected a df of profiles to plot, called df_selected
# from dthis we will take each cell and make an rgb matrix for each using
# Intensity values from ch1, ch2, ch3
# we will then use the rgb matrix to create an image and save it as a .tif file
# in the Output/Plots directory

# selected_rois contains the combinations that we want to generate images for
# we will use these to filter df_selected


for (i in seq_len(NROW(selected_rois))) {
  this_prot <- selected_rois$prot[i]
  this_cond <- selected_rois$cond[i]
  this_organelle <- selected_rois$organelle[i]
  # select the correct data
  df_cell <- df_selected %>%
    filter(prot == this_prot,
           cond == this_cond,
           organelle == this_organelle)
  # reorder the dataframe by rerow
  df_cell <- df_cell[order(df_cell$rerow), ]
  # create the rgb matrices
  black <- matrix(rep(0, nrow(df_cell) * 3),
                  nrow = nrow(df_cell), ncol = 3)
  green <- matrix(rep(df_cell$ch2, 3),
                  nrow = nrow(df_cell), ncol = 3)
  red <- matrix(rep(df_cell$ch1, 3),
                nrow = nrow(df_cell), ncol = 3)
  rgb_matrix <- matrix(c(df_cell$ch1,
                         df_cell$ch2,
                         df_cell$ch1, # rep(0, nrow(df_cell))),
                       nrow = nrow(df_cell), ncol = 3)
  # rgb_matrix is a matrix of integers, we need to convert it to a 3D array
  # so that we can use it to create an image
  black_array <- array(black, dim = c(nrow(df_cell), 1, 3))
  green_array <- array(green, dim = c(nrow(df_cell), 1, 3))
  red_array <- array(red, dim = c(nrow(df_cell), 1, 3))
  rgb_array <- array(rgb_matrix, dim = c(nrow(df_cell), 1, 3))
  rgb_array <- abind(black_array,
                     green_array, green_array, green_array,
                     black_array,
                     red_array, red_array, red_array,
                     black_array,
                     rgb_array, rgb_array, rgb_array,
                     black_array, along = 2)
  img <- Image(rgb_array, colormode = "Color")
  writeImage(img, paste0("Output/Plots/kymo_",
                         this_organelle, "_",
                         this_prot, "_",
                         this_cond, ".tif"))
}
