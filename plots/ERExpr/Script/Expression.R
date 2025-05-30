library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(cowplot)
library(patchwork)

process_data <- function() {
  all_files <- list.files("Data", pattern = "*.csv", full.names = TRUE)
  df_all <- data.frame()
  for (filename in all_files) {
    df_temp <- read.csv(filename)
    fname <- tools::file_path_sans_ext(basename(filename))
    fname <- gsub("^.*?_","",fname)
    df_temp$expt <- fname
    df_temp <- remove_imagej_rows_if_there(df_temp)
    df_all <- rbind(df_all, df_temp)
  }
  return(df_all)
}

remove_imagej_rows_if_there <- function(df) {
  testIf <- names(df)
  if(testIf[1] == "X") {
    df$X <- NULL
  }
  
  return(df)
}


# load in the data
data <- process_data()

# wrangle filenames
data$Label <- gsub("HCT116LBR", "HCT116_LBR", data$Label)
data$Protein <- sapply(strsplit(data$Label,"_"), "[",4)
data$Protein <- gsub("noanchor", "", data$Protein)
data$Protein <- gsub("withanchor", "", data$Protein)
data$Protein <- gsub("dish2", "", data$Protein)
data$Protein <- gsub("StarmChFRB", "", data$Protein)
# because FKBPGFPSec61b and CMVFKBPGFPSec61b are the same protein
data$Protein <- gsub("^CMVFKBPGFPSec61b", "FKBPGFPSec61b", data$Protein)
data$Protein <- factor(data$Protein, levels = c(
  "FKBPGFPSec61b", "PGKFKBPGFPSec61b", "cripCMVFKBPGFPSec61b",
  "LBRFKBPGFP", "LBRFKBPGFPpool5"))


## Plotting
p <- ggplot(data = data, aes(x = Protein, y = IntDen)) + 
  geom_sina(colour = "#00a651", size = 0.5) +
  geom_errorbar(data = data %>%
                  group_by(Protein) %>%
                  summarise(mean = mean(IntDen), sd = sd(IntDen)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black",linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = data %>%
                  group_by(Protein) %>%
                  summarise(mean = mean(IntDen), sd = sd(IntDen)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black",linewidth = 0.4, width = 1) +   #errorbar
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "FKBPGFPSec61b" = "CMV",
    "PGKFKBPGFPSec61b" = "PGK",
    "cripCMVFKBPGFPSec61b" = "crCMV",
    "LBRFKBPGFP" = "CMV",
    "LBRFKBPGFPpool5" = "Knock-in"
  )) +
  labs(x = "", y = "ER Fluorescence") +
  theme_cowplot(8) +
  theme(legend.position = "none")
p

# layout three copies of the plot so that it can be compiled with the other plots
r <- (p | p | p)

ggsave("Output/Plots/expression.pdf", r, width = 176, height = 58, units = "mm")

# find the file with the IntDen value that is closest to the mean so that we
# can use those images to make the figure
closest_to_mean <- function(df, mean_value) {
  df$diff <- abs(df$IntDen - mean_value)
  closest_row <- df[which.min(df$diff), ]
  return(closest_row)
}
# find the mean IntDen value for each protein
mean_values <- data %>%
  group_by(Protein) %>%
  summarise(mean = mean(IntDen))
file_closest_to_mean <- data.frame()
for (i in 1:nrow(mean_values)) {
  protein <- mean_values$Protein[i]
  mean_value <- mean_values$mean[i]
  df_protein <- data[data$Protein == protein, ]
  closest_row <- closest_to_mean(df_protein, mean_value)
  file_closest_to_mean <- rbind(file_closest_to_mean, closest_row)
}

# Stats ----

anova_density <- aov(IntDen ~ Protein, data = data)
summary(anova_density)
TukeyHSD(anova_density)

## Export results ----
# extract Protein and IntDen columns to a new dataframe and save as csv
results <- data %>%
  select(Protein, IntDen) %>%
  distinct() %>%
  arrange(Protein)
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}
write.csv(results, "Output/Data/S12D.csv", row.names = FALSE)
