library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(cowplot)

## Functions ----

scatter_dot_plot <- function(df, x, y, xlab = "", ylab = "", size = 8, filename = NULL) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_sina(colour = "#00a651", size = 0.5) +
    geom_errorbar(data = df %>%
                    group_by(.data[[x]]) %>%
                    summarise(mean = mean(.data[[y]]), sd = sd(.data[[y]])), 
                  aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                  colour = "black",linewidth = 0.2, width = 0) +   #errorbar
    geom_crossbar(data = df %>%
                    group_by(.data[[x]]) %>%
                    summarise(mean = mean(.data[[y]]), sd = sd(.data[[y]])), 
                  aes(y = mean, ymin = mean, ymax = mean),
                  colour = "black",linewidth = 0.4, width = 1) +   #errorbar
    labs(x = xlab, y = ylab) +
    theme_cowplot(size) +
    theme(legend.position = "none")
  if(!is.null(filename)) ggsave(filename, p, width = 3, height = 3, dpi = 300, bg = "white")
  return(p)
}


# in the folder Data/ we have many subfolders with csv files
# recusively search for all csv files and read them into a single dataframe
# the csv files have the same structure
filelist <- list.files(path = "Data/", pattern = ".csv", recursive = TRUE, full.names = TRUE)
df <- do.call(rbind, lapply(filelist, function(x) {
  temp <- read.csv(x); temp$file <- x; temp}))

# the column file contains the path to the csv file
# remove double slashes
df$file <- gsub("//", "/", df$file)
# we can extract the folder name from the path
df$folder <- tolower(dirname(df$file))

# if the folder contains the string "sec61" or "lbr" we will add this to the column prot
df$prot <- ifelse(grepl("sec61", df$folder), "Sec61", 
                  ifelse(grepl("lbr", df$folder), "LBR", NA))
# if the folder contains the string "interphase" we will add this to the column phase, anything else will be "mitosis"
df$phase <- ifelse(grepl("interphase", df$folder), "Interphase", "Mitosis")
# if the folder contains the string "cntl" we will add this to the column treatment, anything else will be "treated"
df$treatment <- ifelse(grepl("cntl", df$folder), "Control", "Rapamycin")
# concatenate prot and treatment
df$prot_treatment <- paste(df$prot, df$treatment, sep = "_")

sdf <- df %>%
  group_by(prot_treatment, prot, treatment, phase) %>%
  summarise(Mean = mean(Pearson.s.Coefficient),
            SD = sd(Pearson.s.Coefficient),
            N = n())

ggplot(data = df, aes(x = prot, y = Pearson.s.Coefficient, group = prot_treatment, colour = treatment, shape = treatment)) +
  geom_sina(size = 1) +
  geom_errorbar(data = sdf, mapping = aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, group = prot_treatment), position = position_dodge(width = 0.9), colour = "black", linewidth = 0.2, width = 0) +
  geom_crossbar(data = sdf, mapping = aes(y = Mean, ymin = Mean, ymax = Mean, group = prot_treatment),  position = position_dodge(width = 0.9), colour = "black", linewidth = 0.4, width = 0.4) +
  scale_colour_manual(values = c("#7f7f7f","#f59331")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "", y = "PCC") +
  facet_wrap(~phase) +
  theme_cowplot(9) +
  lims(y = c(0, 1)) +
  theme(legend.position = "none")
ggsave("Output/Plots/PCC.pdf", width = 3, height = 2, dpi = 300, bg = "white")

# do ANOVA and tukey HSD test on rapamycin treated cells comparing protein and phase
df_rap <- df %>%
  filter(treatment == "Rapamycin")
aov_rap <- aov(Pearson.s.Coefficient ~ prot * phase, data = df_rap)
tukey_rap <- TukeyHSD(aov_rap)
tukey_rap
