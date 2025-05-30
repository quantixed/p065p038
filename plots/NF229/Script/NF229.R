# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("flowCore")
# BiocManager::install("flowWorkspace")
# BiocManager::install("ggcyto")

library(flowCore)
library(ggplot2)
library(ggcyto)
library(flowWorkspace)
library(data.table)
library(dplyr)
library(cowplot)

fcs_files <- list.files(path = "Data",pattern = ".fcs", full.names = TRUE)
print(fcs_files)
fs <- read.flowSet(files = unlist(fcs_files))
colnames(fs)

tflow <- fsApply(fs, exprs, simplify = F)  
tflow <- lapply(tflow, as.data.table)  
tflow <- rbindlist(tflow, idcol = "file")
tflow <- as.data.frame(tflow)

# filter for "clon 5" "pool" and "negative"
tflow <- tflow %>%
  filter(grepl("clon 5", file) | grepl("pool", file) | grepl("negative", file))

# ggplot of "B488-530_30-A" as histogram on a logicle scale
ggplot(tflow, aes(x = `B488-530/30-A`, colour = file)) +
  geom_density() +
  scale_x_logicle() +
  labs(x = "B488-530_30-A", y = "Count") +
  theme_cowplot(9)

# log transform B488-530/30-A to get a plot like FACSAria
tflow$trans <- log10(tflow$`B488-530/30-A` + 1)
# colours from the FACSAria output
my_colours <- c(rgb(2, 85, 170, maxColorValue = 255),
                rgb(255, 51, 52, maxColorValue = 255),
                rgb(131, 130, 1, maxColorValue = 255))
# they look bad though
my_colours <- c("#00a651", "#3f3f3f", "#2276b9")
# plot the log transformed data as a histogram with binwidth of 0.018
ggplot(tflow, aes(x = trans, colour = file, fill = file)) +
  geom_line(stat = "bin", binwidth = 0.018, direction = "mid") +
  geom_area(stat = "bin", binwidth = 0.018, position = "identity", colour = NA, alpha = 0.3) +
  # geom_histogram(binwidth = 0.018, position = "identity", colour = NA, alpha = 0.3) +
  labs(x = "B488-530_30-A", y = "Count") +
  scale_x_continuous(limits = c(0.6,5), labels = scales::label_math(10^.x)) +
  scale_colour_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  theme_cowplot(9) +
  # place legend inside the plot, top left
  theme(legend.position = c(0.1, 0.9),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  annotation_logticks(sides = "b")
# save the plot
ggsave("Output/Plots/NF229.pdf", width = 58, height = 58, units = "mm")
