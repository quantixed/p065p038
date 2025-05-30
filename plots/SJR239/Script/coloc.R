library(tidyverse)
library(ggbeeswarm)
library(viridisLite)
library(patchwork)

# function ----

# Script ----
df <- read.csv("Data/allResults2.csv")
# path column has "foo/bar/buzz" make new columns for each level
df$rep <- sapply(strsplit(df$path, "/"), "[", 1)
df$cond <- sapply(strsplit(df$path, "/"), "[", 2)
df$cell <- sapply(strsplit(df$path, "/"), "[", 3)
df$file <- sapply(strsplit(df$path, "/"), "[", 4)
# now we need to classify the data according to cell column and file column
# first does the cell column start with Int or Mit
df$cell_type <- ifelse(grepl("Int", df$cell), "Interphase", "Metaphase")
# now does the file column contain "after" or not
df$cell_state <- ifelse(grepl("after|rapa", df$file), "Post", "Pre")
# finally parse whether the file column is "ch1" or "ch2"
df$cell_channel <- ifelse(grepl("ch1", df$file), "ch1", "ch2")
# make an id column by concatenating the cond, rep and cell_type columns
df$id <- paste(df$cond, df$rep, df$cell, sep = "_")
# remove rows where the cond name does not begin with LBR
df <- df[grepl("^LBR", df$cond), ]
# each id should have four rows (ch1, ch2, pre and post) remove any instances where this is not the case
summary <- as.data.frame(table(df$id, df$cell_channel, df$cell_state))
duds <- unique(as.character(summary[summary$Freq != 1, 1]))
# df <- df[!df$id %in% duds, ]
# this reveals one id where we have more than one post condition - manually remove this
# if removed from analysis step, this line will become redundant
df <- df[!(df$id == "LBRFKBPGFP6DG5MAPStarFRB_LD347n6_IntCell1" & grepl("600", df$path)),]
# # remove files where the bleach was too strong so detection was poor
# duds <- df$id[df$cell_state == "post" & df$cell_channel == "ch2" & df$obj2 < 10 & grepl("6DG5",df$id)]
# df <- df[!df$id %in% duds, ]
# # remove files where there ER signal was detected as clusters pre in ch2
# duds <- df$id[df$cell_state == "pre" & df$cell_channel == "ch2" & df$obj > 50]
# df <- df[!df$id %in% duds, ]

# calculation then replace NaN with 0
df$prop1 <- df$ch12 / df$obj1
df$prop2 <- df$ch21 / df$obj2
df$prop1[is.nan(df$prop1)] <- 0
df$prop2[is.nan(df$prop2)] <- 0

# we want 6DG5 only
df <- df[grepl("6DG5", df$id), ]

p1 <- ggplot(data = df,
             aes(x = factor(cell_state, levels=c("Pre", "Post")),
                 y = obj1,
                 color = prop1)) +
  geom_quasirandom(size = 0.5) +
  scale_colour_gradientn(colours = viridis(256, option = "D"), limits = c(0,1), oob = scales::squish) +
  facet_grid(. ~ cell_type) +
  labs(x = "MAPPER", y = "Total clusters") +
  theme_bw(9) +
  theme(legend.position = "none")
p2 <- ggplot(data = df,
             aes(x = factor(cell_state, levels=c("Pre", "Post")),
                 y = obj2,
                 color = prop2)) +
  geom_quasirandom(size = 0.5) +
  scale_colour_gradientn(colours = viridis(256, option = "D"), limits = c(0,1), oob = scales::squish) +
  facet_grid(. ~ cell_type) +
  labs(x = "LBR", y = "") +
  theme_bw(9) +
  theme(legend.position = "none")
q1 <- p1 + p2 + theme(legend.position = "right", legend.title = element_blank())
q1 + plot_layout(guides = "collect")
ggsave("Output/Plots/coloc.pdf", q1, width = 98, height = 60, units = "mm")

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}

F2C <- df %>% 
  select(cell_type,cell_state,obj2,prop2)
write.csv(F2C, "Output/Data/F2C.csv", row.names = FALSE)
