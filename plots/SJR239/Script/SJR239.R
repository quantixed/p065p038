library(tidyverse)

# function ----
print_pval <- function(tukey, group1, group2) {
  comp <- paste(group1, group2, sep = "-")
  all <- row.names(tukey$group)
  if (!comp %in% all) {
    comp <- paste(group2, group1, sep = "-")
  }
  pval <- tukey$group[comp, 4]
  cat(paste("p-value for", group1, "vs", group2, "is", pval, "\n"))
}

multi_t_test <- function(data, x, y, val) {
  if(length(x) != length(y)) {
    print("Unequal lengths of x and y")
    return(NULL)
  }
  
  for(i in 1:length(x)) {
    xx <- data[data$group == x[i], val]
    yy <- data[data$group == y[i], val]
    result <- t.test(x = xx, y = yy, paired = TRUE, alternative = "two.sided")
    temp <- data.frame(group1 = x[i], group2 = y[i], pval = result$p.value)
    if(i == 1) {
      out <- temp
    } else {
      out <- rbind(out, temp)
    }
  }
  # order out by pval
  out <- out[order(out$pval), ]
  out$hb <- 0.05 / (nrow(out) - (1:nrow(out) - 1))
  out$res <- ifelse(out$pval < out$hb, "Yes", "No")
  # put back in original order
  out <- out[order(match(out$group1, x)), ]
  return(out)
}

# Script ----
df <- read.csv("Data/allResults.csv")
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
# remove files where the bleach was too strong so detection was poor
duds <- df$id[df$cell_state == "Post" & df$cell_channel == "ch2" & df$obj < 10 & grepl("6DG5",df$id)]
df <- df[!df$id %in% duds, ]
# remove files where there ER signal was detected as clusters pre in ch2
duds <- df$id[df$cell_state == "Pre" & df$cell_channel == "ch2" & df$obj > 50]
df <- df[!df$id %in% duds, ]

facet_names <- as_labeller(c("ch1" = "mScarlet-I3",
                             "ch2" = "GFP",
                             "Interphase" = "Interphase",
                             "Metaphase" = "Metaphase",
                             "LBRFKBPGFP6DG5MAPStarFRB" = "+ MAPPER",
                             "LBRFKBPGFPStarFRB" = "- MAPPER"))

p1 <- ggplot(data = df,
             aes(x = factor(cell_state, levels=c("Pre", "Post")),
                 y = obj,
                 group = id,
                 color = factor(id))) +
  geom_line(linewidth = 1,
            alpha = 0.5) +
  geom_point(size = 2) +
  facet_grid(cell_channel ~ cond + cell_type, labeller = facet_names) +
  labs(x = "", y = "Total clusters") +
  theme_bw() +
  theme(legend.position = "none")
# ggsave("Output/Plots/comparison.png", p1, width = 200, height = 120, dpi = 300, units = "mm")

sp <- df %>% 
  group_by(cell_channel, cond, cell_type, cell_state, rep) %>% 
  summarise(mean = mean(obj),
            sd = sd(obj),
            n = n())

p2 <- ggplot(data = df,
             aes(x = factor(cell_state, levels=c("Pre", "Post")),
                 y = obj,
                 group = id,
                 color = factor(rep))) +
  geom_line(linewidth = 0.2,
            alpha = 0.5) +
  geom_line(data = sp,
            aes(x = factor(cell_state, levels=c("Pre", "Post")),
                y = mean,
                group = rep,
                color = factor(rep)),
            linewidth = 0.5, alpha = 0.5) +
  geom_point(data = sp,
             aes(x = factor(cell_state, levels=c("Pre", "Post")),
                 y = mean,
                 group = rep,
                 color = factor(rep)),
             size = 1, alpha = 0.5, shape = 16) +
  scale_colour_manual(values = c("#4477aa","#228833","#ccbb44","#ee6677")) +
  lims(y = c(0,900)) +
  facet_grid(cell_channel ~ factor(cond, levels = c("LBRFKBPGFPStarFRB","LBRFKBPGFP6DG5MAPStarFRB")) + cell_type, labeller = facet_names) +
  labs(x = "", y = "Total clusters") +
  theme_bw(9) +
  theme(legend.position = "none")

# ggsave("Output/Plots/comparison2.png", p2, width = 200, height = 120, dpi = 300, units = "mm")
ggsave("Output/Plots/comparison2.pdf", p2, width = 90, height = 70, units = "mm")


# statistical test ----

# our questions are:
# Do clusters change pre vs post for ch1 and for ch2? -- see below
# Is there a difference between interphase and mitosis (looks like there are fewer contacts in mitosis)?
# Does MAPPER expression affect LBR labelling of contacts?
# Does LBR labelling affect MAPPER labelling?

# We can use ANOVA with Tukey HSD post hoc tests
# We need to make a new column for to group cond, cell_type and cell_state
sp$group <- factor(paste(sp$cond, sp$cell_type, sp$cell_state, sep = "_"))
df$group <- factor(paste(df$cond, df$cell_type, df$cell_state, sep = "_"))

## compare using replicate means
ch1 <- aov(mean ~ group, data = sp[sp$cell_channel == "ch1", ])
summary(ch1)
tukey <- TukeyHSD(ch1)
tukey
# print p-values of interest
print_pval(tukey, "LBRFKBPGFP6DG5MAPStarFRB_Interphase_Pre","LBRFKBPGFP6DG5MAPStarFRB_Metaphase_Pre")
print_pval(tukey, "LBRFKBPGFP6DG5MAPStarFRB_Interphase_Post","LBRFKBPGFP6DG5MAPStarFRB_Metaphase_Post")

ch2 <- aov(mean ~ group, data = sp[sp$cell_channel == "ch2", ])
summary(ch2)
tukey <- TukeyHSD(ch2)
tukey
# print p-values of interest
print_pval(tukey, "LBRFKBPGFP6DG5MAPStarFRB_Interphase_Post","LBRFKBPGFP6DG5MAPStarFRB_Metaphase_Post")
print_pval(tukey, "LBRFKBPGFPStarFRB_Interphase_Post","LBRFKBPGFPStarFRB_Metaphase_Post")
print_pval(tukey, "LBRFKBPGFP6DG5MAPStarFRB_Interphase_Post","LBRFKBPGFPStarFRB_Interphase_Post")
print_pval(tukey, "LBRFKBPGFP6DG5MAPStarFRB_Metaphase_Post","LBRFKBPGFPStarFRB_Metaphase_Post")

# n.b. similar result using individual data

# Do clusters change pre vs post for ch1 and for ch2?
# we need to do paired t-tests for this question. We can use Holm-Bonferroni correction for multiple testing
# we need to make a new column for to group cond, cell_type and cell_state
df$group <- paste(df$cond, df$cell_type, df$cell_state, sep = "_")
sp$group <- paste(sp$cond, sp$cell_type, sp$cell_state, sep = "_")
# make a list of odd numbered items from groupsToCompare
x <- c("LBRFKBPGFP6DG5MAPStarFRB_Interphase_Pre","LBRFKBPGFP6DG5MAPStarFRB_Metaphase_Pre","LBRFKBPGFPStarFRB_Interphase_Pre","LBRFKBPGFPStarFRB_Metaphase_Pre")
# take x, replace Pre with post, and assign to y
y <- gsub("Pre", "Post", x)
# do the paired comparisons on all data
tt <- multi_t_test(data = df[df$cell_channel == "ch1", ], x = x, y = y, val = "obj")
tt
tt <- multi_t_test(data = df[df$cell_channel == "ch2", ], x = x, y = y, val = "obj")
tt

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}

F2B <- df %>% 
  select(cell_channel, cond, cell_type, cell_state, rep, obj)
write.csv(F2B, "Output/Data/F2B.csv", row.names = FALSE)
