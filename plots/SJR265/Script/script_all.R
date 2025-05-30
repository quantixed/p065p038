## this script can be run from the GitHub repo because results of n1 and n2
## scripts are already in Output/Data

library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(cowplot)

## Variables ----

threshold <- 0.8

## Functions ----

read_in_data <- function(ptrn){
  files <- list.files(path = "Output/Data/", pattern = ptrn, full.names = TRUE)
  # read in the data
  n1 <- read.csv(grep(files, pattern = "n1", value = TRUE))
  n2 <- read.csv(grep(files, pattern = "n2", value = TRUE))
  
  # add a column to designate n1 and n2
  n1$version <- "n1"
  n2$version <- "n2"
  
  # combine the data frames
  combined <- rbind(n1, n2)
  
  return(combined)
}

segment_line <- function(x, filt = c(0, 1, 0)) {
  # is line binarised?
  # check that vector only consists of 0s and 1s
  if (length(unique(x)) > 2) {
    stop("Filter vector must only contain 0s and 1s")
  }
  
  # deal with ends first
  # if first value is 1, insert 0 before it
  if (x[1] == 1) {
    x <- c(0, x)
  }
  # i last value is 1 insert 0 after it
  if (x[length(x)] == 1) {
    x <- c(x, 0)
  }
  # go through the vector and if the sequence filt is found change it to 0
  for (i in 1:(length(x) - length(filt) + 1)) {
    if (all(x[i:(i + length(filt) - 1)] == filt)) {
      x[i:(i + length(filt) - 1)] <- 0
    }
  }
  # we need to find runs of 1s and determine how many runs there are and what is their size
  # we will use rle function
  rle_x <- rle(x)
  # determine how many runs of 1s there are
  n_runs <- sum(rle_x$values == 1)
  # determine the size of each run
  run_size <- rle_x$lengths[rle_x$values == 1]
  
  # return run_size
  return(run_size)
}

## Script ----

# load the dataframes we used for plotting
allResults <- read_in_data("allResults*")

# rename cond so that ctrl = Control and TG = Thapsigargin
allResults$cond <- factor(allResults$cond, levels = c("ctrl", "TG"), labels = c("Control", "Thapsigargin"))

# find the fraction of Norm values for each file that are greater than 1
fractionAboveThresh <- aggregate(Norm ~ file + cond, data = allResults, FUN = function(x) sum(x > threshold) / length(x))


ggplot(data = fractionAboveThresh, aes(x = cond, y = Norm)) + 
  geom_sina(colour = "#00a651", size = 0.5) + #individual points spaced with density
  geom_errorbar(data = fractionAboveThresh %>%
                  group_by(cond) %>%
                  summarise(mean = mean(Norm), sd = sd(Norm)), 
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                colour = "black",linewidth = 0.2, width = 0) +   #errorbar
  geom_crossbar(data = fractionAboveThresh %>%
                  group_by(cond) %>%
                  summarise(mean = mean(Norm), sd = sd(Norm)), 
                aes(y = mean, ymin = mean, ymax = mean),
                colour = "black",linewidth = 0.4, width = 1) +   #errorbar
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "", y = "MCS Fraction") +
  theme_cowplot(8) +
  theme(legend.position = "none")
ggsave("Output/Plots/fractionAboveThresh_all.pdf", width = 2.5, height = 2.5, dpi = 300, bg = "white")

# t test
t.test(Norm ~ cond, data = fractionAboveThresh)

## Export results ----
# if Output/Data directory does not exist, create it
if (!dir.exists("Output/Data")) {
  dir.create("Output/Data", recursive = TRUE)
}

S5B <- fractionAboveThresh %>%
  select(cond, Norm)
write.csv(S5B, "Output/Data/S5B.csv", row.names = FALSE)
