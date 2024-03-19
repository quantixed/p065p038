# devtools::install_github("quantixed/TrackMateR")
library(TrackMateR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(cowplot)
library(patchwork)

## functions ----

process_xml_file <- function(xmlpath, timepath, cell = NULL) {
  
  cat(paste("File is\n", cell, " = ", xmlpath,"\n"))
  tmObj <- readTrackMateXML(XMLpath = xmlpath)
  
  if(is.null(cell)) {
    cat("Could not determine series number\n")
    return()
  }
  # this example is in pixels and frames
  # tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.065, xyunit = "um")
  # we can get a data frame of the correct TrackMate data
  tmDF <- tmObj[[1]]
  # and a data frame of calibration information
  calibrationDF <- tmObj[[2]]
  
  # let's dump any traces that are four frames or less
  # and any that only exist before 100 s
  legal <- tmDF %>% 
    group_by(trace) %>% 
    reframe(trace = trace,
            frames = max(frame) - min(frame) + 1,
            maxframe = max(frame))
  legal <- unique(legal)
  # at this point frame 70 is 100 s
  legal <- legal[legal$frames > 4 & legal$maxframe > 70, ]
  keeps <- tmDF$trace %in% legal$trace
  tmDF <- tmDF[keeps,]
  if(nrow(tmDF) == 0) {
    cat("No valid traces\n")
    return()
  }
  
  # load timestamps
  tDF <- read.csv(file = timepath)
  # we currently have a problem with funky times in the nd2 file
  # check that deltaT increments in an expected way, and correct if not
  if(tDF$deltaT[5] - tDF$deltaT[2] < 2) {
    cat("deltaT is not incrementing as expected. Correcting...\n")
    tDF$deltaT <- c(0,seq(from = 30, by = 1, length.out = nrow(tDF) - 1))
  }
  # t is frames in TrackMate frames also start at 0
  # tDF$t <- tDF$t + 1
  # we only need frames and time
  tDF <- subset(tDF, select = c(t, deltaT))
  
  # to match up the frame numbers with their corresponding time value we can use merge
  tmDF <- merge(x = tmDF, y = tDF, by.x = "frame", by.y = "t", all.x = TRUE, sort = FALSE)
  # set the column called t as the time and delete the extra column, reorder so that we are back with the original order
  tmDF$t <- tmDF$deltaT
  tmDF <- subset(tmDF, select = -deltaT)
  tmDF <- tmDF[order(tmDF$trace, tmDF$t), ]
  
  # the cumulative duration column needs to be fixed
  dur <- numeric()
  dur[1] <- 0
  startdur <- tmDF$t[1]
  
  cat("Recalculating durations...\n")
  
  cat("Number of traces: ", length(unique(tmDF$trace)), "\n", "Number of rows: ", nrow(tmDF), "\n")
  for (i in 2:nrow(tmDF)){
    if(tmDF$trace[i] != tmDF$trace[i-1]) {
      startdur <- tmDF$t[i]
    }
    dur[i] <- tmDF$t[i] - startdur
  }
  tmDF$track_duration <- dur
  
  # in order to average the traces we need to array each into a column according to t
  intensities <- tmDF %>% 
    select(c(trace,t,mean_intensity)) %>% 
    pivot_wider(
      names_from = trace,
      values_from = mean_intensity,
      names_glue = "int_{trace}"
    )
  # reorder by t as the order may be messed up during pivot_wider
  intensities <- intensities[order(intensities$t), ]
  # average all columns except first one (t). We use an.approx from zoo to fill in gaps
  avg <- rowMeans(na.approx(intensities[, 2:ncol(intensities)]), na.rm = TRUE)
  # count valid tracks per time point
  nv <- rowSums(!is.na(intensities[, 2:ncol(intensities)]))
  # assemble simple df for plotting
  avgDF <- data.frame(t = intensities$t,
                      avg = avg,
                      n = nv)
  # at this point we are missing the earlier time points and any laer ones with no spots
  # tDF$deltaT has all time points, we need to add these to avgDF and set the missing ones to NA
  tDF_fill <- data.frame(t = tDF$deltaT,
                         avg = NA,
                         n = 0)
  # drop any rows that are already in avgDF
  tDF_fill <- tDF_fill[!tDF_fill$t %in% avgDF$t, ]
  # now we can rbind the two together
  avgDF <- rbind(avgDF, tDF_fill)
  # and sort by t
  avgDF <- avgDF[order(avgDF$t), ]
  
  p1 <- ggplot(tmDF, aes(x = t, y = mean_intensity)) +
    geom_line(aes(group = trace), colour = "#00a651", linewidth = 0.25, alpha = 0.5) +
    geom_line(avgDF, mapping=aes(x = t, y = avg), colour = "black", alpha = 0.5, linewidth = 0.5) +
    lims(y = c(0,NA), x = c(0,200)) +
    labs(x = "Time (s)", y = "Cluster intensity (AU)") +
    theme_cowplot(8) +
    theme(legend.position = "none")
  
  # in order to average the traces we need to array each into a column according to t
  areas <- tmDF %>% 
    group_by(t) %>% 
    summarise(mean_area = mean(area, na.rm = TRUE),
              sd_area = sd(area,na.rm = TRUE))
  
  p2 <- ggplot(tmDF, aes(x = t, y = area)) +
    geom_line(aes(group = trace), colour = "#00a651", linewidth = 0.25, alpha = 0.5) +
    geom_line(areas, mapping=aes(x = t, y = mean_area), colour = "black", alpha = 0.5, linewidth = 0.5) +
    lims(y = c(0,NA), x = c(0,200)) +
    labs(x = "Time (s)", y = "Area (µm²)") +
    theme_cowplot(8) +
    theme(legend.position = "none")
  
  p3 <- ggplot(avgDF, aes(x = t, y = n)) +
    geom_line(colour = "black", alpha = 0.5, linewidth = 0.5) +
    lims(y = c(0,NA), x = c(0,200)) +
    labs(x = "Time (s)", y = "Total clusters") +
    theme_cowplot(8) +
    theme(legend.position = "none")
  
  r1 <- p1 | p2 | p3
  
  ggsave(paste0("Output/Plots/all_",cell,".pdf"),
         r1, width = 170, height = 40, units = "mm", bg = "white")
  
  # write out summary data
  # check if mean(avgDF$n) is > 5 and if so, write out summary data
  if(mean(avgDF$n) > 5) {
    write.table(avgDF, paste0("Output/Data/spotIntensity_",cell,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(areas, paste0("Output/Data/spotArea_",cell,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    cat("Not enough spots to write out summary data.\n")
  }
}

#' average_waves - a function that is analogous to IgorPro's average waves.
#' 
#' Resamples each "wave" at constant time steps using linear interpolation, to allow averaging.
#'
#' @param df data frame in long format
#' @param tcol string to identify time column
#' @param valcol string to identify value column
#' @param groupcol string to identify grouping column
#' @param tstep time step to resample at
#' @return data frame of summary data per resampled timepoint

average_waves <- function(df, tcol = "t", valcol = "avg", groupcol = "trace", tstep = 1) {
  # in order to average the traces we need to array each into a column according to time
  intensities <- df %>%
    # select only the columns we need
    select(all_of(c(!!tcol,!!valcol,!!groupcol))) %>% 
    # pivot the data so that each trace is a column
    pivot_wider(
      names_from = !!groupcol,
      values_from = !!valcol
    )
  # reorder by t as the order may be messed up during pivot_wider
  intensities <- intensities[order(intensities[[tcol]]), ]
  # make a desired length dataframe with constant time steps
  tmin <- floor(min(intensities[[tcol]]))
  tmax <- ceiling(max(intensities[[tcol]]))
  newdf <- data.frame(t = seq(from = tmin, to = tmax, by = tstep))
  for(i in 2:ncol(intensities)) {
    subdf <- cbind(intensities[,1],intensities[,i])
    subdf <- subdf[complete.cases(subdf),]
    temp <- approx(x = subdf[,1], y = subdf[,2], xout = seq(from = tmin, to = tmax, by = tstep))
    newdf <- cbind(newdf,temp[[2]])
  }
  names(newdf) <- names(intensities)
  # average all columns except first one (t).
  avg <- rowMeans(newdf[,-1], na.rm = TRUE)
  stdev <- apply(newdf[,-1], 1, sd, na.rm = TRUE)
  # count valid tracks per time point
  nv <- rowSums(!is.na(newdf[, -1]))
  # assemble simple df for plotting
  avgDF <- data.frame(t = newdf[[tcol]],
                      avg = avg,
                      sd = stdev,
                      n = nv,
                      sem = stdev / sqrt(nv))
  return(avgDF)
}

## main script ----

# compareDatasets()
# data <- read.csv("Output/Data/allComparison.csv")
# print(mean(data$displacement))

# in Data/ we have 1 xml file per cell and a corresponding deltaT file for timings
xmlFiles <- list.files(path = "Data/all", pattern = "*.xml", full.names = TRUE)
tFiles <- list.files("Data/all", pattern = "*.txt", full.names = TRUE)

# load TrackMate XML file
if(length(xmlFiles) == 0) {
  cat("Caution: Using XML txt file detected as TrackMate XML file.\nIf xmlPath is NULL, TrackMate XML file is missing.\n")
}
i <- 0
for(xmlPath in xmlFiles) {
  stub <- basename(xmlPath)
  stub <- paste0(strsplit(stub,"-")[[1]][1], "_deltaT.txt")
  timePath <- sub(basename(xmlPath), stub, xmlPath)
  process_xml_file(xmlpath = xmlPath, timepath = timePath, cell = i)
  i <- i + 1
}

# now we will read in all the files we just wrote out and make a summary plot
# first we need to get a list of all the files we just wrote out
spotIntensityFiles <- list.files(path = "Output/Data", pattern = "spotIntensity", full.names = TRUE)
spotAreaFiles <- list.files(path = "Output/Data", pattern = "spotArea", full.names = TRUE)

# now we will read in the files, make an extra column to identify them, and then bind them together
spotIntensityDF <- data.frame()
for(spotIntensityFile in spotIntensityFiles) {
  tmp <- read.table(spotIntensityFile, header = TRUE)
  tmp$trace <- basename(spotIntensityFile)
  spotIntensityDF <- rbind(spotIntensityDF, tmp)
}
spotAreaDF <- data.frame()
for(spotAreaFile in spotAreaFiles) {
  tmp <- read.table(spotAreaFile, header = TRUE)
  tmp$trace <- basename(spotAreaFile)
  spotAreaDF <- rbind(spotAreaDF, tmp)
}

# remove any rows with NA values in t column
spotIntensityDF <- spotIntensityDF[!is.na(spotIntensityDF$t),]
spotAreaDF <- spotAreaDF[!is.na(spotAreaDF$t),]

# now we will make a summary table of the summary tables
spotIntensityDF_summary <- average_waves(df = spotIntensityDF, tcol = "t", valcol = "avg", groupcol = "trace", tstep = 1)
spotAreaDF_summary <- average_waves(df = spotAreaDF, tcol = "t", valcol = "mean_area", groupcol = "trace", tstep = 1)
spotCountDF_summary <- average_waves(df = spotIntensityDF, tcol = "t", valcol = "n", groupcol = "trace", tstep = 1)
# now we can make a summary plot for each data frame showing the mean ± sd in ggplot
# we will use the summary tables for this
# first we will make a plot of spot intensity
p1 <- ggplot(data = spotIntensityDF_summary) +
  geom_line(aes(x = t, y = avg)) +
  geom_ribbon(aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Intensity (a.u.)")

# next we will make a plot of spot area
p2 <- ggplot(data = spotAreaDF_summary, aes(x = t, y = avg)) +
  geom_line(aes(x = t, y = avg)) +
  geom_ribbon(aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Area (µm²)")

# now we will make a plot of spot count
p3 <- ggplot(data = spotCountDF_summary, aes(x = t, y = avg)) +
  geom_line(aes(x = t, y = avg)) +
  geom_ribbon(aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Cluster count")

# now we will make a plot of all the summary plots
r1 <- p1 | p2 | p3

ggsave(paste0("Output/Plots/all_summary.pdf"),
       r1, width = 170, height = 40, units = "mm", bg = "white")

# colour palette using the number of traces we have
oneColour <- rep("#00a651", length(unique(spotIntensityDF$trace)))
# first we will make a plot of spot intensity
p5 <- ggplot() +
  geom_path(data = spotIntensityDF, aes(x = t, y = avg, colour = trace), linewidth = 0.25, alpha = 0.4) +
  scale_color_manual(values = oneColour) +
  geom_line(data = spotIntensityDF_summary, aes(x = t, y = avg), linewidth = 1) +
  geom_ribbon(data = spotIntensityDF_summary, aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  lims(y = c(0,NA), x = c(0,200)) +
  theme(legend.position = "none") +
  labs(x = "Time (s)", y = "Intensity (a.u.)")

# next we will make a plot of spot area
p6 <- ggplot() +
  geom_path(data = spotAreaDF, aes(x = t, y = mean_area, colour = trace), linewidth = 0.25, alpha = 0.4) +
  scale_color_manual(values = oneColour) +
  geom_line(data = spotAreaDF_summary, aes(x = t, y = avg), linewidth = 1) +
  geom_ribbon(data = spotAreaDF_summary, aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  theme(legend.position = "none") +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Area (µm²)")

# now we will make a plot of spot count
p7 <- ggplot() + 
  geom_path(data = spotIntensityDF, aes(x = t, y = n, colour = trace), linewidth = 0.25, alpha = 0.4) +
  scale_color_manual(values = oneColour) +
  geom_line(data = spotCountDF_summary, aes(x = t, y = avg), linewidth = 1) +
  geom_ribbon(data = spotCountDF_summary, aes(x = t, ymin = avg - sd, ymax = avg + sd), alpha = 0.2) +
  theme_cowplot(8) +
  theme(legend.position = "none") +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Spot count")

# now we will make a plot of all the summary plots
r2 <- p5 | p6 | p7

ggsave(paste0("Output/Plots/all_traces_summary.pdf"),
       r2, width = 170, height = 40, units = "mm", bg = "white")

## For the figure we will combine specific items from above to save space ----
xmlPath <- "Data/all/230802_LD312_HCT116LBRFKBPGFP_StargazinmChFRB_dish4001-600box.xml"
tmObj <- readTrackMateXML(XMLpath = xmlPath)
# we can get a data frame of the correct TrackMate data
tmDF <- tmObj[[1]]
# and a data frame of calibration information
calibrationDF <- tmObj[[2]]

# let's dump any traces that are four frames or less
# and any that only exist before 100 s
legal <- tmDF %>% 
  group_by(trace) %>% 
  reframe(trace = trace,
          frames = max(frame) - min(frame) + 1,
          maxframe = max(frame))
legal <- unique(legal)
# at this point frame 70 is 100 s
legal <- legal[legal$frames > 4 & legal$maxframe > 70, ]
keeps <- tmDF$trace %in% legal$trace
tmDF <- tmDF[keeps,]
if(nrow(tmDF) == 0) {
  cat("No valid traces\n")
  return()
}

# load timestamps
tDF <- read.csv(file = "Data/all/230802_LD312_HCT116LBRFKBPGFP_StargazinmChFRB_dish4001_deltaT.txt")
# we currently have a problem with funky times in the nd2 file
# check that deltaT increments in an expected way, and correct if not
if(tDF$deltaT[5] - tDF$deltaT[2] < 2) {
  cat("deltaT is not incrementing as expected. Correcting...\n")
  tDF$deltaT <- c(0,seq(from = 30, by = 1, length.out = nrow(tDF) - 1))
}
# t is frames in TrackMate frames also start at 0
# tDF$t <- tDF$t + 1
# we only need frames and time
tDF <- subset(tDF, select = c(t, deltaT))

# to match up the frame numbers with their corresponding time value we can use merge
tmDF <- merge(x = tmDF, y = tDF, by.x = "frame", by.y = "t", all.x = TRUE, sort = FALSE)
# set the column called t as the time and delete the extra column, reorder so that we are back with the original order
tmDF$t <- tmDF$deltaT
tmDF <- subset(tmDF, select = -deltaT)
tmDF <- tmDF[order(tmDF$trace, tmDF$t), ]

# the cumulative duration column needs to be fixed
dur <- numeric()
dur[1] <- 0
startdur <- tmDF$t[1]

cat("Recalculating durations...\n")

cat("Number of traces: ", length(unique(tmDF$trace)), "\n", "Number of rows: ", nrow(tmDF), "\n")
for (i in 2:nrow(tmDF)){
  if(tmDF$trace[i] != tmDF$trace[i-1]) {
    startdur <- tmDF$t[i]
  }
  dur[i] <- tmDF$t[i] - startdur
}
tmDF$track_duration <- dur

# in order to average the traces we need to array each into a column according to t
intensities <- tmDF %>% 
  select(c(trace,t,mean_intensity)) %>% 
  pivot_wider(
    names_from = trace,
    values_from = mean_intensity,
    names_glue = "int_{trace}"
  )
# reorder by t as the order may be messed up during pivot_wider
intensities <- intensities[order(intensities$t), ]
# average all columns except first one (t). We use an.approx from zoo to fill in gaps
avg <- rowMeans(na.approx(intensities[, 2:ncol(intensities)]), na.rm = TRUE)
# count valid tracks per time point
nv <- rowSums(!is.na(intensities[, 2:ncol(intensities)]))
# assemble simple df for plotting
avgDF <- data.frame(t = intensities$t,
                    avg = avg,
                    n = nv)
# at this point we are missing the earlier time points and any laer ones with no spots
# tDF$deltaT has all time points, we need to add these to avgDF and set the missing ones to NA
tDF_fill <- data.frame(t = tDF$deltaT,
                       avg = NA,
                       n = 0)
# drop any rows that are already in avgDF
tDF_fill <- tDF_fill[!tDF_fill$t %in% avgDF$t, ]
# now we can rbind the two together
avgDF <- rbind(avgDF, tDF_fill)
# and sort by t
avgDF <- avgDF[order(avgDF$t), ]

# in order to average the traces we need to array each into a column according to t
areas <- tmDF %>% 
  group_by(t) %>% 
  summarise(mean_area = mean(area, na.rm = TRUE),
            sd_area = sd(area,na.rm = TRUE))

p1 <- ggplot(tmDF, aes(x = t, y = mean_intensity)) +
  geom_line(aes(group = trace), colour = "#00a651", linewidth = 0.25, alpha = 0.5) +
  geom_line(data = avgDF, mapping=aes(x = t, y = avg), colour = "#006733", alpha = 1, linewidth = 0.5) +
  geom_line(data = spotIntensityDF_summary, mapping=aes(x = t, y = avg), colour = "black", linetype = "dotted", linewidth = 0.5) +
  lims(y = c(0,220), x = c(0,200)) +
  labs(x = "Time (s)", y = "Cluster intensity (AU)") +
  theme_cowplot(8) +
  theme(legend.position = "none")

p2 <- ggplot(tmDF, aes(x = t, y = area)) +
  geom_line(aes(group = trace), colour = "#00a651", linewidth = 0.25, alpha = 0.5) +
  geom_line(areas, mapping=aes(x = t, y = mean_area), colour = "#006733", alpha = 1, linewidth = 0.5) +
  geom_line(data = spotAreaDF_summary, mapping=aes(x = t, y = avg), colour = "black", linetype = "dotted", linewidth = 0.5) +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Area (µm²)") +
  theme_cowplot(8) +
  theme(legend.position = "none")

p3 <- ggplot(avgDF, aes(x = t, y = n)) +
  geom_line(colour = "#006733", alpha = 1, linewidth = 0.5) +
  geom_line(data = spotCountDF_summary, mapping=aes(x = t, y = avg), colour = "black", linetype = "dotted", linewidth = 0.5) +
  lims(y = c(0,NA), x = c(0,200)) +
  labs(x = "Time (s)", y = "Total clusters") +
  theme_cowplot(8) +
  theme(legend.position = "none")

r1 <- p3 | p2 | p1

ggsave("Output/Plots/figure_plot.pdf",
       r1, width = 170, height = 40, units = "mm", bg = "white")

