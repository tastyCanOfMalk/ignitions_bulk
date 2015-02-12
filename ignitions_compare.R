if (require("grid") == FALSE) install.packages("grid")
if (require("ggplot2") == FALSE) install.packages("ggplot2")
if (require("gridExtra") == FALSE) install.packages("gridExtra")
if (require("reshape") == FALSE) install.packages("reshape")

# ======================
# TABLE OF CONTENTS
# ==========================
# OPTIONS
# IMPORT & PREPARE
# FUNCTIONS & LIST CREATION
# CREATE SUMMARY & CALCULATIONS TABLE
# REMOVE ROW FUNCTION
#   round(removeRow(1),1)  # removes pyro1
#   round(removeRow(2),1)  # removes pyro2
#   round(removeRow(3),1)  # removes pyro3
# GRAPHIC SETTINGS
# GRAPH CREATION
# GRAPHIC OUTPUT
# EXPORT FUNCTIONS
#   exportPng() 
#   exportPdf() 
# MULTI SAMPLE GRAPHS
#   plotOverlays()  # max of 10 graphs
#   plotTables()  # max of 10 tables
#   plotOverlays()  # compare all data between all samples
#   plotTables()  # compare results from all samples
#   plotCompiledResults()  # means/stdevs of all samples
#   plotBoxplots()  # easy comparison of variables between all samples

#   exportPngSummary()  # exports all multisample graphs

# ======================
# OPTIONS
# ==========================

setwd("C://Users//yue.GLOBAL//Documents//R//ignition")
prefix.range <- c('C', 'D', 'E', 'F', 'G', 'H')
pre.prefix <- "2014-10-28,157"
prefix.list  <- paste0(pre.prefix, 
                       prefix.range[1:length(prefix.range)])

skip.lines <- 13  # default 11, 13 otherwise
#   Needs tweaking if:
#     Error in read.table(file = file, header = header, 
#      sep = sep, quote = quote,  : more columns than column names

read.rows <- 800  # default 840, less if file missing lines
exo.threshold <- 5  # default 10, low exo need around 5

remove1   <- FALSE
remove2   <- FALSE
remove3   <- FALSE  # set to TRUE if removal desired
autosave  <- TRUE  # automatically save sample graphs
autosave1 <- TRUE  # automatically save summary graphics

overlay.list       <- list()  # list for graphs, must be outside loop
results.table.list <- list()  # list for tables, must be outside loop
compiled.results   <- list()  # list for summary of sample data
compiled.results.all   <- list()  # list for summary of sample data

for(i in 1:length(prefix.range)){
  prefix <- prefix.list[[i]]
  
  # ======================
  # IMPORT & PREPARE
  # ==========================
  file.names <- paste0(prefix, '-', seq_len(3), '.log')
  df.list <- lapply(file.names, 
                    function(x) read.csv(x, 
                                         sep = '\t', 
                                         skip = skip.lines, 
                                         nrows = read.rows))
  
  df.total <- do.call(cbind, df.list)  # merge all data
  
  time      <- seq.int(0, 419.5, 0.5)  # create time variable axis
  time.rows <- read.rows - 1  # match length(time) with imported data length
  time      <- time[0:time.rows]  # modify length(time)
  
  keeps    <- c(4, 8, 12)  # keep PYRO (temperature) columns
  df.trim  <- cbind(time, df.total[-nrow(df.total), keeps]) 
  
  colnames(df.trim) <- c("time", "pyro1", "pyro2", "pyro3")
  
  pyro1 <- df.trim$pyro1 * 1000
  pyro2 <- df.trim$pyro2 * 1000
  pyro3 <- df.trim$pyro3 * 1000  # remove decimals from temperature data
  
  df.trim <- as.data.frame(cbind(time, pyro1, pyro2, pyro3)) # relevant data
  
  # ======================
  # FUNCTIONS & LIST CREATION
  # ===========================
  ign.names <- paste0("pyro", seq_len(3))  # for easier function calls
  
  variables.df <- data.frame(pyro1 = as.numeric(),
                             pyro2 = as.numeric(),
                             pyro3 = as.numeric())  # empty variables dataframe
  
  dropTime <- function(x){
    # Computes drop time  coordinates & average baseline temperature
    #   Based on a negative decrease in slope value
    #
    # Args:
    #   x: numeric vector whose droptime & baseline is to be calculated
    #
    # Returns: 
    #   list of values: drop time(x1), drop temp(y1), baseline(y-value)
    for (i in 1:length(x)){
      if (abs(x[i] - x[i + 1] > 5)){
        drop.time <- time[i]
        drop.temp <- x[i]
        baseline <- mean(x[1:i])
        A <- list(drop.time, drop.temp, baseline)
        return(A)
        break
      }
    }
  }
  
  drop.times <- cbind((lapply(mget(ign.names)[1:3], dropTime)))
  drop.times.table <- setNames(do.call(rbind.data.frame,drop.times),
                               dimnames(drop.times)[[2]])
  rownames(drop.times.table) <- c("pyro1", "pyro2", "pyro3")
  colnames(drop.times.table) <- c("drop.time(x1)", 
                                  "drop.temp(y1)",
                                  "baseline.temp")
  
  exoTime <- function(x){
    # Computes exotherm coordinates
    #   Based on a positive increase in slope, "exo.threshold",
    #     may need adjustment based on speed of exotherm reaction
    #     (slower reaction needs a lower threshold value)
    #   
    # Args:
    #   x: numeric vector whose exotherm time is to be calculated
    #
    # Returns: 
    #   list of values: exotherm time(x2), exotherm temp(y2)
    for (i in dropTime(x)[[1]]:length(x)){
      if (abs(x[i+1] - x[i] > exo.threshold && x[i+10] > dropTime(x)[[3]])){
        # value must be within 5-seconds from baseline.temp
        exo.time <- time[i]
        exo.temp <- x[i]
        A <- list(exo.time, exo.temp)
        return(A)
        break
      }
    }
  }
  
  exo.times <- cbind((lapply(mget(ign.names)[1:3], exoTime)))
  exo.times.table <- setNames(do.call(rbind.data.frame,exo.times),
                              dimnames(exo.times)[[2]])
  rownames(exo.times.table) <- c("pyro1", "pyro2", "pyro3")
  colnames(exo.times.table) <- c("exo.time(x2)", 
                                 "exo.temp(y2)")
  
  maxTemp <- function(x){
    # Computes maximum temperature coordinates
    #
    # Args:
    #   x: numeric vector whose max temperature is to be calculated
    #
    # Returns: 
    #   list of values: max temp time(x3), max temp temp(y3)
    max.temp <- max(x)
    max.time <- match(max(x),x)/2
    delta.temp <- max.temp - dropTime(x)[[3]]
    A <- list(max.time, max.temp, delta.temp)
    return(A)    
  }
  
  max.temps <- cbind((lapply(mget(ign.names)[1:3], maxTemp)))
  max.temps.table <- setNames(do.call(rbind.data.frame,max.temps), 
                              dimnames(max.temps)[[2]])
  rownames(max.temps.table) <- c("pyro1", "pyro2", "pyro3")
  colnames(max.temps.table) <- c("max.temp.time(x3)", 
                                 "max.temp(y3)",
                                 "delta.temp")
  
  exoDuration <- function(x){
    # Computes exotherm duration coordinates
    #
    # Args:
    #   x: numeric vector whose exotherm duration is to be calculated
    #
    # Returns: 
    #   list of values: exo duration(x4), exo duration temp temp(y4)
    for (i in maxTemp(x)[[1]]:length(x)){
      if (x[i*2] <= dropTime(x)[3]){
        baseline1.time <- (i)
        baseline1.temp <- (x[i*2])
        A <- list(baseline1.time, baseline1.temp)
        return(A)
        break
      }
    }
  }
  
  exo.durations <- cbind((lapply(mget(ign.names)[1:3], exoDuration)))
  exo.durations.table <- setNames(do.call(rbind.data.frame, exo.durations),
                                  dimnames(exo.durations)[[2]])
  rownames(exo.durations.table) <- c("pyro1", "pyro2", "pyro3")
  colnames(exo.durations.table) <- c("exo.duration.time(x4)", 
                                     "exo.duration.temp(y4)")
  
  # ======================
  # CREATE SUMMARY & CALCULATIONS TABLE
  # ===========================
  summary.table <- cbind(drop.times.table, 
                         exo.times.table, 
                         max.temps.table, 
                         exo.durations.table)  # all results in one table
  
  calculations.table <- cbind(exo.times.table[[1]] - drop.times.table[[1]],
                              round(max.temps.table[[3]], 0),
                              exo.durations.table[[1]] - exo.times.table[[1]])
  
  rownames(calculations.table) <- c("pyro1", "pyro2", "pyro3")
  colnames(calculations.table) <- c("ignite time", 
                                     "delta temp",
                                     "exo duration")
  
  ign.time.mean     <- round(mean(calculations.table[1:3, 1]), 0)
  max.temp.mean     <- round(mean(calculations.table[1:3, 2]), 0)
  exo.duration.mean <- round(mean(calculations.table[1:3, 3]), 0)
  
  means           <- cbind(ign.time.mean, max.temp.mean, exo.duration.mean)
  rownames(means) <- c("mean")
  
  
  ign.time.sd       <- round(sd(calculations.table[1:3, 1]), 0)
  max.temp.sd       <- round(sd(calculations.table[1:3, 2]), 0)
  exo.duration.sd   <- round(sd(calculations.table[1:3, 3]), 0)
  
  sds           <-  cbind(ign.time.sd, max.temp.sd, exo.duration.sd)
  rownames(sds) <- c("std.dev")
  
  calculations.table.total <- rbind (calculations.table, means, sds)
  
  # ======================
  # REMOVE ROW FUNCTION
  # ===========================
  removeRow <- function(x){
    # Allows removal of suspected outlier sample row from calculations.table
    # 
    # Args:
    #   x: row number to be removed (1, 2, or 3)
    # 
    # Returns: 
    #   a new table with specified row values replaced with NA
    temprow <- matrix(c(rep.int(NA, length(data))), nrow = 1, ncol = 3)
    calc.table.new <- rbind(calculations.table[-x,], temprow)
    
    ign.time.mean     <- mean(calc.table.new[1:3, 1], na.rm = TRUE)
    max.temp.mean     <- mean(calc.table.new[1:3, 2], na.rm = TRUE)
    exo.duration.mean <- mean(calc.table.new[1:3, 3], na.rm = TRUE)
    means           <- cbind(ign.time.mean, max.temp.mean, exo.duration.mean)
    rownames(means) <- c("mean")
    
    ign.time.sd     <- sd(calc.table.new[1:3, 1], na.rm = TRUE)
    max.temp.sd     <- sd(calc.table.new[1:3, 2], na.rm = TRUE)
    exo.duration.sd <- sd(calc.table.new[1:3, 3], na.rm = TRUE)
    sds           <- cbind(ign.time.sd, max.temp.sd, exo.duration.sd)
    rownames(sds) <- c("std.dev")
    
    calc.table.new <- rbind(calculations.table[-x,], temprow, means, sds)
    return(calc.table.new)
  }
  
  # Function call removes sample data from table based on conditional
  #   does not allow for removal of more than 1 line
  
  if (remove1 == TRUE){
    calculations.table.total <- round(removeRow(1),1)  # removes pyro1
    prefix <- paste0(prefix, "-R1")  # for filename
  }
  
  if (remove2 == TRUE){
    calculations.table.total <- round(removeRow(2),1)  # removes pyro2
    prefix <- paste0(prefix, "-R2")  # for filename
  }
  
  if (remove3 == TRUE){
    calculations.table.total <- round(removeRow(3),1)  # removes pyro3
    prefix <- paste0(prefix, "-R3") # for filename
  }
  
  # ======================
  # GRAPHIC SETTINGS
  # ===========================
  color.list <- list("black",     #1  pyro values
                     "red",       #2  drop time
                     "red",       #3  baseline avg
                     "blue",      #4  ignition time
                     "orange",    #5  max temp
                     "green",     #6  duration
                     "gray92",    #7  graph panel bg
                     "gray80",    #8  graph panel outline
                     "white",     #9  major gridlines
                     "gray90",    #10 minor gridlines
                     "lightblue", #11 p5 results fill
                     "white"      #12 p6 fill 
                     )
  
  size.list <- list(1,     #1 pyro values
                    3,     #2 coordinate values
                    0.25,  #3 baseline avg
                    .7,    #4 major gridlines
                    .5,    #5 minor gridlines
                    20     #6 p6 font size
                    )
  
  alpha.list <- list(.25,   #1 pyro values
                     1,     #2 coordinate values
                     0.10,  #3 baseline avg
                     1,     #4 combined graph pyro values
                     .08    #5 combined graph baseline alpha
  )
  
  xy.scale.list <- list(0,    #1 x-axis min
                        420,  #2 x-axis max
                        30,   #3 x-axis scale/breaks
                        950,  #4 y-axis min
                        1350, #5 y-axis max
                        #max(summary.table[[7]])+30, #5 y-axis max automatic
                        100,  #6 y-axis scale/breaks
                        0.5,  #7 x-axis vjust
                        1.0   #8 y-axis vjust
                        )
  
  # ======================
  # GRAPH CREATION
  # ===========================
  p1 <- ggplot(data = df.trim, aes(x = time, y  = pyro1)) + 
    
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_text(vjust = xy.scale.list[[7]]),
          axis.title.y = element_text(vjust = xy.scale.list[[8]])) + 
    scale_x_continuous("time (sec)",
                       limits = c(xy.scale.list[[1]],
                                  xy.scale.list[[2]]), 
                       breaks = seq.int(xy.scale.list[[1]],
                                        xy.scale.list[[2]],
                                        xy.scale.list[[3]])) +
    scale_y_continuous("temperature (°C)",
                       limits = c(xy.scale.list[[4]],
                                  xy.scale.list[[5]]), 
                       breaks = seq.int(xy.scale.list[[4]],
                                        xy.scale.list[[5]],
                                        xy.scale.list[[6]])) + 
    ggtitle("pyro1") +
    geom_point(color = color.list[[1]],  # pyro values 
               size  = size.list[[1]],
               alpha = alpha.list[[1]]) +
    geom_point(x = summary.table[[1, 1]],  # drop time x1
               y = summary.table[[1, 2]],  # drop time y1
               color = color.list[[2]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(y = summary.table[[1, 3]],  # baseline average y-value
               color = color.list[[3]],
               size  = size.list[[3]],
               alpha = alpha.list[[3]]) +
    geom_point(x = summary.table[[1, 4]],  # ignition time x2
               y = summary.table[[1, 5]],  # ignition time y2
               color = color.list[[4]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[1, 6]],  # max temp x3
               y = summary.table[[1, 7]],  # max temp y3
               color = color.list[[5]], 
               size  = size.list[[2]], 
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[1, 9]],  # exo duration x4
               y = summary.table[[1, 10]],  # exo duration y4
               color = color.list[[6]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]])
  
  p2 <- ggplot(data = df.trim, aes(x = time, y  = pyro2)) + 
    
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_text(vjust = xy.scale.list[[7]]),
          axis.title.y = element_text(vjust = xy.scale.list[[8]])) + 
    scale_x_continuous("time (sec)",
                       limits = c(xy.scale.list[[1]],
                                  xy.scale.list[[2]]), 
                       breaks = seq.int(xy.scale.list[[1]],
                                        xy.scale.list[[2]],
                                        xy.scale.list[[3]])) +
    scale_y_continuous("temperature (°C)",
                       limits = c(xy.scale.list[[4]],
                                  xy.scale.list[[5]]), 
                       breaks = seq.int(xy.scale.list[[4]],
                                        xy.scale.list[[5]],
                                        xy.scale.list[[6]])) + 
    ggtitle("pyro2") +
    geom_point(color = color.list[[1]],  # pyro values 
               size  = size.list[[1]],
               alpha = alpha.list[[1]]) +
    geom_point(x = summary.table[[2, 1]],  # drop time x1
               y = summary.table[[2, 2]],  # drop time y1
               color = color.list[[2]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(y = summary.table[[2, 3]],  # baseline average y-value
               color = color.list[[3]],
               size  = size.list[[3]],
               alpha = alpha.list[[3]]) +
    geom_point(x = summary.table[[2, 4]],  # ignition time x2
               y = summary.table[[2, 5]],  # ignition time y2
               color = color.list[[4]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[2, 6]],  # max temp x3
               y = summary.table[[2, 7]],  # max temp y3
               color = color.list[[5]], 
               size  = size.list[[2]], 
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[2, 9]],  # exo duration x4
               y = summary.table[[2, 10]],  # exo duration y4
               color = color.list[[6]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]])
  
  p3 <- ggplot(data = df.trim, aes(x = time, y  = pyro3)) + 
  
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_text(vjust = xy.scale.list[[7]]),
          axis.title.y = element_text(vjust = xy.scale.list[[8]])) + 
    scale_x_continuous("time (sec)",
                       limits = c(xy.scale.list[[1]],
                                  xy.scale.list[[2]]), 
                       breaks = seq.int(xy.scale.list[[1]],
                                        xy.scale.list[[2]],
                                        xy.scale.list[[3]])) +
    scale_y_continuous("temperature (°C)",
                       limits = c(xy.scale.list[[4]],
                                  xy.scale.list[[5]]), 
                       breaks = seq.int(xy.scale.list[[4]],
                                        xy.scale.list[[5]],
                                        xy.scale.list[[6]])) + 
    ggtitle("pyro3") +
    geom_point(color = color.list[[1]],  # pyro values 
               size  = size.list[[1]],
               alpha = alpha.list[[1]]) +
    geom_point(x = summary.table[[3, 1]],  # drop time x1
               y = summary.table[[3, 2]],  # drop time y1
               color = color.list[[2]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(y = summary.table[[3, 3]],  # baseline average y-value
               color = color.list[[3]],
               size  = size.list[[3]],
               alpha = alpha.list[[3]]) +
    geom_point(x = summary.table[[3, 4]],  # ignition time x2
               y = summary.table[[3, 5]],  # ignition time y2
               color = color.list[[4]],
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[3, 6]],  # max temp x3
               y = summary.table[[3, 7]],  # max temp y3
               color = color.list[[5]], 
               size  = size.list[[2]], 
               alpha = alpha.list[[2]]) +
    geom_point(x = summary.table[[3, 9]],  # exo duration x4
               y = summary.table[[3, 10]],  # exo duration y4
               color = color.list[[6]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]])
  
  df.trim.melt <- melt(df.trim, id = c("time"))
  
  p4 <- ggplot(data = df.trim.melt, aes(x = time,
                                        y = value, 
                                        group = variable, 
                                        color = variable)) + 
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_text(vjust = xy.scale.list[[7]]),
          axis.title.y = element_text(vjust = xy.scale.list[[8]])) + 
    scale_x_continuous("time (sec)",
                       limits = c(xy.scale.list[[1]],
                                  xy.scale.list[[2]]), 
                       breaks = seq.int(xy.scale.list[[1]],
                                        xy.scale.list[[2]],
                                        xy.scale.list[[3]])) +
    scale_y_continuous("temperature (°C)",
                       limits = c(xy.scale.list[[4]],
                                  xy.scale.list[[5]]), 
                       breaks = seq.int(xy.scale.list[[4]],
                                        xy.scale.list[[5]],
                                        xy.scale.list[[6]])) + 
    ggtitle(paste0("pyro1-3 overlay: ", prefix)) +
    scale_colour_discrete(name = "Sample") +
    geom_line(alpha = alpha.list[[4]]) +
    geom_line(y     = summary.table[[1, 3]], 
              color = "red",
              alpha = alpha.list[[5]]) +
    geom_line(y     = summary.table[[2, 3]], 
              color = "green", 
              alpha = alpha.list[[5]]) +
    geom_line(y     = summary.table[[3, 3]], 
              color = "blue", 
              alpha = alpha.list[[5]]) +
    geom_point(x = mean(summary.table[[1]]),  # mean drop time x1
               y = mean(summary.table[[2]]),  # mean drop temp y1
               color = color.list[[3]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) + 
    geom_point(x = mean(summary.table[[4]]),  # mean exo time x1
               y = mean(summary.table[[5]]),  # mean exo temp y1
               color = color.list[[4]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) + 
    geom_point(x = mean(summary.table[[6]]),  # mean max time x1
               y = mean(summary.table[[7]]),  # mean max temp y1
               color = color.list[[5]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]]) + 
    geom_point(x = mean(summary.table[[9]]),  # mean duration time x1
               y = mean(summary.table[[10]]),  # mean duration temp y1
               color = color.list[[6]], 
               size  = size.list[[2]],
               alpha = alpha.list[[2]])
  
  p5 <- qplot(1:10, 1:10, geom = "blank", main = prefix) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(color = "white")) + 
    scale_x_discrete("",breaks = NULL) +
    scale_y_discrete("",breaks = NULL) +
    annotation_custom(grob = tableGrob(calculations.table.total,
                                     gpar.corefill = gpar(fill = color.list[[11]],
                                                          alpha=0.5, 
                                                          col = NA),
                                     h.even.alpha = 0.5))
  
  p6.sample.title <- c("Sample: ", 
                    paste0(prefix),
                    "Report created: ", 
                    date())  # appears in last graph
  
  p6 <- qplot(1:10, 1:10, geom = "blank") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(color = "white")) + 
    scale_x_discrete("",breaks = NULL) +
    scale_y_discrete("",breaks = NULL) +
    annotation_custom(grob = tableGrob(p6.sample.title,
                                 gpar.corefill = gpar(fill = color.list[[12]],
                                                      alpha=0.5, 
                                                      col = NA),
                                 h.even.alpha = 0.5,
                                 gpar.coretext = gpar(fontsize = size.list[[6]])))
  
  # ======================
  # GRAPHIC OUTPUT
  # ===========================
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # function taken from http://goo.gl/6PJFSB
    # comments removed
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
      print(plots[[1]])
      } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
  }
  
  multiplot(p1, p2, p3, p4, p5, p6, cols=2)
  
  # ======================
  # EXPORT FUNCTIONS
  # ===========================
  exportPdf <- function(){
    file.name.pdf <- paste0(prefix, ".pdf")
    pdf(file = file.name.pdf,
        width = 11,
        height = 8.5)
  multiplot(p1, p2, p3, p4, p5, p6, cols=2)
  dev.off()
  }
  
  exportPng <- function(){
    file.name.png <- paste0(prefix, ".png")
    dev.copy(png, file = file.name.png,
             width = 1200,
             height = 900)
    multiplot(p1, p2, p3, p4, p5, p6, cols=2)
    dev.off()
  }
  
  if(autosave == TRUE) {
    exportPng()
  }
  
  overlay.list[[i]]         <- p4  # stores individual p4 graphs to list
  results.table.list[[i]]   <- p5  # stores individual p5 charts to list
  compiled.results.all[[i]] <- calculations.table.total[1:3,]  # pyro1-3
  compiled.results[[i]]     <- calculations.table.total[4:5,]  # mean & std.dev
}

# ======================
# MULTI SAMPLE GRAPHS
# ===========================
multiOverlay <- function(x, y = 2){
  # Calls multiplot function for plotting graphs from different samples
  # 
  # Args:
  #   x: number of graphs to be plotted
  #   y: number of columns
  # 
  # Returns: multiplot of individual sample overlay graphs
  #   
  if(length(x) == 1){
    multiplot(overlay.list[[1]], 
              cols = y)
  }   
  if(length(x) == 2){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              cols = y)
  } 
  if(length(x) == 3){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              cols = y)
  } 
  if(length(x) == 4){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              cols = y)
  } 
  if(length(x) == 5){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              cols = y)
  } 
  if(length(x) == 6){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              overlay.list[[6]],
              cols = y)
  } 
  if(length(x) == 7){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              overlay.list[[6]],
              overlay.list[[7]],
              cols = y)
  } 
  if(length(x) == 8){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              overlay.list[[6]],
              overlay.list[[7]],
              overlay.list[[8]],
              cols = y)
  } 
  if(length(x) == 9){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              overlay.list[[6]],
              overlay.list[[7]],
              overlay.list[[8]],
              overlay.list[[9]],
              cols = y)
  } 
  if(length(x) == 10){
    multiplot(overlay.list[[1]], 
              overlay.list[[2]],
              overlay.list[[3]],
              overlay.list[[4]],
              overlay.list[[5]],
              overlay.list[[6]],
              overlay.list[[7]],
              overlay.list[[8]],
              overlay.list[[9]],
              overlay.list[[10]],
              cols = y)
  } 
}
multiTable <- function(x, y = 2){
  # Calls multiplot function for plotting results table from different samples
  # 
  # Args:
  #   x: number of graphs to be plotted
  #   y: number of columns
  # 
  # Returns: multiplot of individual results tables
  #   
  if(length(x) == 1){
    multiplot(results.table.list[[1]], 
              cols = y)
  }   
  if(length(x) == 2){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              cols = y)
  } 
  if(length(x) == 3){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              cols = y)
  } 
  if(length(x) == 4){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              cols = y)
  } 
  if(length(x) == 5){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              cols = y)
  } 
  if(length(x) == 6){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              results.table.list[[6]],
              cols = y)
  } 
  if(length(x) == 7){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              results.table.list[[6]],
              results.table.list[[7]],
              cols = y)
  } 
  if(length(x) == 8){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              results.table.list[[6]],
              results.table.list[[7]],
              results.table.list[[8]],
              cols = y)
  } 
  if(length(x) == 9){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              results.table.list[[6]],
              results.table.list[[7]],
              results.table.list[[8]],
              results.table.list[[9]],
              cols = y)
  } 
  if(length(x) == 10){
    multiplot(results.table.list[[1]], 
              results.table.list[[2]],
              results.table.list[[3]],
              results.table.list[[4]],
              results.table.list[[5]],
              results.table.list[[6]],
              results.table.list[[7]],
              results.table.list[[8]],
              results.table.list[[9]],
              results.table.list[[10]],
              cols = y)
  } 
  
}

plotOverlays <- function(){
  # Calls multiOverlay function
  #
  # Returns: multiplot of individual sample overlay graphs
  #   
  multiOverlay(overlay.list)
}
plotTables <- function(){
  # Calls multiTable function
  #
  # Returns: multiplot of individual sample overlay graphs
  #   
  multiTable(results.table.list)  
}

# Below code is pretty crappy and needs revisited in the future.
# only supports 10 total samples per run due to the if statements
if(length(prefix.list) == 1){
  #for boxplot
  bplot.total <- cbind(compiled.results.all[1])
  sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1])
}
if(length(prefix.list) == 2){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"))
    
    #  For boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2])
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2])
  }
if(length(prefix.list) == 3){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"))
    #  For boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3])
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3])
  }
if(length(prefix.list) == 4){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"))
    # for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4])
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4])
  }
if(length(prefix.list) == 5){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"))
    #  for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5])
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5])
}
if(length(prefix.list) == 6){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5],
                                    compiled.results[6])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total),
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"),
                                          paste0(prefix.list[[6]], "-mean"),
                                          paste0(prefix.list[[6]], "-stdev"))
    #  for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5],
                         compiled.results.all[6]
    )
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5],
                         prefix.list[6], prefix.list[6], prefix.list[6])
  }
if(length(prefix.list) == 7){
    compile
    d.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5],
                                    compiled.results[6],
                                    compiled.results[7])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"),
                                          paste0(prefix.list[[6]], "-mean"),
                                          paste0(prefix.list[[6]], "-stdev"),
                                          paste0(prefix.list[[7]], "-mean"),
                                          paste0(prefix.list[[7]], "-stdev"))
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5],
                         compiled.results.all[6],
                         compiled.results.all[7])
    
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5],
                         prefix.list[6], prefix.list[6], prefix.list[6],
                         prefix.list[7], prefix.list[7], prefix.list[7])
  }
if(length(prefix.list) == 8){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5],
                                    compiled.results[6],
                                    compiled.results[7],
                                    compiled.results[8])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"),
                                          paste0(prefix.list[[6]], "-mean"),
                                          paste0(prefix.list[[6]], "-stdev"),
                                          paste0(prefix.list[[7]], "-mean"),
                                          paste0(prefix.list[[7]], "-stdev"),
                                          paste0(prefix.list[[8]], "-mean"),
                                          paste0(prefix.list[[8]], "-stdev"))
    #  for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5],
                         compiled.results.all[6],
                         compiled.results.all[7],
                         compiled.results.all[8])
    
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5],
                         prefix.list[6], prefix.list[6], prefix.list[6],
                         prefix.list[7], prefix.list[7], prefix.list[7],
                         prefix.list[8], prefix.list[8], prefix.list[8])
  }
if(length(prefix.list) == 9){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5],
                                    compiled.results[6],
                                    compiled.results[7],
                                    compiled.results[8],
                                    compiled.results[9])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"),
                                          paste0(prefix.list[[6]], "-mean"),
                                  
                                          paste0(prefix.list[[6]], "-stdev"),
                                          paste0(prefix.list[[7]], "-mean"),
                                          paste0(prefix.list[[7]], "-stdev"),
                                          paste0(prefix.list[[8]], "-mean"),
                                          paste0(prefix.list[[8]], "-stdev"),
                                          paste0(prefix.list[[9]], "-mean"),
                                          paste0(prefix.list[[9]], "-stdev"))
    #  for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5],
                         compiled.results.all[6],
                         compiled.results.all[7],
                         compiled.results.all[8])
    
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5],
                         prefix.list[6], prefix.list[6], prefix.list[6],
                         prefix.list[7], prefix.list[7], prefix.list[7],
                         prefix.list[8], prefix.list[8], prefix.list[8])
  }
if(length(prefix.list) == 10){
    compiled.results.total <- cbind(compiled.results[1], 
                                    compiled.results[2],
                                    compiled.results[3],
                                    compiled.results[4],
                                    compiled.results[5],
                                    compiled.results[6],
                                    compiled.results[7],
                                    compiled.results[8],
                                    compiled.results[9],
                                    compiled.results[10])
    
    compiled.results.total <- setNames(do.call(rbind.data.frame, 
                                               compiled.results.total), 
                                       dimnames(compiled.results.total)[[2]])
    
    colnames(compiled.results.total) <- c("ignite time", 
                                          "delta temp", 
                                          "exo duration")
    
    rownames(compiled.results.total) <- c(paste0(prefix.list[[1]], "-mean"),
                                          paste0(prefix.list[[1]], "-stdev"),
                                          paste0(prefix.list[[2]], "-mean"),
                                          paste0(prefix.list[[2]], "-stdev"),
                                          paste0(prefix.list[[3]], "-mean"),
                                          paste0(prefix.list[[3]], "-stdev"),
                                          paste0(prefix.list[[4]], "-mean"),
                                          paste0(prefix.list[[4]], "-stdev"),
                                          paste0(prefix.list[[5]], "-mean"),
                                          paste0(prefix.list[[5]], "-stdev"),
                                          paste0(prefix.list[[6]], "-mean"),
                                          paste0(prefix.list[[6]], "-stdev"),
                                          paste0(prefix.list[[7]], "-mean"),
                                          paste0(prefix.list[[7]], "-stdev"),
                                          paste0(prefix.list[[8]], "-mean"),
                                          paste0(prefix.list[[8]], "-stdev"),
                                          paste0(prefix.list[[9]], "-mean"),
                                          paste0(prefix.list[[9]], "-stdev"),
                                          paste0(prefix.list[[10]], "-mean"),
                                          paste0(prefix.list[[10]], "-stdev"))
    #  for boxplot
    bplot.total <- cbind(compiled.results.all[1],
                         compiled.results.all[2],
                         compiled.results.all[3],
                         compiled.results.all[4],
                         compiled.results.all[5],
                         compiled.results.all[6],
                         compiled.results.all[7],
                         compiled.results.all[8],
                         compiled.results.all[9],
                         compiled.results.all[10]
    )
    sample.name <- rbind(prefix.list[1], prefix.list[1], prefix.list[1],
                         prefix.list[2], prefix.list[2], prefix.list[2],
                         prefix.list[3], prefix.list[3], prefix.list[3],
                         prefix.list[4], prefix.list[4], prefix.list[4],
                         prefix.list[5], prefix.list[5], prefix.list[5],
                         prefix.list[6], prefix.list[6], prefix.list[6],
                         prefix.list[7], prefix.list[7], prefix.list[7],
                         prefix.list[8], prefix.list[8], prefix.list[8],
                         prefix.list[9], prefix.list[9], prefix.list[9],
                         prefix.list[10], prefix.list[10], prefix.list[10])
  }

plotCompiledResults <- function(){
  # Prints compiled results in a table for biewing
  # 
  # Args:
  #   
  # Returns: tableGrob of means/std.devs of all sample results
  #   
  qplot(1:10, 1:10, geom = "blank", main = "All samples") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(color = "white")) + 
    scale_x_discrete("",breaks = NULL) +
    scale_y_discrete("",breaks = NULL) +
    annotation_custom(grob = tableGrob(compiled.results.total,
                                       gpar.corefill = gpar(fill = color.list[[11]],
                                                            alpha=0.5, 
                                                            col = NA),
                                       h.even.alpha = 0.5))
}

plotBoxplots <- function(){
  # Formats data & creates individual boxplots, then calls multiplot()
  # 
  # Args:
  #   
  # Returns: boxplots of important data for comparison
  #   
  bplot.data <- setNames(do.call(rbind.data.frame, bplot.total), 
                         dimnames(bplot.total)[[2]])
  
  colnames(bplot.data) <- c("ignite.time",
                            "delta.temp", 
                            "exo.duration")
  bplot <- cbind(bplot.data, sample.name)
  bplot.melt <- melt(bplot, id = c("sample.name",
                                   "exo.duration",
                                   "delta.temp",
                                   "ignite.time"))
  
  g1 <- ggplot(data = bplot.melt, aes(x = sample.name,
                                      y = ignite.time, 
                                      fill = sample.name)) + 
    geom_boxplot() + 
    ggtitle("Ignite time") + 
    ylab("seconds") +
    guides(fill = FALSE) + 
    stat_summary(fun.y = mean, geom = "point", shape = 5, size = 4) + 
    geom_jitter(position = position_jitter(width = .2)) +
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = xy.scale.list[[8]]),
          axis.text.x = element_text(angle = 25, hjust = 1),
          text = element_text(size = 15))
  
  g2 <- ggplot(data = bplot.melt, aes(x = sample.name,
                                      y = delta.temp, 
                                      fill = sample.name)) + 
    geom_boxplot() + 
    ggtitle("Delta temp") + 
    ylab("°C") +    
    guides(fill = FALSE) + 
    stat_summary(fun.y = mean, geom = "point", shape = 5, size = 4) + 
    geom_jitter(position = position_jitter(width = .2)) +
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = xy.scale.list[[8]]),
          axis.text.x = element_text(angle = 25, hjust = 1),
          text = element_text(size = 15))
  
  g3 <- ggplot(data = bplot.melt, aes(x = sample.name,
                                      y = exo.duration, 
                                      fill = sample.name)) + 
    geom_boxplot() + 
    ggtitle("Exotherm duration") + 
    ylab("seconds") +  
    guides(fill = FALSE) + 
    stat_summary(fun.y = mean, geom = "point", shape = 5, size = 4) + 
    geom_jitter(position = position_jitter(width = .2)) +
    theme(panel.background = element_rect(fill  = color.list[[7]],
                                          color = color.list[[8]]),
          panel.grid.major = element_line(color = color.list[[9]],
                                          size  = size.list[[4]]),
          panel.grid.minor = element_line(color = color.list[[10]], 
                                          size  = size.list[[5]]),
          panel.grid.major.y = element_line(color = color.list[[9]],
                                            size  = size.list[[4]]),
          panel.grid.minor.y = element_line(color = color.list[[10]],
                                            size  = size.list[[5]]),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = xy.scale.list[[8]]),
          axis.text.x = element_text(angle = 25, hjust = 1),
          text = element_text(size = 15))
  
  multiplot(g1,g2,g3, cols=3)
}

plotOverlays()  
  # Plots individual sample graphs
  # useful for comparing deviations between all samples tested
plotTables()    
  # Prints results tables from all samples
  # useful if need to copy all datapoints
plotCompiledResults()  
  # Plots means/stdevs for all samples
  # useful table summary when individual datapoints are less
  # important than mean & std.dev
plotBoxplots()
  # Plots boxplots comparing samples & variables
  # useful as quick summation of data and pattern comparison between samples

exportPngSummary <- function(){
  # Export all summary style graphics
  # 
  # Args:
  #   
  # Returns: individual .png files comparing all samples & data
  #
  file.name.png <- paste0(pre.prefix, "-Overlays.png")
  dev.copy(png, file = file.name.png,
           width = 1200,
           height = 900)
  plotOverlays()  
  dev.off()
  
  file.name.png <- paste0(pre.prefix, "-Tables.png")
  dev.copy(png, file = file.name.png,
           width = 1200,
           height = 900)
  plotTables()    
  dev.off()
  
  file.name.png <- paste0(pre.prefix, "-SummaryTable.png")
  dev.copy(png, file = file.name.png,
           width = 500,
           height = 500)
  plotCompiledResults()  
  dev.off()
  
  file.name.png <- paste0(pre.prefix, "-Boxplots.png")
  dev.copy(png, file = file.name.png,
           width = 1600,
           height = 800)
  plotBoxplots()
  dev.off()
}

if(autosave1 == TRUE){
  exportPngSummary()
}
