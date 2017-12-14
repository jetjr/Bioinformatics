#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)

#SETWD: Location of centrifuge_report.tsv files. Should all be in same directory
setwd("/Users/jetjr/neutropenicfever/BMT/new/clark/")

#OUTPUT Directory: Location to store bubble plot and summary data file
out.dir <- "/Users/jetjr/neutropenicfever/"

temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)
temp = sapply(strsplit(temp, "_"), "[", 3)
sample_names <- as.list(sub(".fasta.abundance.csv", "", temp))
myfiles = Map(cbind, myfiles, sample = sample_names)

#Proportion calculations: Each species "Number of Unique Reads" is divided by total "Unique Reads"
props1 = lapply(myfiles, function(x) { 
  names(x) <- c("name", "Lineage", "Count", "Proportion_All", "Percent", "sample")
  x$Percent <- as.numeric(as.character(x$Percent))
  x$method <- "CLARK"
  return(x[,c("name", "Percent", "sample", "method")])
})

## START CENTRIFUGE PARSING ##
setwd("/Users/jetjr/neutropenicfever/centrifuge/post-qc/")

temp = list.files(pattern="*report.tsv")
myfiles = lapply(temp, read.delim)
temp = sapply(strsplit(temp, "_"), "[", 3)
sample_names <- as.list(sub("-qc-centrifuge", "", temp))
myfiles = Map(cbind, myfiles, sample = sample_names)

#Filter settings, default is to remove human and synthetic constructs
filter <- llply(myfiles, subset, name != "Homo sapiens")
filter2 <- llply(filter, subset, name != "synthetic construct")

#Proportion calculations: Each species "Number of Unique Reads" is divided by total "Unique Reads"
props2 = lapply(filter2, function(x) { 
  x$Percent <- (x$numUniqueReads / sum(x$numUniqueReads)) * 100
  x$method <- "Centrifuge"
  return(x[,c("name","Percent","sample", "method", "numUniqueReads", "numReads")])
})

#Final dataframe created for plotting, can change proportion value (Default 1%)
final1 <- llply(props1, subset, Percent > 1)
df1 <- ldply(final1, data.frame)

#Final dataframe created for plotting, can change proportion value (Default 1%)
final2 <- llply(props2, subset, Percent > 2)
df2 <- ldply(final2, data.frame)

combine <- rbind(df1, df2)
combine$name <- as.character(combine$name)
combine$name[2] <- "Human parvovirus B19"
combine$name[8] <- "Torque teno virus"
combine$name[9] <- "Torque teno virus"
combine$name[36] <- "Torque teno virus"
combine$name[37] <- "Torque teno virus"
combine$name[38] <- "Torque teno virus"

combine$sample <- factor(combine$sample, levels(combine$sample)[c(5,6,1,2,3,4)])

names(df) <- c("x", "Proportion", "z")

#SCATTER PLOT WITH POINT SIZE
df <- read.csv("/Users/jetjr/neutropenicfever/BMT/centrifugeBMT_bubble.csv")
df["group"] <- c(rep("1", 7), rep("2", 8), rep("3", 13), rep("4", 11), rep("5", 5), rep("6", 8), rep("7", 5), rep("8", 7), rep("9", 6), rep("10",3))
df["color"] <- c(rep("red", 10), rep("blue", 34), rep("yellow", 23), rep("black", 9), rep("orange", 3), rep("purple", 5), rep("tomato", 4), rep("deeppink", 2), rep("turquoise2", 2), rep("coral", 6))
#Set file name and bubble plot title. Stored in out.dir
file_name <- "Clark_Dot"
plot_title <- "Figure 2: Centrifuge vs CLARK Classification"

png(filename=paste0(out.dir, paste0(file_name,".png")), width = 1000, height = 800)
p2 <- ggplot(df1, aes(as.factor(sample), as.factor(name))) + geom_point(aes(size = Percent))
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=15))
p2 <- p2 + labs(y = "Organism", x = "Method")
p2 <- p2 + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
p2 <- p2 + guides(colour=F)
#p2 <- p2 + facet_grid(. ~ sample)
print(p2)
dev.off()

write.csv(df, file = paste0(out.dir, file_name, ".csv"))

