#!/usr/bin/env Rscript

library(optparse)
library(plyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(
    c("-d", "--dir"),
    default = "",
    type = "character",
    help = "Centrifuge outdir",
    metavar = "character"
  ),
  make_option(
    c("-o", "--outdir"),
    default = file.path(getwd(), "plots"),
    type = "character",
    help = "out directory",
    metavar = "character"
  )
);

opt_parser = OptionParser(option_list = option_list);
opt        = parse_args(opt_parser);
cent.dir   = opt$dir
out.dir    = opt$outdir

#SETWD: Location of centrifuge_report.tsv files. Should all be in same directory
setwd(cent.dir)

temp = list.files(pattern="*report.tsv")
myfiles = lapply(temp, read.delim)
sample_names <- as.list(sub("-centrifuge_report.tsv", "", temp))
myfiles = Map(cbind, myfiles, sample = sample_names)

#Proportion calculations: Each species "Number of Unique Reads" is divided by total "Unique Reads". Get hit ratios. 

props = lapply(myfiles, function(x) { 
  x$proportion <- ((x$numUniqueReads / sum(x$numUniqueReads)) * 100)
  x$abundance <- x$abundance * 100
  x$hitratio <- x$numUniqueReads / x$numReads
  return(x[,c("name","proportion", "abundance", "genomeSize", "sample", "numReads", "numUniqueReads", "taxID", "hitratio")])
})

#Final dataframe created for plotting, can change proportion value (Default 1%)
final <- llply(props, subset, abundance > 0.01)
final2 <- llply(final, subset, proportion > 0.01)
df <- ldply(final2, data.frame)

error <- qt(0.975,df=length(df$hitratio)-1)*sd(df$hitratio)/sqrt(length(df$hitratio))
low <- mean(df$hitratio)-error


exclude = lapply(final2, function(x){
  ids <- c(x$taxID[which(x$hitratio < low)], "9606", "32630")
  cat(ids, file=paste0(x$sample[1],"_exclude"), sep = ",")
})

prefer = lapply(final2, function(x){
  ids <- c(x$taxID[which(x$hitratio > low)])
  cat(ids, file=paste0(x$sample[1],"_prefer"), sep = ",")
})

write.csv(df, file = paste0(out.dir, "original_report.csv"))
