#!/usr/bin/env Rscript

#-------------EDIT HERE----------
cent.dir <- "/rsgrps/bh_class/username/taxonomy/"
out.dir <- "/rsgrps/bh_class/username/taxonomy/barplots/"
#--------------------------------

file.names <- dir(cent.dir, pattern="-centrifuge_report.tsv")

gen_barplot <- function (data) {
  data_title <- gsub("-centrifuge_report.tsv", "", data)
  data <- read.delim(paste0(i, data))
  total_reads <- sum(data$numReads)
  proportion_classified <- data$numReads / total_reads
  data["proportion_classified"] <- proportion_classified
  read_subset <- subset(data, proportion_classified > 0.005, select = c("name", "numReads", "proportion_classified"))
  read_subset$numReads <- as.numeric(read_subset$numReads)
  png(filename=paste0(out.dir,data_title,"_taxonomy.png"), width = 600, height = 600)
  op <- par(mar=c(15, 8, 4, 2) + 0.1, mgp = c(10, 1, 0))
  p1 <- barplot(read_subset$proportion_classified, main=paste0("Read Proportional Classification: ",data_title), names.arg = read_subset$name, las=2, cex.names = 1, cex.axis = 1, ylab="Proportion Classified", ylim = c(0, 0.90))
  grid(nx=NA, ny=NULL)
  print(p1)
  dev.off()
}

for (i in cent.dir) {
  lapply(file.names, gen_barplot)
}
