#CENTRIFUGE TAXONOMIC DATA
#PART 1: GENERATE BARPLOTS OF PROPORTIONAL AND SCORED ABUNDANCES FOR EACH SAMPLE AGAINST BOTH NT AND BHV INDEX

library(ggplot2)
library(lattice)
library(ggrepel)

centnt.dir <- "/Users/jetjr/neutropenicfever/centrifuge-nt/"
centbhv.dir <-"/Users/jetjr/neutropenicfever/centrifuge-bhv/"
clark.dir <-"/Users/jetjr/neutropenicfever/CLARK-workflow/"

file.names <- dir(centbhv.dir, pattern="fa.tsv")

gen_abundance <- function (data) {
  data_title <- gsub(".fa.tsv", "", data)
  data <- read.delim(paste0(i, data))
  total_reads <- sum(data$numReads)
  proportion_classified <- data$numReads / total_reads
  data["proportion_classified"] <- proportion_classified 
  read_subset <- subset(data, proportion_classified > 0.01, select = c("name", "numReads", "proportion_classified"))
  read_subset$numReads <- as.numeric(read_subset$numReads)
  read_p <- ggplot(read_subset, aes(x=name, y=proportion_classified)) + geom_bar(stat="identity") + ylab("Proportion Classified") + xlab("Taxonomic Classification") + ggtitle(paste0("Centrifuge Proportion Classified: ", data_title))
  read_p <- read_p + theme(axis.text.x=element_text(angle =-65, hjust = 0))
  ggsave(filename=paste0(data_title,"_proportional.jpeg"), path=i)
  abund_subset <- subset(data, abundance > 0.01, select =c("name", "abundance"))
  abund_p <- ggplot(abund_subset, aes(x=name, y=abundance)) + geom_bar(stat="identity") + ylab("Abundance") + xlab("Classification") + ggtitle(paste0("Centrifuge Abundance Value: ", data_title))
  abund_p <- abund_p + theme(axis.text.x=element_text(angle =-65, hjust = 0))
  ggsave(filename=paste0(data_title,"_abundance.jpeg"), path=i)
}

for (i in c(centnt.dir, centbhv.dir)) {
  lapply(file.names, gen_abundance)
}

#-----------------------------------------------------------------------------------------------
#PART 2: DOUBLE BAR PLOT COMPARISONS OF DIFFERENT CENTRIFUGE INDEXES 
comp_abundance <- function (data) {
  data_title <- gsub(".fa.tsv", "", data)
  data1 <- read.delim(paste0(centnt.dir, data))
  total_reads <- sum(data1$numReads)
  proportion_classified1 <- data1$numReads / total_reads
  data1["proportion_classified"] <- proportion_classified1 
  data_subset1 <- subset(data1, proportion_classified1 > 0.01, select = c("name", "proportion_classified"))
  data_subset1.1 <- subset(data1, abundance > 0.01, select = c("name", "abundance"))
  data_subset1["index"] <- "centnt"
  data_subset1.1["index"] <- "centnt"

  data2 <- read.delim(paste0(centbhv.dir, data))
  total_reads <- sum(data2$numReads)
  proportion_classified2 <- data2$numReads / total_reads
  data2["proportion_classified"] <- proportion_classified2 
  data_subset2 <- subset(data2, proportion_classified2 > 0.01, select = c("name", "proportion_classified"))
  data_subset2.1 <- subset(data2, abundance > 0.01, select = c("name", "abundance"))
  data_subset2["index"] <- "centbhv"
  data_subset2.1["index"] <- "centbhv"

  bind <- rbind(data_subset1, data_subset2)
  bind2 <- rbind(data_subset1.1, data_subset2.1)
  colors=c("ivory1", "deepskyblue3")
  jpeg(filename=paste0("/Users/jetjr/neutropenicfever/cent_comparison/",data_title,"_prop_nt_v_bhv", ".jpeg"), width = 750, height = 750)
  prop_p <- barchart(proportion_classified~name, data=bind, main=paste0("Read Proporitional Classification Comparision: bhv vs nt index (", data_title,")"),
         groups=index, scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)), auto.key=list(space="top", columns=2, cex.title=1), 
         par.settings=list(superpose.polygon=list(col=colors)))
  print(prop_p)
  dev.off()
  jpeg(filename=paste0("/Users/jetjr/neutropenicfever/cent_comparison/",data_title,"_abund_nt_v_bhv", ".jpeg"), width = 750, height = 750)
  abund_p <- barchart(abundance~name, data=bind2, main=paste0("Abundance Score Classification Comparision: bhv vs nt index (", data_title,")"),
                     groups=index, scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)), auto.key=list(space="top", columns=2, cex.title=1), 
                     par.settings=list(superpose.polygon=list(col=colors)))
  print(abund_p)
  dev.off()
  
}

lapply(file.names, comp_abundance)

#-----------------------------------------------------------------------------------------
#PART 3: Centfrifuge bhv vs CLARK BARCHARTS
#INPUT RESPECTIVE FILE BASENAMES OF INTEREST

file.base <- gsub(".fa.tsv", "", file.names)

clarkcent <- function (data) {

clark <- read.csv(paste0(clark.dir, data, ".abundance.csv"))
names(clark) <- c("name", "lineage", "count", "proportion_all", "proportion_classified")
clark$proportion_classified <- as.numeric(as.character(clark$proportion_classified))
clark_subset <- subset(clark, proportion_classified > 1, select = c("name", "proportion_classified"))
clark_subset["index"] <- "CLARK"


centrifuge <- read.delim(paste0(centbhv.dir, data, ".fa.tsv"))
total_reads <- sum(centrifuge$numReads)
proportion_classified <- ( centrifuge$numReads / total_reads ) * 100
centrifuge["proportion_classified"] <- proportion_classified
cent_subset <- subset(centrifuge, proportion_classified > 1, select = c("name", "proportion_classified"))
cent_subset["index"] <- "Centrifuge"

combine <- rbind(clark_subset, cent_subset)
combine$name <- reorder(combine$name, combine$proportion_classified)

jpeg(filename=paste0("/Users/jetjr/neutropenicfever/clarkvcent_comparison/pub/",data,"_clark_v_cent", ".jpeg"), width = 300, height = 300)
clark_p <- barchart(proportion_classified~name, data=combine, main=paste0("Read Proporitional Classification Comparision: CLARK vs Centrifuge (", data,")"),
                   groups=index, scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)), ylab = "Proportion Classified", xlab="Taxonomy", auto.key=list(space="top", columns=2, cex.title=1), 
                   par.settings=my.theme)
print(clark_p)
dev.off()
}

lapply(file.base, clarkcent)

#---------------------------------------------------------------------------
#PART 4: COMPARATIVE SCATTERPLOTS
#PLOT ABUNDANCE WITH NUMREADS AS WEIGHT

data <- "IonXpress_058_NF001"
readvabund <- function (data) {
  data_title <- gsub(".fa.tsv", "", data)
  data1 <- read.delim(paste0(centbhv.dir, data))
  #remove abundance of 0
  data1_subset <- subset(data1, select = c("name", "numReads", "abundance"))
  jpeg(filename=paste0("/Users/jetjr/neutropenicfever/read_v_abund/",data_title,"_read_v_abund", ".jpeg"), width = 700, height = 600)
  p1 <- ggplot(data1_subset, aes(x = abundance, y = numReads, label = as.character(name))) + geom_text_repel(data=subset(data1_subset, abundance > 0.02 ),aes(label=as.character(name)), 
          size=4, box.padding = unit(1.4, 'lines'), point.padding = unit(2.0, 'lines'), color="black") + geom_point(color="dimgrey") + ggtitle(paste0("Centrifuge Reads vs Abundance Score: ", data_title)) + ylab("Number of Reads") + xlab ("Abundance Score")
  p1 <- p1 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
  p1 <- p1 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))
  p1 <- p1 + theme(axis.text.x = element_text(angle = 00, size=10,color="grey12"))
  p1 <- p1 + theme(axis.text.y = element_text(angle = 00, size=10,color="grey12"))
  p1 <- p1 + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color ="black")) + theme_bw()
  print(p1)
  dev.off()
}
lapply(file.names, readvabund)

#THEME

font.settings <- list(font = 3, fontfamily = "serif")
my.theme <- list(axis.text = font.settings,
                 superpose.polygon=list(col=colors))


