#SCRIPT GENERATES JPEG IMAGES OF BAR ABUNDANCE GRAPHS (CLARK DATA)

library(ggplot2)

dat.dir <- "/Users/jetjr/neutropenicfever/CLARK_Results/"

file.names <- dir(dat.dir, pattern="abundance.csv")

rel_abundance <- function (data) {
  data_title <- gsub(".csv", "", data)
  data <- read.csv(paste0(dat.dir, data))
  names(data) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
  data$Proportion_Classified <- as.numeric(as.character(data$Proportion_Classified))
  data_subset <- subset(data, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
  data_p <- ggplot(data_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle(data_title) + ylab("Kmer Count")
  data_p <- data_p + theme(axis.text.x=element_text(angle =-65, hjust = 0))
  ggsave(filename=paste0(data_title,".jpeg"), path=dat.dir)
}

#Preferred method to generate jpeg images of bar abundance graphs
lapply(file.names, rel_abundance)

#Alternate method
for (i in file.names) {
  data_name <- gsub(".csv", ".jpeg", i)
  rel_abundance(i)
  ggsave(filename=data_name, path=dat.dir, device = "jpeg")
}


#SCRATCH ----------------------------------------------------------------------------------------
ion001_abundance <- read.csv(paste0(dat.dir, "IonXpress_001_full_abundance.csv"))
colnames(ion001_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion001_abundance$Proportion_Classified <- as.numeric(as.character(ion001_abundance$Proportion_Classified))
ion001_subset <- subset(ion001_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion001_p <- ggplot(ion001_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress001 Abundance Plot") + ylab("Kmer Count")
ion001_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))

ion002_abundance <- read.csv(paste0(dat.dir, "IonXpress_002_full_abundance.csv"))
colnames(ion002_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion002_abundance$Proportion_Classified <- as.numeric(as.character(ion002_abundance$Proportion_Classified))
ion002_subset <- subset(ion002_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion001_p <- ggplot(ion002_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress002 Abundance Plot") + ylab("Kmer Count")
ion001_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))

ion017_abundance <- read.csv(paste0(dat.dir, "IonXpress_017_full_abund2.csv"))
colnames(ion017_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion017_abundance$Proportion_Classified <- as.numeric(as.character(ion017_abundance$Proportion_Classified))
ion017_subset <- subset(ion017_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion017_p <- ggplot(ion017_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress017 Abundance Plot") + ylab("Kmer Count")
ion017_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))

ion018_abundance <- read.csv(paste0(dat.dir, "IonXpress_018_full_abund.csv"))
colnames(ion018_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion018_abundance$Proportion_Classified <- as.numeric(as.character(ion018_abundance$Proportion_Classified))
ion018_subset <- subset(ion018_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion018_p <- ggplot(ion018_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress018 Abundance Plot") + ylab("Kmer Count")
ion018_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))

ion019_abundance <- read.csv(paste0(dat.dir, "IonXpress_019_full_abund.csv"))
colnames(ion019_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion019_abundance$Proportion_Classified <- as.numeric(as.character(ion019_abundance$Proportion_Classified))
ion019_subset <- subset(ion019_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion019_p <- ggplot(ion019_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress019 Abundance Plot") + ylab("Kmer Count")
ion019_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))

ion058_abundance <- read.csv(paste0(dat.dir, "IonXpress_058_full_abundance.csv"))
colnames(ion058_abundance) <- c("TaxaName", "Lineage", "Count", "Proportion_All", "Proportion_Classified")
ion058_abundance$Proportion_Classified <- as.numeric(as.character(ion058_abundance$Proportion_Classified))
ion058_subset <- subset(ion058_abundance, Proportion_Classified > 1, select = c("TaxaName", "Count", "Proportion_Classified"))
ion058_p <- ggplot(ion058_subset, aes(x=TaxaName, y=Count)) + geom_bar(stat="identity") + ggtitle("IonXpress058 Abundance Plot") + ylab("Kmer Count")
ion058_p + theme(axis.text.x=element_text(angle = -65, hjust = 0))