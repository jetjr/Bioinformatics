#!/usr/bin/env Rscript

library(ggplot2)

setwd("/rsgrps/bhurwitz/jetjr/scripts/bioinformatics/")

df <- read.delim("s_aureus_ref.coverage")
names(df) <- c("seq", "Position", "Depth")

cbPalette <- c("#605856", "#1C6E8C", "#274156")
p1 <- ggplot(df, aes(x=Position, y=Depth)) + geom_area()
#p1 <- p1 + facet_wrap(~ name, ncol=1, scales = "free")
p1 <- p1 + ggtitle("Coverage") + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"))
p1 <- p1 + scale_fill_manual(values=cbPalette) + guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"), text = element_text(size=14))
p1
ggsave("coverage_graph.png")

