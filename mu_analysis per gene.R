getwd()
setwd("C:\\Users\\emoraa\\OneDrive - Kemri Wellcome Trust\\Desktop\\delta")
#required library
require(ggplot2)
require(dplyr)
require(tidyverse)
require(Biostrings)
require(stringr)
###loading data
mutations <- read.csv("myresults.csv")
head(mutations)

Lineages <- read.csv("delta_lineages.csv")
plots <- merge (mutations, Lineages, by = "seqname")
write.csv(plots, "plots.csv")
head(mutations)
View(plots)
lineage <- c( "AY.10","AY.11", "AY.12","AY.16","AY.28","AY.3","AY.33" ,"AY.39","AY.39.1","AY.4","AY.4.2" ,"AY.40","AY.41" ,"AY.5","AY.7","B.1.617.2")
plot1 <-mutations %>% filter(protein== "NSP3")
ggplot(plot1, aes(x= variant)) +
  geom_bar(position="dodge", stat= "count", fill=lineage)



+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


