'22-04-2021

Edidah Moraa

This script does takes GISAID metadata, cleans and maps the various sequencing platforms to respective countries globally using the rworldmap package.
The input files include 
  a) GISAID metadata.csv
  b) GISAID sequencing platform ,
  c) World co-ordinates ISO-2 codes.tsv'


#install.packages("pacman")
library(pacman)
#loading libraries
p_load(tidyverse, patchwork,stringr,rworldmap,rnaturalearthdata,rnaturalearth,classInt,RColorBrewer,ggplot)


##Generating data for downloading using the gsaid metadata file
#loading the gisaid metadata file
my_dta <- as_tibble(read_tsv("gsaid_datasets.csv"))

## slicing the data to sections
groups <- (split(my_dta, (seq(nrow(my_dta))-1) %/% 9500)) #here I want 9500 rows per file until last row is reached

for (i in seq_along(groups)) {
  write.csv(groups[[i]], paste0("sample_output_file", i, ".csv")) #iterate and write file
}


#loading data
seq_platform <- as_tibble(read_csv("seqtech_dataset.csv"))

#data preparation

##data subsetting based on relevant column
seq_platform <-seq_platform %>% select ("Virus name","Accession ID","Location","Sequencing technology")

##column renaming
nms <-c("virus_name","accession_id","location","seq_technology")
seq_platform <- set_names(seq_platform, nms)

## duplicates removal
seq_platform <- seq_platform[!duplicated(seq_platform$accession_id), ]

##split the location column to 
seq_platform_det <- seq_platform %>%
  separate(location, into=c("continent", "country", "region", "subregion1", "subregion2"), sep="\\/", remove=FALSE, extra="warn", fill="right")%>%
  select("virus_name","accession_id","seq_technology", "continent", "country", "region", "subregion1")

##writeback a tsv file(for storage)
write_tsv(seq_platform_det, "seq_platform.tsv")
rm(seq_platform)

#load the tsv file to r and clean to ensure column uniformity on the seq_technology platform,country. Select for necessary columns.
seq_platform <- as_tibble(read_tsv("seq_platform.tsv"))
seq_platform <- seq_platform %>% mutate(tech = case_when( str_detect(string=seq_technology, pattern = "novaseq|noveseq|NovaSeq|Miniseq|Miseq|MiSeq Dx|NextSeq|Nextseq|Nextseq 500|MiSeq|iSeq") ~ "illumina",
                                                          str_detect(string=seq_technology, pattern = "Illumina|ILLUMINA|illumina|lllumina|Ilumina") ~ "illumina",
                                                          str_detect(string=seq_technology, pattern = "covidseq|DRAGEN") ~ "illumina",
                                                          str_detect(string=seq_technology, pattern = "Oxford|oxford|ONT")  ~ "ONT",
                                                          str_detect(string=seq_technology, pattern = "Nanopore|NANOPORE")  ~ "ONT",
                                                          str_detect(string=seq_technology, pattern = "minION|MinION|MinIon|Minion")  ~ "ONT",
                                                          str_detect(string=seq_technology, pattern = "Nanopore")  ~ "ONT",
                                                          str_detect(string=seq_technology, pattern = "Nanopore MinION")  ~ "ONT",
                                                          str_detect(string = seq_technology, pattern = "GridION|GridIon|minIon") ~ "ONT",
                                                          str_detect(string = seq_technology, pattern = "promethion/Promethion") ~ "ONT",
                                                          str_detect(string = seq_technology, pattern = "ion torrent|IonTorrent|Ion Torrent|Ion torrent|ION_TORRENT") ~ "ion torrent",
                                                          str_detect(string = seq_technology, pattern = "s5|S5") ~ "ion torrent",
                                                          str_detect(string = seq_technology, pattern = "bgi|mgi|BGI|MGI|DNBSEQ") ~ "mgi",
                                                          str_detect(string = seq_technology, pattern = "pacbio|Pacbio|PacBio|Sequel|PACBIO") ~ "pacbio",
                                                          str_detect(string = seq_technology, pattern = "clearLabs|ClearLabs|Clear Labs Clear Dx|Clear Labs|Clear Dx System") ~ "clearlabs",
                                                          str_detect(string = seq_technology, pattern = "sanger|Sanger|consensus|Consensus") ~ "unknown",
                                                          str_detect(string = seq_technology, pattern = "unknown|Sequencing technology|") ~ "unknown",
                                                        
                                                          ))
seq_platform$country<- tolower(seq_platform$country)
seq_platform <- seq_platform %>% select(-c(subregion1,virus_name))

#Write back a tsv file of the final data frame.
write_tsv(cleaned_seqplatform, "cleaned_seqplatform.tsv")
cleaned_seqplatform <- read_tsv("cleaned_seqplatform.tsv")

##harmonize erroneous datasets
seq_platform <- cleaned_seqplatform %>% mutate(country= if_else(accession_id == "EPI_ISL_5066147", "india", country),
                                         country = if_else(accession_id == "EPI_ISL_5066158", "india", country),
                                         country = if_else(accession_id == "EPI_ISL_6885322", "belgium", country),
                                         country = if_else(accession_id == "EPI_ISL_4252028", "kenya", country),
                                         country = if_else(accession_id == "EPI_ISL_3152245", "belgium", country),
                                         country = if_else(accession_id == "EPI_ISL_8142215", "india", country))


#import the coordinate data set, cleaning, selecting and write back tsv file
coordinates <- read_tsv("coordinates_1.tsv")
coordinates$country = tolower(coordinates$countrylower)
coordinates <- coordinates%>%select(c(latitude,longitude,country, country_code))

##writeback a tsv
write_tsv(coordinates, "coordinates_1.tsv")


# join the two dataframes: coordinates and sequencing technology dataset
seq_platform_sf<- full_join(seq_platform, coordinates, by = "country")

##filter for unsplit data; lack countrycode resignation followed by dropping the datasets
seq_platform_sf_na <- seq_platform_sf %>% filter(is.na(longitude))
seq_platform_sf <- seq_platform_sf %>% drop_na(country)

#write back a tsv
write_tsv(seq_platform_sf, "seq_platform_sf")
seq_platform_sf <- read_tsv("seq_platform_sf")

####data visualization
 # barplot
totalseq <- read_tsv("seq_platform_sf.tsv")
totalseq1 <- seq_platform_sf %>% mutate(tech = case_when( str_detect(string= tech, pattern = "illumina") ~ "Illumina",
                                                          str_detect(string = tech, pattern = "pacbio") ~ "PacBio",
                                                          str_detect(string=tech, pattern = "clearlabs") ~ "Clearlabs",
                                                          str_detect(string = tech, pattern = "ion torrent") ~ "Ion Torrent",
                                                          str_detect(string=tech, pattern = "mgi") ~ "MGI DNBSEQ",
                                                          str_detect(string = tech, pattern = "ONT") ~ "Oxford Nanopore"
))

# remove Clearlabs row
totalseq1 = subset(totalseq1, totalseq1$tech != "Clearlabs")
library(scales)

# convert no. of sequences to percentage
per_totalseq <- totalseq1 %>% 
  group_by(tech) %>% summarise(total = n())%>% 
  mutate(perc = total / sum(total)*100)

# generate bar plot for distribution of GISAID sequences according to sequencing platforms
plot1 <- ggplot(per_totalseq, aes(x=reorder(tech, -perc), y= perc, fill =tech))+ geom_bar(stat="identity") +
  scale_fill_manual(values = c("coral","red","lightseagreen","skyblue2","mediumvioletred")) +
  xlab("Sequencing platforms") +
  ylab("Proportion of sequences (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0,0),limits=c(0,100)) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour ="black"),
        text = element_text(size = 16))

plot1



## plotting the maps
####illumina plot
Illumina  <- seq_platform_sf %>% filter(tech == "illumina")
Illumina_1 <- Illumina %>% group_by(country) %>% summarise(total = n())
#writeback a tsv file
write_tsv(Illumina_1, "Illumina_1.tsv")
Illumina1 <- read_tsv("Illumina_1.tsv")
Illumina2 <- full_join(coordinates, Illumina1, by = "country")

#defining intervals
###generate categories

Illumina2$totalgrp<-cut(Illumina2$total, breaks = c(-Inf,100,1000,10000,100000,Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
Illumina2 <- Illumina2 %>% filter(!is.na(Illumina2$totalgrp))

#color scheme
colorpalette <- brewer.pal(6,"Oranges")

#plot a map
Illumina <- Illumina2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
Illuminaplatform_map <- mapCountryData(Illumina, nameColumnToPlot="totalgrp", addLegend = TRUE, mapTitle = "Illumina Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(Illuminaplatform_map
                        ,legendLabels="limits"
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))

###ont plot
ONT <- seq_platform_sf %>% filter(tech == "ONT")
ONT_1 <- ONT %>% group_by(country) %>% summarise(total = n())
##writeback a tsv
write_tsv(ONT_1, "ONT_1.tsv")
ONT_1 <- read_tsv("ONT_1.tsv")
ONT2 <- full_join(coordinates, ONT_1, by = "country")

#defining categories
ONT2$totalgrp<-cut(ONT2$total, breaks = c(-Inf,100,1000,10000,100000, Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
ONT2 <- ONT2 %>% filter(!is.na(ONT2$totalgrp))

#color scheme
colorpalette <- brewer.pal(5,"Blues")

#plot a map
ONT <- ONT2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
ONTplatform_map <- mapCountryData(ONT, nameColumnToPlot="totalgrp", addLegend = TRUE, mapTitle = "ONT Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(ONTplatform_map
                        ,legendLabels="limits"
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))
## Ion Torrent platform
Ion_Torrent <- seq_platform_sf %>% filter(tech == "ion torrent")
Ion_Torrent_1 <- Ion_Torrent %>% group_by(country) %>% summarise(total = n())
#writeback a tsv file
write_tsv(Ion_Torrent_1, "Ion_Torrent_1.tsv")
Ion_Torrent1 <- read_tsv("Ion_Torrent_1.tsv")
Ion_Torrent2 <- full_join(coordinates, Ion_Torrent1, by = "country")

##ion torrent map
Ion_Torrent2$totalgrp<-cut(Ion_Torrent2$total, breaks = c(-Inf,100,1000,10000,100000, Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
Ion_Torrent2  <- Ion_Torrent2 %>% filter(!is.na(Ion_Torrent2$totalgrp))

#color scheme
colorpalette <- brewer.pal(5,"Reds")

#plot a map
Ion_Torrent <- Ion_Torrent2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
Ion_Torrentplatform_map <- mapCountryData(Ion_Torrent, nameColumnToPlot="totalgrp", addLegend = TRUE, mapTitle = "Ion_Torrent Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(Ion_Torrentplatform_map
                        ,legendLabels="limits"
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))
##pacbio map

PacBio <- seq_platform_sf %>% filter(tech == "pacbio")
PacBio1 <- PacBio %>% group_by(country) %>% summarise(total = n())
#writeback a tsv file
write_tsv(PacBio1, "PacBio1.tsv")
PacBio1 <- read_tsv("PacBio1.tsv")
PacBio2 <- full_join(coordinates, PacBio1, by = "country")

## categorize 
PacBio2$totalgrp<-cut(PacBio2$total, breaks = c(-Inf,100,1000,10000,100000, Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
PacBio2 <- PacBio2 %>% filter(!is.na(PacBio2$totalgrp))

#color scheme
colorpalette <- brewer.pal(5,"RdPu")

#plot a map
PacBio <- PacBio2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
PacBioplatform_map <- mapCountryData(PacBio, nameColumnToPlot="totalgrp", addLegend = TRUE, mapTitle = "PacBio Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(PacBioplatform_map
                        ,legendLabels="limits"
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))
##mgi
DNBSeq <- seq_platform_sf %>% filter(tech == "mgi")
DNBSeq1 <- DNBSeq %>% group_by(country) %>% summarise(total = n())
#writeback a tsv file
write_tsv(DNBSeq1, "DNBSeq1.tsv")
DNBSeq1 <- read_tsv("DNBSeq1.tsv")
DNBSeq2 <- full_join(coordinates, DNBSeq1, by = "country")

## categorize data
DNBSeq2$totalgrp<-cut(DNBSeq2$total, breaks = c(-Inf,100,1000,10000,100000, Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
DNBSeq2 <- DNBSeq2 %>% filter(!is.na(DNBSeq2$totalgrp))

#color scheme
colorpalette <- brewer.pal(5,"GnBu")

#plot a map
DNBSeq <- DNBSeq2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
DNBSeqplatform_map <- mapCountryData(DNBSeq, nameColumnToPlot="totalgrp", addLegend = TRUE, mapTitle = "DNBSeq Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(DNBSeqplatform_map
                        ,legendLabels="limits"
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))
## clearlabs
clearlabs <- seq_platform_sf %>% filter(tech == "clearlabs")
clearlabs1 <- clearlabs %>% group_by(countrylower) %>% summarise(total = n())
#writeback a tsv file
write_tsv(clearlabs1, "clearlabs1.tsv")
clearlabs <- read_tsv("clearlabs1.tsv")
clearlabs2 <- full_join(coordinates, clearlabs1, by = "countrylower")
#categorize
clearlabs2$totalgrp<-cut(clearlabs2$total, breaks = c(-Inf,100,1000,10000,100000, Inf),labels = c("<100","100-1000","1001-10000","10001-100000",">100000"))
clearlabs2 <- clearlabs2 %>% filter(!is.na(clearlabs2$totalgrp))

#color scheme
colorpalette <- brewer.pal(5,"RdPu")

#plot a map
clearlabs <- clearlabs2 %>%
  joinCountryData2Map(., joinCode = "ISO2" , nameJoinColumn = "country_code")
clearlabsplatform_map <- mapCountryData(clearlabs, nameColumnToPlot="totalgrp", addLegend = TRUE, xlim = NA, ylim = NA,mapTitle = "Clearlabs Sequencing Platform Distribution", colourPalette=colorpalette)
do.call(addMapLegend, c(clearlabsplatform_map
                        ,labelFontSize=0.1
                        ,legendWidth=0.1
                        ,legendShrink=0.1 
                        ,legendMar=3
                        ,horizontal=TRUE
                        ,legendIntervals='page'
))

