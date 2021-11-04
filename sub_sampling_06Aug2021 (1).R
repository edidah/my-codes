
dev.off()

library(tidyverse); library(stringr); library(lubridate); library(scales); library(janitor)

setwd("~/Documents/Policy_Briefs_Sequencing_Reportts/Country/Seychelles/ML/GLOBAL _SEQUENCES")

my_dta <- as_tibble(read_tsv("~/Documents/Policy_Briefs_Sequencing_Reportts/Country/Seychelles/ML/GLOBAL _SEQUENCES/metadata.tsv"))


gisaid_dta <- my_dta%>%
  rename("Pango_Lineage"="Pango lineage", id="Virus name")%>%
  select(-c("Type", "Additional location information", "Patient age"))%>%
  separate(Location, into=c("continent", "country", "region", "subregion1", "subregion2", "subregion3"), sep="\\/", remove=FALSE, extra="warn", fill="right")%>%
  mutate(datecollection=as.Date(`Collection date`, format="%Y-%m-%d"))%>%
  filter(!is.na(datecollection))%>%
  mutate(monthcollection=as.Date(cut(datecollection, breaks="month")))%>%### Modify to month year
  mutate(yearcollection=as.Date(cut(datecollection, breaks="year")))%>%### 
  mutate(sample_id=paste(id, `Collection date`,`Submission date`, sep="|"))

head(gisaid_dta)
names(gisaid_dta)
str(gisaid_dta)

subsampling_VOC <-gisaid_dta%>%
  filter(datecollection > as.Date("2020-06-11") & datecollection < as.Date("2021-06-04"))%>%
  filter(Pango_Lineage=="B.1.617.2"|Pango_Lineage=="B.1.1.7"|Pango_Lineage=="B.1.351"|Pango_Lineage=="B.1"|
           Pango_Lineage=="B.1.333"|Pango_Lineage=="B.1.596"|Pango_Lineage=="B.1.1.50"|
           Pango_Lineage=="B.1.1"|Pango_Lineage=="B.1.535"|Pango_Lineage=="B.1.612")%>%
  group_by(continent, yearcollection, monthcollection, Pango_Lineage)%>%
  sample_n(5, replace =T)%>%
  distinct(id, .keep_all=TRUE)%>%
  ungroup()%>%
  select(sample_id)




my_lineages <-c("B.1.617.2", "B.1.617.1","B.1.1.7","B.1.351","B.1", "B.1.333","B.1.596","B.1.1.50","B.1.1","B.1.535","B.1.612")
#tabyl(subsampling_VOC, Pango_Lineage)

#write.table(subsampling_VOC, file = "subsampled_seychelles2.txt", row.names = F, na = "", quote = FALSE, col.names = F)



sequence_dta <- as_tibble(read.csv("~/Dropbox/COVID-19/writeups/Seychelles/subsample1/subsampledout_seychelles_info_details.csv"))%>%
  distinct(seq_id, .keep_all=T)%>%
  rename(sample_id="seq_id")



lineages_dta <- as_tibble(read.csv("~/Dropbox/COVID-19/writeups/Seychelles/subsample1/lineage_report.csv"))%>%
  distinct(taxon, .keep_all=T)%>%
  rename(sample_id="taxon")

seq_lineage <- merge(lineages_dta,sequence_dta, by=c("sample_id"), all = T )%>%
  filter(actualseq!="")%>%
  separate(sample_id, into=c("id", "datecollection", "datesubmission"), sep="\\|", remove=FALSE, extra="warn", fill="right")%>%
  separate(id, into=c("virus", "country", "seqid", "year", "other"), sep="\\/", remove=FALSE, extra="warn", fill="right")%>%
  filter(lineage%in%my_lineages)%>%
  mutate(taxon_id=paste(lineage, country, seqid, datecollection, sep="|"))%>%
  arrange(lineage,country, datecollection)%>%
  select(taxon_id,actualseq)

#tabyl(seq_lineage,country)
names(seq_lineage)

write.csv(seq_lineage, file="seq_lineage_seyschelles_06Aug2021.csv", row.names = F, na="")
  