library(tidyverse)
library(sunburstR)
library(ggplot2)
##loading data
totalmutations <- read.csv("plots.csv")
View(totalmutations)
###spike region mutations
s_mutations <- totalmutations %>%
  filter(protein =='S') %>%
  mutate(
    path = paste(protein, variant, lineage, sep = "-")
  ) %>%
  slice(2:100) %>%
  mutate(
    V2 = 1
  )
sunburst(data = data.frame(xtabs(V2~path, s_mutations)), legend = FALSE,
         colors = c("D99527", "6F7239", "CE4B3C", "C8AC70", "018A9D"))
###ORF3a mutations****
ORF3a <- totalmutations %>%
  filter(protein == "ORF3a") %>%
  mutate(
    path = paste(protein, variant, lineage, sep = "-")
  ) %>%
  slice(2:100) %>%
  mutate(
    V2 = 1)
sunburst(data = data.frame(xtabs(V2~path, ORF3a)), legend = FALSE,
         colors = c("D99527", "6F7239", "CE4B3C", "C8AC70", "018A9D"))
m_mutations <-  totalmutations %>%
  filter(protein == "M") %>%
  mutate(
    path = paste(protein, variant, lineage, sep = "-")
  ) %>%
  slice(2:100) %>%
  mutate(
    V2 = 1)
sunburst(data = data.frame(xtabs(V2~path, m_mutations)), legend = FALSE,
         colors = c("D99527", "6F7239", "CE4B3C", "C8AC70", "018A9D"))
n_mutations <- totalmutations %>%
  filter(protein == "N") %>%
  mutate(
    path = paste(protein, variant, lineage, sep = "-")
  ) %>%
  slice(2:100) %>%
  mutate(
    V2 = 1)
sunburst(data = data.frame(xtabs(V2~path, n_mutations)), legend = FALSE,
         colors = c("D99527", "6F7239", "CE4B3C", "C8AC70", "018A9D"))
