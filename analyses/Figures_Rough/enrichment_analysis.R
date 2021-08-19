# Need to do enrichment analysis regarding the R-loops discovered in each technology
library(tidyverse)


load("analyses/rmapfftsmall.rda")

rmapfftsmall
load("../RMapDB-shiny/data/rltab.rda")

rltab2 <- rltab %>%
  mutate(modeLst = map(Modes, function(x) {strsplit(x, split = "\n") %>% unlist()})) %>%
  unnest(cols = modeLst)

rlLst <- rltab2 %>%
  select(id, modeLst) %>%
  group_by(modeLst) %>%
  {setNames(group_split(.), group_keys(.)[[1]])} %>%
  lapply(pull, var = id)

library(SuperExactTest)  
  
res=supertest(rlLst, n=length(unique(unlist(rlLst))))
plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))





