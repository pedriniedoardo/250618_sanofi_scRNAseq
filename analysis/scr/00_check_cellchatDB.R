# AIM ---------------------------------------------------------------------
# check if the expression of CD40 and its ligand is reported in cellchat


# libraries ---------------------------------------------------------------
library(tidyverse)
library(CellChat)


# load data ---------------------------------------------------------------
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# wrangling ---------------------------------------------------------------
CellChatDB$interaction %>%
  filter(receptor %in% c("CD40")) %>%
  as.list()



