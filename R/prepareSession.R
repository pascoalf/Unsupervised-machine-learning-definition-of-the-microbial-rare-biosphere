#
library(ulrb)
#
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(purrr)
library(vegan)
library(readxl)
library(gghalves)
library(gridExtra)
library(RColorBrewer)
library(cluster)
#library(dada2)
library(microbenchmark)
library(FuzzyQ)
##
qualitative_colors <- 
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## 
##
years_colors <- RColorBrewer::brewer.pal(n = 6, "Blues")
##
set.seed(123)

