library(shiny)
library(shiny)
library(tidyverse)
library(cowplot)
library(corrplot)

df <- read_csv("data.csv") 
gene_names <- names(df[8:20])