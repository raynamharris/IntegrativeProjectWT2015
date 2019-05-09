pkgs = c("tidyverse", "knitr", "rmarkdown",
         "plyr", "reshape2", "cowplot", "factoextra", "car",  
         "pheatmap", "viridis", "magrittr", "genefilter",  
         "ggrepel","colorblindr", "xtable", "Hmisc", "corrplot",
         "magick", "ggrepel", "stringr", "kableExtra",
         "devtools", "UpSetR","BiocManager")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
devtools::install_github("clauswilke/ggtextures")