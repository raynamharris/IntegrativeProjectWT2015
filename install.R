pkgs = c("tidyverse", "knitr", "rmarkdown", "kableExtra",
         "forcats", 
         "cowplot", "factoextra", "FactoMineR", "apaTables", "png", "grid", "scales",
         "Rtsne", 
         "corrr", "Hmisc", "ggrepel",
         "stringr", 
         "devtools", "BiocManager",
         "magick", 
         "gsl", "OpenMx", "MBESS")
install.packages(pkgs)

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("BiocParallel")
