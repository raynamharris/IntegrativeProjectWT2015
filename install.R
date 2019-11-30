pkgs = c("tidyverse", "knitr", "rmarkdown", "kableExtra",
         "forcats", 
         "cowplot", "factoextra", "FactoMineR", "apaTables", "png", "grid",
         "Rtsne", 
         "corrr",
         "devtools", "BiocManager",
         "magick", 
         "gsl", "OpenMx", "MBESS")
install.packages(pkgs)

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("BiocParallel")
