pkgs = c("tidyverse", "knitr", "rmarkdown", "kableExtra",
         "forcats", 
         "cowplot", "factoextra", "FactoMineR", "apaTables", "png", "grid",
         "Rtsne", 
         "corrr",
         "devtools", "BiocManager")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("BiocParallel")

