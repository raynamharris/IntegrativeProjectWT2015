Rationale
=========

This script uses the "colorblindr" library to see how my figures look
through color blind simulators.

Setup
-----

    #devtools::install_github("wilkelab/cowplot")
    #install.packages("colorspace", repos = "http://R-Forge.R-project.org")
    #devtools::install_github("clauswilke/colorblindr")
    #install.packages("magick")

    library(ggplot2)
    library(colorblindr)
    library(cowplot)
    library(magick)

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/99_colorblindr/')

Hippocampus color scheme
------------------------

    imgDG <- magick::image_read("../figures/00_schematics/DG.png")
    p <- ggdraw() + draw_image(imgDG) # turn png into ggplot object
    p2 <- edit_colors(p, tritan, sev = 1.0)
    p3 <- edit_colors(p, deutan, sev = 1.0)
    p4 <- edit_colors(p, protan, sev = 1.0)

    imgCA3 <- magick::image_read("../figures/00_schematics/CA3.png")
    q <- ggdraw() + draw_image(imgCA3) # turn png into ggplot object
    q2 <- edit_colors(q, tritan, sev = 1.0)
    q3 <- edit_colors(q, deutan, sev = 1.0)
    q4 <- edit_colors(q, protan, sev = 1.0)

    imgCA1 <- magick::image_read("../figures/00_schematics/CA1.png")
    r <- ggdraw() + draw_image(imgCA1) # turn png into ggplot object
    r2 <- edit_colors(r, tritan, sev = 1.0)
    r3 <- edit_colors(r, deutan, sev = 1.0)
    r4 <- edit_colors(r, protan, sev = 1.0)

    plot_grid(p, p2, p3, p4, q, q2, q3, q4, r, r2, r3, r4, nrow=3)

![](../figures/99_colorblindr/colorfulhippocampi-1.png)

Figure 1: Behavior
------------------

    fig1 <- magick::image_read("../figures/figures-01.png")
    f <- ggdraw() + draw_image(fig1) # turn png into ggplot object
    f2 <- edit_colors(f, tritan, sev = 1.0)
    f4 <- edit_colors(f, protan, sev = 1.0)
    plot_grid(f2, f4, nrow=1)

![](../figures/99_colorblindr/fig1-1.png)

Figure 2: RNA-seq all data
--------------------------

    fig2 <- magick::image_read("../figures/figures-02.png")
    f <- ggdraw() + draw_image(fig2) # turn png into ggplot object
    f2 <- edit_colors(f, tritan, sev = 1.0)
    f4 <- edit_colors(f, protan, sev = 1.0)
    plot_grid(f2, f4, nrow=1)

![](../figures/99_colorblindr/fig2-1.png)

Figure 3: RNA-seq data subset by subfiel
----------------------------------------

    fig3 <- magick::image_read("../figures/figures-03.png")
    f <- ggdraw() + draw_image(fig2) # turn png into ggplot object
    f2 <- edit_colors(f, tritan, sev = 1.0)
    f4 <- edit_colors(f, protan, sev = 1.0)
    plot_grid(f2, f4, nrow=1)

![](../figures/99_colorblindr/fig3-1.png)
