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
    knitr::opts_chunk$set(fig.path = '../figures/09_colorblindr/')

Hippocampus Color scheme
------------------------

    imgDG <- magick::image_read("../figures/00_schematics/DG.png")
    imgCA3 <- magick::image_read("../figures/00_schematics/CA3.png")
    imgCA1 <- magick::image_read("../figures/00_schematics/CA1.png")

    p <- ggdraw() + draw_image(imgDG) # turn png into ggplot object
    p2 <- edit_colors(p, tritan, sev = 1.0)
    p3 <- edit_colors(p, deutan, sev = 1.0)
    p4 <- edit_colors(p, protan, sev = 1.0)

    q <- ggdraw() + draw_image(imgCA3) # turn png into ggplot object
    q2 <- edit_colors(q, tritan, sev = 1.0)
    q3 <- edit_colors(q, deutan, sev = 1.0)
    q4 <- edit_colors(q, protan, sev = 1.0)

    r <- ggdraw() + draw_image(imgCA1) # turn png into ggplot object
    r2 <- edit_colors(r, tritan, sev = 1.0)
    r3 <- edit_colors(r, deutan, sev = 1.0)
    r4 <- edit_colors(r, protan, sev = 1.0)
    #plot_grid(p, p2, p3, p4, q, q2, q3, q4, r, r2, r3, r4, nrow=3)

Figure 1
--------

    fig1 <- magick::image_read("../figures/02_RNAseq_ConsistentYoked/02_RNAseq_ConsistentYoked-01.png")
    fig1

    ##   format width height colorspace filesize
    ## 1    PNG  2072   1389       sRGB   247148

    r <- ggdraw() + draw_image(fig1) # turn png into ggplot object
    r2 <- edit_colors(r, tritan, sev = 1.0)
    r4 <- edit_colors(r, protan, sev = 1.0)
    plot_grid(r2, r4, nrow=2)

![](../figures/09_colorblindr/fig1-1.png)

Figure 2
--------

    fig2 <- magick::image_read("../figures/02_RNAseq_ConsistentConflict/02_RNAseq_Conflict-01.png")
    fig2

    ##   format width height colorspace filesize
    ## 1    PNG  1252   1987       sRGB   237491

    r <- ggdraw() + draw_image(fig2) # turn png into ggplot object
    r2 <- edit_colors(r, tritan, sev = 1.0)
    r4 <- edit_colors(r, protan, sev = 1.0)
    plot_grid(r2, r4, nrow=1)

![](../figures/09_colorblindr/fig2-1.png)

Supplementary behavior
----------------------

    behavior <- magick::image_read("../figures/01_behavior/01_avoidancebehvaior-01.png")
    behavior

    ##   format width height colorspace filesize
    ## 1    PNG  2245   1190       sRGB   226683

    b <- ggdraw() + draw_image(behavior) # turn png into ggplot object
    b2 <- edit_colors(b, tritan, sev = 1.0)
    b4 <- edit_colors(b, protan, sev = 1.0)
    plot_grid(b2, b4, nrow=2)

![](../figures/09_colorblindr/behavior-1.png)

    CA1conflictconsistent <- magick::image_read("../figures/02e_RNAseq_GO/CA1conflictconsistent-1.png")
    CA1consistentyoked <- magick::image_read("../figures/02e_RNAseq_GO/CA1consistentyoked-1.png")
    CA1yokedyoked <- magick::image_read("../figures/02e_RNAseq_GO/CA1yokedyoked-1.png")
    CA3conflictconsistent <- magick::image_read("../figures/02e_RNAseq_GO/CA3conflictconsistent-1.png")
    CA3consistentyoked <- magick::image_read("../figures/02e_RNAseq_GO/CA3consistentyoked-1.png")
    CA3yokedyoked <- magick::image_read("../figures/02e_RNAseq_GO/CA3yokedyoked-1.png")
    DGconflictconsistent <- magick::image_read("../figures/02e_RNAseq_GO/DGconflictconsistent-1.png")
    DGconsistentyoked <- magick::image_read("../figures/02e_RNAseq_GO/DGconsistentyoked-1.png")
    DGyokedyoked <- magick::image_read("../figures/02e_RNAseq_GO/DGyokedyoked-1.png")

    r1 <- ggdraw() + draw_image(CA1consistentyoked) # turn png into ggplot object
    r2 <- ggdraw() + draw_image(CA1conflictconsistent) # turn png into ggplot object
    r3 <- ggdraw() + draw_image(CA1yokedyoked) # turn png into ggplot object

    q1 <- ggdraw() + draw_image(CA3consistentyoked) # turn png into ggplot object
    q2 <- ggdraw() + draw_image(CA3conflictconsistent) # turn png into ggplot object
    q3 <- ggdraw() + draw_image(CA3yokedyoked) # turn png into ggplot object

    p1 <- ggdraw() + draw_image(DGconsistentyoked) # turn png into ggplot object
    p2 <- ggdraw() + draw_image(DGconflictconsistent) # turn png into ggplot object
    p3 <- ggdraw() + draw_image(DGyokedyoked) # turn png into ggplot object


    numentrance2 <- magick::image_read("../figures/01_behavior/numentrance-2.png")
    numentrance3 <- magick::image_read("../figures/01_behavior/numentrance-3.png")
    numentrance4 <- magick::image_read("../figures/01_behavior/numentrance-4.png")

    n2 <- ggdraw() + draw_image(numentrance2) # turn png into ggplot object
    n3 <- ggdraw() + draw_image(numentrance3) # turn png into ggplot object
    n4 <- ggdraw() + draw_image(numentrance4) # turn png into ggplot object

    plot_grid(n2, p1, q1, r1, nrow=2) #consistent

![](../figures/09_colorblindr/GO-1.png)

    plot_grid(n3, p2, q2, r2, nrow=2) #conflcit

![](../figures/09_colorblindr/GO-2.png)

    plot_grid(n4, p3, q3, r3, nrow=2) #yoked

![](../figures/09_colorblindr/GO-3.png)

    plot_grid(p, p1, p2, p3, nrow=2) #DG

![](../figures/09_colorblindr/GO-4.png)

    plot_grid(q, q1, q2, q3, nrow=2) #CA3

![](../figures/09_colorblindr/GO-5.png)

    plot_grid(r, r1, r2, r3, nrow=2) #CA1

![](../figures/09_colorblindr/GO-6.png)
