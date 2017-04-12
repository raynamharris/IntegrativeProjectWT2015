colorvalAPA <- c("#404040", "#f4a582", "#ca0020")
#ca0020 ## red - conflict 
#f4a582 ## pink - same
#404040 ## darkgrey - yoked

colorvalAPA2 <- c("##bababa", "#404040","#f4a582", "#ca0020")
#ca0020 ## red - conflict 
#f4a582 ## pink - same
#bababa ## grey - yoked to same
#404040 ## darkgrey - yoked to conflict

colorvalPunch <- c("#d95f02","#7570b3", "#1b9e77")



## pheatmap
colorpalette <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )

ann_colors1 = list(
  APA =  c(Yoked = (values=c("#404040")), Same = (values=c("#f4a582")), Conflict = (values=c("#ca0020"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))

