colorvalPunch <- c("#d95f02","#1b9e77", "#7570b3")
## DG  "#d95f02"
## CA3 "#1b9e77"
## CA1 "#7570b3")

colorvalAPA <- c("#404040", "#f4a582", "#ca0020")
#404040 ## darkgrey - yoked
#f4a582 ## pink - conflict
#ca0020 ## red - consistent 

colorvalAPA2 <-  c( "#404040", "#bababa",  "#ca0020", "#f4a582")

#404040 ## darkgrey - yoked-consistent
#bababa ## light grey - yoked-conflict
#ca0020 ## red - consistent 
#f4a582 ## pink - conflict

colorvalAPA3 <-  c( "#bababa","#404040", "#f4a582", "#ca0020")
#bababa ## light grey - Yoked_Conflict
#404040 ## darkgrey - Yoked_NoConflict
#f4a582 ## pink - Trained_Conflict
#ca0020 ## red - Trained_NoConflict 

colorvalAPA4<-  c( "#bababa","#f4a582")
#bababa ## light grey - yoked-conflict
#f4a582 ## pink - conflict

colorvalAPA5<-  c( "#404040","#ca0020")
#404040 ## darkgrey - yoked-consistent
#ca0020 ## red - consistent 

colorvalAPA6<-  c( "#f4a582","#ca0020")
#f4a582 ## pink - conflict
#ca0020 ## red - consistent 

colorvalAPA7<-  c( "#bababa","#404040")
#bababa ## light grey - yoked-conflict
#404040 ## darkgrey - yoked-consistent


DGConflictControl <- c("#bababa","#f4a582")
DGConsistentControl <- c("#bababa","#ca0020")
DGConflictConsistent <- c("#f4a582","#ca0020")


volcano1 <-  c("consistent" = "#ca0020",
           "yoked_consistent" = "#404040", 
           "none" = "#d9d9d9")

volcano2 <-  c("consistent" = "#ca0020",
                
               "none" = "#d9d9d9")

volcano3 <-  c("yoked_consistent" = "#404040",
               "yoked_conflict" = "#bababa", 
               "none" = "#d9d9d9")

volcano4 <-  c("conflict" = "#f4a582",
               "yoked_conflict" = "#bababa", 
               "none" = "#d9d9d9")





## pheatmap
colorpalette <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )


ann_colors1 = list(
  Group =  c(control = (values=c("#404040")), 
             conflict = (values=c("#f4a582")), 
             consistent = (values=c("#ca0020"))),
  Region = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))

ann_colors2 = list(
  APA =  c(Control = (values=c("#404040")), 
           Conflict = (values=c("#f4a582")), 
             Consistent = (values=c("#ca0020"))),
  Punch = c(CA1 = (values=c("#7570b3")),
             CA3 = (values=c("#1b9e77")), 
             DG = (values=c("#d95f02"))))

ann_colors3 = list(
  APA3 =  c(Yoked_Conflict = (values=c("#bababa")), 
           Yoked_NoConflict = (values=c("#404040")), 
           Trained_Conflict = (values=c("#f4a582")),
           Trained_NoConflict = (values=c("#ca0020"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))

ann_colors4 = list(
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))),
  APA2 =  c('yoked_conflict' = (values=c("#bababa")),
            'yoked_consistent' = (values=c("#404040")), 
            'conflict' = (values=c("#f4a582")),
            'consistent' = (values=c("#ca0020"))))

ann_colors5 = list(
  APA2 =  c('yoked_consistent' = (values=c("#404040")), 
            'consistent' = (values=c("#ca0020"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))

ann_colors6 = list(
  APA2 =  c('conflict' = (values=c("#f4a582")), 
            'consistent' = (values=c("#ca0020"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))

ann_colors7 = list(
  APA2 =  c('yoked_consistent' = (values=c("#404040")), 
            'yoked_conflict' = (values=c("#bababa"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))






APAsession = list(
  APA =  c(control = (values=c("#404040")), 
           conflict = (values=c("#f4a582")),
           consistent = (values=c("#ca0020"))),
  Session =  c(Hab = (values=c("#eff3ff")), 
               T1 = (values=c("#dfc27d")),
               T2 = (values=c("#bf812d")),
               T3 = (values=c("#8c510a")),
               Retest = (values=c("#543005")),
               T4_C1 = (values=c("#cbc9e2")),
               T5_C2 = (values=c("#9e9ac8")),
               T6_C3 = (values=c("#756bb1")),
               Retention = (values=c("#54278f"))))

APAsession2 = list(
  APA =  c('yoked-conflict' = (values=c("#bababa")),
           'yoked-consistent' = (values=c("#404040")), 
           'conflict' = (values=c("#f4a582")),
           'consistent' = (values=c("#ca0020"))),
  Session =  c(Hab = (values=c("#eff3ff")), 
               T1 = (values=c("#dfc27d")),
               T2 = (values=c("#bf812d")),
               T3 = (values=c("#8c510a")),
               Retest = (values=c("#543005")),
               T4_C1 = (values=c("#cbc9e2")),
               T5_C2 = (values=c("#9e9ac8")),
               T6_C3 = (values=c("#756bb1")),
               Retention = (values=c("#54278f"))))


imgDG <- magick::image_read("../figures/00_schematics/DG.png")
imgCA3 <- magick::image_read("../figures/00_schematics/CA3.png")
imgCA1 <- magick::image_read("../figures/00_schematics/CA1.png")