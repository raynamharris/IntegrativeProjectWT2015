# Behavior colors 

# for mean and sd plots
ann_colors_APA2 = list(
  APA2 =  c('yoked-consistent' = (values=c("#404040")), 
            'consistent' = (values=c("#ca0020")),
            'yoked-conflict' = (values=c("#bababa")),
            'conflict' = (values=c("#f4a582"))))

# for beahvir heatmap
APA2session = list(
   APA2 =  c('yoked-consistent' = (values=c("#404040")), 
            'consistent' = (values=c("#ca0020")),
            'yoked-conflict' = (values=c("#bababa")),
            'conflict' = (values=c("#f4a582"))),
   Session =  c(Hab = (values=c("#ffffff")), 
               T1 = (values=c("#f0f0f0")),
               T2 = (values=c("#d9d9d9")),
               T3 = (values=c("#bdbdbd")),
               Retest = (values=c("#969696")),
               T4_C1 = (values=c("#737373")),
               T5_C2 = (values=c("#525252")),
               T6_C3 = (values=c("#252525")),
               Retention = (values=c("#000000"))))

colorvalPunch <- c("#d95f02","#1b9e77", "#7570b3")
## DG  "#d95f02"
## CA3 "#1b9e77"
## CA1 "#7570b3")

colorvalAPA00 <-  c( "#404040","#ca0020", "#bababa", "#f4a582")
#404040 ## darkgrey - Yoked_NoConflict
#ca0020 ## red - Trained_NoConflict 
#bababa ## light grey - Yoked_Conflict
#f4a582 ## pink - Trained_Conflict

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
           "neither" = "#d9d9d9")

volcano2 <-  c("consistent" = "#ca0020",
               "conflict" = "#f4a582",
               "neither" = "#d9d9d9")

volcano3 <-  c("yoked_consistent" = "#404040",
               "yoked_conflict" = "#bababa", 
               "neither" = "#d9d9d9")

volcano4 <-  c("conflict" = "#f4a582",
               "yoked_conflict" = "#bababa", 
               "neither" = "#d9d9d9")

volcano5 <-  c("yes" = "red",
               "no" = "black", 
               "none" = "#d9d9d9")
               
volcanoDGvCA3 <- c("CA3" = "#1b9e77",
                   "DG" = "#d95f02", 
                    "neither" = "#d9d9d9")  
                    
volcanoDGvCA1 <-  c("CA1" = "#7570b3",
                    "DG" = "#d95f02", 
                     "neither" = "#d9d9d9")            

volcanoCA3vCA1 <-  c("CA1" = "#7570b3",
                     "CA3" = "#1b9e77", 
                     "neither" = "#d9d9d9")


colorvalavoidance <- c("black","red")


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

ann_colors8 = list(
  avoidance =  c('yes' = (values=c("red")), 
            'no' = (values=c("black"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))


ann_colorsLevel = list(
  Level = c(
    Behavior = (values=c("#e6ab02")),
    Physiology = (values=c("#a6761d")),
    DG = (values=c("#d95f02")),
    CA3 = (values=c("#1b9e77")), 
    CA1 = (values=c("#7570b3"))))
