# this file stores color palettes

# for mean and sd plots of behavior 
ann_colors_APA2 = list(
  APA2 =  c('yoked-consistent' = (values=c("#404040")), 
            'consistent' = (values=c("#ca0020")),
            'yoked-conflict' = (values=c("#bababa")),
            'conflict' = (values=c("#f4a582"))))


# for beahvior and physiology

colorvalAPA00 <-  c( "#404040","#ca0020", "#bababa", "#f4a582")
#404040 ## darkgrey - Yoked_NoConflict
#ca0020 ## red - Trained_NoConflict 
#bababa ## light grey - Yoked_Conflict
#f4a582 ## pink - Trained_Conflict


# for all data included RNA-seq plots

colorvalPunch <- c("#d95f02","#1b9e77", "#7570b3")
## DG  "#d95f02"
## CA3 "#1b9e77"
## CA1 "#7570b3")

colorvalsubfield <- c("#d95f02","#1b9e77", "#7570b3")

volcanoDGvCA3 <- c("CA3" = "#1b9e77",
                   "DG" = "#d95f02", 
                   "neither" = "#d9d9d9")  

volcanoDGvCA1 <-  c("CA1" = "#7570b3",
                    "DG" = "#d95f02", 
                    "neither" = "#d9d9d9")            

volcanoCA3vCA1 <-  c("CA1" = "#7570b3",
                     "CA3" = "#1b9e77", 
                     "neither" = "#d9d9d9")


pheatmapcolors = list(
  subfield = c(DG = (values=c("#d95f02")),
              CA3 = (values=c("#1b9e77")),
              CA1 = (values=c("#7570b3"))),
  treatment = c('yoked consistent' = (values=c("#404040")), 
            'consistent' = (values=c("#ca0020")),
              'yoked conflict' = (values=c("#bababa")),
            'conflict' = (values=c("#f4a582"))))

pheatmapcolors2 = list(
  Punch = c(DG = (values=c("#d95f02")),
               CA3 = (values=c("#1b9e77")),
               CA1 = (values=c("#7570b3"))),
  Treatment = c('yoked consistent' = (values=c("#404040")), 
                'consistent' = (values=c("#ca0020")),
                'yoked conflict' = (values=c("#bababa")),
                'conflict' = (values=c("#f4a582"))))



# for subfield volcano plots

volcano1 <-  c("consistent" = "#ca0020",
           "yoked_consistent" = "#404040", 
           "NS" = "#d9d9d9")

volcano2 <-  c("consistent" = "#ca0020",
               "conflict" = "#f4a582", 
               "NS" = "#d9d9d9")

volcano3 <-  c("yoked_consistent" = "#404040",
               "yoked_conflict" = "#bababa", 
               "NS" = "#d9d9d9")

# for avoidance volcano plots and heatmaps

volcano5 <-  c("yes" = "red",
               "no" = "black", 
               "none" = "#d9d9d9")
               
colorvalavoidance <- c("black","red")

ann_colors8 = list(
  avoidance =  c('yes' = (values=c("red")), 
                 'no' = (values=c("black"))),
  Punch = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))


# for integrative correlation anaylsis

colorvalAPA6<-  c( "#f4a582","#ca0020")
#f4a582 ## pink - conflict
#ca0020 ## red - consistent 

ann_colorsLevel = list(
  Level = c(
    Behavior = (values=c("#e6ab02")),
    Physiology = (values=c("#a6761d")),
    DG = (values=c("#d95f02")),
    CA3 = (values=c("#1b9e77")), 
    CA1 = (values=c("#7570b3"))))
