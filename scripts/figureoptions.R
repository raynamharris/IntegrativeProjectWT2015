# this file stores color palettes and themes

theme_ms <- function () { 
  theme_classic(base_size = 8) +
    theme(
      panel.grid.major  = element_blank(),  # remove major gridlines
      panel.grid.minor  = element_blank(),  # remove minor gridlines
      plot.title = element_text(hjust = 0.5, face = "bold") # center & bold 
    )
}



# for mean and sd plots of behavior 
ann_colors_treatment = list(
  treatment =  c('standard yoked' = "#404040", 
                 'standard trained' = "#ca0020",
                 'conflict yoked' = "#969696",
                 'conflict trained' = "#f4a582"),
  Session = c(  "T1"      = "#fff7f3",
                "T2"      = "#fde0dd",
                "T3"      = "#fcc5c0",
                "Retest"  = "#fa9fb5",
                "T4_C1"   = "#f768a1",
                "T5_C2"   = "#dd3497",
                "T6_C3"   = "#ae017e",
                "Retention" = "#7a0177"   )
)

levelstreatment =  c('standard.yoked' , 'standard.trained' ,
                     'conflict.yoked' ,  'conflict.trained' )

levelsTrials <- c("Hab", "T1" , "T2", "T3" ,
                  "Retest" , "T4_C1" , "T5_C2" ,
                  "T6_C3" , "Retention" )

levelssubfield <- c("DG", "CA3", "CA1")

# for beahvior and physiology

colorvalAPA00 <-  c( "#404040","#ca0020", "#969696", "#f4a582")
#404040 ## darkgrey - Yoked_NoConflict
#ca0020 ## red - Trained_NoConflict 
#969696 ## light grey - Yoked_Conflict
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
  treatment = c( 'standard yoked' = (values=c("#404040")), 
                 'standard trained' = (values=c("#ca0020")),
                 'conflict yoked' = (values=c("#969696")),
                 'conflict trained' = (values=c("#f4a582"))))

pheatmapcolors2 = list(
  combinedgroups = c(yoked = "black",
                     trained = "red"),
  subfield = c(DG = "#d95f02"))

pheatmapcolors3 = list(
  combinedgroups = c(yoked = "black",
                     trained = "red"),
  subfield = c(CA1 = "#1b9e77"))

# for subfield volcano plots

trainedcolors <-  c( "standard.yoked" = "#404040", "standard.trained" = "#ca0020")
conflictcolors <- c("conflict.yoked" = "#969696","conflict.trained" = "#f4a582")
treatmentcolors <- c( "standard.yoked" = "#404040", "standard.trained" = "#ca0020",
                      "conflict.yoked" = "#969696","conflict.trained" = "#f4a582")
treatmentcolors2 <- c( "standard\nyoked" = "#404040", "standard\ntrained" = "#ca0020",
                      "conflict\nyoked" = "#969696","conflict\ntrained" = "#f4a582")



volcano1 <-  c("standard.trained" = "#ca0020",
               "standard.yoked" = "#404040", 
               "NS" = "#d9d9d9")

volcano2 <-  c("standard.trained" = "#ca0020",
               "conflict.trained" = "#f4a582", 
               "NS" = "#d9d9d9")

volcano3 <-  c("standard.yoked" = "#404040",
               "home.cage" = "#e0e0e0", 
               "NS" = "#d9d9d9")

volcano4 <-  c("conflict.yoked" = "#969696",
               "home.cage" = "#e0e0e0", 
               "NS" = "#d9d9d9")

volcano5 <-  c("conflict.trained" = "#f4a582",
               "conflict.yoked" = "#969696", 
               "NS" = "#d9d9d9")

volcano6 <-  c("standard.yoked" = "#404040",
               "conflict.yoked" = "#969696", 
               "NS" = "#d9d9d9")

volcano7 <-  c("trained" = "darkred",
               "yoked" = "black", 
               "NS" = "#d9d9d9")




# subfield box plots

fivegroups <- c( "#e0e0e0", "#404040",  "#ca0020", "#969696", "#f4a582")
fourgroups <- c(           "#404040",  "#ca0020", "#969696", "#f4a582")


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
