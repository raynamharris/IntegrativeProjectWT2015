# this file stores color palettes and themes

# the for all figures going in manuscript
theme_ms <- function () { 
  theme_classic(base_size = 8) +
    theme(
      panel.grid.major  = element_blank(),  # remove major gridlines
      panel.grid.minor  = element_blank(),  # remove minor gridlines
      plot.title = element_text(hjust = 0.5, face = "bold") # center & bold 
    )
}

# facor levels
levelstreatment =  c('standard.yoked' , 'standard.trained' ,'conflict.yoked' ,  'conflict.trained' )
levelstraining = c("yoked", "trained")
levelssubfield = c("DG", "CA3", "CA1")


levelstreatmentlegend = c('standard yoked' , 'standard trained' ,'conflict yoked' ,  'conflict trained' )


treatment_names <- list(
  'standard.yoked' = "standard yoked",
  'standard.trained' = "standard trained",
  'conflict.yoked' = "conflict yoked",
  'conflict.trained' = "conflict trained"
)

treatment_labeller <- function(variable,value){
    return(treatment_names[value])
  }



# colors schemes

treatmentcolors <- c( "standard.yoked" = "#404040", "standard.trained" = "#ca0020",
                      "conflict.yoked" = "#969696","conflict.trained" = "#f4a582")

colorvalsubfield <- c("DG" = "#d95f02", "CA3" = "#1b9e77", "CA1" = "#7570b3")

volcano1 <-  c("trained" = "darkred", "yoked" = "black",  "NS" = "#d9d9d9")
volcano2 <-  c("standard.trained" = "#ca0020", "conflict.trained" = "#f4a582", "NS" = "#d9d9d9")
volcano3 <-  c("standard.yoked" = "#404040", "conflict.yoked" = "#969696",  "NS" = "#d9d9d9")


allcolors <- c(treatmentcolors, colorvalsubfield, volcano1)
