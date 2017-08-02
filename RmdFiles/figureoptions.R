colorvalAPA <- c("#404040", "#f4a582", "#ca0020")
#404040 ## darkgrey - yoked
#ca0020 ## red - conflict 
#f4a582 ## pink - same

colorvalPunch <- c("#d95f02","#1b9e77", "#7570b3")

## pheatmap
colorpalette <-  colorRampPalette(c("Deep Sky Blue 3", "white", "red"))( 30 )

ann_colors1 = list(
  Group =  c(control = (values=c("#404040")), 
             consistent = (values=c("#f4a582")), 
             conflict = (values=c("#ca0020"))),
  Region = c(CA1 = (values=c("#7570b3")),
            CA3 = (values=c("#1b9e77")), 
            DG = (values=c("#d95f02"))))


session_colors1 = list(
  APA =  c(Control = (values=c("#404040")), 
           Consistent = (values=c("#f4a582")),
           Conflict = (values=c("#ca0020"))),
  Session =  c(Hab = (values=c("#eff3ff")), 
               T1 = (values=c("#bdd7e7")),
               T2 = (values=c("#6baed6")),
               T3 = (values=c("#3182bd")),
               Retest = (values=c("#08519c")),
               T4_C1 = (values=c("#cbc9e2")),
               T5_C2 = (values=c("#9e9ac8")),
               T6_C3 = (values=c("#756bb1")),
               Retention = (values=c("#54278f"))))