# WGCNA ----
# Install and Load WGCNA package
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#install.packages("flashClust")
library(WGCNA)
library(flashClust)
library(magrittr) ## to use the weird pipe
library(dplyr)
library(DESeq2) ## for gene expression analysis


options(stringsAsFactors=FALSE)
enableWGCNAThreads()

########################################################    
#        Prep Variance Stabilized Data
########################################################    

countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
colData <- read.csv("../data/02a_colData.csv", header = T)
colData <- colData[c(1,2,4,8)]

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Punch + APA2 + Punch*APA2)

dds # view the DESeq object - note numnber of genes
dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
dds # view number of genes afternormalization and the number of samples
dds <- DESeq(dds) # Differential expression analysis

vsd=getVarianceStabilizedData(dds) 
datExpr0 <- vsd
datExpr0 <- datExpr0[rowMeans(datExpr0[, -1])>1, ] ## remove rows with rowsum > some value
datExpr0 <- t(datExpr0) ## transpose data
datExpr0 <- as.data.frame(datExpr0)



# test that all samples good to go
gsg=goodSamplesGenes(datExpr0, verbose = 1)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts
#head(gsg)

#-----Make a trait data frame from just sample info without beahvior
datTraits <- colData

datTraits$Mouse <- as.integer(factor(datTraits$Mouse))
datTraits$RNAseqID <- as.integer(factor(datTraits$RNAseqID))
datTraits$APA2 <- as.integer(factor(datTraits$APA2))
datTraits$Punch <- as.integer(factor(datTraits$Punch))

datTraits$Mouse <- as.numeric(factor(datTraits$Mouse))
datTraits$RNAseqID <- as.numeric(factor(datTraits$RNAseqID))
datTraits$APA2 <- as.numeric(factor(datTraits$APA2))
datTraits$Punch <- as.numeric(factor(datTraits$Punch))

str(datTraits)


#######   #################    ################   #######    
#                 Call sample outliers
#######   #################    ################   #######   

#-----Sample dendrogram and traits
A=adjacency(t(datExpr0),type="signed")
#-----Calculate whole network connectivity
k=as.numeric(apply(A,2,sum))-1

#######   #################    ################   ####### 
#      Standardized connectivity
#######   #################    ################   ####### 

Z.k=scale(k)
thresholdZ.k=-2.5 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
#-----Convert traits to colors
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
str(traitColors)
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlier=outlierColor,traitColors)

#######   #################    ################   ####### 
#      Plot the sample dendrogram
#######   #################    ################   ####### 

pdf(file="../figures/02f_RNAseq_WGCNA/SampleDendro.pdf", width=6, height=6)
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()


# Plot a line to show the cut

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.45, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#-----Remove outlying samples 
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr0=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
A=adjacency(t(datExpr0),type="distance")
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)
dim(datExpr0)
dim(datTraits)


#######   #################    ################   #######    
#                     Choose soft threshold
#######   #################    ################   #######     

powers= c(seq(1,10,by=1), seq(from =12, to=20, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose =5,networkType="signed") #call network topology analysis function

pdf(file="../figures/02f_RNAseq_WGCNA/softthresholdo.pdf", width=6, height=5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

softPower=16



#######   #################    ################   #######    
#                    Construct network
#######   #################    ################   #######     


adjacency=adjacency(datExpr0, power=softPower, type="signed" )  #add  if adjusting 
TOM= TOMsimilarity(adjacency, TOMType="signed")
dissTOM= 1-TOM
geneTree= flashClust(as.dist(dissTOM), method="average")

pdf(file="../figures/02f_RNAseq_WGCNA/geneTree1.pdf", width=6, height=5)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()

#######   #################    ################   #######    
#                    Make modules
#######   #################    ################   ####### 

minModuleSize=100
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)
dynamicColors= labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#-----Merge modules whose expression profiles are very similar
MEList= moduleEigengenes(datExpr0, colors= dynamicColors)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs, use = 'pairwise.complete.obs')
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.0006
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors= merge$colors
mergedMEs= merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#######   #################    ################   #######    
#                Relate modules to traits
#######   #################    ################   ####### 

datt=datExpr0

#-----Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
#-----Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#-----Correlations of genes with eigengenes
moduleGeneCor=cor(MEs, datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, datTraits);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#---------------------Module-trait heatmap

#quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
dev.off()
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
######--------------------end--------------------#######



#---------------------Eigengene heatmap
which.module="blue" #replace with module of interest
datME=MEs
datExpr=datt
#quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=(datTraits$genoAPA), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample")

## saved as 4-<color>-alldata 



# Define variable weight containing the weight column of datTrait
newdatTraits <- datTraits 
newdatTraits$names <- rownames(datTraits)
newdatTraits$names <- NULL

punch = as.data.frame(newdatTraits$Punch);
names(punch) = "punch"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, punch, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(punch), sep="");
names(GSPvalue) = paste("p.GS.", names(punch), sep="");



#### output "gene" list ----
## see https://github.com/ClaireMGreen/TDP-43_Code/blob/afb43cddb8ec1a940fbcfa106a1cc3cf77568b7e/WGCNA2.R

#find out what the IDs are of the genes that are contained within a module. 

blue <- as.data.frame(colnames(datExpr0)[moduleColors=='blue'])
blue$module <- "blue"
colnames(blue)[1] <- "trainscript_length"
red <- as.data.frame(colnames(datExpr0)[moduleColors=='red'])
red$module <- "red"
colnames(red)[1] <- "trainscript_length"
green <- as.data.frame(colnames(datExpr0)[moduleColors=='green'])
green$module <- "green"
colnames(green)[1] <- "trainscript_length"
yellow <- as.data.frame(colnames(datExpr0)[moduleColors=='yellow'])
yellow$module <- "yellow"
colnames(yellow)[1] <- "trainscript_length"
brown <- as.data.frame(colnames(datExpr0)[moduleColors=='brown'])
brown$module <- "brown"
colnames(brown)[1] <- "trainscript_length"
turquoise <- as.data.frame(colnames(datExpr0)[moduleColors=='turquoise'])
turquoise$module <- "turquoise"
colnames(turquoise)[1] <- "trainscript_length"
black <- as.data.frame(colnames(datExpr0)[moduleColors=='black'])
black$module <- "black"
colnames(black)[1] <- "trainscript_length"
magenta <- as.data.frame(colnames(datExpr0)[moduleColors=='magenta'])
magenta$module <- "magenta"
colnames(magenta)[1] <- "trainscript_length"
pink <- as.data.frame(colnames(datExpr0)[moduleColors=='pink'])
pink$module <- "pink"
colnames(pink)[1] <- "trainscript_length"
purple <- as.data.frame(colnames(datExpr0)[moduleColors=='purple'])
purple$module <- "purple"
colnames(purple)[1] <- "trainscript_length"
