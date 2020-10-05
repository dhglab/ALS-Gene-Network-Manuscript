# Script to run weighted gene co-expression analysis (WGCNA)
# Jerry C. Wang

# Set up R
rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(88)

# Load libraries
library(WGCNA)

# Define directories
dataDir <- "../data/"
metaDataDir <- "../metadata/"
resultsDir <- "../results/"

load(paste0(dataDir,"GTEX.expr4.RData"))
load(paste0(metaDataDir,"GTEX.meta4.RData"))
GTEX.expr.reg <- t(GTEX.expr.reg)

###########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(GTEX.expr.reg, powerVector = powers, verbose = 5)

# Plot the results:
pdf(paste0(resultsDir,"Soft_threshold_selection",".pdf"),height=9,width=5)
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),ylim = c(-1,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.65,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 10;
adjacency = adjacency(GTEX.expr.reg, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

datTraits <- GTEX.meta.total[,c("PC1","PC2","PC3","PC4","PC5","GTEX.meta.SMRIN","GTEX.meta.DTHHRDY"
                                ,"GTEX.meta.AGE","GTEX.meta.SMCENTER.D1.A1","GTEX.meta.ETHNCTY.98","GTEX.meta.TRISCHD")]
traitmat=as.data.frame(cbind(as.numeric(datTraits[,1]),as.numeric(datTraits[,2]),as.numeric(datTraits[,3]),as.numeric(datTraits[,4]),
                             as.numeric(datTraits[,5]),as.factor(datTraits[,6]),as.factor(datTraits[,7]),as.numeric(datTraits[,8]),
                             as.numeric(datTraits[,9]),as.numeric(datTraits[,10]),as.numeric(datTraits[,11])))
# convert categorical variables in factor and numeric as numeric
rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("seqPC1","seqPC2","seqPC3","seqPC4", "seqPC5", "RIN","Hardy","Age","SMCENTER","Ethnicity","TRISCHD")

geneSigs=matrix(NA,nrow=ncol(traitmat),ncol=ncol(GTEX.expr.reg)) # create a vector to hold the data
for(i in 1:ncol(geneSigs)) {
  # calculate r correlation value for numeric variables
  # calculate adjusted R^2s square-root for categorical variables (factors)
  exprvec=as.numeric(GTEX.expr.reg[,i]) # get the expression vector for ith gene
  seqPC1r=bicor(traitmat[,1],exprvec,use="pairwise.complete.obs")
  seqPC2r=bicor(traitmat[,2],exprvec,use="pairwise.complete.obs")
  seqPC3r=bicor(traitmat[,3],exprvec,use="pairwise.complete.obs")
  seqPC4r=bicor(traitmat[,4],exprvec,use="pairwise.complete.obs")
  seqPC5r=bicor(traitmat[,5],exprvec,use="pairwise.complete.obs")
  RINr=bicor(traitmat[,6],exprvec,use="pairwise.complete.obs")
  Hardyr=cor(traitmat[,7],exprvec,method = "spearman",use="pairwise.complete.obs")
  Ager=bicor(traitmat[,8],exprvec,use="pairwise.complete.obs")
  SMCENTERr=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,9])))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(traitmat[,9])))$coefficients[2,1])
  Ethnicityr=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,10])))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(traitmat[,10])))$coefficients[2,1])
  TRISCHDr=bicor(traitmat[,11],exprvec,use="pairwise.complete.obs")
  
  geneSigs[,i]=c(seqPC1r,seqPC2r,seqPC3r,seqPC4r,seqPC5r,RINr,Hardyr,Ager,SMCENTERr,Ethnicityr,TRISCHDr)
  print(paste("Loading... ",i))
}
colnames(geneSigs)=colnames(GTEX.expr.reg)
rownames(geneSigs)=colnames(traitmat)
# convert to colors
geneSigsColor=matrix(NA,nrow=nrow(geneSigs),ncol=ncol(GTEX.expr.reg)) # create a vector to hold the data
for ( i in 1:nrow(geneSigsColor)) {
  geneSigsColor[i,] =numbers2colors(as.numeric(geneSigs[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 
}

rownames(geneSigsColor)=rownames(geneSigs)
colnames(geneSigsColor)=colnames(geneSigs)
# Try out tree cutting parameters
mColorh <- mLabelh <- colorLabels <- NULL  
for (minModSize in c(50,100,150)) {
  for (dthresh in c(0.1,0.2,0.25)) {
    for (ds in c(2,4)) {
      print("Trying parameters:")
      print(c(minModSize,dthresh,ds))
      tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE, minClusterSize = minModSize, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
      merged <- mergeCloseModules(exprData = GTEX.expr.reg,colors = tree$labels,
                                  cutHeight = dthresh)
      mColorh <- cbind(mColorh,labels2colors(merged$colors))
      mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
    }
  }
} 
mColorh1=cbind(mColorh,t(geneSigsColor))
mLabelh1=c(mLabelh,rownames(geneSigsColor))

pdf(paste0(resultsDir,"Signed_Dendro_all_param_SP",softPower,".pdf"),height=25,width=20)
plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()

# make final cut based on parameters chosen
mms=100
ds=2
dthresh=0.1
tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))

merged <- mergeCloseModules(exprData = GTEX.expr.reg,colors = tree$labels, cutHeight = dthresh)
mColorh <- cbind(labels2colors(merged$colors),t(geneSigsColor))
mLabelh <- c("Merged Colors",rownames(geneSigsColor))

pdf(paste0(resultsDir,"Signed_Dendro_selected_param_SP",softPower,".pdf"),height=25,width=20)
plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power=",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh));
dev.off()

mergedColors = labels2colors(merged$colors);

# Eigengenes of the new merged modules:
MEList=moduleEigengenes(GTEX.expr.reg, colors = mergedColors,softPower= softPower, nPC=1)
MEs=MEList$eigengenes
MEs=orderMEs(MEs)
rownames(MEs) = rownames(GTEX.expr.reg)
moduleColors = mergedColors
names(moduleColors) <- colnames(GTEX.expr.reg)
save(geneTree,moduleColors,MEs,file=paste0(resultsDir,"Modules.RData"))


#####################  SESSION INFO ###########################

#R version 3.6.3 (2020-02-29)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Catalina 10.15.6

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#ttached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] WGCNA_1.67            fastcluster_1.1.25    dynamicTreeCut_1.63-1

#loaded via a namespace (and not attached):
#[1] Biobase_2.46.0        bit64_0.9-7           splines_3.6.3        
#[4] foreach_1.5.0         Formula_1.2-3         stats4_3.6.3         
#[7] latticeExtra_0.6-29   blob_1.2.1            fit.models_0.63      
#[10] yaml_2.2.1            robustbase_0.93-6     impute_1.60.0        
#[13] pillar_1.4.4          RSQLite_2.2.0         backports_1.1.8      
#[16] lattice_0.20-41       glue_1.4.1            digest_0.6.25        
#[19] RColorBrewer_1.1-2    checkmate_2.0.0       colorspace_1.4-1     
#[22] htmltools_0.5.0       preprocessCore_1.48.0 Matrix_1.2-18        
#[25] pcaPP_1.9-73          pkgconfig_2.0.3       mvtnorm_1.1-1        
#[28] purrr_0.3.4           GO.db_3.10.0          scales_1.1.1         
#[31] jpeg_0.1-8.1          htmlTable_2.0.1       tibble_3.0.1         
#[34] generics_0.0.2        IRanges_2.20.2        ggplot2_3.3.2        
#[37] ellipsis_0.3.1        nnet_7.3-14           BiocGenerics_0.32.0  
#[40] survival_3.2-3        magrittr_1.5          crayon_1.3.4         
#[43] memoise_1.1.0         doParallel_1.0.15     MASS_7.3-51.6        
#[46] foreign_0.8-76        tools_3.6.3           data.table_1.12.8    
#[49] lifecycle_0.2.0       matrixStats_0.56.0    stringr_1.4.0        
#[52] S4Vectors_0.24.4      munsell_0.5.0         cluster_2.1.0        
#[55] AnnotationDbi_1.48.0  compiler_3.6.3        rlang_0.4.6          
#[58] grid_3.6.3            iterators_1.0.12      rstudioapi_0.11      
#[61] htmlwidgets_1.5.1     robust_0.5-0.0        base64enc_0.1-3      
#[64] gtable_0.3.0          codetools_0.2-16      DBI_1.1.0            
#[67] rrcov_1.5-2           R6_2.4.1              gridExtra_2.3        
#[70] knitr_1.29            dplyr_1.0.0           bit_1.1-15.2         
#[73] Hmisc_4.4-0           stringi_1.4.6         parallel_3.6.3       
#[76] Rcpp_1.0.4            vctrs_0.3.1           rpart_4.1-15         
#[79] acepack_1.4.1         png_0.1-7             DEoptimR_1.0-8       
#[82] tidyselect_1.1.0      xfun_0.15            
