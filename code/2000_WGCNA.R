library(tidyverse)
library(corrplot)
library(WGCNA)
library(readxl)
library(data.table)

library(flashClust)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Loading the data: WGCNA requires genes be given in the columns

RawData <- read.csv("Bruneck2000_scaled_data.csv")
gene.names=as.character(RawData$Node)
rownames(RawData) <- gene.names
RawData <- select(RawData,-1,-2)
dim(RawData)
t_Rawdata <- t(RawData)

colnames(t_Rawdata)

colSums(is.na(t_Rawdata))

Remove <- c("VEGF.A.l2",
            "TRANCE.l2",
            "TNFSF14.l2",
            "TF",
            "SIRT2.l2",
            "SCF.l2",
            "OPG.l2",
            "MMP.10.l2",
            "MMP.1.l2",
            "MIP.1.alpha",
            "MCP.1.l2",
            "IL.8.l2",
            "IL.6RA",
            "IL.6.l2",
            "IL.4.l2",
            "IL.18R1",
            "IL.18.l2",
            "HGF.l2",
            "FGF.23.l2",
            "EN.RAGE.l2",
            "CXCL6.l2",
            "CXCL1.l2",
            "CX3CL1.l2",
            "CSF.1.l2",
            "CD40.l2",
            "CCL20.l2",
            "CASP.8.l2",
            "Beta.NGF.l2",
            "U6...24"
)

rows_to_remove <- which(row.names(RawData) %in% Remove)
RawData <- RawData[-rows_to_remove,]
gene.names <- rownames(RawData)

std.gene.names <- c("ORM1",
                    "ORM2",
                    "SERPINA1",
                    "A1BG",
                    "SERPINF2",
                    "LRG1",
                    "SERPINA3",
                    "ABCF1",
                    "AFM",
                    "ALB",
                    "AMBP",
                    "AGT",
                    "SERPINC1",
                    "LPA",
                    "APOA1",
                    "APOA2",
                    "APOA4",
                    "APOB",
                    "APOC1",
                    "APOC2",
                    "APOC3",
                    "APOD",
                    "APOE",
                    "APOH",
                    "APOL1",
                    "APOM",
                    "SLC4A1",
                    "C1QB",
                    "C1QC",
                    "C1R",
                    "C1S",
                    "C4BPA",
                    "SERPINA6",
                    "CD5L",
                    "CP",
                    "CFB",
                    "CFH",
                    "CFI",
                    "CLU",
                    "C2",
                    "C3",
                    "C5",
                    "C6",
                    "C7",
                    "C8A",
                    "PSPHP1",
                    "CPN2",
                    "F13A1",
                    "FBLN1",
                    "FCN3",
                    "AHSG",
                    "FGA",
                    "FGG",
                    "FN1",
                    "GSN",
                    "GPX3",
                    "HBA1",
                    "HBD",
                    "HBE1",
                    "ERVMER34-1",
                    "SERPIND1",
                    "MRS2",
                    "HPR",
                    "HRG",
                    "SERPING1",
                    "IGHA1",
                    "IGHA2",
                    "IGHG1",
                    "IGHG2",
                    "IGHG3",
                    "IGHG4",
                    "IGHM",
                    "JCHAIN",
                    "ITIH1",
                    "ITIH2",
                    "ITIH4",
                    "KLKB1",
                    "KNG1",
                    "MBL2",
                    "SERPINF1",
                    "PGLYRP2",
                    "PF4",
                    "PLG",
                    "RBP4",
                    "SAA4",
                    "SELENOP",
                    "SHBG",
                    "STOM",
                    "CLEC3B",
                    "SERPINA7",
                    "THRB",
                    "Total-IGHG",
                    "TF",
                    "TTR",
                    "GC",
                    "VTN",
                    "AZGP1",
                    "CXCL8",
                    "VEGFA",
                    "ADM",
                    "CD40LG",
                    "GDF15",
                    "PGF",
                    "SELE",
                    "EGF",
                    "TNFRSF11B",
                    "SRC",
                    "IL1R1",
                    "IL6",
                    "CSTB",
                    "CCL2",
                    "KLK6",
                    "IL6R",
                    "F2R",
                    "CMA1",
                    "KLK11",
                    "TEK",
                    "TNFRSF1A",
                    "CD274",
                    "IL27",
                    "CSF1",
                    "CXCL1",
                    "OLR1",
                    "TNFRSF10B",
                    "FGF23",
                    "KITLG",
                    "IL18R1",
                    "TNFRSF1B",
                    "MMP3",
                    "HSPB1",
                    "TNFSF14",
                    "PRL",
                    "MPO",
                    "GH1",
                    "MMP1",
                    "RETN",
                    "FASN",
                    "PAPPA",
                    "PTX3",
                    "REN",
                    "CHI3L1",
                    "IL1RL1",
                    "TIMELESS",
                    "NGF",
                    "XPNPEP2",
                    "TNFSF11",
                    "HGF",
                    "SELPLG",
                    "MB",
                    "THBD",
                    "IL16",
                    "MMP10",
                    "PRAP1",
                    "CCL4",
                    "CTSD",
                    "AGER",
                    "CCL3",
                    "MMP7",
                    "CXCL6",
                    "ITGB1BP2",
                    "CXCL16",
                    "DKK1",
                    "SIRT2",
                    "GAL",
                    "AGRP",
                    "S100A12",
                    "CD40",
                    "PLAT",
                    "HBEGF",
                    "ESM1",
                    "IL4",
                    "VEGFD",
                    "MMP12",
                    "SPON1",
                    "CASP8",
                    "CTSL",
                    "CX3CL1",
                    "FABP4",
                    "NPPB",
                    "LEP",
                    "CCL20",
                    "MUC16",
                    "IKBKG",
                    "FST",
                    "PECAM1",
                    "NT-pro-BNP",
                    "RNASE3",
                    "BDNF",
                    "CCL7",
                    "GDNF",
                    "CDCP1",
                    "CD244",
                    "IL7",
                    "TGFB1",
                    "PLAUR",
                    "IL17C",
                    "IL17A",
                    "CXCL11",
                    "AXIN1",
                    "TNFSF10",
                    "IL20RA",
                    "CXCL9",
                    "CST5",
                    "IL2RB",
                    "IL1A",
                    "OSM",
                    "IL2",
                    "TSLP",
                    "CCL4L1",
                    "CD6",
                    "SLAMF1",
                    "TGFA",
                    "CCL13",
                    "CCL11",
                    "IL10RA",
                    "FGF5",
                    "LIFR",
                    "FGF21",
                    "CCL19",
                    "IL15RA",
                    "IL10RB",
                    "IL22RA1",
                    "PDGFB",
                    "CXCL5",
                    "IL12B",
                    "IL24",
                    "IL13",
                    "ARTN",
                    "IL10",
                    "TNF",
                    "CCL23",
                    "CD5",
                    "FLT3LG",
                    "CXCL10",
                    "EIF4EBP1",
                    "IL20",
                    "CCL28",
                    "DNER",
                    "IL33",
                    "IFNG",
                    "FGF19",
                    "LIF",
                    "NRTN",
                    "CCL8",
                    "CCL25",
                    "TNFRSF9",
                    "NTF3",
                    "TNFSF12",
                    "SULT1A1",
                    "STAMBP",
                    "IL5",
                    "ADA",
                    "LTA",
                    "miR-122",
                    "miR-126",
                    "miR-150",
                    "miR-155",
                    "miR-191",
                    "miR-195",
                    "miR-197",
                    "miR-21",
                    "miR-223",
                    "miR-24",
                    "miR-28-3p",
                    "U6",
                    "miR-27b",
                    "miR-92a",
                    "miR-146a",
                    "miR-148a",
                    "miR-320",
                    "miR-335",
                    "miR-126*",
                    "miR-223*",
                    "miR-378"
)

rownames(RawData) <- std.gene.names
t_Rawdata <- t(RawData)

# Choosing a soft-threshold to fit a scale-free topology to the network

powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t_Rawdata,dataIsExpr = TRUE,powerVector = powers,
                      corFnc = cor,corOptions = list(use = 'p', method = 'spearman'), 
                      networkType = "unsigned")

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 3;

#calculate the adjacency matrix
adj= adjacency(t_Rawdata,type = "unsigned", power = softPower, corFnc = "cor", corOptions ="use = 'p', method = 'spearman'");

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(t_Rawdata, networkType = "unsigned", TOMType = "unsigned", power = softPower);
colnames(TOM) =rownames(TOM) =std.gene.names
dissTOM=1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);
# Set the minimum module size
minModuleSize = 2;
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(t_Rawdata[,restGenes], power = softPower)
colnames(diss1) =rownames(diss1) =std.gene.names[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

# Extract modules

module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=std.gene.names[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

# Look at expression patterns of these genes, as they're clustered

module.order <- unlist(tapply(1:ncol(t_Rawdata),as.factor(dynamicColors),I))
m <- t(t(t_Rawdata[,module.order])/apply(t_Rawdata[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])

# Quantify module similarity by eigengene correlation (Eigengenes: Module representatives)

MEList = moduleEigengenes(t_Rawdata, colors = dynamicColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

# Corrplot messing about

sizeGrWindow(18, 9)
par(mfrow = c(1,1));
cex1 = 0.9;

corrplot(cor(MEs, method = "spearman", use = "pairwise.complete.obs"), 
         method = "circle", order = "hclust", 
         tl.cex = 0.8)

# Determine sign on correlation (positive/negative)

sizeGrWindow(18, 9)
par(mfrow = c(1,1));
cex1 = 0.9;

Corr <- (cor(t_Rawdata, method = "spearman", use = "pairwise.complete.obs"))

corrplot(adj, method = "color", order = "hclust", tl.cex = 0.5)
corrplot(Corr, method = "color", order = "hclust", tl.cex = 0.5)

# RNA/Protein lists

miRNAs <- grep("^miR*|^let*",gene.names, value = TRUE)
proteins <- grep("^miR*|^let*",gene.names, value = TRUE, invert = TRUE)

miR_F <- grep("^miR*|^let*",std.gene.names, value = TRUE)
proteins_F <- grep("^miR*|^let*",std.gene.names, value = TRUE, invert = TRUE)

write(as.character(miR_F), file = "miRNAs.txt", sep = "\n")
write(as.character(proteins_F), file = "proteins.txt", sep = "\n")

colnames(RawData)

# Processing for Modulation by Target Prediction

TarPred <- read_csv("miRNA data/Consolidated_Data_FLAG.csv")
TarPred <- TarPred[,-1]
TarPred$miRNA <- gsub('hsa-', '', TarPred$miRNA)

miRNA_list.A <- c('miR-122-5p',
                  'miR-126-5p',
                  'miR-150-5p',
                  'miR-155-5p',
                  'miR-191-5p',
                  'miR-195-5p',
                  'miR-197-5p',
                  'miR-21-5p',
                  'miR-223-5p',
                  'miR-24-5p',
                  'miR-28-3p',
                  'U6',
                  'miR-27b-5p',
                  'miR-92a-5p',
                  'miR-146a-5p',
                  'miR-148a-5p',
                  'miR-320-5p',
                  'miR-335-5p',
                  'miR-378-5p'
)

miRNA_list.B <- c('miR-122-3p',
                  'miR-126-3p',
                  'miR-150-3p',
                  'miR-155-3p',
                  'miR-191-3p',
                  'miR-195-3p',
                  'miR-197-3p',
                  'miR-21-3p',
                  'miR-223-3p',
                  'miR-24-3p',
                  'miR-27b-3p',
                  'miR-92a-3p',
                  'miR-146a-3p',
                  'miR-148a-3p',
                  'miR-320-3p',
                  'miR-335-3p',
                  'miR-378-3p'
)

Unweight.data <- cbind(t_Rawdata, t_Rawdata[,254:274])
Unweight.data <- Unweight.data[,c(-271, -272, -285, -286, -293, -294)]

View(head(Unweight.data))

colnames(Unweight.data)[254:272] <- miRNA_list.A
colnames(Unweight.data)[273:289] <- miRNA_list.B

adj.unweight <- adjacency(Unweight.data,type = "unsigned", power = softPower, corFnc = "cor", corOptions ="use = 'p', method = 'spearman'")
Corr.unweight <- cor(Unweight.data, method = "spearman", use = "pairwise.complete.obs")

adj.weight <- data.frame(matrix(nrow=289,ncol=289))
row.names(adj.weight) = row.names(adj.unweight)
colnames(adj.weight) = colnames(adj.unweight)

as.integer(TarPred %>% 
             filter(miRNA=='miR-122-5p',Gene=='PGF') %>% 
             select('No. Sources'))

row = 'miR-122-5p'
col = 'PGF'
as.integer(TarPred %>% 
             filter(miRNA==row,Gene==col) %>% 
             select('No. Sources'))

for(i in 1:dim(adj.unweight)[1]) {
  for(j in 1:dim(adj.unweight)[2]){
    row = rownames(adj.weight)[i]
    col = colnames(adj.weight)[j]
    count = as.integer(TarPred %>% filter(miRNA==row,Gene==col) %>% select('No. Sources'))
    count2 = as.integer(TarPred %>% filter(miRNA==col,Gene==row) %>% select('No. Sources'))
    if(grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      adj.weight[i,j] = adj.unweight[i,j]
    } else if(!grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      adj.weight[i,j] = adj.unweight[i,j]
    } else if(grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      if(is.na(count) || count == '') {
        adj.weight[i,j] = 0
      } else if(count >= 2) {
        adj.weight[i,j] = adj.unweight[i,j]*1
      } else if(count == 1) {
        adj.weight[i,j] = adj.unweight[i,j]*0.75
      } else {
        adj.weight[i,j] = 0
      }
    } else if(!grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      if(is.na(count2) || count2 == '') {
        adj.weight[i,j] = 0
      } else if(count2 >= 2) {
        adj.weight[i,j] = adj.unweight[i,j]*1
      } else if(count2 == 1) {
        adj.weight[i,j] = adj.unweight[i,j]*0.75
      } else {
        adj.weight[i,j] = 0
      }
    }
  }
}

rows_to_remove.adj <- which(row.names(adj.weight) %in% proteins_F)
cols_to_remove.adj <- which(colnames(adj.weight) %in% miR_F)

adj.weight.F <- adj.weight[-rows_to_remove.adj,-cols_to_remove.adj]

write.csv(adj.weight.F,"Weighted_Adjacencies.csv")

Corr.weight <- data.frame(matrix(nrow=289,ncol=289))
row.names(Corr.weight) = row.names(Corr.unweight)
colnames(Corr.weight) = colnames(Corr.unweight)

for(i in 1:dim(Corr.unweight)[1]) {
  for(j in 1:dim(Corr.unweight)[2]){
    row = rownames(Corr.weight)[i]
    col = colnames(Corr.weight)[j]
    count = as.integer(TarPred %>% filter(miRNA==row,Gene==col) %>% select('No. Sources'))
    count2 = as.integer(TarPred %>% filter(miRNA==col,Gene==row) %>% select('No. Sources'))
    if(grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      Corr.weight[i,j] = Corr.unweight[i,j]
    } else if(!grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      Corr.weight[i,j] = Corr.unweight[i,j]
    } else if(grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      if(is.na(count) || count == '') {
        Corr.weight[i,j] = 0
      } else if(count >= 2) {
        Corr.weight[i,j] = Corr.unweight[i,j]*1
      } else if(count == 1) {
        Corr.weight[i,j] = Corr.unweight[i,j]*0.75
      } else {
        Corr.weight[i,j] = 0
      }
    } else if(!grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      if(is.na(count2) || count2 == '') {
        Corr.weight[i,j] = 0
      } else if(count2 >= 2) {
        Corr.weight[i,j] = Corr.unweight[i,j]*1
      } else if(count2 == 1) {
        Corr.weight[i,j] = Corr.unweight[i,j]*0.75
      } else {
        Corr.weight[i,j] = 0
      }
    }
  }
}

rows_to_remove.corr <- which(row.names(Corr.weight) %in% proteins_F)
cols_to_remove.corr <- which(colnames(Corr.weight) %in% miR_F)

Corr.weight.F <- Corr.weight[-rows_to_remove.adj,-cols_to_remove.adj]

write.csv(Corr.weight.F,"Weighted_Spearman.csv")

# Plot weighted matrices

corrplot(adj, method = "color", order = "hclust", tl.cex = 0.5)
corrplot(Corr, method = "color", order = "hclust", tl.cex = 0.5)

adj.weight2 <- as.matrix(adj.weight)
Corr.weight2 <- as.matrix(Corr.weight)

corrplot(adj.weight2, method = "color", order = "hclust", tl.cex = 0.5)
corrplot(Corr.weight2, method = "color", order = "hclust", tl.cex = 0.5)

# Writing files for ARACNe-AP
# Protein vs. Protein data
# miRNA vs. miRNA/Protein data

write(as.character(std.gene.names), file = "GeneNames.txt", sep = "\n")
RawData2 <- RawData %>% rownames_to_column(' ') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ') %>% 
  t()
RawData2[1,1] <- ''
write.table(RawData2, file = "GeneData.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

rows_to_remove2 <- which(row.names(RawData) %in% miR_F)
Protein_Data <- RawData[-rows_to_remove2,]
Protein_Data <- Protein_Data %>% rownames_to_column(' ') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ') %>% 
  t()
Protein_Data[1,1] <- ''
write.table(Protein_Data, file = "Protein_Data.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

# Spectral Plot pre-processing

write_csv(as.data.frame(Corr),"Corr_Table.csv")
SpecData_PreProcessing <- read_excel("SpecData.xlsx", 
                                     sheet = "miR-Pro Neg=Abs")
SpecData <- as.data.frame(SpecData_PreProcessing)
row.names(SpecData) <- SpecData[,1]
SpecData <- SpecData[,-1]

# Plotting weighted matrices

install.packages("plot.matrix")
library(plot.matrix)

png(height=1200, width=5000, pointsize=1, file="WeightedCorr2000.png")
plot(Corr.weight.F)
dev.off()

Corr.WGCNA <- c('PF4',
                'STOM',
                'CD40LG',
                'EGF',
                'SRC',
                'F2R',
                'CD274',
                'CXCL1',
                'OLR1',
                'IL18R1',
                'HSPB1',
                'TNFSF14',
                'MPO',
                'MMP1',
                'CXCL6',
                'ITGB1BP2',
                'DKK1',
                'SIRT2',
                'S100A12',
                'CD40',
                'HBEGF',
                'CASP8',
                'IKBKG',
                'PECAM1',
                'RNASE3',
                'CD244',
                'IL7',
                'TGFB1',
                'AXIN1',
                'OSM',
                'CCL13',
                'CXCL5',
                'EIF4EBP1',
                'CCL8',
                'SULT1A1',
                'STAMBP')

rows_to_keep.CorF <- which(colnames(Corr.weight.F) %in% Corr.WGCNA)

Corr.weight.Mat <- Corr.weight.F[,rows_to_keep.CorF]

write.csv(Corr.weight.Mat,"Weighted_Spearman_WGCNA.csv")

data <- read.csv("Weighted_Spearman_WGCNA.csv", row.names=1,header=TRUE)
par(cex=0.8)
corrplot(as.matrix(data), tl.col="black",tl.srt=45,insig = "blank", order="hclust")
