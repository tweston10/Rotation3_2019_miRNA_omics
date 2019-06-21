install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
source("http://bioconductor.org/biocLite.R") 
biocLite("impute")

library(tidyverse)
library(corrplot)
library(WGCNA)
library(readxl)
library(data.table)

library(flashClust)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Loading the data: WGCNA requires genes be given in the columns

RawData <- read.csv("Bruneck2015_scaled_data.csv")
gene.names=as.list(RawData$Node)
rownames(RawData) <- gene.names
RawData <- select(RawData,-1,-2)
dim(RawData)
t_Rawdata <- t(RawData)

# Removing duplicate CVD/INF panel proteins

unique(gene.names)
y <- list()
y <- gsub("[0-9]+_(.*)_\\(...\\)","\\1",gene.names)
unique(y)
z <- as.data.frame(table(y))
Dups <- z[!z$Freq<2,]

NAs <- as.data.frame(colSums(is.na(t_Rawdata))) %>% rownames_to_column('gene')
grep(paste(as.character(Dups$y),collapse = "|"),NAs$gene,value = TRUE)

Dup_all <- paste(as.vector(Dups$y))
Dup_CVD <- paste(as.vector(Dups$y),"(CVD)", sep = "_")
Dup_INF <- paste(as.vector(Dups$y),"(INF)", sep = "_")

NAs %>% filter(grepl(paste(Dup_all, collapse="|"), gene))

Remove <- c("101_IL-8_(CVD)",
"102_VEGF-A_(CVD)",
"110_OPG_(CVD)",
"113_IL-6_(INF)",
"115_MCP-1_(CVD)",
"120_TRAIL_(CVD)",
"127_CSF-1_(CVD)",
"128_CXCL1_(CVD)",
"163_CCL4_(CVD)",
"131_FGF-23_(CVD)",
"132_SCF_(CVD)",
"133_IL-18_(CVD)",
"138_TNFSF14_(CVD)",
"142_MMP-1_(CVD)",
"153_Beta-NGF_(CVD)",
"155_TRANCE_(CVD)",
"156_HGF_(CVD)",
"161_MMP-10_(CVD)",
"166_CCL3_(CVD)",
"168_CXCL6_(CVD)",
"172_SIRT2_(INF)",
"175_EN-RAGE_(CVD)",
"176_CD40_(CVD)",
"184_CASP-8_(CVD)",
"186_CX3CL1_(CVD)",
"190_CCL20_(CVD)")

rows_to_remove <- which(row.names(RawData) %in% Remove)
RawData <- RawData[-rows_to_remove,]
gene.names <- rownames(RawData)

std.gene.names <- c("ADM",
                    "CD40LG",
                    "GDF15",
                    "PGF",
                    "SELE",
                    "EGF",
                    "SRC",
                    "IL1RN",
                    "IL6",
                    "CSTB",
                    "KLK6",
                    "LGALS3",
                    "F2R",
                    "KLK11",
                    "TEK",
                    "F3",
                    "TNFRSF1A",
                    "PDGFB",
                    "IL27",
                    "OLR1",
                    "TNFRSF10B",
                    "IL6R",
                    "TNFRSF1B",
                    "MMP3",
                    "HSPB1",
                    "PRL",
                    "MPO",
                    "GH1",
                    "RETN",
                    "FAS",
                    "PAPPA",
                    "PTX3",
                    "REN",
                    "CHI3L1",
                    "IL1RL1",
                    "HAVCR1",
                    "XPNPEP2",
                    "SELPLG",
                    "MB",
                    "THBD",
                    "IL16",
                    "PLAUR",
                    "CTSD",
                    "AGER",
                    "MMP7",
                    "CXCL16",
                    "DKK1",
                    "SIRT2",
                    "GAL",
                    "AGRP",
                    "PLAT",
                    "HBEGF",
                    "ESM1",
                    "VEGFD",
                    "MMP12",
                    "SPON1",
                    "CTSL",
                    "FABP4",
                    "NPPB",
                    "LEP",
                    "MUC16",
                    "IKBKG",
                    "FST",
                    "PECAM1",
                    "NT-pro-BNP",
                    "RNASE3",
                    "CXCL8",
                    "VEGFA",
                    "BDNF",
                    "CCL7",
                    "GDNF",
                    "CDCP1",
                    "CD244",
                    "IL7",
                    "TNFRSF11B",
                    "TGFB1",
                    "PLAU",
                    "IL17C",
                    "CCL2",
                    "IL17A",
                    "CXCL11",
                    "AXIN1",
                    "TNFSF10",
                    "CXCL9",
                    "CST5",
                    "OSM",
                    "CXCL1",
                    "CCL4",
                    "CD6",
                    "KITLG",
                    "IL18",
                    "SLAMF1",
                    "TGFA",
                    "CCL13",
                    "CCL11",
                    "TNFSF14",
                    "FGF23",
                    "IL10RA",
                    "FGF5",
                    "MMP1",
                    "LIFR",
                    "FGF21",
                    "CCL19",
                    "IL15RA",
                    "IL10RB",
                    "IL18R1",
                    "CD274",
                    "NGF",
                    "CXCL5",
                    "TNFSF11",
                    "HGF",
                    "IL12B",
                    "MMP10",
                    "IL10",
                    "CCL23",
                    "CD5",
                    "CCL3",
                    "FLT3LG",
                    "CXCL6",
                    "CXCL10",
                    "EIF4EBP1",
                    "CCL28",
                    "DNER",
                    "S100A12",
                    "CD40",
                    "FGF19",
                    "CCL8",
                    "CASP8",
                    "CCL25",
                    "CX3CL1",
                    "TNFRSF9",
                    "NTF3",
                    "TNFSF12",
                    "CCL20",
                    "SULT1A1",
                    "STAMBP",
                    "ADA",
                    "LTA",
                    "CSF1",
                    "miR-16-5p",
                    "miR-19b-3p",
                    "miR-20b-5p",
                    "miR-21-5p",
                    "miR-24-3p",
                    "miR-26a-5p",
                    "miR-27b-3p",
                    "miR-28-5p",
                    "miR-29a-3p",
                    "miR-30b-5p",
                    "miR-92a-3p",
                    "miR-93-5p",
                    "miR-99b-5p",
                    "miR-122-5p",
                    "miR-126-3p",
                    "miR-126-5p",
                    "miR-130a-3p",
                    "miR-140-5p",
                    "miR-143-3p",
                    "miR-145-5p",
                    "miR-150-5p",
                    "miR-191-5p",
                    "miR-192-5p",
                    "miR-197-3p",
                    "miR-221-3p",
                    "miR-222-3p",
                    "miR-223-3p",
                    "miR-223-5p",
                    "miR-320a",
                    "miR-324-5p",
                    "miR-375",
                    "miR-486-5p",
                    "miR-574-3p",
                    "miR-885-5p",
                    "let-7b-5p",
                    "let-7d-5p",
                    "RNY4-3p",
                    "RNY4-5p",
                    "RNY5-5p"
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

softPower = 5;

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

corrplot(adj, diag = FALSE, method = "circle", order = "hclust", tl.cex = 0.5)

corrplot(cor(MEs, method = "spearman", use = "pairwise.complete.obs"), 
         method = "circle", order = "hclust", 
         tl.cex = 0.8)

# Writing files for ARACNe-AP

write(as.character(std.gene.names), file = "GeneNames.txt", sep = "\n")
RawData2 <- RawData %>% rownames_to_column(' ')
write.table(RawData2, file = "GeneData.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)

miRNAs <- grep("[0-9]+_(.*)_\\(...\\)|RNY*",gene.names, value = TRUE, invert = TRUE)
proteins <- grep("[0-9]+_(.*)_\\(...\\)|RNY*",gene.names, value = TRUE)

proteins2 <- sub("[0-9]+_(.*)_\\(...\\)", '\\1', proteins)

miR.F <- grep("miR*|let*",std.gene.names, value = TRUE)
prot.F <- grep("miR*|let*",std.gene.names, value = TRUE, invert = TRUE)

write(as.character(miR.F), file = "miRNAs.txt", sep = "\n")
write(as.character(prot.F), file = "proteins.txt", sep = "\n")

rows_to_remove2 <- which(row.names(RawData) %in% miR.F)
Protein_Data <- RawData[-rows_to_remove2,]
Protein_Data <- Protein_Data %>% rownames_to_column(' ')

write.table(Protein_Data, file = "Protein_Data.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)

# Protein vs. Protein data
# miRNA vs. miRNA/Protein data

# Determine sign on correlation (positive/negative)

sizeGrWindow(18, 9)
par(mfrow = c(1,1));
cex1 = 0.9;

Corr <- (cor(t_Rawdata, method = "spearman", use = "pairwise.complete.obs"))
corrplot(Corr, method = "circle", order = "hclust", tl.cex = 0.5)

# Spectral Plot pre-processing

write_csv(as.data.frame(Corr),"Corr_Table.csv")
SpecData_PreProcessing <- read_excel("SpecData_PreProcessing.xlsx", 
                                     sheet = "Rem. CCL3")
SpecData <- as.data.frame(SpecData_PreProcessing)
row.names(SpecData) <- SpecData[,1]
SpecData <- SpecData[,-1]

# Processing for Modulation by Target Prediction

TarPred <- read_csv("~/King's BRC DTP Work/Rotation 3/data/miRNA_data/Consolidated_Data_FLAG.csv")
TarPred <- TarPred[,-1]
TarPred$miRNA <- gsub('hsa-', '', TarPred$miRNA)

adj.weight <- data.frame(matrix(nrow=178,ncol=178))
row.names(adj.weight) = row.names(adj)
colnames(adj.weight) = colnames(adj)

as.integer(TarPred %>% 
             filter(miRNA=='miR-122-5p',Gene=='PGF') %>% 
             select('No. Sources'))

row = 'miR-122-5p'
col = 'PGF'
as.integer(TarPred %>% 
             filter(miRNA==row,Gene==col) %>% 
             select('No. Sources'))

for(i in 1:dim(adj)[1]) {
  for(j in 1:dim(adj)[2]){
    row = rownames(adj.weight)[i]
    col = colnames(adj.weight)[j]
    count = as.integer(TarPred %>% filter(miRNA==row,Gene==col) %>% select('No. Sources'))
    count2 = as.integer(TarPred %>% filter(miRNA==col,Gene==row) %>% select('No. Sources'))
    if(grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      adj.weight[i,j] = adj[i,j]
    } else if(!grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      adj.weight[i,j] = adj[i,j]
    } else if(grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      if(is.na(count) || count == '') {
        adj.weight[i,j] = 0
      } else if(count >= 2) {
        adj.weight[i,j] = adj[i,j]*1
      } else if(count == 1) {
        adj.weight[i,j] = adj[i,j]*0.75
      } else {
        adj.weight[i,j] = 0
      }
    } else if(!grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      if(is.na(count2) || count2 == '') {
        adj.weight[i,j] = 0
      } else if(count2 >= 2) {
        adj.weight[i,j] = adj[i,j]*1
      } else if(count2 == 1) {
        adj.weight[i,j] = adj[i,j]*0.75
      } else {
        adj.weight[i,j] = 0
      }
    }
  }
}

rows_to_remove.adj <- which(row.names(adj.weight) %in% prot.F)
cols_to_remove.adj <- which(colnames(adj.weight) %in% miR.F)

adj.weight.F <- adj.weight[-rows_to_remove.adj,-cols_to_remove.adj]

write.csv(adj.weight.F,"Weighted_Adjacencies.csv")

Corr.weight <- data.frame(matrix(nrow=178,ncol=178))
row.names(Corr.weight) = row.names(Corr)
colnames(Corr.weight) = colnames(Corr)

for(i in 1:dim(Corr)[1]) {
  for(j in 1:dim(Corr)[2]){
    row = rownames(Corr.weight)[i]
    col = colnames(Corr.weight)[j]
    count = as.integer(TarPred %>% filter(miRNA==row,Gene==col) %>% select('No. Sources'))
    count2 = as.integer(TarPred %>% filter(miRNA==col,Gene==row) %>% select('No. Sources'))
    if(grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      Corr.weight[i,j] = Corr[i,j]
    } else if(!grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      Corr.weight[i,j] = Corr[i,j]
    } else if(grepl("^miR*|^let*",row) && !grepl("^miR*|^let*",col)) {
      if(is.na(count) || count == '') {
        Corr.weight[i,j] = 0
        } else if(count >= 2) {
          Corr.weight[i,j] = Corr[i,j]*1
        } else if(count == 1) {
          Corr.weight[i,j] = Corr[i,j]*0.75
        } else {
          Corr.weight[i,j] = 0
        }
    } else if(!grepl("^miR*|^let*",row) && grepl("^miR*|^let*",col)) {
      if(is.na(count2) || count2 == '') {
        Corr.weight[i,j] = 0
      } else if(count2 >= 2) {
        Corr.weight[i,j] = Corr[i,j]*1
      } else if(count2 == 1) {
        Corr.weight[i,j] = Corr[i,j]*0.75
      } else {
        Corr.weight[i,j] = 0
      }
    }
  }
}

rows_to_remove.corr <- which(row.names(Corr.weight) %in% prot.F)
cols_to_remove.corr <- which(colnames(Corr.weight) %in% miR.F)

Corr.weight.F <- Corr.weight[-rows_to_remove.adj,-cols_to_remove.adj]

write.csv(Corr.weight.F,"Weighted_Spearman.csv")

# Plot weighted matrices

corrplot(adj, method = "color", order = "hclust", tl.cex = 0.5)
corrplot(Corr, method = "color", order = "hclust", tl.cex = 0.5)

adj.weight2 <- as.matrix(adj.weight)
Corr.weight2 <- as.matrix(Corr.weight)

corrplot(adj.weight2, method = "color", order = "hclust", tl.cex = 0.5)
corrplot(Corr.weight2, method = "color", order = "hclust", tl.cex = 0.5)

