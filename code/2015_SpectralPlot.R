# heatmap 
# source - http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html
# https://rpubs.com/nurakawa/spectral-clustering

install.packages(c("kernlab","cluster","anocva","fpc"))

library(kernlab)
library(cluster)
library (anocva)
library (fpc) # for ch score

setwd("C:\\Bhawana\\Mayr\\ObesityNature\\Spanish cohort\\Control_ObeseT0_Protein")

# read corr data with protein + rna 
numericData <- 1 - SpecData

# rownames(distData ) <- SpecData[,1:1] # change rownames of matrix

############# FINAL#############
sc4 <- specc(numericData, centers=4)
withinss(sc4)
size(sc4)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc4[x]))
}
calinhara(numericData,y,cn=4) 

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_4.png")
plot(sc4,xaxt = 'n')
axis(side=1,at= seq(sc4),labels=rownames(numericData),las=2)
dev.off()

############################################
sc3 <- specc(numericData, centers=3)
withinss(sc3)
size(sc3)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc3[x]))
}
calinhara(numericData,y,cn=3) 
png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_3.png")
plot(sc3,xaxt = 'n')
axis(side=1,at= seq(sc3),labels=rownames(numericData),las=2)
dev.off()


sc5 <- specc(numericData, centers=5)
withinss(sc5)
size(sc5)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc5[x]))
}
calinhara(numericData,y,cn=5) 
png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_5.png")
plot(sc5,xaxt = 'n')
axis(side=1,at= seq(sc5),labels=rownames(numericData),las=2)
dev.off()

sc7 <- specc(numericData, centers=7)
withinss(sc7)
size(sc7)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc7[x]))
}
calinhara(numericData,y,cn=7) 
png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_7.png")
plot(sc7,xaxt = 'n')
axis(side=1,at= seq(sc7),labels=rownames(numericData),las=2)
dev.off()

sc10 <- specc(numericData, centers=10)
withinss(sc10)
#centers(sc)
#size(sc)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc10[x]))
}
calinhara(numericData,y,cn=10) 

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_10.png")
plot(sc10,xaxt = 'n')
axis(side=1,at= seq(sc10),labels=rownames(numericData),las=2)
dev.off()


sc12 <- specc(numericData, centers=12)
withinss(sc12)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc12[x]))
}
calinhara(numericData,y) # 26.303

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_12.png")
plot(sc12,xaxt = 'n')
axis(side=1,at= seq(sc12),labels=rownames(numericData),las=2)
dev.off()

sc15 <- specc(numericData, centers=15)
withinss(sc15)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc15[x]))
}
calinhara(numericData,y) # 26.303

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_15.png")
plot(sc15,xaxt = 'n')
axis(side=1,at= seq(sc15),labels=rownames(numericData),las=2)
dev.off()

sc18 <- specc(numericData, centers=18)
withinss(sc18)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc18[x]))
}
calinhara(numericData,y) # 26.303

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_18.png")
plot(sc18,xaxt = 'n')
axis(side=1,at= seq(sc18),labels=rownames(numericData),las=2)
dev.off()

sc20 <- specc(numericData, centers=20)
withinss(sc20)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
	y <- as.integer(c(y,sc20[x]))
}
calinhara(numericData,y) # 26.303

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_20.png")
plot(sc20,xaxt = 'n')
axis(side=1,at= seq(sc20),labels=rownames(numericData),las=2)
dev.off()


# check margin
#graphics.off()
#par("mar")
#par(mar=c(1,1,1,1))

#png(height=1700, width=600, pointsize=15, file="Spectral_RNA_AGI.png")
#plot(corrData , col=sc)
#dev.off()

#https://www.rdocumentation.org/packages/fpc/versions/2.1-11.1/topics/calinhara
#calinhara(distData,size(sc20))
#sc10$get(`Cluster memberships`) # access white space attributes
# attributes(sc10) # get attribite names

#https://www.rdocumentation.org/packages/fpc/versions/2.1-11.1/topics/cluster.stats
#cluster.stats(distData,size(sc20),silhouette=TRUE)

#https://www.rdocumentation.org/packages/fpc/versions/2.1-11.2/topics/pamk
#pamk(distData ,krange=2:20,criterion="asw", usepam=TRUE,
#     scaling=FALSE, alpha=0.001, diss=inherits(distData , "dist"),
#     critout=FALSE, seed=NULL)
# x <- numericData[,2:ncol(numericData)]
#silhouette( size(sc20),x)
#write.csv(sc.'Cluster memberships', file = "RNA_AGI_SpecClusterMember.csv"
#head(sc)
#sc

sc6 <- specc(numericData, centers=6)
withinss(sc6)
size(sc6)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
  y <- as.integer(c(y,sc6[x]))
}
calinhara(numericData,y,cn=6) 

sc8 <- specc(numericData, centers=8)
withinss(sc8)
size(sc8)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
  y <- as.integer(c(y,sc8[x]))
}
calinhara(numericData,y,cn=8) 

png(height=600, width=1700, pointsize=15, file="Spectral_RNA_AGI_8.png")
plot(sc8,xaxt = 'n')
axis(side=1,at= seq(sc12),labels=rownames(numericData),las=2)
dev.off()

sc9 <- specc(numericData, centers=9)
withinss(sc9)
size(sc9)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
  y <- as.integer(c(y,sc9[x]))
}
calinhara(numericData,y,cn=9) 

sc3 <- specc(numericData, centers=3)
withinss(sc3)
size(sc3)
# loop to extract Cluster membership
y=list()
for (x in c(1:ncol(numericData))){
  y <- as.integer(c(y,sc3[x]))
}
calinhara(numericData,y,cn=3) 
png(height=600, width=1700, pointsize=15, file="Spectral_Final_3.png")
plot(sc3,xaxt = 'n')
axis(side=1,at= seq(sc3),labels=rownames(numericData),las=2)
dev.off()