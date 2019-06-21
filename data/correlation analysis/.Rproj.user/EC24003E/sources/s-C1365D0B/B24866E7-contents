install.packages("BiocManager") 
BiocManager::install("WGCNA")

install.packages("corrplot")

library(tidyverse)
library(corrplot)
library(WGCNA)
library(readxl)
library(data.table)

Bru_miR <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2015/Bruneck 2015_platelet_poor_plasma_RNAs.xlsx", 
                      sheet = "PPP_celnorm_RQ", skip = 1)
Bru_Olink_CVD <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2015/Bruneck2015_platelet_poor_plasma_OlinkCVD1-INF.xlsx", 
                            sheet = "CVD1_NaN", na = "NaN", skip = 1)
Bru_Olink_INF <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2015/Bruneck2015_platelet_poor_plasma_OlinkCVD1-INF.xlsx", 
                            sheet = "INF_NaN", na = "NaN", skip = 1)

Bru_Olink_CVD$`Study Number` <- gsub('\\s+', '', Bru_Olink_CVD$`Study Number`)
Bru_Olink_INF$`Study Number` <- gsub('\\s+', '', Bru_Olink_INF$`Study Number`)

colnames(Bru_Olink_CVD)[4:95] <- paste(colnames(Bru_Olink_CVD)[4:95], "(CVD)", sep = "_")
colnames(Bru_Olink_INF)[4:95] <- paste(colnames(Bru_Olink_INF)[4:95], "(INF)", sep = "_")

names(Bru_Olink_CVD)[2] <- "BruneckID"
names(Bru_Olink_INF)[2] <- "BruneckID"

df <- merge(Bru_Olink_CVD, Bru_Olink_INF, by = 'BruneckID', all = TRUE)
df2 <- merge(df, Bru_miR, by = 'BruneckID', all = TRUE)

names(df2)
drops = c("Date of Birth.x", "NPX.x", "#Flagged.x", "#Chip name.x", "Date of Birth.y", 
          "NPX.y", "#Flagged.y", "#Chip name.y", "DoB", "London ID")
df3 <- df2[ , !(names(df2) %in% drops)]
df4 <- df3[1:338,]

df4 <- df4[,colSums(is.na(df4)) <= 170]
colSums(is.na(df4))

Bru_Scaled <- scale(df4[2:205])
Bru_Coll <- cbind(df4[1],Bru_Scaled)
apply(Bru_Coll,2,sd)

Bru_Coll2 <- as.data.frame(t(Bru_Coll[,-1]))
colnames(Bru_Coll2) <- Bru_Coll$BruneckID
Bru_Coll2$Node <- factor(row.names(Bru_Coll2))

Bru_Coll3 <- Bru_Coll2 %>% 
  select('Node', 1:338)
rownames(Bru_Coll3) <- NULL

write.csv(Bru_Coll3, file = "Bruneck2015_scaled_data.csv")
