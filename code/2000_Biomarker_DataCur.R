library(tidyverse)
library(readxl)
library(data.table)

Bru_miR <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2000/Bruneck 2000 Plasma cel normalisation.xlsx", skip = 10)
Bru_Agilent <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2000/Bruneck 2000 Agilent.xlsx")
Bru_Olink <- read_excel("~/King's BRC DTP Work/Rotation 3/data/Bruneck_data_for_timir/bruneck2000/OLINK_BRUNECK_2000 (1) (2).xlsx")

A_proteins <- c("ID", Bru_Agilent$Proteins)
A_ID <- colnames(Bru_Agilent)

Bru_Agilent <- as.data.frame(t(Bru_Agilent))
Bru_Agilent <- rownames_to_column(Bru_Agilent, 'ID')
colnames(Bru_Agilent) <- A_proteins
Bru_Agilent <- Bru_Agilent[-1,]

Bru_miR$`Bruneck code` <- gsub('\\s+', '', Bru_miR$`Bruneck code`)
Bru_Olink$`ID` <- gsub('\\s+', '', Bru_Olink$`ID`)

Bru_miR <- Bru_miR[,-1]
names(Bru_miR)[1] <- 'ID'

df <- merge(Bru_Agilent, Bru_Olink, by = 'ID', all = TRUE)
df2 <- merge(df, Bru_miR, by = 'ID', all = TRUE)

df3 <- df2[rowSums(is.na(df2)) <= 152,]
df4 <- df3[,colSums(is.na(df3)) <= 366]

colSums(is.na(df4))

tmp <- mutate_all(df4[2:304], function(x) as.numeric(as.character(x)))

colSums(is.na(tmp))

Bru_Scaled <- scale(tmp)
Bru_Coll <- cbind(df4[1],Bru_Scaled)
apply(Bru_Coll,2,sd)

Bru_Coll2 <- as.data.frame(t(Bru_Coll[,-1]))
colnames(Bru_Coll2) <- Bru_Coll$ID

Bru_Coll3 <- rownames_to_column(Bru_Coll2, 'Node')

write.csv(Bru_Coll3, file = "Bruneck2000_scaled_data.csv")

# Statistical Comparison of Outcomes - Processing

Bru_Outcome <- read_excel("Bruneck2000_Outcome_ForTimir.xlsx")
names(Bru_Outcome)[2] <- 'ID'

StatComp.All <- scale(Bru_miR[2:23]) %>% 
  cbind(Bru_miR[1]) %>% 
  merge(Bru_Outcome, by = 'ID', all = TRUE) %>% 
  filter(!is.na(CVD0010)) %>% 
  filter(!is.na(`miR-122`))

StatComp.Data <- StatComp.All[1:23] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ')

StatComp.Data[1,1] <- ''

write.table(StatComp.Data, file = "miR_database.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

StatComp.Lab <- StatComp.All['CVD0010']
StatComp.Lab$CVD0010 <- as.character(StatComp.Lab$CVD0010)
StatComp.Lab$CVD0010[StatComp.Lab$CVD0010 == "1"] <- "Yes"
StatComp.Lab$CVD0010[StatComp.Lab$CVD0010 == "0"] <- "No"

write.table(t(StatComp.Lab), file = "miR_labels.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

StatComp.CoMor <- StatComp.All[c('AGE','SEX')]

write.table(StatComp.Cont, file = "miR_comorbidities.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

StatComp.Cont <- cbind(StatComp.All[1],StatComp.All[25:32])

write.table(StatComp.Cont, file = "miR_continuous.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)

Stat.P1 <- c('PLF4',
             'STOM',
             'CD40.L',
             'EGF',
             'SRC',
             'PAR.1',
             'PDGF.subunit.B',
             'CXCL1.l1',
             'LOX.1',
             'IL.18.l1',
             'HSP.27',
             'TNFSF14.l1',
             'MPO',
             'MMP.1.l1',
             'CXCL6.l1',
             'ITGB1BP2',
             'Dkk.1',
             'SIRT2.l1',
             'EN.RAGE.l1',
             'CD40.l1',
             'HB.EGF',
             'CASP.8.l1',
             'NEMO',
             'PECAM.1',
             'ECP',
             'CD244',
             'IL.7',
             'LAP.TGF.beta.1',
             'AXIN1',
             'OSM',
             'MCP.4',
             'CXCL5',
             'X4E.BP1',
             'MCP.2',
             'ST1A1',
             'STAMPB')
  
Stat.P2 <- c('ABCF1_HUMAN',
             'AFAM_HUMAN',
             'RAGE',
             'FETUA',
             'ALBU_HUMAN',
             'APOA2_HUMAN',
             'APOA4_HUMAN',
             'APOC2',
             'APOC3',
             'APOD',
             'APOL1',
             'ZA2G',
             'C1QB',
             'CCL19',
             'CCL23',
             'CCL28',
             'CCL3',
             'CCL4.l1',
             'MCP.3',
             'CFAB',
             'CFAI',
             'CHI3L1',
             'CERU',
             'CST5',
             'CSTB',
             'CTSD',
             'CX3CL1.l1',
             'IL.8.l1',
             'CXCL9',
             'ESM.1',
             'FBLN1',
             'FCN3',
             'FIBA',
             'FGF.21',
             'FGF.5',
             'FIBG',
             'Flt3L',
             'GDF.15',
             'GH',
             'GELS',
             'IGHA2',
             'IGHG3',
             'IGHG4',
             'IL.10',
             'IL.10RB',
             'IL.15RA',
             'IL.6.l1',
             'SCF.l1',
             'hK11',
             'A2GL_HUMAN',
             'BNP',
             't.PA',
             'U.PAR',
             'CO9',
             'EN.RAGE.l1',
             'SELE',
             'PSGL.1',
             'CBG',
             'SHBG',
             'TRFE',
             'TRAIL.R2',
             'TNFRSF9',
             'TWEAK',
             'VTNC',
             'mAmP')

Stat.P3 <- c('CCL23',
             'CCL3',
             'CCL4',
             'CCL7',
             'CCL8',
             'CD244',
             'CD274',
             'CD40',
             'CSTB',
             'CXCL1',
             'F3',
             'FABP4',
             'HGF',
             'IL10RA',
             'MMP1',
             'PLAT',
             'RNASE3',
             'SPON1',
             'TNFRSF10B',
             'TNFRSF1A')

StatComp.P <- merge(df, Bru_Outcome, by = 'ID', all = TRUE) %>% 
  filter(!is.na(CVD0010))

StatComp.P <- StatComp.P[rowSums(!is.na(StatComp.P)) > 200,]

View(head(StatComp.P))
rowSums(!is.na(StatComp.P))

rows_to_keep.S1 <- which(colnames(StatComp.P) %in% Stat.P1)
rows_to_keep.S2 <- which(colnames(StatComp.P) %in% Stat.P2)
rows_to_keep.S3 <- which(colnames(StatComp.P) %in% Stat.P3)

Stat_data.WGCNA <- StatComp.P[,rows_to_keep.S1] %>% 
  mutate_all(function(x) as.numeric(as.character(x)))
Stat_data.ARACNe <- StatComp.P[,rows_to_keep.S2] %>% 
  mutate_all(function(x) as.numeric(as.character(x)))
Stat_data.ARACNe2 <- StatComp.P[,rows_to_keep.S3] %>% 
  mutate_all(function(x) as.numeric(as.character(x)))

View(head(Stat_data.WGCNA))
View(head(Stat_data.ARACNe))

Stat_data.WGCNA <- cbind(StatComp.P$ID,Stat_data.WGCNA) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ')
Stat_data.ARACNe <- cbind(StatComp.P$ID,Stat_data.ARACNe) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ')
Stat_data.ARACNe2 <- cbind(StatComp.P$ID,Stat_data.ARACNe2) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(' ')

Stat_data.WGCNA[1,1] <- ''
Stat_data.ARACNe[1,1] <- ''

write.table(Stat_data.WGCNA, file = "Stat_Dataset_WGCNA.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)
write.table(Stat_data.ARACNe, file = "Stat_Dataset_ARACNe.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

Stat_label.WGCNA <- StatComp.P['CVD0010']
Stat_label.WGCNA$CVD0010 <- as.character(Stat_label.WGCNA$CVD0010)
Stat_label.WGCNA$CVD0010[Stat_label.WGCNA$CVD0010 == "1"] <- "Yes"
Stat_label.WGCNA$CVD0010[Stat_label.WGCNA$CVD0010 == "0"] <- "No"

write.table(t(Stat_label.WGCNA), file = "Stat_Labels.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

write.table(StatComp.P[c('AGE','SEX')], file = "Stat_CoMorbids_AgeSex.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

write.table(StatComp.P[c('AGE','SEX','Triglycerides','Total_cholesterol','HDL_C','LDL_C')], 
            file = "Stat_CoMorbids_AgeSexLipids.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE)

View(head(StatComp.P))

StatComp.Cont <- cbind(StatComp.P[1],StatComp.P[284:291])

write.table(StatComp.Cont, file = "Stat_Cont.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)

# Cox Regression Analysis - Processing

Cox_data <- merge(Bru_miR, Bru_Outcome, by = 'ID', all = TRUE) %>% 
  filter(!is.na(CVD0010)) %>% 
  filter(!is.na(`miR-122`))

write_csv(Cox_data, path = "Cox_dataset.csv")
write_csv(Cox_data[1,2:23], path = "Cox_labels.csv")

Cox.P1 <- c('AXIN1',
            'CASP.8.l1',
            'CD244',
            'PDGF.subunit.B',
            'CD40.l1',
            'CD40.L',
            'CXCL1.l1',
            'CXCL5',
            'CXCL6.l1',
            'Dkk.1',
            'EGF',
            'X4E.BP1',
            'PAR.1',
            'HB.EGF',
            'HSP.27',
            'NEMO',
            'IL.7',
            'ITGB1BP2',
            'MMP.1.l1',
            'PECAM.1',
            'PLF4',
            'SIRT2.l1',
            'SRC',
            'STAMPB',
            'STOM',
            'ST1A1',
            'LAP.TGF.beta.1',
            'TNFSF14.l1')

Cox.P2 <- c('ABCF1_HUMAN',
            'FETUA',
            'APOA4_HUMAN',
            'APOC2',
            'CCL23',
            'CCL3',
            'CCL4.l1',
            'MCP.3',
            'MCP.2',
            'PDGF.subunit.B',
            'CSTB',
            'CXCL1.l1',
            'CXCL16',
            'FGF.21',
            'IGHG3',
            'IL.10RA',
            'MMP.1.l1',
            'CA.125',
            't.PA',
            'PRL',
            'SELE',
            'PSGL.1',
            'CBG',
            'TRFE',
            'TRAIL.R2')

Cox_data.P <- merge(df, Bru_Outcome, by = 'ID', all = TRUE) %>% 
  filter(!is.na(CVD0010))

rows_to_remove.P1 <- which(colnames(Cox_data.P) %in% Stat.P1)
rows_to_remove.P2 <- which(colnames(Cox_data.P) %in% Stat.P2)

Cox_data.P1 <- Cox_data.P[,rows_to_remove.P1] %>% 
  mutate_all(function(x) as.numeric(as.character(x)))
Cox_data.P2 <- Cox_data.P[,rows_to_remove.P2] %>% 
  mutate_all(function(x) as.numeric(as.character(x)))

rowSums(is.na(Cox_data.P1))

View(head(Cox_data.P1))
View(head(Cox_data.P2))

Cox_data.P1 <- cbind(Cox_data.P$ID,Cox_data.P1,Cox_data.P[,284:291])
Cox_data.P2 <- cbind(Cox_data.P$ID,Cox_data.P2,Cox_data.P[,284:291])

Cox_data.P1 <- Cox_data.P1[rowSums(is.na(Cox_data.P1)) <= 10,]
Cox_data.P2 <- Cox_data.P2[rowSums(is.na(Cox_data.P2)) <= 10,]

write_csv(Cox_data.P1, path = "Cox_dataset_p1.csv")
write_csv(Cox_data.P1[1,2:29], path = "Cox_labels_p1.csv")

write_csv(Cox_data.P2, path = "Cox_dataset_p2.csv")
write_csv(Cox_data.P2[1,2:26], path = "Cox_labels_p2.csv")
