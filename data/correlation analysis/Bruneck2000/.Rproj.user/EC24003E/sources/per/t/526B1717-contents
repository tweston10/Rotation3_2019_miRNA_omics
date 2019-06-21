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

# Cox Regression Analysis - Processing

Cox_data <- merge(Bru_miR, Bru_Outcome, by = 'ID', all = TRUE) %>% 
  filter(!is.na(CVD0010)) %>% 
  filter(!is.na(`miR-122`))

write_csv(Cox_data, path = "Cox_dataset.csv")
write_csv(Cox_data[1,2:23], path = "Cox_labels.csv")
