library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

miRNA_list <- list('hsa-miR-122-3p',
                   'hsa-miR-122-5p',
                   'hsa-miR-126-3p',
                   'hsa-miR-126-3p',
                   'hsa-miR-126-5p',
                   'hsa-miR-126-5p',
                   'hsa-miR-146a-3p',
                   'hsa-miR-146a-5p',
                   'hsa-miR-148a-3p',
                   'hsa-miR-148a-5p',
                   'hsa-miR-150-3p',
                   'hsa-miR-150-5p',
                   'hsa-miR-155-3p',
                   'hsa-miR-155-5p',
                   'hsa-miR-191-3p',
                   'hsa-miR-191-5p',
                   'hsa-miR-195-3p',
                   'hsa-miR-195-5p',
                   'hsa-miR-197-3p',
                   'hsa-miR-197-5p',
                   'hsa-miR-21-3p',
                   'hsa-miR-21-5p',
                   'hsa-miR-223-3p',
                   'hsa-miR-223-3p',
                   'hsa-miR-223-5p',
                   'hsa-miR-223-5p',
                   'hsa-miR-24-3p',
                   'hsa-miR-24-5p',
                   'hsa-miR-27b-3p',
                   'hsa-miR-27b-5p',
                   'hsa-miR-28-3p',
                   'hsa-miR-320-3p',
                   'hsa-miR-320-5p',
                   'hsa-miR-335-3p',
                   'hsa-miR-335-5p',
                   'hsa-miR-378-3p',
                   'hsa-miR-378-5p',
                   'hsa-miR-92a-3p',
                   'hsa-miR-92a-5p'
)

miRWalk.list <- list.files(pattern = 'miRWalk*')
df.list <- lapply(miRWalk.list, read_csv)

df <- rbindlist(df.list, use.names = TRUE, fill = TRUE, idcol = NULL)
miRWalk <- df %>% select(mirnaid, genesymbol, bindingp, accessibility, position)
names(miRWalk)
names(miRWalk)[names(miRWalk)=="mirnaid"] <- "miRNA"
names(miRWalk)[names(miRWalk)=="genesymbol"] <- "Gene"
names(miRWalk)[names(miRWalk)=="bindingp"] <- "Score"
names(miRWalk)[names(miRWalk)=="accessibility"] <- "Accessibility"
names(miRWalk)[names(miRWalk)=="position"] <- "Position"

miRDB_Targets <- read_csv("miRDB_target_prediction.csv")
mirDIP_Targets <- read_csv("mirDIP_Targets.csv")
PolymiRTS_Targets <- read_excel("PolymiRTS_3.0_Targets.xlsx")
TargetScan_Targets <- read_csv("TargetScan/Collated_data_FIXED.csv", na = "N/A")
DIANA_Targets <- read_csv("DIANA/Collated_data_FIX.csv")

PolymiRTS_Targets <- fill(PolymiRTS_Targets, 'Gene Symbol', 'Description')

miRDB <- miRDB_Targets %>% select('miRNA Name', 'Gene Symbol', 'Target Score')
names(miRDB)[1:3] = c("miRNA", "Gene", "Score")

mirDIP <- mirDIP_Targets %>% select('MicroRNA', 'Gene Symbol', 'Integrated Score')
names(mirDIP)[1:3] = c("miRNA", "Gene", "Score")

PolymiRTS <- PolymiRTS_Targets %>% select('miR ID', 'Gene Symbol')
names(PolymiRTS)[1:2] = c("miRNA", "Gene")
PolymiRTS$Score <- 1

TargetScan <- TargetScan_Targets %>% 
  filter(`Representative miRNA` %in% miRNA_list) %>% 
  select('Representative miRNA', 'Target gene', 'Cumulative weighted context++ score')
names(TargetScan)[1:3] = c("miRNA", "Gene", "Score")

DIANA <- DIANA_Targets %>% select('Mirna Name', 'Gene name', 'miTG score')
names(DIANA)[1:3] = c("miRNA", "Gene", "Score")

Consol.List <- list(miRWalk, miRDB, mirDIP, PolymiRTS, TargetScan, DIANA)
df <- rbindlist(Consol.List, use.names = TRUE, fill = TRUE, idcol = "Algorithm")
df$Algorithm[df$Algorithm %in% 1] <- "miRWalk"
df$Algorithm[df$Algorithm %in% 2] <- "miRDB"
df$Algorithm[df$Algorithm %in% 3] <- "mirDIP"
df$Algorithm[df$Algorithm %in% 4] <- "PolymiRTS"
df$Algorithm[df$Algorithm %in% 5] <- "TargetScan"
df$Algorithm[df$Algorithm %in% 6] <- "DIANA"

write.csv(df, file = "Consolidated_Data.csv")

df2 <- df %>% 
  select('Algorithm', 'miRNA', 'Gene', 'Score', 'Accessibility') %>% 
  distinct(Algorithm, miRNA, Gene, .keep_all = TRUE)

Table1 <- df2 %>% 
  select(-Accessibility) %>% 
  group_by(miRNA, Gene) %>% 
  spread(Algorithm, Score) %>% 
  arrange(miRNA, Gene)
View(Table1)

Access <- df2 %>% 
  select(miRNA, Gene, Accessibility) %>% 
  filter(!is.na(Accessibility))

Table2 <- df2 %>% 
  select(-Accessibility) %>% 
  group_by(Gene, miRNA) %>% 
  spread(Algorithm, Score) %>% 
  arrange(Gene, miRNA)
View(Table2)

write.csv(Table1, file = "Consolidated_Data_ByMiRNA.csv")
write.csv(Table2, file = "Consolidated_Data_ByGene.csv")

NormTable <- Table1 %>% 
  mutate(DIANA_N = 0.5+((DIANA-0)/(2*(1-0))),
         miRDB_N = 0.5+((miRDB-0)/(2*(100-0))),
         mirDIP_N = 0.5+((mirDIP-0)/(2*(1-0))),
         miRWalk_N = 0.5+((miRWalk-0)/(2*(1-0))),
         TargetScan_N = 0.5+(((1-(2^TargetScan))-0)/(2*(1-0)))) %>% 
  select('miRNA', 'Gene', 'DIANA_N', 'miRDB_N', 'mirDIP_N', 
         'miRWalk_N', 'PolymiRTS', 'TargetScan_N')
View(NormTable)

NormTable2 <- NormTable %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(Av_Score = (DIANA_N + miRDB_N + mirDIP_N + miRWalk_N + TargetScan_N + PolymiRTS)/6)
View(NormTable2)

FinalTable <- NormTable2 %>% 
  mutate('No. Sources' = 
           if_else(DIANA_N > 0,1,0) +
           if_else(miRDB_N > 0,1,0) +
           if_else(mirDIP_N > 0,1,0) +
           if_else(miRWalk_N > 0,1,0) +
           if_else(PolymiRTS > 0,1,0) +
           if_else(TargetScan_N > 0,1,0))
View(FinalTable)

FinalTable2 <- FinalTable %>% 
  mutate('2+Tools' = if_else(`No. Sources`>1,1,0),
         '3+Tools' = if_else(`No. Sources`>2,1,0))
View(FinalTable2)

nrow(FinalTable2)
sum(FinalTable2$`2+Tools`)
sum(FinalTable2$`3+Tools`)

range(miRDB$Score)
range(mirDIP$Score)
range(miRWalk$Score)
range(PolymiRTS$Score)
range(TargetScan$Score)
range(DIANA$Score)

write.csv(FinalTable2, file = "Consolidated_Data_FLAG.csv")
