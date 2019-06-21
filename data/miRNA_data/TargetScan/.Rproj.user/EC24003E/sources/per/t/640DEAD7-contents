install.packages("data.table")

library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

file.list <- list.files(pattern='*.xlsx')
df.list <- lapply(file.list, read_excel)

df <- rbindlist(df.list, use.names = TRUE, fill = TRUE, idcol = NULL)

miRNA_list <- list('hsa-miR-16-5p',
                   'hsa-miR-19b-3p',
                   'hsa-miR-20b-5p',
                   'hsa-miR-21-5p',
                   'hsa-miR-24-3p',
                   'hsa-miR-26a-5p',
                   'hsa-miR-27b-3p',
                   'hsa-miR-28-5p',
                   'hsa-miR-29a-3p',
                   'hsa-miR-30b-5p',
                   'hsa-miR-92a-3p',
                   'hsa-miR-93-5p',
                   'hsa-miR-99b-5p',
                   'hsa-miR-122-5p',
                   'hsa-miR-126-3p',
                   'hsa-miR-126-5p',
                   'hsa-miR-130a-3p',
                   'hsa-miR-140-5p',
                   'hsa-miR-143-3p',
                   'hsa-miR-145-5p',
                   'hsa-miR-150-5p',
                   'hsa-miR-191-5p',
                   'hsa-miR-192-5p',
                   'hsa-miR-197-3p',
                   'hsa-miR-221-3p',
                   'hsa-miR-222-3p',
                   'hsa-miR-223-3p',
                   'hsa-miR-223-5p',
                   'hsa-miR-320a',
                   'hsa-miR-324-5p',
                   'hsa-miR-375',
                   'hsa-miR-486-5p',
                   'hsa-miR-574-3p',
                   'hsa-miR-885-5p',
                   'hsa-let-7b-5p',
                   'hsa-let-7d-5p')

miRNA_data <- df %>% filter(`Representative miRNA` %in% miRNA_list)

miRNA_data2 <- miRNA_data %>% 
  mutate_all(~replace(., is.na(.), 0))

write.csv(miRNA_data2, file = "Collated_data.csv")
