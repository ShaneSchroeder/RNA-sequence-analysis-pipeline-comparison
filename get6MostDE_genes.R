library(dplyr)
DEG_Data <- read.table("~/rsem/exp/geneMat.de.txt")
DEG_Data_desc <- DEG_Data %>% 
  arrange(desc(RealFC))
data <- DEG_Data_desc[1:6,]
write.csv(data, "~/rsem/exp/6mostSigGenes.csv")
