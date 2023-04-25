args <- commandArgs(trailingOnly = TRUE)
tcr_file <- args[1]
bystander_file <- args[2] 

#ro_bystander <- read.csv("~/bystande.csv", header = T)
#bystander_filter <- read.table("/picb/lilab5/liuzipei/vdjdb/vdjdb.filter2.txt", header = F)
options(stringsAsFactors = F)
bystander <- read.table("/picb/lilab5/liuzipei/vdjdb/vdjdb_full.txt", sep = "\t", header = T, fill = T, quote = "") 
bystander_filter <- bystander[bystander$species == "HomoSapiens" & bystander$antigen.species %in% c("CMV","EBV","InfluenzaA") & bystander$vdjdb.score %in% c(NA,"2","3"), ]

tcr <- read.csv(tcr_file, header = T)
bystander <- tcr[tcr$cdr3 %in% bystander_filter$cdr3.beta,]

write.csv(bystander, file = bystander_file, row.names = F, quote = F)
# #table(ro_bystander$CDR3B %in% bystander_filter$V1)
# #table(ro_bystander$CDR3B %in% bystander$cdr3.beta)
# table(ro_bystander$CDR3B %in% bystander_filter$cdr3.beta)
# 
# #table(tcr1$cdr3 %in% bystander$cdr3.beta)
# table(tcr1$cdr3 %in% bystander_filter$cdr3.beta)
# table(tcr1$cdr3[tcr1$cdr3 %in% bystander_filter$cdr3.beta] %in% ro_bystander$CDR3B)
# 
# table(tcr2$cdr3 %in% bystander_filter$cdr3.beta)


