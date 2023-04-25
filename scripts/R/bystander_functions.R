options(stringsAsFactors = F)
get_virus_cdr3b <- function(df){
  bystander <- read.table("/picb/lilab5/liuzipei/vdjdb/vdjdb_full.txt", sep = "\t", header = T, fill = T, quote = "") 
  bystander_filter <- bystander[bystander$species == "HomoSapiens" & bystander$antigen.species %in% c("CMV","EBV","InfluenzaA") & bystander$vdjdb.score %in% c(NA,"1","2","3"), ]
  return(bystander_filter$cdr3.beta)
}

# #table(ro_bystander$CDR3B %in% bystander_filter$V1)
# #table(ro_bystander$CDR3B %in% bystander$cdr3.beta)
# table(ro_bystander$CDR3B %in% bystander_filter$cdr3.beta)
# 
# #table(tcr1$cdr3 %in% bystander$cdr3.beta)
# table(tcr1$cdr3 %in% bystander_filter$cdr3.beta)
# table(tcr1$cdr3[tcr1$cdr3 %in% bystander_filter$cdr3.beta] %in% ro_bystander$CDR3B)
# 
# table(tcr2$cdr3 %in% bystander_filter$cdr3.beta)


