
library(data.table)

varTraitData = fread("dnanexus/regenie/conditional_per_var/NC.single.cases.FDR0.05.tsv")
varTraitData = unique(varTraitData[,.(TRAIT,VAR)])

ldData = fread("dnanexus/regenie/conditional_per_var/ld_results_parsed_reformat.txt")
colnames(ldData) = c("topVar","condVar","R2","D_prime")
ldData$condVar = data.table(unlist(lapply(ldData$condVar, function(x) unlist(strsplit(x,"[:]"))[1])))$V1
cutoff = 0.5
ldData = ldData[D_prime > cutoff]

cojoData = fread("dnanexus/regenie/conditional_per_var/all_trait.conditional.cojo", header = F)
cojoData$TRAIT = data.table(unlist(lapply(cojoData$V2, function(x) unlist(strsplit(x,"[.]"))[1])))$V1



for (i in seq(nrow(varTraitData))){
  entry = varTraitData[i,]
  var = entry$VAR
  trait = entry$TRAIT
  trait = strsplit(trait,"[p]")[[1]][[2]]
  print(trait)
  
  if (var %in% ldData$topVar){
    ldVars = unique(ldData[topVar == var]$condVar)
    condVars = cojoData[TRAIT == trait]
    
    if (nrow(condVars) == 0){    print("PROBLEM with CONDITIONAL vars")  }
    
    if (length(ldVars) > 0){  
      subtract = condVars[!(condVars$V1 %in% ldVars)]
      if (nrow(condVars) != nrow(subtract)){ write.table(subtract$V1, paste0("dnanexus/",gsub(":", "_", var),".",trait,".subtract.txt"), row.names=F, col.names=F, quote=F ) }
      
      overlap = condVars[V1 %in% ldVars]
      if (nrow(condVars) != nrow(overlap)){ write.table(overlap$V1, paste0("dnanexus/",gsub(":", "_", var),".",trait,".intersect.txt"), row.names=F, col.names=F, quote=F ) }
      
    }
  }
}
