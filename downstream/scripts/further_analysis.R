#! R

##what questions do we have regarding overlap of msVIPER analysis, KM based on multiple factors?
viperSet <- c("UCS", "UCEC", "UCEC.epi", "UCEC.str")
viperList <- lapply(viperSet, function(f){
  uset <- get(paste0("viper.", f))
  uset <- uset[[2]][uset[[2]]$q.value < 0.05,]
  return(uset)})
names(viperList) <- viperSet

rownames(viper.UCEC.epi[[2]][viper.UCEC.epi[[2]]$q.value<0.01,]) %in% rownames(viper.UCEC.str[[2]][viper.UCEC.str[[2]]$q.value<0.05,])
rownames(viper.UCEC.epi[[2]][viper.UCEC.epi[[2]]$q.value<0.05,]) %in% rownames(viper.UCEC.str[[2]][viper.UCEC.str[[2]]$q.value<0.05,])
rownames(viper.UCEC.epi[[2]][viper.UCEC.epi[[2]]$q.value<0.05,])[rownames(viper.UCEC.epi[[2]][viper.UCEC.epi[[2]]$q.value<0.05,]) %in% rownames(viper.UCEC.str[[2]][viper.UCEC.str[[2]]$q.value<0.05,])]
viper.UCS[rownames(viper.UCS) %in% "ILK",]
viper.UCS[[2]][rownames(viper.UCS[[2]]) %in% "ILK",]
viper.UCEC[[2]][rownames(viper.UCEC[[2]]) %in% "ILK",]
rn.viper.UCEC <- rownames(viper.UCEC[[2]])
rn.viper.UCS <- rownames(viper.UCS[[2]])
rn.viper.UCEC.epi <- rownames(viper.UCEC.epi[[2]])
rn.viper.UCEC.str <- rownames(viper.UCEC.str[[2]])
rn.UCEC.UCS <- rn.viper.UCEC[rn.viper.UCEC %in% rn.viper.UCS]
rn.UCEC.UCS
rn.viper.UCEC.str <- rownames(viper.UCEC.str[[2]][viper.UCEC.str[[2]]$q.value < 0.05,])
rn.viper.UCEC.str
rn.viper.UCEC.epi <- rownames(viper.UCEC.epi[[2]][viper.UCEC.epi[[2]]$q.value < 0.05,])
rn.viper.UCEC <- rownames(viper.UCEC[[2]][viper.UCEC[[2]]$q.value < 0.05,])
rn.viper.UCS <- rownames(viper.UCS[[2]][viper.UCS[[2]]$q.value < 0.05,])
rn.UCEC.UCS <- rn.viper.UCEC[rn.viper.UCEC %in% rn.viper.UCS]
rn.UCEC.UCS
rn.UCEC.UCS.epi <- rn.UCEC.UCS[rn.UCEC.UCS %in% rn.UCEC.epi]
rn.UCEC.UCS.epi <- rn.UCEC.UCS[rn.UCEC.UCS %in% rn.viper.UCEC.epi]
rn.UCEC.UCS.epi
rn.UCEC.UCS.epi.str <- rn.UCEC.UCS.epi[rn.UCEC.UCS.epi %in% rn.viper.UCEC.str]
rn.UCEC.UCS.epi.str
