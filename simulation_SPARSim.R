sim_fc = function(sce_mat,ct1 = Chu_param_preset$Chu_C6,ngene = 300,ncell = 100,seed = 2022,OS = "Windows"){
  set.seed(seed)
  nct = nrow(sce_mat)
  ndat = ncol(sce_mat)
  len = length(ct1[["intensity"]])
  
  simdat = list()
  SPARSim_result = list()
  if(OS == "windows"){
    cl <- makeCluster(ndat)
    registerDoSNOW(cl)
  }else if (OS == "Linux"){
    registerDoMC(cores=ndat)
  }
  
  
  SPARSim_result = foreach(j = 1:ndat, .packages = c("MASS","RMTL","psych","corpcor","fields","glmnet","SPARSim"), .errorhandling='pass') %dopar% {
    simdat[[j]] = list()
    simdat[[j]][[1]] = ct1
    for (i in 2:nct) {
      ind = sce_mat[i,j]
      
      fc_mul <-c(rep(ind,ngene),rep(1,len-ngene))
      
      simdat[[j]][[i]] <- SPARSim_create_DE_genes_parameter(
        sim_param = ct1,
        fc_multiplier = fc_mul,
        N_cells = ncell,
        condition_name = paste0("cell_",LETTERS[i]))
    }
    res = SPARSim_simulation(simdat[[j]])$count_matrix
    return(res)
  }
  
  if(OS == "windows"){
    stopCluster(cl)
  }
  
  simdata_sc_fc = lapply(SPARSim_result, function(x) return(log1p(x)))
  
  simdata_sc_fc_sub = lapply(simdata_sc_fc, function(x){
    x = x[c(1:ngene),]
  })
  
  return(simdata_sc_fc_sub)
}