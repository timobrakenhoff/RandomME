### SIMULATION FUNCTION for SMART data selection to create Heatmaps

me_mod <- function(df,out,expo,confs,ME_conf,ME_exp_sd,ME_conf_sd,sims=1,plt=F,sav=F,seed=1){
  
  #make fu val from outcome
  fu_out <- paste0("fu_",out)
  
  #setseed
  set.seed(seed)
  
  #Make formula for surival model 
  surv_mod <- as.formula(paste("surv.ob~",paste(c(expo,confs),collapse="+")))
  
  #Make grid of ME (ME as % of total variance)
  ME.list <- list(ME_exp=ME_exp_sd,ME_con=ME_conf_sd)
  
  ME.grid <- expand.grid(ME.list)
  
  #Make results array
  res.mat <- array(NA,dim=c(nrow(ME.grid),3,sims))
  colnames(res.mat) <- c("exp_coef","exp_ci_l","exp_ci_u")
  
  #Loop over all scenarios (rows of grid)
  for(i in 1:nrow(ME.grid)){
    #counter
    cat(paste("\n",i,": ME_exp_sd=",ME.grid[i,1],"; ME_conf_sd=",ME.grid[i,2],"\n"))
    
    ## Run sims per scenario
    for(j in 1:sims){
      
      #counter
      cat(j)
      #Reset DF
      df.l <- df
      
      #Amount of ME in SD to add (# SD = ratio of total/(1-ratio of total))
      me.sd.exp  <- ME.grid[i,1]/(1-ME.grid[i,1])
      me.sd.conf <- ME.grid[i,2]/(1-ME.grid[i,2])
      
      #Add ME to exp
      df.l[,expo] <- df[,expo] + me.sd.exp*rnorm(length(df[,expo]))
      
      #Add ME to conf
      df.l[,ME_conf] <- df[,ME_conf] + me.sd.conf*rnorm(length(df[,ME_conf]))
      
      #Run model and extract coef of exp
      surv.ob      <- Surv(time = df.l[,fu_out], event = as.numeric(df.l[,out]))
      mod.1        <- coxph(surv_mod,data=df.l)
      mod.1.coef   <- summary(mod.1)$coefficients[1,1]
      mod.1.ci     <- log(summary(mod.1)$conf.int[1,3:4])
      
      #Collect results in res.mat
      res.mat[i,,j] <- c(mod.1.coef,mod.1.ci)
      
    }}
  
  #Summarize info over array layers
  res.mat.sum <- apply(res.mat,c(1,2),mean)  
  
  #Make final matrix with all info
  tot.mat <- cbind(ME.grid*100,res.mat.sum)
  
  #Add column which is difference between coef and reference
  tot.mat$exp_dif     <- tot.mat$exp_coef-tot.mat$exp_coef[1]
  
  #Add column with relative difference in percentage
  tot.mat$exp_dif_rel <- (tot.mat$exp_coef/tot.mat$exp_coef[1]-1)*100
  
  #Add rounded off column for text in plot
  tot.mat$exp_dif_rel_t <- round(tot.mat$exp_dif_rel)
  
  #Round of coef
  tot.mat$exp_coef <- round(tot.mat$exp_coef,3)
  
  #Results lists
  in.list <- list(out=out,expo=expo,confs=confs,ME_conf=ME_conf,ME_exp_sd=ME_exp_sd,
                  ME_conf_sd=ME_conf_sd,sims=sims,plt=plt,sav=sav,seed=seed)
  
  
  out.list <- list(in.list=in.list,res.mat=res.mat,tot.mat=tot.mat)
  
  if(plt==T){
    print("Making Plot")
    
    gheat.mod <- plot_mod(expo=expo,confs=confs,ME_conf=ME_conf,tot.mat=tot.mat,sims=sims,sav=sav)
    
    out.list <- c(out.list,gheat.mod=gheat.mod)
  }
  
 out.list
}