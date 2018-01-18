## Plotting function for SMART ME Simulation (sim_smart_25.R)

plot_mod <- function(expo=expo,confs=confs,ME_conf=ME_conf,tot.mat=tot.mat,sims=sims,sav=sav,lims,mdpt){

  #Give names for the axes (FIX FOR BPD)
  if(expo=="bp_s"){
    exp.name = "SBP"
  } else {
    exp.name = "CIMT"  
  }
  
  if(ME_conf=="bp_s"){
    conf.name= "SBP"
  } else {
    if(ME_conf=="abi"){
      conf.name=  "ABI"
    } else {
      conf.name = "DBP"
    }
  }
  
  #Make percentage labels for each cell 
  cell.lab <- sprintf("%s%d%%", ifelse(tot.mat$exp_dif_rel_t>0,"",""), tot.mat$exp_dif_rel_t)
  
  #Put reference effect in first cell (REMOVED)
  #cell.lab[1] <- paste0(cell.lab[1],"\n(",tot.mat$exp_coef[1],")")
  
  #Make heatmap with ggplot
  gheat.mod <- ggplot(tot.mat, aes_string(x="ME_con", y="ME_exp", fill="exp_dif_rel")) +
                geom_tile(color="white") +
                geom_text(aes(label=cell.lab),parse=F, fontface="bold", size=3.5) +
                scale_fill_gradient2(low="#D55E00",mid="white",high = "#56B4E9",midpoint=mdpt,
                                     limit=lims,guide=FALSE) +
                coord_equal()+
                labs(y=paste("Percentage of total variance of",exp.name,"due to measurement error"), 
                     x=paste("Percentage of total variance of",conf.name,"due to measurement error"))+
                scale_y_continuous(breaks=unique(tot.mat[,1]))+
                scale_x_continuous(breaks=unique(tot.mat[,1]))+
                theme(panel.background = element_rect(fill='white', colour='grey'),
                      plot.title=element_text(hjust=0),
                      axis.ticks=element_blank(),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14),
                      legend.title=element_text(size=12),
                      legend.text=element_text(size=10))
  
  #Save plot
  if(sav==T){
    sav.name <- paste0("Analysis-R\\R_Plots\\","heat-",expo,"_",ME_conf,"-",
                       paste(confs[-which(confs==ME_conf)],collapse=""),
                       "_s",sims)
    sav.nam.pdf <- paste0(sav.name,".pdf")
    sav.nam.png <- paste0(sav.name,".png")
    
    ggsave(gheat.mod,filename=sav.nam.pdf,width=7,height=7)
    ggsave(gheat.mod,filename=sav.nam.png,width=7,height=7)
  }
  
gheat.mod  
}