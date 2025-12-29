#This is the script that should be run inside the Analysis fof SCR Bayes script to generate the 
#and save it in the result folders. 
#loading the required packages for ease in coding 

trace_plot_ckp<-function(mcmc_chains,trace_plot_parameters,save_dir){
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  require(bayesplot)
  require(rstanarm)
  setwd(save_dir)
  if(dir.exists("trace_plots")==F){
    dir.create("trace_plots")
  }
  setwd("./trace_plots/")
  phrase<-"The Trace plot for: "
  title<-paste(phrase,trace_plot_parameters,sep = " ")
  plot<-mcmc_trace(mcmc_chains,pars = trace_plot_parameters)+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle(title)
  ggsave(paste(trace_plot_parameters,".png",sep = ""),width = 40,height = 40,units = "cm")
}

