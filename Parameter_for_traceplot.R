##important parameters for which the trace.plot might be needed
trace_plot_par_generator<-function(data)
  {
    imp_parameters<-c("sigma","sigma2","lam0","psi","psi.sex","Nsuper","D.adj")
    if("gelman_estimate"%in%names(data)==T)
    {
    data<-subset(data,is.na(match(data$parameter,imp_parameters))==F)
    traceplot_parameter<-data$paramenter[which(data$gelman_estimate>1.1)]
    } else if("mcmcse"%in%names(data)==T)
      {
    data<-subset(data,is.na(match(data$output_parameter,imp_parameters))==F)
    data$estimate_conv<-(data$mcmcse/data$mean)*100
    traceplot_parameter<-data$output_parameter[which(data$estimate_conv >1 )]
    } else
    traceplot_parameter =c("warning; Gelman and estimate not present or labels are diferent")
    return(traceplot_parameter)
  }

