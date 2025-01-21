
#### generate plot
library(ggplot2)

generate_one_plot<-function(type_of_plot="MSE",kappa1=1,kappa2=1,kesee=0.9){
  tau_0_seq=seq(0.1,2,0.1)
  design_type="Gaussian"
  beta_true_gen="t3"
  delta=4
  m=20/delta
  load(paste0("TheorySolution_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData"))
  load(paste0("empirical_result_matrix_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData"))
  
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=rep(0,length(tau_0_seq)))
  empir_index=1:3
  if(type_of_plot=="MSE"){
    empir_index=1:3
    plot_data_set$theory_limit=TheorySolution[,1]^2+kappa1^2*(TheorySolution[,3]-1)^2+kappa2^2*TheorySolution[,4]^2
  }else{
    empir_index=4:6
    plot_data_set$theory_limit=kappa1^2*TheorySolution[,3]/(kappa1*sqrt(TheorySolution[,1]^2+kappa1^2*TheorySolution[,3]^2+kappa2^2*TheorySolution[,4]^2))
  }
  plot_data_set[,c(2,3,4)]=empirical_result_matrix[,empir_index]
  plot_data_set$delta=delta
  dataset1=plot_data_set
  
  
  delta=2
  m=20/delta
  load(paste0("TheorySolution_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData"))
  load(paste0("empirical_result_matrix_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData"))
  
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=rep(0,length(tau_0_seq)))
  empir_index=1:3
  if(type_of_plot=="MSE"){
    empir_index=1:3
    plot_data_set$theory_limit=TheorySolution[,1]^2+kappa1^2*(TheorySolution[,3]-1)^2+kappa2^2*TheorySolution[,4]^2
  }else{
    empir_index=4:6
    plot_data_set$theory_limit=kappa1^2*TheorySolution[,3]/(kappa1*sqrt(TheorySolution[,1]^2+kappa1^2*TheorySolution[,3]^2+kappa2^2*TheorySolution[,4]^2))
  }
  plot_data_set[,c(2,3,4)]=empirical_result_matrix[,empir_index]
  plot_data_set$delta=delta
  dataset2=plot_data_set
  merged_data=rbind(dataset1,dataset2)
  if(type_of_plot=="MSE"){
    figg=ggplot(merged_data, aes(x = x, group = delta)) +
      geom_line(aes(y = theory_limit, color = as.factor(delta))) +
      #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
      geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
      scale_color_manual(values = c("red", "blue"), 
                         labels = c(expression(delta==2), expression(delta==4))) +
      scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                         labels = c(expression(delta == 2), expression(delta == 4))) +
      labs(x = expression(tau[0]), y = "Square Error", color = "" ) +
      ylim(0, 4.2) +
      theme(panel.background = element_blank(),  # Blank background
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            panel.border = element_rect(colour = "black", fill = NA) # Add the frame
            ,legend.position="none",plot.title = element_text(hjust = 0.5))+ggtitle(bquote(kappa[1] == .(kappa1))) 
    
  }else{
    figg=ggplot(merged_data, aes(x = x, group = delta)) +
      geom_line(aes(y = theory_limit, color = as.factor(delta))) +
      #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
      geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
      scale_color_manual(values = c("red", "blue"), 
                         labels = c(expression(delta==2), expression(delta==4))) +
      scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                         labels = c(expression(delta == 2), expression(delta == 4))) +
      labs(x = expression(tau[0]), y = "Cosine Similarity", color = "" ) +
      ylim(0, 1) +
      theme(panel.background = element_blank(),  # Blank background
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            panel.border = element_rect(colour = "black", fill = NA ) # Add the frame
            , legend.position = "none",plot.title = element_text(hjust = 0.5))+ggtitle(bquote(kappa[1] == .(kappa1))) 
    
  }
  only_for_legend_plot=ggplot(merged_data, aes(x = x, group = delta)) +
    geom_line(aes(y = theory_limit, color = as.factor(delta))) +
    #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
    geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
    scale_color_manual(values = c("red", "blue"), 
                       labels = c(expression(delta==2), expression(delta==4))) +
    scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                       labels = c(expression(delta == 2), expression(delta == 4))) +
    labs(x = expression(tau[0]), y = "Cosine Similarity", color =  " ",shape= " " ) +
    ylim(0, 10) +
    theme(legend.position = "right",legend.text = element_text(size = 18),
          legend.key.size = unit(2.5, "lines"),panel.background = element_blank())
    
  g_legend <- ggplotGrob(only_for_legend_plot)
  legend <- g_legend$grobs[[which(g_legend$layout$name == "guide-box-right")]]
  return(list(fig=figg,legend=legend))
}
generate_one_plot("MSE",kappa1 = 1.5)$fig
legend=generate_one_plot("MSE",kappa1 = 1.5)$legend



library(gridExtra)
grid.arrange(generate_one_plot("MSE",kappa1 = 0.5)$fig, generate_one_plot("MSE",kappa1 = 1.5)$fig,
             generate_one_plot("Cor",kappa1 = 0.5)$fig, generate_one_plot("Cor",kappa1 = 1.5)$fig, nrow = 2,ncol=2,
             right=legend)

grid.arrange(generate_one_plot("MSE",kappa1 = 1)$fig, generate_one_plot("MSE",kappa1 = 2)$fig,
             generate_one_plot("Cor",kappa1 = 1)$fig, generate_one_plot("Cor",kappa1 = 2)$fig, nrow = 2,ncol=2,
             right=legend)

