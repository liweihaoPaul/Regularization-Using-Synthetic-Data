# plot code for may3


p = 200
n = 500
s = 40
q = 0.1
M=20*p
delta=n/p
signal_strengthscala = sqrt(n)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
fdrcor_dataframe <- data.frame(cor=c(0,0.1,0.2,0.3,0.4), m1=c(0,0.1,0.2,0.3,0.4),
                               m2=c(0,0.1,0.2,0.3,0.4),m3=c(0,0.1,0.2,0.3,0.4))
powercor_dataframe <- data.frame(cor=c(0,0.1,0.2,0.3,0.4), m1=c(0,0.1,0.2,0.3,0.4),
                                 m2=c(0,0.1,0.2,0.3,0.4),m3=c(0,0.1,0.2,0.3,0.4))
for(rho in c(0,0.1,0.2,0.3,0.4)){
  load(paste0('MAY3_p',p,'MAPvary_cor_',rho,'s',s,'sig_sqrt_n.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=15,byrow = T)
  fdr_col_index=c(7,10,13)
  power_col_index=c(8,11,14)
  fdrcor_dataframe[which(rho==c(0,0.1,0.2,0.3,0.4)),c(2,3,4)]=apply(result_matrix[,fdr_col_index],2,mean)
  powercor_dataframe[which(rho==c(0,0.1,0.2,0.3,0.4)),c(2,3,4)]=apply(result_matrix[,power_col_index],2,mean)
}
long_dfpower <- melt(powercor_dataframe, id.vars = 'cor', variable.name = 'Method', value.name = 'Error')
long_dfpower$Method <- rep(c('MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)
long_dffdr <- melt(fdrcor_dataframe, id.vars = 'cor', variable.name = 'Method', value.name = 'Error')
long_dffdr$Method <- rep(c('MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)

plotcorup=ggplot(long_dffdr, aes(x = cor, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  labs(y = "FDR",x="") +
  theme_minimal()+ylim(0.0,0.15)
plotcordown=ggplot(long_dfpower, aes(x = cor, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  labs(x = "Toeplitz correlation",
       y = "Power") +
  theme_minimal()+ylim(0,0.9)


######### now for signal strength part
rho=0.2
sig_sequence=sqrt(n)*c(1,1.25,1.5,1.75,2)
fdrsig_dataframe <- data.frame(Signal=c(1,1.25,1.5,1.75,2), m1=c(0,0.1,0.2,0.3,0.4),
                               m2=c(0,0.1,0.2,0.3,0.4),m3=c(0,0.1,0.2,0.3,0.4))
powersig_dataframe <- data.frame(Signal=c(1,1.25,1.5,1.75,2), m1=c(0,0.1,0.2,0.3,0.4),
                                 m2=c(0,0.1,0.2,0.3,0.4),m3=c(0,0.1,0.2,0.3,0.4))

for(signal_strengthscalr in sig_sequence){
  load(paste0('May3_p',p,'MAPvary_signal_',signal_strengthscalr,'s',s,'cor',rho,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=15,byrow = T)
  
  fdr_col_index=c(7,10,13 )
  power_col_index=c(8,11,14 )
  fdrsig_dataframe[which(signal_strengthscalr== sig_sequence),c(2,3,4)]=apply(result_matrix[,fdr_col_index],2,mean)
  powersig_dataframe[which(signal_strengthscalr == sig_sequence),c(2,3,4)]=apply(result_matrix[,power_col_index],2,mean)
  
}
long_dfpower <- melt(powersig_dataframe, id.vars = 'Signal', variable.name = 'Method', value.name = 'Error')
long_dfpower$Method <-  rep(c('MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)
long_dffdr <- melt(fdrsig_dataframe, id.vars = 'Signal', variable.name = 'Method', value.name = 'Error')
long_dffdr$Method <-  rep(c('MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)

plotsigup=ggplot(long_dffdr, aes(x = Signal, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  labs(y = "FDR",x="") +
  theme_minimal()+ylim(0.0,0.15)+
  scale_x_continuous(breaks=c(1,1.25,1.5,1.75,2))
plotsigdown=ggplot(long_dfpower, aes(x = Signal, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  labs(x = "Signal strength",
       y = "Power") +
  theme_minimal()+ylim(0,0.9)+
  scale_x_continuous(breaks=c(1,1.25,1.5,1.75,2))


grid.arrange(plotcorup, plotsigup, plotcordown,plotsigdown, nrow = 2,ncol=2)



