# create dict plot for M=20p, tau_0=1/4, 


library(ggplot2)

# Load the first dataset
load("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.25_delta_4.RData")
fit_smooth_4 = smooth.spline(correspondence_eta_kappa1[,1], correspondence_eta_kappa1[,2])
plot_index_4 = which(abs(fit_smooth_4$y - correspondence_eta_kappa1[,2]) < 0.02)
data_4 = data.frame(Kappa1 = correspondence_eta_kappa1[plot_index_4,1], 
                    Eta2 = correspondence_eta_kappa1[plot_index_4,2],
                    Group = 'delta=4, tau_0=0.25')

# Load the second dataset
load("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.25_delta_2.RData")
fit_smooth_2 = smooth.spline(correspondence_eta_kappa1[,1], correspondence_eta_kappa1[,2])
plot_index_2 = which(abs(fit_smooth_2$y - correspondence_eta_kappa1[,2]) < 0.015)
data_2 = data.frame(Kappa1 = correspondence_eta_kappa1[plot_index_2,1], 
                    Eta2 = correspondence_eta_kappa1[plot_index_2,2],
                    Group = 'delta=2, tau_0=0.25')

# Combine the data
data_combined = rbind(data_4, data_2)
data_combined$Group <- factor(data_combined$Group, 
                              levels = c("delta=4, tau_0=0.25", "delta=2, tau_0=0.25"))

# Plot
ggplot(data_combined, aes(x = Kappa1, y = Eta2, color = Group)) +
  geom_line() +
  scale_color_manual(values = c('blue', 'red'),
                     labels = c(expression( delta == 4   ), 
                                expression( delta == 2  ))) +
  labs(x = expression(kappa[1]), 
       y = expression(eta[M]^2),
       color = "Case") +
  theme_minimal()
