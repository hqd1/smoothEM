rm(list=ls())
require(dplyr)
require(tidyr)
library(readr)
require(plotly)
require(RColorBrewer)
require(gridExtra)
library(ggplot2)
library(cowplot)
jBrewColors <- brewer.pal(n = 8, name = "Dark2")
jBrewColors2 <- brewer.pal(n = 9, name = "YlGnBu")
#jBrewCombined <- c(jBrewColors[2],jBrewColors2[5])
jBrewCombined <- c("#afafaf",jBrewColors[2])
theme_set(theme_bw()+
            theme(axis.title = element_text(size = 14, face = "bold",colour = jBrewColors[8])))

source("/Users/huydang/Box/My_Backpack/Archived_Phd/Research/Smooth_spikes_project/RCode/ACI/Functions.R")

##### generate data #####
mean.poisson = function(x,arg.others){
  n = arg.others[[1]]
  dom = arg.others[[2]]
  if (length(n)!= dim(dom)[1]) stop('incompatible dimensions between arguments.')
  ret = -1
  for (i in 1:length(n)){
    if (dom[i,1] <= x & x <= dom[i,2]) return(n[i]*dnorm(x,mean(dom[i,]), 0.05))
  }
  return(ret)
  # n*sin(8*x*pi)+1
  # if(0<x & x<0.25) return(dnorm(x,1/8,n*1/64))
  # if(0.25<x & x<0.5) return(dnorm(x,3/8,n*1/64))
}
size = seq(2,0,length.out = 11)
i = 1
nspikes = seq(1,0,length.out = 11)
fit6 = matrix(NA, nrow = length(nspikes)-1, ncol = 7)
colnames(fit6) = c("l1","l2","linf","TPR","FPR","FNR","time")
set.seed(6)
arg.others1 = c(130,70,30)
#arg.others2 = matrix(c(0.1,0.25,0.4,0.45,0.8,0.95), ncol = 2, byrow = T)
arg.others2 = matrix(c(0.2,0.25,0.4,0.45,0.8,0.85), ncol = 2, byrow = T)
spikes.data.main = generate.data(range = c(0,1),by = 0.002,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
                                 spikes.dist = "normal",spikes.param = c(size[i]*2*3, 0.01), #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = 1,plot = "yes")
##### data and traditional smoothing #####
plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc, smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes))
plot1 = ggplot(plot_data, aes(x = grid, y = y))+
  geom_point(size = 0.1, color = "grey69")+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 0.6)+
  theme(legend.position = "none")+
  labs(x = "Grid", y = "Spiky Data")

spline.basis = create.bspline.basis(rangeval=c(0,1), nbasis=300, norder = 4)  #create basis
Lfd = vec2Lfd(c(0,1), range = c(0,1))

#smooth fit for different lambda values
smooth_fit = vapply(0:5, function(i){
  plot_data %>%
    with(smooth.basis(argvals = grid,y = y, fdParobj = fdPar(spline.basis, Lfd, lambda = 10^(-i)))$fd %>%
           eval.fd(grid,.))
}, plot_data$grid) %>%
  as_tibble() %>%
  setNames(c("lambda = 1e-0", "lambda = 1e-1", "lambda = 1e-2", "lambda = 1e-3", "lambda = 1e-4", "lambda = 1e-5"))

smooth_data = smooth_fit %>%
  tibble(plot_data,.) %>%
  gather(key = "fit_id", value = "fit_value", starts_with("lambda"))

plot2 = ggplot(smooth_data, aes(x = grid, y = fit_value))+
  geom_line(aes(color = fit_id), size = 0.6)+
  geom_point(aes(x = grid, y = y), color = "grey69", size = 0.1)+
  labs(x = "Grid", y = "Spiky Data")+
  theme(legend.title = element_blank())
plot_grid(plot1,plot2, ncol = 2,rel_widths = c(0.41,0.59))
##### Fit the mu*1(spikes) vector #####
mu_fit = vapply(0:5, function(i){
  plot_data %>%
    with(smooth.basis(argvals = grid,y = mu, fdParobj = fdPar(spline.basis, Lfd, lambda = 10^(-i)))$fd %>%
           eval.fd(grid,.))
}, plot_data$grid) %>%
  as_tibble() %>%
  setNames(c("lambda = 1e-0", "lambda = 1e-1", "lambda = 1e-2", "lambda = 1e-3", "lambda = 1e-4", "lambda = 1e-5")) %>%
  tibble(plot_data,.) %>%
  gather(key = "fit_id", value = "mu_fit_value", starts_with("lambda"))

plot3 = ggplot(mu_fit, aes(x = grid, y = mu_fit_value))+
  geom_line(aes(color = fit_id, linetype = fit_id), size = 0.6)+
  geom_point(aes(x = grid, y = mu), color = "grey69", size = 0.1)+
  labs(x = "Grid", y = "Spikes")+
  theme(legend.title = element_blank())

# Fit the true smooth
smooth_fit = vapply(0:5, function(i){
  plot_data %>%
    with(smooth.basis(argvals = grid,y = y - mu, fdParobj = fdPar(spline.basis, Lfd, lambda = 10^(-i)))$fd %>%
           eval.fd(grid,.))
}, plot_data$grid) %>%
  as_tibble() %>%
  setNames(c("lambda = 1e-0", "lambda = 1e-1", "lambda = 1e-2", "lambda = 1e-3", "lambda = 1e-4", "lambda = 1e-5")) %>%
  tibble(plot_data,.) %>%
  gather(key = "fit_id", value = "smooth_fit_value", starts_with("lambda"))

plot4 = ggplot(smooth_fit, aes(x = grid, y = smooth_fit_value))+
  geom_line(aes(color = fit_id, linetype = fit_id), size = 0.6)+
  geom_point(aes(x = grid, y = y - mu), color = "grey69", size = 0.1)+
  labs(x = "Grid", y = "Noisy smooth")+
  theme(legend.position = "none")

plot_grid(plot4,plot3, ncol = 2, rel_widths = c(0.41,0.59))

##### Plot uniform spikes case for simulation 1 #####
size = seq(2,0,length.out = 11)
i = 6
nspikes = seq(1,0,length.out = 11)
fit6 = matrix(NA, nrow = length(nspikes)-1, ncol = 7)
colnames(fit6) = c("l1","l2","linf","TPR","FPR","FNR","time")
noise.sd = 1
set.seed(6)
spikes.data.main = generate.data(range = c(0,1),by = 0.0005,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "uniform", spikes.loc.param = NA,nspikes = 400,
                                 spikes.dist = "fixed",spikes.param = size[5]*2*3*noise.sd, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = noise.sd,plot = "no")
plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc,smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes)) %>%
  mutate(spikes = ifelse(spikes == TRUE, as.logical(rbinom(length(spikes),1,prob = nspikes[5])), spikes), #change some spikes to nonspikes
         y = ifelse(spikes == FALSE,y - mu, y)) #adjust y of these new nonspikes
         
plot5 = ggplot(plot_data, aes(x = grid, y = y))+
  geom_point(size = 0.1, color = "grey69")+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 0.6)+
  theme(legend.position = "none")+
  labs(x = "Grid", y = "Spiky Data")
#change i and j (size and nspikes)
i = 5;j = 3; by = 0.002
i = 1; j = 8; by = 0.002
set.seed(6)
spikes.data.main = generate.data(range = c(0,1),by = by,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "uniform", spikes.loc.param = NA,nspikes = 0.22/by,
                                 spikes.dist = "fixed",spikes.param = size[i]*2*3*noise.sd, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = noise.sd,plot = "no")
plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc,smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes)) %>%
  mutate(spikes = ifelse(spikes == TRUE, as.logical(rbinom(length(spikes),1,prob = nspikes[j])), spikes), #change some spikes to nonspikes
         y = ifelse(spikes == FALSE,y - mu, y)) #adjust y of these new nonspikes
mean(plot_data$spikes)
assign(paste0("plot",i*10), ggplot(plot_data, aes(x = grid, y = y))+
  geom_point(size = 0.1, col = "grey69")+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 0.6)+
  theme(legend.position = "none")+
  labs(x = "Grid", y = "Spiky Data"))

plot_grid(plot50,plot10, ncol = 2)

sim_data = plot_data %>% 
  dplyr::select(y,grid)
EM.result = classify(sim_data, range = c(0,1), nbasis = 300, norder = 5, Lcoef = c(0,0,1), smooth.true = NA, MSOM = TRUE, VIOM = FALSE, alpha.prior = 0.35)
smooth.data = fit(sim_data, range = c(0,1),nbasis = 300, norder = 5, Lcoef = c(0,0,1), spikes.flag = EM.result$EM.result$spikes.EM, plot = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
sim_data_smooth = sim_data %>%
  mutate(spikes = grid %in% grid[spikes_flag],
         smooth = as.vector(eval.fd(grid,smooth.data$fd)),
         true_smooth = plot_data$smooth,
         .keep = "all")%>%
  gather(key = "smooth_id", value = "smooth_value", smooth:true_smooth)

plot1 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, fill = spikes), size = 0.8, shape = 21, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1)+
  scale_fill_manual(labels = c("flagged non-spike", "flagged spikes" ), values=jBrewCombined) +
  scale_colour_manual(labels = c("true smooth", "fitted curve"),values=c(jBrewColors[5], "#000099")) +
  scale_linetype_manual(labels = c("true smooth", "fitted curve"), values = c(1,3))+
  labs(x = "Grid", y = "Spiky Data")+
  theme(legend.position = "none")

plot2 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, fill = spikes), size = 0.8, shape = 21, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1)+
  scale_fill_manual(labels = c("flagged non-spike", "flagged spikes" ), values=jBrewCombined) +
  scale_colour_manual(labels = c("true smooth", "fitted curve"),values=c(jBrewColors[5], "#000099")) +
  scale_linetype_manual(labels = c("true smooth", "fitted curve"), values = c(1,3))+
  labs(x = "Grid", y = "Spiky Data")+
  theme(legend.title = element_blank())
plot_grid(plot1,plot2, ncol = 2,rel_widths = c(0.41,0.59)) #width 850 height 400




##### Plot non hmpp spikes case for simulation 2 #####
mean.poisson = function(x,arg.others){
  n = arg.others[[1]]
  dom = arg.others[[2]]
  if (length(n)!= dim(dom)[1]) stop('incompatible dimensions between arguments.')
  ret = -1
  for (i in 1:length(n)){
    if (dom[i,1] <= x & x <= dom[i,2]) return(n[i]*dnorm(x,mean(dom[i,]), 0.05))
  }
  return(ret)
}
size = seq(2,0,length.out = 11)
nspikes = seq(1,0,length.out = 11)
i = 1
fit6 = matrix(NA, nrow = length(nspikes)-1, ncol = 7)
colnames(fit6) = c("l1","l2","linf","TPR","FPR","FNR","time")

# arg.others1 = c(150,150,150)
# arg.others2 = matrix(c(0.05,0.25,0.4, 0.6,0.8,0.9), ncol = 2, byrow = T)

arg.others1 = sample(c(70,30,50, 50, 60), size = 5, replace = F)
arg.others2 = matrix(c(0.02,0.12,0.2, 0.3,0.45,0.55, 0.65, 0.75, 0.85, 1), ncol = 2, byrow = T)

i = 5; j = 1
i = 1; j = 6
set.seed(6)
spikes.data.main = generate.data(range = c(0,1),by = 0.001,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
                                 spikes.dist = "normal",spikes.param = c(size[i]*2*3, 0.01), #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = 1,plot = "yes")

plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc,smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes)) %>%
  mutate(spikes = ifelse(spikes == TRUE, as.logical(rbinom(length(spikes),1,prob = nspikes[j])), spikes), #change some spikes to nonspikes
         y = ifelse(spikes == FALSE,y - mu, y)) #adjust y of these new nonspikes
mean(plot_data$spikes)
assign(paste0("plot",i*10), ggplot(plot_data, aes(x = grid, y = y))+
  geom_point(size = 0.1, color = "grey69")+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 0.6)+
  theme(legend.position = "none")+
  labs(x = "Grid", y = "Spiky Data"))

sim_data = plot_data %>% 
  dplyr::select(y,grid)
EM.result = classify(sim_data, range = c(0,1), nbasis = 300, norder = 5, Lcoef = c(0,0,1), smooth.true = NA, MSOM = TRUE, VIOM = FALSE)
smooth.data = fit(sim_data, range = c(0,1),nbasis = 300, norder = 5, Lcoef = c(0,0,1), spikes.flag = EM.result$EM.result$spikes.EM, plot = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
sim_data_smooth = sim_data %>%
  mutate(spikes = grid %in% grid[spikes_flag],
         smooth = as.vector(eval.fd(grid,smooth.data$fd)),
         true_smooth = plot_data$smooth,
         .keep = "all")%>%
  gather(key = "smooth_id", value = "smooth_value", smooth:true_smooth)

plot1 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, fill = spikes), size = 0.8, shape = 21, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1)+
  scale_fill_manual(labels = c("flagged non-spike", "flagged spikes" ), values=jBrewCombined) +
  scale_colour_manual(labels = c("true smooth", "fitted curve"),values=c(jBrewColors[5], "#000099")) +
  scale_linetype_manual(labels = c("true smooth", "fitted curve"), values = c(1,3))+
  labs(x = "Grid", y = "Spiky Data")+
  theme(legend.position = "none")

plot2 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, fill = spikes), size = 0.8, shape = 21, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1)+
  scale_fill_manual(labels = c("flagged non-spike", "flagged spikes" ), values=jBrewCombined) +
  scale_colour_manual(labels = c("true smooth", "fitted curve"),values=c(jBrewColors[5], "#000099")) +
  scale_linetype_manual(labels = c("true smooth", "fitted curve"), values = c(1,3))+
  labs(x = "Grid", y = "Spiky Data")+
  theme(legend.title = element_blank())

plot_grid(plot50,plot10, ncol = 2) # 800 by 350
plot_grid(plot1,plot2, ncol = 2,rel_widths = c(0.41,0.59)) #width 850 height 400




##### Application Plots: heatwave #####
heatwave = read_csv("//Users/huydang/Box/My_Backpack/Archived_Phd/Research/Smooth_spikes_project/Data/heat-wave-index-usa.csv") %>%
  rename(hwi = `Heat Wave Index (NOAA)`) %>%
  mutate(y = hwi, grid = Year, .keep = "none") %>%
  filter(grid >= 1910)
ggplot(heatwave, aes(x = grid, y = y))+
  geom_point(col = "grey69")+
  labs(x = "time", y = "heatwave index")

set.seed(6)
EM.result = classify(heatwave, range = c(1910,2015), norder = 4, Lcoef = c(0,1),lambda = c(50,40,30,20,10), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
smooth.data = fit(heatwave, range = c(1910,2015), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = c(50,40,30,20,10), plot = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
heatwave_smooth = heatwave %>%
  mutate(smooth = as.vector(eval.fd(grid,smooth.data$fd)),
         spikes = grid %in% grid[spikes_flag],
         .keep = "all")
library(ggforce)
plot7 = ggplot(heatwave_smooth, aes(x = grid, y = y, col = spikes))+
  scale_color_manual(values=jBrewCombined)+
  geom_point()+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
  labs(y = "heatwave index")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  facet_zoom(ylim = c(0,12), zoom.size = 0.5)

#plot7 + facet_zoom(y = y >= 0 & y <= 10)
#area with high summer temperature
summer_temp = read_csv("/Users/huydang/Box/My_Backpack/Archived_Phd/Research/Smooth_spikes_project/Data/high-summer-temp-usa.csv") %>%
  rename(Pct_area_daily_high = `Daily Highs (% of land area)`) %>%
  mutate(y= Pct_area_daily_high, grid = Year, .keep = "none")
set.seed(6)
EM.result = classify(summer_temp, range = c(1910,2015), norder = 4, Lcoef = c(0,1),lambda = c(50,40,30,20,10), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
smooth.data = fit(summer_temp, range = c(1910,2015), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = c(50,40,30,20,10), plot = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
temp_smooth = summer_temp %>%
  mutate(smooth = as.vector(eval.fd(grid,smooth.data$fd)),
         spikes = grid %in% grid[spikes_flag],
         .keep = "all")
plot8 = ggplot(temp_smooth, aes(x = grid, y = y, col = spikes))+
  scale_color_manual(values=jBrewCombined)+
  geom_point()+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
  labs(y = "% area")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  facet_zoom(ylim = c(0,12), zoom.size = 0.5)
plot9 = plot_grid(plot7,plot8, nrow = 2)
grid.arrange(
  arrangeGrob(plot9,  
              bottom=grid::textGrob(label= "time",
                                    gp= grid::gpar(fontsize=14, fontface="bold", col=jBrewColors[8]))))

##### Application Plots: power consumption #####
##### Format of each data file: 3 columns corresponding to  
########## Meter ID
########## Five digit code:  
############### Day code: digits 1-3 (day 1 = 1st January 2009)
############### Time code: digits 4-5 (1-48 for each 30 minutes with 1= 00:00:00 – 00:29:59)
########## Electricity consumed during 30 minute interval (in kWh)

# power_1392 = read.table("/Users/huydang/Box/My_Backpack/Archived_Phd/Research/Smooth_spikes_project/Data/38_CER Electricity_Gas/CER Electricity Revised March 2012/File1.txt", header = F) %>%
#   rename(meter_id = V1, time = V2, power_consumption = V3) %>%
#   filter(meter_id == 1392) %>%
#   arrange(time)
power_1976 = read.table("/Users/huydang/Box/My_Backpack/Archived_Phd/Research/Smooth_spikes_project/Data/38_CER Electricity_Gas/CER Electricity Revised March 2012/File1.txt", header = F) %>%
  rename(meter_id = V1, time = V2, power_consumption = V3) %>%
  filter(meter_id == 1976) %>%
  arrange(time)

##### EDA 
ggplot(power_1392, aes(x = time, y = power_consumption))+
  geom_point(col = "grey69")+
  labs(x = "time", y = "power consumption (in kWh)")

# one_day_data = lapply(196:276, function(i){
#   #195 corresponds to Tuesday, July 14, 2009
#   #213 to 243 correspons to the month of August, 2009
#   power_1392 %>%
#     filter(between(time,i*100, (i+1)*100)) %>%
#     mutate(y = power_consumption,grid = seq(0.5,24,by = 0.5),.keep ="none")
# })
one_day_data_jan2010_id1976 = lapply(366:396, function(i){
  #195 corresponds to Tuesday, July 14, 2009
  #366 to 396 correspons to the month of Jan, 2010: coldest month in ireland
  power_1976 %>%
    filter(between(time,i*100, (i+1)*100)) %>%
    mutate(y = power_consumption,grid = seq(0.5,24,by = 0.5),.keep ="none")
})

one_day_data_jul2010_id1976 = lapply(547:577, function(i){
  #195 corresponds to Tuesday, July 14, 2009
  #547 to 577 correspons to the month of july, 2010: hottest month in ireland
  power_1976 %>%
    filter(between(time,i*100, (i+1)*100)) %>%
    mutate(y = power_consumption,grid = seq(0.5,24,by = 0.5),.keep ="none")
})


ggplot(one_day_data_jan2010_id1976[[4]], aes(x = grid, y = y))+
  geom_point(col = "grey69")+
  labs(x = "time", y = "power consumption (in kWh)")
ggplot(one_day_data_jul2010_id1976[[4]], aes(x = grid, y = y))+
  geom_point(col = "grey69")+
  labs(x = "time", y = "power consumption (in kWh)")

##### Run Algorithm 
# Lcoef starts with differential order 0, i.e. 
# Lx(t) = b0(t)x(t) + b1(t)Dx(t) + b2(t)(D^2)x(t) + ... + bm(t)(D^m)x(t)
# then  Lcoef = c(0,0,1) means b0 = b1 = 0, b2 = 1
EM.result = classify(one_day_data_jan2010_id1976[[4]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(50,40,30,20, 10,1, 0.1, 0.01, 0.001,0.000001), sigma.jiggle = 1, smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(one_day_data_jan2010_id1976[[4]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
grid = seq(0.5,24,by = 0.5)
daily_smooth = one_day_data_jan2010_id1976[[4]] %>%
  mutate(smooth = as.vector(eval.fd(grid,smooth.data$fd)),
         spikes = grid %in% grid[spikes_flag],.keep = "all")
ggplot(daily_smooth, aes(x = grid, y = y, color = spikes))+
  scale_color_manual(values=jBrewCombined)+
  geom_point()+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())#+facet_zoom(ylim = c(0,0.3), zoom.size = 0.5)+geom_hline(yintercept=0.5, linetype="dashed", color = "red")
##### Plot Results 

daily_smooth = list()
for (i in (1:2)){
  EM.result = classify(one_day_data[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1),lambda = c(50,40,30,20,10), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
  smooth.data = fit(one_day_data[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = 0.1, plot = FALSE)
  spikes_flag = EM.result$EM.result$spikes.EM
  daily_smooth[[i]] = one_day_data[[i]] %>%
    mutate(smooth = as.vector(eval.fd(grid,smooth.data$fd)),
           spikes = grid %in% grid[spikes_flag],.keep = "all")
  assign(paste0("plot",i), ggplot(daily_smooth[[i]], aes(x = grid, y = y, color = spikes))+
           scale_color_manual(values=jBrewCombined)+
           geom_point()+
           geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
           geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
           theme(legend.position = "none")+
           theme(axis.title.x = element_blank())+
           theme(axis.title.y = element_blank())+
           facet_zoom(ylim = c(0,0.3), zoom.size = 0.5))
}

plot3 = plot_grid(plot1,plot2, nrow = 2)
grid.arrange(
  arrangeGrob(plot3,  
              bottom=grid::textGrob(label= "time",
                                    gp= grid::gpar(fontsize=14, fontface="bold", col=jBrewColors[8])),
              left=grid::textGrob(label= "power consumption (in kWh)", rot=90, 
                                  gp= grid::gpar(fontsize=14, fontface="bold", col=jBrewColors[8]))))
