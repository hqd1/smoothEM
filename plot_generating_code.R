rm(list=ls())
require(dplyr)
require(tidyr)
library(readr)
require(plotly)
require(RColorBrewer)
require(gridExtra)
library(ggplot2);library(ggforce);library(ggnewscale) #for new_color_scale()
library(cowplot)
jBrewColors <- brewer.pal(n = 8, name = "Dark2")
jBrewColors2 <- brewer.pal(n = 9, name = "YlGnBu")
#jBrewCombined <- c(jBrewColors[2],jBrewColors2[5])
jBrewCombined <- c("#afafaf",jBrewColors[2])
theme_set(theme_bw()+
            theme(axis.title = element_text(size = 14, face = "bold",colour = jBrewColors[8])))

source("Functions.R")
##### generate data #####
# the following function generates the mean function for homo poisson case
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

# size = seq(2,0,length.out = 11)
# nspikes = seq(1,0,length.out = 11)
# i = 1
# set.seed(4)
# arg.others1 = c(130,90,30)
# arg.others2 = matrix(c(0.1,0.3,0.45,0.55,0.8,0.95), ncol = 2, byrow = T)
# #arg.others2 = matrix(c(0.2,0.25,0.4,0.45,0.8,0.85), ncol = 2, byrow = T)
# spikes.data.main = generate.data(range = c(0,1),by = 0.002,poly.coefs = c(10,90,-250,180,-40),
#                                  spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
#                                  spikes.dist = "fixed",spikes.param = size[i]*2*3, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
#                                  noise.sd = 1,plot = "yes")
##### data and traditional smoothing #####
size = seq(2,0,length.out = 11)
nspikes = seq(1,0,length.out = 11)
arg.others1 = sample(c(70,30,50, 50, 40), size = 5, replace = F)
arg.others2 = matrix(c(0.02,0.12,0.2, 0.3,0.45,0.55, 0.65, 0.75, 0.85, 0.95), ncol = 2, byrow = T)
i = 1; j = 1
set.seed(6)
spikes.data.main = generate.data(range = c(0,1),by = 0.002,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
                                 spikes.dist = "normal",spikes.param = c(size[i]*2*3, 0.01), #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = 1,plot = "no")

plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc,smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes)) %>%
  mutate(spikes = ifelse(spikes == TRUE, as.logical(rbinom(length(spikes),1,prob = nspikes[j])), spikes), #change some spikes to nonspikes; nspikes == 0 means all spikes is removed. This code is used to automate simulation, but here we are not removing any spikes
         y = ifelse(spikes == FALSE,y - mu, y)) #adjust y of these new nonspikes
mean(plot_data$spikes)
plot1 = ggplot(plot_data, aes(x = grid, y = y))+
         geom_point(size = 0.8, fill = "grey69", color = "grey69", shape = 21)+
         geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 1.2)+
         theme(legend.position = "none")+
         scale_y_continuous(limits = c(-15, 35))+
         labs(x = "x", y = "y")
#here return and run the same thing with i = 1; j = 6
plot1
###### Plot the analysis of nonhomo pp case for the plot_data above ######
sim_data = plot_data %>%
  dplyr::select(y,grid)
EM.result = classify(sim_data, range = c(0,1), nbasis = 300, norder = 5, Lcoef = c(0,1), smooth.true = NA, MSOM = TRUE, VIOM = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(sim_data, range = c(0,1),nbasis = 300, norder = 5, Lcoef = c(0,1), spikes.flag = spikes_flag, plot = FALSE)
sim_data_smooth = sim_data %>%
  mutate(spikes = grid %in% grid[spikes_flag],
         smooth = as.vector(smooth.data),
         true_smooth = as.vector(plot_data$smooth),
         .keep = "all")%>%
  gather(key = "smooth_id", value = "smooth_value", smooth:true_smooth)

sim_data_smooth$smooth_id <- factor(sim_data_smooth$smooth_id, c("smooth", "true_smooth"))
plot10 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, fill = spikes, shape = spikes), size = 1.5, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1.2)+
  scale_y_continuous(limits = c(-15, 35))+
  scale_fill_manual(labels = c("flagged non-spikes", "flagged spikes" ), values=jBrewCombined) +
  #scale_fill_manual(labels = c("flagged non-spikes", "flagged spikes" ), values= c("gray", "light red")) +
  scale_shape_manual(labels = c("flagged non-spikes", "flagged spikes" ),values = c(21, 24))+
  scale_colour_manual(labels = c("fitted smooth","true smooth"), values=c(true_smooth = jBrewColors[5], smooth = "#000099")) +
  scale_linetype_manual(labels = c("fitted smooth","true smooth"), values = c(true_smooth = 1, smooth = 4))+
  aes(group=rev(smooth_id))+
  labs(x = "x", y = "y")+
  theme(legend.position="none")

plot_grid(plot1,plot10, ncol = 2,rel_widths = c(0.41,0.59)) #width 850 height 400


#####Figure 1: smooth fit for different lambda values######
spline.basis = create.bspline.basis(rangeval=c(0,1), nbasis=300, norder = 4)  #create basis
Lfd = vec2Lfd(c(0,1), range = c(0,1))

smooth_fit = vapply(0:5, function(i){
  plot_data %>%
    with(smooth.basis(argvals = grid,y = y, fdParobj = fdPar(spline.basis, Lfd, lambda = 10^(-i)))$fd %>%
           eval.fd(grid,.))
}, plot_data$grid) %>%
  as_tibble() %>%
  setNames(c("lambda = 1e-0", "lambda = 1e-1", "lambda = 1e-2", "lambda = 1e-3", "lambda = 1e-4", "lambda = 1e-5"))

smooth_data = smooth_fit %>%
  cbind(plot_data,.) %>%
  mutate(!!paste0('true f') := smooth,.keep = "unused") %>%
  gather(key = "fit_id", value = "fit_value", starts_with("lambda") | starts_with("true"))

plot100 = ggplot(smooth_data, aes(x = grid, y = fit_value))+
  geom_point(aes(x = grid, y = y), fill = "grey69", color = "grey69", shape = 21, size = 0.8)+
  geom_line(aes(color = fit_id, linetype = fit_id), size = 1)+
  scale_color_manual(values=c(jBrewColors[1:6],"black"))+
  scale_y_continuous(limits = c(-15, 35))+
  labs(x = "x", y = "y")+
  theme(legend.title = element_blank(),
        legend.position = c(0.25, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"))


#####Fit the mu*1(spikes) vector #####
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
plot3
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

###### Do it again for uniform case ########
size = seq(2,0,length.out = 11)
nspikes = seq(1,0,length.out = 11)
noise.sd = 1
set.seed(12)
i = 1;j = 1
spikes.data.main = generate.data(range = c(0,1),by = 0.002,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "uniform", spikes.loc.param = NA,nspikes = 120,
                                 spikes.dist = "fixed",spikes.param = size[i]*2*3*noise.sd, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = noise.sd,plot = "no")
#####
plot_data = spikes.data.main$data %>%
  mutate(spikes = grid %in% spikes.data.main$spikes.loc,smooth = spikes.data.main$smooth, mu = size[i]*2*3*as.numeric(spikes)) %>%
  mutate(spikes = ifelse(spikes == TRUE, as.logical(rbinom(length(spikes),1,prob = nspikes[j])), spikes), #change some spikes to nonspikes
         y = ifelse(spikes == FALSE,y - mu, y)) #adjust y of these new nonspikes
mean(plot_data$spikes)
plot2 = ggplot(plot_data, aes(x = grid, y = y))+
         geom_point(size = 0.8, fill = "grey69", color = "grey69", shape = 21)+
         geom_line(aes(x = grid, y = smooth), color = jBrewColors[5], size = 1.2)+
         theme(legend.position = "none")+
         scale_y_continuous(limits = c(-15, 35))+
         labs(x = "x", y = "y")

###### Plot the analysis of uniform case for the plot_data above ######
sim_data = plot_data %>%
  dplyr::select(y,grid)
EM.result = classify(sim_data, range = c(0,1), nbasis = 300, norder = 5, Lcoef = c(0,1), smooth.true = NA, MSOM = TRUE, VIOM = FALSE, alpha.prior = 0.35)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(sim_data, range = c(0,1),nbasis = 300, norder = 5, Lcoef = c(0,1), spikes.flag = spikes_flag, plot = FALSE)
sim_data_smooth = sim_data %>%
  mutate(spikes = grid %in% grid[spikes_flag],
         smooth = as.vector(smooth.data),
         true_smooth = as.vector(plot_data$smooth),
         .keep = "all")%>%
  gather(key = "smooth_id", value = "smooth_value", smooth:true_smooth)

sim_data_smooth$smooth_id <- factor(sim_data_smooth$smooth_id, c("smooth", "true_smooth"))
plot20 = ggplot(sim_data_smooth, aes(x = grid, y = smooth_value))+
  geom_point(aes(x = grid, y = y, shape = spikes, fill = spikes), size = 1.5, stroke = 0)+
  geom_line(aes(color = smooth_id, linetype = smooth_id), size = 1.2)+
  scale_y_continuous(limits = c(-15, 35))+
  scale_fill_manual(labels = c("flagged non-spikes", "flagged spikes" ), values=jBrewCombined) +
  scale_shape_manual(labels = c("flagged non-spikes", "flagged spikes" ),values = c(21, 24))+
  scale_colour_manual(labels = c("fitted smooth","true smooth"), values=c(true_smooth = jBrewColors[5], smooth = "#000099")) +
  scale_linetype_manual(labels = c("fitted smooth","true smooth"), values = c(true_smooth = 1, smooth = 4))+
  aes(group=rev(smooth_id))+
  labs(x = "x", y = "y")+
  theme(legend.position="none")
  #guides(fill = guide_legend(override.aes = list(size = 2)))+
  # theme(legend.title = element_blank(),
  #       legend.text = element_text(size=16),
  #       legend.position="bottom")


plot20
#####Figure: smooth fit for different lambda values######
spline.basis = create.bspline.basis(rangeval=c(0,1), nbasis=300, norder = 4)  #create basis
Lfd = vec2Lfd(c(0,1), range = c(0,1))
smooth_fit = vapply(0:5, function(i){
  plot_data %>%
    with(smooth.basis(argvals = grid,y = y, fdParobj = fdPar(spline.basis, Lfd, lambda = 10^(-i)))$fd %>%
           eval.fd(grid,.))
}, plot_data$grid) %>%
  as_tibble() %>%
  setNames(c("lambda = 1e-0", "lambda = 1e-1", "lambda = 1e-2", "lambda = 1e-3", "lambda = 1e-4", "lambda = 1e-5"))

smooth_data = smooth_fit %>%
  tibble(plot_data,.) %>%
  mutate(!!paste0('true f') := smooth,.keep = "unused") %>%
  gather(key = "fit_id", value = "fit_value", starts_with("lambda") | starts_with("true"))

plot200 = ggplot(smooth_data, aes(x = grid, y = fit_value))+
  geom_point(aes(x = grid, y = y), fill = "grey69", color = "grey69", shape = 21, size = 0.8)+
  geom_line(aes(color = fit_id, linetype = fit_id), size = 1)+
  scale_color_manual(values=c(jBrewColors[1:6],"black"))+
  scale_y_continuous(limits = c(-15, 35))+
  labs(x = "x", y = "y")+
  theme(legend.title = element_blank(),
        legend.position = c(0.25, 0.25), 
        legend.background = element_rect(fill = "white", color = "black"))

plot_grid(plot200,plot100, ncol = 2) #999 by 500
# grobs = ggplotGrob(plot20)$grobs
# legend = grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
legend = get_legend(plot20)
plot20 = plot20+theme(legend.position="none")
#pgrid = plot_grid(plot2,plot1, plot20, plot10, ncol = 2,rel_widths = c(0.5,0.5))
pgrid = plot_grid(plot20, plot10, ncol = 2,rel_widths = c(0.5,0.5))
plot_grid(pgrid, legend, nrow = 2, rel_heights = c(1, .1)) #999 by 500

##### Application Plots: heatwave #####
heatwave = read_csv("Data/heat-wave-index-usa.csv") %>%
  rename(hwi = `Heat Wave Index (NOAA)`) %>%
  mutate(y = hwi, grid = Year, .keep = "none") %>%
  filter(grid >= 1910)
ggplot(heatwave, aes(x = grid, y = y))+
  geom_point(col = "grey69")+
  labs(x = "time", y = "heatwave index")

set.seed(6)
EM.result = classify(heatwave, range = c(1910,2015), norder = 5, lambda = c(1000,500,100,50,10,5,1,1e-01,1e-02, 1e-03, 1e-04), Lcoef = c(0,1), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
smooth.data = fit(heatwave, range = c(1910,2015), norder = 5, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
heatwave_smooth = heatwave %>%
  mutate(smooth = as.vector(smooth.data),
         spikes = grid %in% grid[spikes_flag],
         .keep = "all")
plot7 = ggplot(heatwave_smooth)+
  geom_point(aes(x = grid, y = y, color= spikes, fill = spikes, shape= spikes))+
  scale_fill_manual(values=jBrewCombined)+
  scale_color_manual(values=jBrewCombined)+
  scale_shape_manual(values = c(21, 24))+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
  labs(y = "heatwave index")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  facet_zoom(ylim = c(0,12), zoom.size = 0.5)
plot7
#area with high summer temperature
summer_temp = read_csv("Data/high-summer-temp-usa.csv") %>%
  rename(Pct_area_daily_high = `Daily Highs (% of land area)`) %>%
  mutate(y= Pct_area_daily_high, grid = Year, .keep = "none")
EM.result = classify(summer_temp, range = c(1910,2015), norder = 4, Lcoef = c(0,1),lambda = c(50,10,1,1e-01,1e-02, 1e-03, 1e-05, 1e-08), smooth.true = NA, MSOM = TRUE, VIOM = TRUE, beta = 1)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(summer_temp, range = c(1910,2015), norder = 4, Lcoef = c(0,1), spikes.flag = spikes_flag, lambda = EM.result$lambda, plot = FALSE)
temp_smooth = summer_temp %>%
  mutate(smooth = as.vector(smooth.data),
         spikes = grid %in% grid[spikes_flag],
         .keep = "all")
plot8 = ggplot(temp_smooth, aes(x = grid, y = y, col = spikes))+
  geom_point(aes(x = grid, y = y, color= spikes, fill = spikes, shape= spikes))+
  scale_fill_manual(values=jBrewCombined)+
  scale_color_manual(values=jBrewCombined)+
  scale_shape_manual(values = c(21, 24))+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5])+
  labs(y = "% area")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  facet_zoom(ylim = c(0,12), zoom.size = 0.5)
plot8
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
power_1976 = read.table("Data/File1.txt", header = F) %>%
  rename(meter_id = V1, time = V2, power_consumption = V3) %>%
  filter(meter_id == 1976) %>%
  arrange(time)
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
daily_smooth_jan2010_id1976 = vector("list",31)
daily_smooth_jul2010_id1976 = vector("list",31)
for (i in 1:31){
  EM.result = classify(one_day_data_jan2010_id1976[[i]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1, 0.1, 0.01, 0.001,0.000001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE, alpha.prior = 0.45)
  spikes_flag = EM.result$EM.result$spikes.EM
  smooth.data = fit(one_day_data_jan2010_id1976[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
  grid = seq(0.5,24,by = 0.5)
  daily_smooth_jan2010_id1976[[i]] = one_day_data_jan2010_id1976[[i]] %>%
    mutate(smooth = as.vector(smooth.data),
           spikes = grid %in% grid[spikes_flag],.keep = "all")
  EM.result = classify(one_day_data_jul2010_id1976[[i]], range = c(0.5,24),norder = 5, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1, 0.1, 0.01, 0.001,0.000001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE, alpha.prior = 0.45)
  spikes_flag = EM.result$EM.result$spikes.EM
  smooth.data = fit(one_day_data_jul2010_id1976[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
  grid = seq(0.5,24,by = 0.5)
  daily_smooth_jul2010_id1976[[i]] = one_day_data_jul2010_id1976[[i]] %>%
    mutate(smooth = as.vector(smooth.data),
           spikes = grid %in% grid[spikes_flag],.keep = "all")
}
##### Plot Results #####
smooth_jan2020_id1976 = daily_smooth_jan2010_id1976[[1]] %>% mutate(grid = grid, smooth.1 = smooth, .keep = "none")
smooth_jul2020_id1976 = daily_smooth_jul2010_id1976[[1]] %>% mutate(grid = grid, smooth.1 = smooth, .keep = "none")
for (i in 2:31){
  smooth_jan2020_id1976[,paste0("smooth.",i)] = daily_smooth_jan2010_id1976[[i]]$smooth
  smooth_jul2020_id1976[,paste0("smooth.",i)] = daily_smooth_jul2010_id1976[[i]]$smooth
}
#now do the same to id 1977
######id 1977 ######
power_1977 = read.table("/Users/huydang/Library/Mobile Documents/com~apple~CloudDocs/Archived-PhD/Research/Smooth_spikes_project/Data/38_CER Electricity_Gas/CER Electricity Revised March 2012/File1.txt", header = F) %>%
  rename(meter_id = V1, time = V2, power_consumption = V3) %>%
  filter(meter_id == 1977) %>%
  arrange(time)

one_day_data_jan2010_id1977 = lapply(366:396, function(i){
  #195 corresponds to Tuesday, July 14, 2009
  #366 to 396 correspons to the month of Jan, 2010: coldest month in ireland
  power_1977 %>%
    filter(between(time,i*100, (i+1)*100)) %>%
    mutate(y = power_consumption,grid = seq(0.5,24,by = 0.5),.keep ="none")
})

one_day_data_jul2010_id1977 = lapply(547:577, function(i){
  #195 corresponds to Tuesday, July 14, 2009
  #547 to 577 correspons to the month of july, 2010: hottest month in ireland
  power_1977 %>%
    filter(between(time,i*100, (i+1)*100)) %>%
    mutate(y = power_consumption,grid = seq(0.5,24,by = 0.5),.keep ="none")
})

daily_smooth_jan2010_id1977 = vector("list",31)
daily_smooth_jul2010_id1977 = vector("list",31)
for (i in 1:31){
  EM.result = classify(one_day_data_jan2010_id1977[[i]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1, 0.1, 0.01, 0.001,0.000001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE, alpha.prior = 0.45,sigma.jiggle = robust.sig)
  spikes_flag = EM.result$EM.result$spikes.EM
  smooth.data = fit(one_day_data_jan2010_id1977[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
  grid = seq(0.5,24,by = 0.5)
  daily_smooth_jan2010_id1977[[i]] = one_day_data_jan2010_id1977[[i]] %>%
    mutate(smooth = smooth.data,
           spikes = grid %in% grid[spikes_flag],.keep = "all")
  EM.result = classify(one_day_data_jul2010_id1977[[i]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1, 0.1, 0.01, 0.001,0.000001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE, alpha.prior = 0.45,sigma.jiggle = robust.sig)
  spikes_flag = EM.result$EM.result$spikes.EM
  smooth.data = fit(one_day_data_jul2010_id1977[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
  grid = seq(0.5,24,by = 0.5)
  daily_smooth_jul2010_id1977[[i]] = one_day_data_jul2010_id1977[[i]] %>%
    mutate(smooth = smooth.data,
           spikes = grid %in% grid[spikes_flag],.keep = "all")
}
smooth_jan2020_id1977 = daily_result_jan2010_id1977[[1]]$smooth.spikes %>% mutate(grid = grid, smooth.jan.1 = smooth, .keep = "none")
smooth_jul2020_id1977 = daily_result_jul2010_id1977[[1]]$smooth.spikes %>% mutate(grid = grid, smooth.jul.1 = smooth, .keep = "none")
for (i in 2:31){
  smooth_jan2020_id1977[,paste0("smooth.jan.",i)] = daily_result_jan2010_id1977[[i]]$smooth.spikes$smooth
  smooth_jul2020_id1977[,paste0("smooth.jul.",i)] = daily_result_jul2010_id1977[[i]]$smooth.spikes$smooth
}

###### plot all curves in both jan and jul for both sme and household#######
smooth_jan2020_id1976 %>%
  reshape2::melt(id = "grid") %>%
  ggplot(aes(x = grid, y = value, color = variable))+
  geom_line()+ylim(0,2)+theme(legend.position = "none")

smooth_jan2020_id1976[,c("grid","smooth.11")] %>%
  reshape2::melt(id = "grid") %>%
  ggplot(aes(x = grid, y = value, color = variable))+
  geom_line()

Greens = brewer.pal(n = 9, name = "Greens")
Reds = brewer.pal(n = 9, name = "Reds")
temp=smooth_jan2020_id1976 %>%
  full_join(smooth_jul2020_id1976)%>%
  mutate(mean.jan = rowMeans(dplyr::select(smooth_jan2020_id1976, starts_with("smooth.jan."))),
         mean.jul = rowMeans(dplyr::select(smooth_jul2020_id1976, starts_with("smooth.jul.")))) %>%
  reshape2::melt(id = "grid")

plot1976 = ggplot(temp, aes(x = grid, y = value))+
  geom_line(data = dplyr::filter(temp,grepl("smooth",variable)), aes(group = variable, color = variable),show.legend = FALSE)+
  scale_color_manual(name = "", values = c(rep(Greens[2], 31), rep(Reds[2], 31)))+
  new_scale_color()+
  geom_line(data = dplyr::filter(temp,grepl("mean",variable)),aes(group = variable, color = variable))+
  scale_color_manual(name = "",values = c("mean.jan" = Greens[7], "mean.jul" = Reds[7]), labels = c("Mean Jan","Mean Jul"))+ylim(0,2)+
  theme(axis.line = element_line(color = "black"),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(colour="black", size = 10, face = "bold"),
        legend.position = c(.5,0.9),
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold",vjust = -7, hjust = 0.5)) +background_grid()+
  labs(title = "SME",x = "Time",y = "Smooth Power Consumption")

temp=smooth_jan2020_id1977 %>%
  full_join(smooth_jul2020_id1977)%>%
  mutate(mean.jan = rowMeans(dplyr::select(smooth_jan2020_id1977, starts_with("smooth.jan."))),
         mean.jul = rowMeans(dplyr::select(smooth_jul2020_id1977, starts_with("smooth.jul.")))) %>%
  reshape2::melt(id = "grid")
plot1977 = ggplot(temp, aes(x = grid, y = value))+
  geom_line(data = dplyr::filter(temp,grepl("smooth",variable)), aes(group = variable, color = variable),show.legend = FALSE)+
  scale_color_manual(name = "", values = c(rep(Greens[2], 31), rep(Reds[2], 31)))+
  new_scale_color()+
  geom_line(data = dplyr::filter(temp,grepl("mean",variable)),aes(group = variable, color = variable))+
  scale_color_manual(name = "",values = c("mean.jan" = Greens[7], "mean.jul" = Reds[7]), labels = c("Mean Jan","Mean Jul"))+ylim(0,0.5)+
  theme(axis.line = element_line(color = "black"),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(colour="black", size = 10, face = "bold"),
        legend.position = c(.5,0.9),
        legend.direction = "horizontal",
        plot.title = element_text(#colour="black", size = 20,
          face = "bold",vjust = -7, hjust = 0.5)) +background_grid()+
  labs(title = "Household",x = "Time",y = "")

plot_grid(plot1976,plot1977, nrow = 1) #999 by 500

##### plot 1 individual household and sme usage ######
i = 5
EM.result = classify(one_day_data_jan2010_id1977[[i]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1,0.1,0.01,0.001,0.0001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(one_day_data_jan2010_id1977[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
grid = seq(0.5,24,by = 0.5)
daily_smooth = one_day_data_jan2010_id1977[[i]] %>%
  mutate(smooth = smooth.data,
         spikes = grid %in% grid[spikes_flag],.keep = "all")
plot1977 = ggplot(daily_smooth)+
  geom_point(aes(x = grid, y = y,color = spikes, fill = spikes, shape = spikes),show.legend = FALSE)+
  scale_color_manual(values=jBrewCombined)+
  scale_fill_manual(values=jBrewCombined)+
  scale_shape_manual(values=c(21,24))+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5],show.legend = FALSE)+
  geom_hline(yintercept=0.325, linetype="dashed", color = "red")+
  theme(axis.line = element_line(color = "black"),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(colour="black", size = 10, face = "bold"),
        legend.position = c(.5,0.9),
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold",vjust = -7, hjust = 0.5))+background_grid()+
  labs(title = "Household",x = "Time",y = "")

EM.result = classify(one_day_data_jan2010_id1976[[i]], range = c(0.5,24),norder = 4, Lcoef = c(0,1),lambda = c(1000,100,50,40,30,20, 10,1,0.1,0.01,0.001,0.0001), smooth.true = NA, MSOM = TRUE, VIOM = TRUE)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(one_day_data_jan2010_id1976[[i]], range = c(0.5,24), norder = 4, Lcoef = c(0,1), spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
grid = seq(0.5,24,by = 0.5)
daily_smooth = one_day_data_jan2010_id1976[[i]] %>%
  mutate(smooth = smooth.data,
         spikes = grid %in% grid[spikes_flag],.keep = "all")
plot1976 = ggplot(daily_smooth)+
  geom_point(aes(x = grid, y = y,color = spikes, fill = spikes, shape = spikes),show.legend = FALSE)+
  scale_color_manual(values=jBrewCombined)+
  scale_fill_manual(values=jBrewCombined)+
  scale_shape_manual(values=c(21,24))+
  geom_line(aes(x = grid, y = smooth), color = jBrewColors[5],show.legend = FALSE)+
  geom_hline(yintercept=0.325, linetype="dashed", color = "red")+
  theme(axis.line = element_line(color = "black"),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(colour="black", size = 10, face = "bold"),
        legend.position = c(.5,0.9),
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold",vjust = -7, hjust = 0.5)) + background_grid()+
  labs(title = "SME",x = "Time",y = "Power Consumption (in kWh)")
plot_grid(plot1976,plot1977, nrow = 1) #999 by 500

####histogram of param estimates#####
mu_hat.jan2010.id1976 = c();mu_hat.jul2010.id1976 = c();mu_hat.jan2010.id1977 = c();mu_hat.jul2010.id1977 = c()
for (i in 1:31){
  mu_hat.jan2010.id1976[i] = daily_result_jan2010_id1976[[i]]$param.est[2]
  mu_hat.jul2010.id1976[i] = daily_result_jul2010_id1976[[i]]$param.est[2]
  mu_hat.jan2010.id1977[i] = daily_result_jan2010_id1977[[i]]$param.est[2]
  mu_hat.jul2010.id1977[i] = daily_result_jul2010_id1977[[i]]$param.est[2]
}

mu_hat.jan2010.id1976;mu_hat.jul2010.id1976;mu_hat.jan2010.id1977;mu_hat.jul2010.id1977
plot(1:31,mu_hat.jan2010.id1976, type = "l") #jan 4th 2010 is a monday
plot(1:31,mu_hat.jul2010.id1976, type = "l") #jul 1st 2010 is a thursday
plot(1:31,mu_hat.jan2010.id1977, type = "l") #jan 4th 2010 is a monday
plot(1:31,mu_hat.jul2010.id1977, type = "l") #jul 1st 2010 is a thursday

alpha_hat.jan2010.id1976 = c();alpha_hat.jul2010.id1976 = c();alpha_hat.jan2010.id1977 = c();alpha_hat.jul2010.id1977 = c()
for (i in 1:31){
  alpha_hat.jan2010.id1976[i] = daily_result_jan2010_id1976[[i]]$param.est[1]
  alpha_hat.jul2010.id1976[i] = daily_result_jul2010_id1976[[i]]$param.est[1]
  alpha_hat.jan2010.id1977[i] = daily_result_jan2010_id1977[[i]]$param.est[1]
  alpha_hat.jul2010.id1977[i] = daily_result_jul2010_id1977[[i]]$param.est[1]
}

alpha_hat.jan2010.id1976;alpha_hat.jul2010.id1976;alpha_hat.jan2010.id1977;alpha_hat.jul2010.id1977
plot(1:31,alpha_hat.jan2010.id1976, type = "l") #jan 4th 2010 is a monday
plot(1:31,alpha_hat.jul2010.id1976, type = "l") #jul 1st 2010 is a thursday
plot(1:31,alpha_hat.jan2010.id1977, type = "l") #jan 4th 2010 is a monday
plot(1:31,alpha_hat.jul2010.id1977, type = "l") #jul 1st 2010 is a thursday

spikes_time_stamp.jan2010.id1976 = c()
spikes_time_stamp.jul2010.id1976 = c()
spikes_time_stamp.jan2010.id1977 = c()
spikes_time_stamp.jul2010.id1977 = c()

for (i in 1:31){
  spikes_time_stamp.jan2010.id1976 = c(spikes_time_stamp.jan2010.id1976, daily_result_jan2010_id1976[[i]]$smooth.spikes %>%
                                         filter(spikes == T) %>% pull(grid))
  spikes_time_stamp.jul2010.id1976 = c(spikes_time_stamp.jul2010.id1976, daily_result_jul2010_id1976[[i]]$smooth.spikes %>%
                                         filter(spikes == T) %>% pull(grid))
  spikes_time_stamp.jan2010.id1977 = c(spikes_time_stamp.jan2010.id1977, daily_result_jan2010_id1977[[i]]$smooth.spikes %>%
                                         filter(spikes == T) %>% pull(grid))
  spikes_time_stamp.jul2010.id1977 = c(spikes_time_stamp.jul2010.id1977, daily_result_jul2010_id1977[[i]]$smooth.spikes %>%
                                         filter(spikes == T) %>% pull(grid))
}

par(mfrow = c(2,2))
hist(spikes_time_stamp.jan2010.id1976, breaks = 24, xlim = c(0,24),
     xlab = "Time",ylab = "Frequency of spikes", main  = "SME/January", cex.lab = 1.3)
hist(spikes_time_stamp.jan2010.id1977, breaks = 24, xlim = c(0,24),
     xlab = "Time",ylab = "Frequency of spikes", main  = "Household/January", cex.lab = 1.3)
hist(spikes_time_stamp.jul2010.id1976, breaks = 24, xlim = c(0,24),
     xlab = "Time",ylab = "Frequency of spikes", main  = "SME/July", cex.lab = 1.3)
hist(spikes_time_stamp.jul2010.id1977, breaks = 24, xlim = c(0,24),
     xlab = "Time",ylab = "Frequency of spikes", main  = "Household/July", cex.lab = 1.3)

############### benchmark
##### generate data #####
# the following function generates the mean function for homo poisson case
##### hompp spikes #####
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
RWSS_GCV <- function(x, lambda1, criterion1 = 0.01, b, lambda2 = lambda1, criterion2 = 0.05){
  m <- length(x)
  wi <- rep(1, m)
  ss <- smooth.spline(x, lambda = lambda1, all.knots = TRUE, w = wi)
  t.1 <- ss$y
  r_org <- x - t.1
  sigma_MAV <- median(abs(r_org)) / 0.6745
  r <- r_org / sigma_MAV
  flag <- sum(abs(r[r<0]))
  while ( 1 ){
    wi <- ifelse(r < 0, 1, ifelse(1 - (r / b) ^ 2 > 0, 1 - (r / b) ^ 2, 0))
    ss <- smooth.spline(t.1, all.knots = TRUE, w = wi, cv = FALSE)
    t.2 <- ss$y
    r_org <- t.1 - t.2
    sigma_MAV <- median(abs(r_org[r_org < b * sigma_MAV])) / 0.6745
    r <- r_org / sigma_MAV
    if (abs(flag/sum(abs(r[r<0])) - 1) < criterion1){
      break
    }
    flag <- sum(abs(r[r<0]))
    t.1 = t.2
  }
  
  # WSS parts
  # pre-baseline
  b <- t.1
  # initialize the weights
  weights <- rep(1, m)
  # redo the smoothing.spline fitting
  ss_re <- smooth.spline(b, lambda = lambda2, all.knots = TRUE, w = weights)
  t_1 <- ss_re$y
  
  residuals <- log((b - t_1)^2)
  variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
  flag <- sum(abs(b - t_1))
  index = which((b-t_1) < 0)
  
  while(1){
    weights <- 1 / exp(variance$y)
    weights <- weights/sum(weights) # scale
    weights[index] = 1-weights[index]
    
    ss_re <- smooth.spline(t_1, lambda = lambda2, all.knots = TRUE, w = weights, cv=FALSE)
    res2 <- (t_1 - ss_re$y)^2
    residuals <- log(ifelse(res2 == 0, min(res2[res2 != 0]), res2))
    variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
    
    if( abs(sum(abs(t_1 - ss_re$y)) / flag - 1) < criterion2 ){break}
    flag <- sum(abs(t_1 - ss_re$y))
    
    index = which((t_1-ss_re$y) < 0)
    t_1 <- ss_re$y
    
    
    if (length(index)<0.1*length(t_1)) {break}
  }
  return(t_1)
}

size = seq(2,0,length.out = 11)
nspikes = seq(1,0,length.out = 11)

i = 6; j = 6; by = 0.002; 
noise.sd = 1
set.seed(6)
arg.others1.2000 = sample(c(150,130,150, 150, 130), size = 5, replace = F)
arg.others1.1000 = sample(c(100,70,100, 70, 70), size = 5, replace = F)
arg.others1.500 = sample(c(50,30,50, 30, 30), size = 5, replace = F)
arg.others1.200 = sample(c(20,10,20, 20, 10), size = 5, replace = F)
arg.others1 = eval(parse(text = paste0('arg.others1.',1/by)))
arg.others2 = matrix(c(0.02,0.12,0.2, 0.3,0.45,0.55, 0.65, 0.75, 0.85, 1), ncol = 2, byrow = T)

#arg.others2 = matrix(c(0.2,0.25,0.4,0.45,0.8,0.85), ncol = 2, byrow = T)
spikes.data.main = generate.data(range = c(0,1),by = by,smooth = dbeta(seq(0,1, by = by),4,1),#3*sin(seq(0,1, by = by)*3*pi),#
                                 spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
                                 spikes.dist = "fixed",spikes.param = size[i]*2*3, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = noise.sd,plot = "no")
spikes.data = spikes.data.main
index.remove = (rbinom(length(spikes.data$spikes.index),size = 1,prob = nspikes[j]) == 0)
spikes.data$data$y[spikes.data$spikes.index[index.remove]] = spikes.data$data$y[spikes.data$spikes.index[index.remove]]-spikes.data$height
spikes.data$spikes.index = spikes.data$spikes.index[!index.remove]
spikes.data$spikes.loc = spikes.data$spikes.loc[!index.remove]
cat("alpha = ",length(spikes.data$spikes.index)*by, "\n")
# rwss-gcv
lambda1 <- lambda2 <- 1e-5; b = 4
yhat_rwss_gcv <- RWSS_GCV(spikes.data$data$y, lambda1 = lambda1, b=b,lambda2 = lambda2)
# smoothEM
EM.result = classify(spikes.data$data,range = c(0,1), Lcoef = c(0,1), norder = 5, nbasis = 300,smooth.true = spikes.data$smooth, alpha.prior = 0.4, plot = FALSE, MSOM = TRUE, VIOM = FALSE)
spikes_flag = EM.result$EM.result$spikes.EM
smooth.data = fit(spikes.data$data, range = c(0,1),nbasis = 300, norder = 5, Lcoef = c(0,1),spikes.flag = spikes_flag, plot = FALSE)
sim_data_smooth = spikes.data$data %>%
  mutate(spikes = grid %in% grid[spikes.data$spikes.index], #grid %in% grid[spikes_flag],
         smoothEM = as.vector(smooth.data),
         true_smooth = as.vector(spikes.data$smooth),
         .keep = "all")%>%
  gather(key = "smooth_id", value = "smooth_value", smoothEM:true_smooth)
# gam's adaptive spline
b <- mgcv::gam(spikes.data$data$y~s(spikes.data$data$grid,bs="ad", m = 5)) ## adaptive
y_as <- b$fitted.values
# plot(b,shade=TRUE,main=round(cor(fitted(b),spikes.data$smooth),digits=4))
# lines(spikes.data$data$grid,spikes.data$smooth-mean(spikes.data$smooth),col=2)

# plot-------
rwss_dat <- sim_data_smooth %>%
  filter(smooth_id == "smoothEM") %>%
  dplyr::select(y, grid, spikes) %>%
  mutate(smooth_id = "RWSS-GCV", smooth_value = yhat_rwss_gcv) %>%
  dplyr::select(y, grid, spikes, smooth_id, smooth_value)
adaptive_dat <- sim_data_smooth %>%
  filter(smooth_id == "smoothEM") %>%
  dplyr::select(y, grid, spikes) %>%
  mutate(smooth_id = "mgcv-AS", smooth_value = y_as) %>%
  dplyr::select(y, grid, spikes, smooth_id, smooth_value)
benchmark_dat <- rbind(sim_data_smooth, rwss_dat, adaptive_dat) %>%
  tibble() %>%
  mutate(smooth_id = ifelse(smooth_id == "true_smooth", "truth", smooth_id)) %>%
  mutate(smooth_id = factor(smooth_id, levels = c("truth","smoothEM","RWSS-GCV","mgcv-AS"))) %>%
  mutate(spikes = ifelse(spikes == TRUE, "spikes", "non-spikes")) %>%
  arrange(smooth_id, grid)

fill <- c("#bdbdbd",jBrewColors2[5])
shape <- c(21,24)
names(fill) <- names(shape) <- c("non-spikes", "spikes")
col <- c(jBrewColors[5], "#000099", "#FF33CC", "#E69F00")
lt <- c(1, 2, 3, 4)
names(col) <- names(lt) <- levels(benchmark_dat$smooth_id)

plot10 = ggplot(benchmark_dat)+
  geom_point(aes(x = grid, y = y, shape = factor(spikes)), fill = "#bdbdbd", size = 1.5, stroke = 0)+
  geom_line(aes(x = grid, y = smooth_value, color = smooth_id, linetype = smooth_id), size = 1.2)+
  #scale_y_continuous(limits = c(-15, 35))+
  #scale_fill_manual(values=fill) +
  scale_shape_manual(name = "", values = shape)+
  scale_colour_manual(name = "", values=col) +
  scale_linetype_manual(name = "",values = lt)+
  aes(group=rev(smooth_id))+
  labs(x = "x", y = "y")+
  theme(legend.position="bottom",
        legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size = 13))+
  guides(shape = guide_legend(order = 1))


plot10

legend <- get_legend(plot1)
plot10 <- plot10 + theme(legend.position="none")
plot1 <- plot1 + theme(legend.position="none")
pgrid <- plot_grid(plotlist = list(plot10, plot1), ncol = 2)
plot_grid(pgrid, legend, nrow = 2,rel_heights = c(0.9,0.1))
