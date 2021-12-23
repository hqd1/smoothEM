#####DO NOT DELETE-------
require(dplyr)
require(plotly)
require(RColorBrewer)


size = seq(2,0,length.out = 11)
#x = size[1:10]
x = size[1:10]
y = seq(1,0,length.out = 11)[1:10]*0.2 #this is percentage of spikes
l1 = c();l2 = c();linf = c();TPR = c();FPR = c();FNR = c();time = c();l2_param = c()

tmp = expand.grid(x=x, y = y)
for (i in 31:40){
  load(paste("/Users/huydang/Library/Mobile Documents/com~apple~CloudDocs/Archived-PhD/Research/Smooth_spikes_project/RCode/ACI/dec21-uniform/simulationResult",i,".rdata", sep =""))
  l1 = c(l1, fit6[1:10,1])
  l2 = c(l2, fit6[1:10,2])
  linf = c(linf, fit6[1:10,3])
  TPR = c(TPR, fit6[1:10,4])
  FPR = c(FPR, fit6[1:10,5])
  FNR = c(FNR, fit6[1:10,6])
  time= c(time,fit6[1:10,7])
  l2_param = c(l2_param,fit6[1:10,8])
}
df = tibble(size = tmp[,1], spikes.perc = tmp[,2]) %>%
  arrange(desc(size))%>%
  mutate(l1 = l1, l2 = l2, linf = linf,TPR = TPR, FPR = FPR, FNR = FNR, time = time, l2_param = l2_param)

#l2 contour
summary(df$l2)
quantile(df$l2, probs = seq(0,1,by = 0.1))
p <- plot_ly(df,
             x = ~size, 
             y = ~spikes.perc, 
             #z = ~matrix(FPR, nrow = 10, ncol = 9), 
             z = ~l2,
             type = "contour",
             #colorscale = 'Jet',
             autocontour = F,
             colorbar=list(titlefont = list(size = 25),tickfont=list(size=20), lenmode='fraction', len=1),
             contours = list(
               start = 0.1,
               end = 0.9,
               size = 0.2
             )
)%>%
  colorbar(title = "L2") %>%
  layout(
    xaxis = list(title = TeX("\\mu^*/6\\sigma^*"),titlefont = list(size = 25),tickfont = list(size = 20)), 
    yaxis = list(titlefont = list(size = 25),tickfont = list(size = 20)))%>%
    config(mathjax = "cdn")

p #save: width 800 height 450
#FNR CONTOUR
p <- plot_ly(df,
             x = ~size, 
             y = ~spikes.perc, 
             #z = ~matrix(FPR, nrow = 10, ncol = 9), 
             z = ~FNR,
             type = "contour",
             #colorscale = 'Jet',
             autocontour = F,
             colorbar=list(titlefont = list(size = 25),tickfont=list(size=20), lenmode='fraction', len=1),
             contours = list(
               start = 0,
               end = 0.9,
               size = 0.1
             )
)%>%
  colorbar(title = "FNR")%>%
  layout(
    xaxis = list(title = TeX("\\mu^*/6\\sigma^*"),titlefont = list(size = 25),tickfont = list(size = 20)),
    yaxis = list(titlefont = list(size = 25),tickfont = list(size = 20)))%>%
  config(mathjax = "cdn")

p
#l2.param contour plot
p <- plot_ly(df,
             x = ~size, 
             y = ~spikes.perc, 
             #z = ~matrix(FPR, nrow = 10, ncol = 9), 
             z = ~l2_param,
             type = "contour",
             #colorscale = 'Jet',
             autocontour = F,
             colorbar=list(titlefont = list(size = 25),tickfont=list(size=20), lenmode='fraction', len=1),
             contours = list(
               start = 0,
               end = 3,
               size = 0.1
             )
)%>%
  colorbar(title = "SSE")%>%
  layout(
    xaxis = list(title = TeX("\\mu^*/6\\sigma^*"),titlefont = list(size = 25),tickfont = list(size = 20)),
    yaxis = list(titlefont = list(size = 25),tickfont = list(size = 20)))%>%
  config(mathjax = "cdn")
p

#FPR CONTOUR
p <- plot_ly(df,
             x = ~size, 
             y = ~spikes.perc, 
             #z = ~matrix(FPR, nrow = 10, ncol = 9), 
             z = ~FPR,
             type = "contour",
             #colorscale = 'Jet',
             autocontour = F,
             colorbar=list(titlefont = list(size = 25),tickfont=list(size=20), lenmode='fraction', len=1),
             contours = list(
               start = 0,
               end = 0.9,
               size = 0.1
             )
)%>%
  colorbar(title = "FPR")%>%
  layout(
    xaxis = list(title = TeX("\\mu^*/6\\sigma^*"),titlefont = list(size = 25),tickfont = list(size = 20)),
    yaxis = list(titlefont = list(size = 25),tickfont = list(size = 20)))%>%
  config(mathjax = "cdn")

p

# # 3d plots
# jBrewColors <- brewer.pal(n = 9, name = "PuRd")
# p <- plot_ly(df,x = ~size, y = ~spikes.perc, z = ~l2, 
#              text = ~paste("Size: ", size, 'spikes.perc:', spikes.perc),
#              marker = list( color = ~l2,colorscale = jBrewColors, showscale = TRUE
#                             #,size = ~nspikes,sizeref = 8, sizemode = 'diameter'
#              )) %>%
#   add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'size'),
#                       yaxis = list(title = 'spikes.perc'),
#                       zaxis = list(title = 'l2'))
#   )
# p
# # contour with custom levels
# jBrewColors <- brewer.pal(n = 9, name = "Reds")
# l1 = c()
# l2 = c()
# linf = c()
# TPR = c()
# FPR = c()
# FNR = c()
# time = c()
# x = size[1:10]
# y = round(seq(1,0,length.out = 10)*552) #this is nspikes
# for (i in 10:1){
#   load(paste("/storage/home/h/hqd1/SpikyDataProject/simulationResult",i,".rdata", sep =""))
#   l1 = rbind(l1, rev(fit6[1:10,1]))
#   l2 = rbind(l2, rev(fit6[1:10,2]))
#   linf = rbind(linf, rev(fit6[1:10,3]))
#   TPR = rbind(TPR, rev(fit6[1:10,4]))
#   FPR = rbind(FPR, rev(fit6[1:10,5]))
#   FNR = rbind(FNR, rev(fit6[1:10,6]))
#   time= rbind(time,rev(fit6[1:10,7]))
# }
# contour(x= size[10:1], y = round(seq(0,1,length.out = 10)*552), z = l2, 
#         levels = quantile(as.vector(l2), prob = seq(0.1,0.9, length.out = 9)),
#         col = jBrewColors,
#         main = "l2 contour plot")
# 
# #Extract l2 for certain combinations of spikes.perc and size
# l2.some = df %>% 
#   filter(spikes.perc > 0.16) %>%
#   select(spikes.perc, size, l2)
# 
# l2.mat = matrix(l2.some$l2, nrow = 4, ncol = 10, byrow = FALSE)
# rownames(l2.mat) = c(".276", ".246",".214",".184")
# colnames(l2.mat) = seq(0.2,2, by = 0.2)
# l2.mat
# #Extract time for certain combinations of spikes.perc and size
# time.some = df %>% 
#   dplyr::select(spikes.perc, size, time)
# 
# time.mat = matrix(time.some$time, nrow = 10, ncol = 10, byrow = FALSE)
# # rownames(time.mat) = c(".276", ".246",".214",".184")
# colnames(time.mat) = seq(0.2,2, by = 0.2)
# time.mat
# 
# #END DO NOT DELETE
