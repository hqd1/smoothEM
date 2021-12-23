require(doParallel)
require(foreach)
require(dplyr)
require(plotly)
require(RColorBrewer)
setwd("/storage/home/h/hqd1/SpikyDataProject")
source("Functions.R")
registerDoParallel(cores=20)
##### uniform spikes #####
size = seq(2,0,length.out = 11)
arg_pbs = commandArgs(trailingOnly=TRUE)
id = as.numeric(arg_pbs[1]); i = as.numeric(arg_pbs[2]); by = as.numeric(arg_pbs[3]); cat(id, "\t",i, "\t",by, "\n")
noise.sd = 1
nspikes = seq(1,0,length.out = 11)
fit6 = matrix(NA, nrow = length(nspikes)-1, ncol = 8)
colnames(fit6) = c("l1","l2","linf","TPR","FPR","FNR","time","l2.param")
set.seed(6)
spikes.data.main = generate.data(range = c(0,1),by = by,poly.coefs = c(10,90,-250,180,-40),
                                 spikes.loc.dist = "uniform", spikes.loc.param = NA,nspikes = (1/by)*0.22,
                                 spikes.dist = "fixed",spikes.param = size[i]*2*3*noise.sd, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
                                 noise.sd = noise.sd,plot = "no")

##### hompp spikes #####
# mean.poisson = function(x,arg.others){
#   n = arg.others[[1]]
#   dom = arg.others[[2]]
#   if (length(n)!= dim(dom)[1]) stop('incompatible dimensions between arguments.')
#   ret = -1
#   for (i in 1:length(n)){
#     if (dom[i,1] <= x & x <= dom[i,2]) return(n[i]*dnorm(x,mean(dom[i,]), 0.05))
#   }
#   return(ret)
# }
# size = seq(2,0,length.out = 11)
# arg_pbs = commandArgs(trailingOnly=TRUE)
# id = as.numeric(arg_pbs[1]); i = as.numeric(arg_pbs[2]); by = as.numeric(arg_pbs[3]); cat(id, "\t",i, "\t",by, "\n")
# noise.sd = 1
# nspikes = seq(1,0,length.out = 11)
# fit6 = matrix(NA, nrow = length(nspikes)-1, ncol = 8)
# colnames(fit6) = c("l1","l2","linf","TPR","FPR","FNR","time","l2.param")
# set.seed(6)
# arg.others1.2000 = sample(c(150,130,150, 150, 130), size = 5, replace = F)
# arg.others1.1000 = sample(c(100,70,100, 70, 70), size = 5, replace = F)
# arg.others1.500 = sample(c(50,30,50, 30, 30), size = 5, replace = F)
# arg.others1.200 = sample(c(20,10,20, 20, 10), size = 5, replace = F)
# arg.others1 = eval(parse(text = paste0('arg.others1.',1/by)))
# arg.others2 = matrix(c(0.02,0.12,0.2, 0.3,0.45,0.55, 0.65, 0.75, 0.85, 1), ncol = 2, byrow = T)
# 
# #arg.others2 = matrix(c(0.2,0.25,0.4,0.45,0.8,0.85), ncol = 2, byrow = T)
# spikes.data.main = generate.data(range = c(0,1),by = by,poly.coefs = c(10,90,-250,180,-40),
#                                  spikes.loc.dist = "nonhomopp", spikes.loc.param = mean.poisson, arg.others = list(arg.others1,arg.others2),
#                                  spikes.dist = "fixed",spikes.param = size[i]*2*3, #spikes.param: 3*sd*2sides -> when size[i] = 1, the bottom smooth can climb up to the top
#                                  noise.sd = noise.sd,plot = "no")
# cat("alpha = ",length(spikes.data.main$spikes.index)*by, "\n")
##### common part: analysis #####
for(j in (1:(length(nspikes)-1))){
  cat("inner it = ", j, "\n")
  TPR = rep(NA,20)
  FPR = TPR
  FNR = TPR
  l1  = TPR
  l2  = TPR
  linf  = TPR
  tmp = TPR
  l2.param = TPR
  boot = foreach(k = 1:20, .combine = 'rbind',.inorder=FALSE) %dopar% {
    source("Functions.R")
    cat("rep = ", k, "\n")
    ptm = proc.time()
    spikes.data = spikes.data.main
    set.seed(k)
    index.remove = (rbinom(length(spikes.data$spikes.index),size = 1,prob = nspikes[j]) == 0)
    spikes.data$data$y[spikes.data$spikes.index[index.remove]] = spikes.data$data$y[spikes.data$spikes.index[index.remove]]-spikes.data$height
    spikes.data$spikes.index = spikes.data$spikes.index[!index.remove]
    spikes.data$spikes.loc = spikes.data$spikes.loc[!index.remove]
    EM.result = classify(spikes.data$data,range = c(0,1), norder = 5, nbasis = 300,smooth.true = spikes.data$smooth, alpha.prior = 0.35, plot = FALSE, MSOM = TRUE, VIOM = FALSE)
    al.s = length(spikes.data$spikes.index)/length(spikes.data$data$y)
    mu.s = ifelse(length(spikes.data$spikes.index) == 0, 0, size[i]*2*3*noise.sd)
    if(any(is.na(EM.result$EM.result$param))||length(EM.result$EM.result$param) == 0){
      l2.param[k] = NA
    }else{l2.param[k] = sum((EM.result$EM.result$param - c(al.s,mu.s,noise.sd))^2)}
    if(length(EM.result$EM.result$spikes.EM) == 0){
      true.pos = 0
      false.neg = length(spikes.data$spikes.index)
      false.pos = 0
      #fit = fit(spikes.data$data,range = c(0,1), norder = 4, nbasis = 300,smooth.true = spikes.data$smooth,spikes.ind = spikes.data$spikes.index, spikes.flag = EM.result$EM.result$spikes.EM, plot = FALSE)
    } else {
      true.pos = sum(EM.result$EM.result$spikes.EM %in% spikes.data$spikes.index)
      false.neg = sum(!(spikes.data$spikes.index %in% EM.result$EM.result$spikes.EM))
      false.pos = sum(!(EM.result$EM.result$spikes.EM %in% spikes.data$spikes.index))
      #fit = fit(spikes.data$data,range = c(0,1), norder = 4, nbasis = 300,smooth.true = spikes.data$smooth,spikes.ind = spikes.data$spikes.index, spikes.flag = EM.result$EM.result$spikes.EM, lambda = EM.result$lambda, plot = FALSE)
    }
    TPR[k] = ifelse(length(spikes.data$spikes.index) == 0 && length(EM.result$EM.result$spikes.EM) == 0, 1, true.pos/length(spikes.data$spikes.index))
    FPR[k] = false.pos/(length(spikes.data$smooth)-length(spikes.data$spikes.index))
    FNR[k] = ifelse(length(spikes.data$spikes.index) == 0 && false.neg == 0, 0,false.neg/length(spikes.data$spikes.index))
    fit1 = fit(spikes.data$data,range = c(0,1), norder = 5, nbasis = 300,smooth.true = spikes.data$smooth, spikes.flag = EM.result$EM.result$spikes.EM, plot = FALSE)
    l1[k] = fit1$l1
    l2[k] = fit1$l2
    linf[k] = fit1$linf
    tmp[k] = (proc.time() - ptm)[3]
    c(l1[k],l2[k],linf[k],TPR[k],FPR[k],FNR[k],tmp[k],l2.param[k])
  }
  fit6[j,] = apply(boot, 2, mean)
}
save(fit6, file = paste0("/storage/home/hqd1/SpikyDataProject/simulationResult",id,".rdata"))



