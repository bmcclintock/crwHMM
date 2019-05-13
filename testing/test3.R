library(momentuHMM)
set.seed(1)
n <- 250
simDat <- simData(nbStates=2,obsPerAnimal=n,dist=list(step="weibull",angle="wrpcauchy"),Par=list(step=c(2,2,0.1,2),angle=c(0,0,0.01,0.8)),beta=matrix(c(-2.5,-1.5),1,2),lambda=0.5,states=TRUE)

d <- simDat[which(!(simDat$time %in% seq(2,n-1))),]
names(d)[8:10] <- c("smaj","smin","eor")
names(d)[3:4] <- c("lon","lat")
names(d)[2] <- "date"
names(d)[1] <- "id"
d$lc <- "3"

d <- d[,c("id","date","lc","lon","lat")]
d$date <- as.POSIXct(d$date*3600,origin="2019-05-10")

fit1 <- fit_ssm(d, time.step = 1, nbStates=1)
fit2 <- fit_ssm(d, time.step = 1, nbStates=2,verbose=2,parameters=list(log_beta=rep(fit1$ssm[[1]]$opt$par[["log_beta"]],2),log_sigma=rep(fit1$ssm[[1]]$opt$par[["log_sigma"]],2),l_tau=fit1$ssm[[1]]$opt$par[3:4]))
