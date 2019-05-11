library(momentuHMM)
set.seed(1)
simDat <- simData(nbStates=2,obsPerAnimal=1000,dist=list(step="weibull",angle="wrpcauchy"),Par=list(step=c(2,2,0.1,2),angle=c(0,0,0.01,0.8)),beta=matrix(c(-2.5,-1.5),1,2),lambda=0.5,states=TRUE)

d <- simDat[which(!(simDat$time %in% seq(2,999))),]
names(d)[8:10] <- c("smaj","smin","eor")
names(d)[3:4] <- c("lon","lat")
names(d)[2] <- "date"
names(d)[1] <- "id"
d$lc <- "3"

d <- d[,c("id","date","lc","lon","lat")]
d$date <- as.POSIXct(d$date*3600,origin="2019-05-10")

# this test crashes
fit1 <- fit_ssm(d, time.step = 1, nbStates=1)
