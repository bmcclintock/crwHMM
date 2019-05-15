library(momentuHMM)
library(sp)

URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/",
              "10255/move.373/Elliptical%20Time-Density%20Model%20%28Wall%",
              "20et%20al.%202014%29%20African%20Elephant%20Dataset%20%",
              "28Source-Save%20the%20Elephants%29.csv")
rawData <- read.csv(url(URL))

# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])

# convert times from factors to POSIX
rawData$time <- as.POSIXct(rawData$time,tz="GMT")

# project to UTM coordinates using package rgdal
library(rgdal)
llcoord <- SpatialPoints(rawData[,3:4], 
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))

# add UTM locations to data frame
rawData$x <- attr(utmcoord,"coords")[,1]
rawData$y <- attr(utmcoord,"coords")[,2]

rawData$ln.sd.x <- rawData$ln.sd.y <- 3.565449 
rawData$error.corr <- 0
crwOut<-crawlWrap(rawData,timeStep="hour",theta=c(1,1,6.855, -0.007),fixPar=c(NA,NA,NA,NA),err.model=list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr))
plot(crwOut,ask=FALSE)

rawDat <- as.data.frame(llcoord)
rawDat$date <- rawData$time
rawDat$id <- rawData$ID
rawDat$lc <- "3"
rawDat$smaj <- rawDat$smin <- 50
rawDat$eor <- 0
rawDat <- rawDat[,c("id","date","lc","lon","lat","smaj","smin","eor")]

t1 <- proc.time()
fG_fit <- foieGras::fit_ssm(rawDat,time.step=1,verbose=2)
e1 <- proc.time() - t1

swim_data <- rawDat[,c("date","lon","lat")]
t2 <- proc.time()
swim_fit <- swim::fitSwim(data=swim_data,ts=1)
e2 <- proc.time() - t2

t3 <- proc.time()
fit <- fit_ssm(rawDat, time.step = 1, verbose=2)
e3 <- proc.time() - t3

t4 <- proc.time()
fit2 <- fit_ssm(rawDat, time.step = 1, nbStates=2,verbose=2)
e4 <- proc.time() - t4
