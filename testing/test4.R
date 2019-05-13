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

crwOut<-crawlWrap(rawData,timeStep="hour",theta=c(6.855, -0.007),fixPar=c(NA,NA))
plot(crwOut,ask=FALSE)

rawDat <- as.data.frame(llcoord)
rawDat$date <- rawData$time
rawDat$id <- rawData$ID
rawDat$lc <- "3"
rawDat$smaj <- rawDat$smin <- 50
rawDat$eor <- 0
rawDat <- rawDat[,c("id","date","lc","lon","lat","smaj","smin","eor")]

fG_fit1 <- foieGras::fit_ssm(rawDat,model="crw",time.step=1)
fit1 <- fit_ssm(rawDat,model="crw", time.step = 1)


# create momentuHMMData object from crwData object
elephantData <- prepData(data=crwOut, covNames="temp")

# label states
stateNames <- c("encamped","exploratory")
# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")

# initial parameters
Par0_m1 <- list(step=c(100,500,100,200),angle=c(0.3,0.7))

# fit model
m1 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m1, 
             estAngleMean = list(angle=FALSE), stateNames = stateNames)

simDat <- simData(model=m1,obsPerAnimal=1000,lambda=0.5,errorEllipse=list(M=50,m=50,r=0),initialPosition=coordinates(utmcoord)[1,])
                  
crwHMM_data <- simDat[which(!(simDat$time %in% seq(2,999))),c("ID","time","x","y","error_semimajor_axis","error_semiminor_axis","error_ellipse_orientation")]                  
colnames(crwHMM_data) <- c("id","date","lon","lat","smaj","smin","eor")
coordinates(crwHMM_data) <- c("lon","lat")
proj4string(crwHMM_data) <- CRS("+proj=utm +zone=30 ellps=WGS84")
crwHMM_data <- as.data.frame(spTransform(crwHMM_data,CRS("+proj=longlat +datum=WGS84")))
crwHMM_data$date <- as.POSIXct(crwHMM_data$date*3600,origin="2019-05-10")
crwHMM_data$lc <- "3"
crwHMM_data <- crwHMM_data[,c("id","date","lc","lon","lat","smaj","smin","eor")]


fit2 <- fit_ssm(crwHMM_data, time.step = 1, nbStates=1,verbose=2,parameters=list(log_beta=crwOut$crwFits$`Salif Keita`$estPar[2],log_sigma=crwOut$crwFits$`Salif Keita`$estPar[1]))
