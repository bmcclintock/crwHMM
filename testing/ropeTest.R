set.seed(1,kind = "Mersenne-Twister",normal.kind = "Inversion")
data(rope)
rope2 <- as.data.frame(rope)
newrope <- NULL
# add fake "resting" locations
for(i in unique(rope2$id)){
  irope <- subset(rope2,id==i)
  resting <- data.frame(id=i,date=seq(tail(irope$date,1)+24*3600,tail(irope$date,1)+7*24*60*60,length=30),lc=2,lon=tail(irope$lon,1)+rnorm(30,0,.2),lat=tail(irope$lat,1)+rnorm(30,0,.2))#,smaj=50,smin=50,eor=90)
  irope <- rbind(irope,resting)
  newrope <- rbind(newrope,irope)
}

fit1 <- fit_ssm(newrope, time.step = 24, nbStates=1)
fit2 <- fit_ssm(newrope, time.step = 24, verbose=2,nbStates=2)
