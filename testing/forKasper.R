remotes::install_github('bmcclintock/crwHMM@develop')
library(crwHMM)
library(SparseM)

data(rope)

## fit two (multistate) continuous-time correlated random walk models with identical starting values
## number of fixed parameters estimated = 2 + 2*nbStates + nbStates*(nbStates-1) + nbStates-1

# fit equivalent crawl model in TMB
# very fast when nbStates=1 (4 fixed parameters estimated)
fit1 <- crwHMM::fit_ssm(subset(rope,id=="r11"), time.step = 24, nbStates=1)
fit1$ssm[[1]]$tmb$env$spHess(random=TRUE) # has mostly '.' in off-diagonal
SparseM::image(fit1$ssm[[1]]$tmb$env$spHess(random=TRUE)) #look at covariance structure
table(as.logical(fit1$ssm[[1]]$tmb$env$spHess(random=TRUE))==0) # 2360 non-zero elements out of 156816

# slows to a crawl when nbStates>1 (9 fixed parameters estimated)
fit2 <- crwHMM::fit_ssm(subset(rope,id=="r11"), time.step = 24, nbStates=2, verbose=2)
fit2$ssm[[1]]$tmb$env$spHess(random=TRUE) # has 0's instead of '.' in off-diagonal
SparseM::image(fit2$ssm[[1]]$tmb$env$spHess(random=TRUE)) #look at covariance structure
table(as.logical(fit2$ssm[[1]]$tmb$env$spHess(random=TRUE))==0) # 2634 non-zero elements out of 156816