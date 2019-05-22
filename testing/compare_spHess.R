library(tidyverse)
library(crawl)
devtools::load_all(".")

data(ellie)
ellie <- ellie %>% select(-smaj,-smin,-eor) %>% prefilter(ellie)

fit_crwHMM <- 
  sfilter(
    x = ellie,
    time.step = 24,
    fit.to.subset = FALSE,
    optim="nlminb",
    verbose = FALSE
  )
fit_crwHMM$tmb$env$spHess(random=TRUE) # has mostly '.' in off-diagonal
SparseM::image(fit_crwHMM$tmb$env$spHess(random=TRUE)) #look at covariance structure


# 2 states
fit2_crwHMM <- 
  sfilter(
    x = ellie,
    time.step = 24,
    fit.to.subset = FALSE,
    optim="nlminb",
    verbose = 2,
    nbStates=2
  )
fit2_crwHMM$tmb$env$spHess(random=TRUE)  # has 0's instead of '.' in off-diagonal
SparseM::image(fit2_crwHMM$tmb$env$spHess(random=TRUE)) #look at covariance structure

