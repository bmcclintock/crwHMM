library(tidyverse)
library(crawl)
devtools::load_all(".")

data(ellie)
ellie <- ellie %>% select(-smaj,-smin,-eor) %>% prefilter(ellie)

fit_crwHMM <- 
  sfilter(
    x = ellie,
    time.step = 12,
    fit.to.subset = TRUE,
    optim="optim",
    verbose = FALSE
  )

fit_crwHMM$opt$par
plot(y~x, fit_crwHMM$predicted, asp=1, cex=0.25, pch=16)

############
# crawl fit

# same start as crwHMM fit
theta <- c(fit_crwHMM$inits$l_tau, fit_crwHMM$inits$log_sigma, fit_crwHMM$inits$log_beta)

fit_crawl <- 
  crwMLE(
    data = ellie, 
    err.model = list(x=~log(amf_x), y=~log(amf_y)), 
    Time.name = "date", time.scale = "day",
    fixPar = c(NA,1,NA,1,NA,NA),
    theta = theta
    )

### Print table comparison
t(fit_crwHMM$opt$par) %>% as.data.frame() %>% 
rbind(
  .,
  fit_crawl$par[c(6,5,1,3)]
) %>% `rownames<-`(c("crwHMM","crawl"))



data("northernFurSeal")
nfs <- northernFurSeal %>% st_as_sf(coords=c("long","lat")) %>% st_set_crs(4326) %>% 
  st_transform(3832)
nfs <- nfs %>% mutate(id="1") %>% rename(date=GMT, lc=loc_class) %>% 
  select(id, date, lc, geometry)
nfs <- prefilter(nfs)
fit_crwHMM <- 
  sfilter(
    x = nfs,
    time.step = 1,
    fit.to.subset = FALSE,
    optim="nlminb",
    verbose = FALSE
  )


