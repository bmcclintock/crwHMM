library(tidyverse)
library(sf)
library(crawl)
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
    verbose = 2
  )
