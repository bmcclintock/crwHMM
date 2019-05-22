##' @title fit the state-space model to \code{prefilter}-ed data
##'
##' @description generates initial values for model parameters and unobserved states;
##' structures data and initial values for C++ \code{TMB} template;
##' fits state-space model; minimises the joint log-likelihood via the selected
##' optimizer (\code{nlminb} or \code{optim}); structures and passes output
##' object to \code{fit_ssm}
##'
##' @details called by \code{fit_ssm}. \code{sfilter} can only fit to an
##' individual track, use \code{fit_ssm} to fit to multiple tracks (see ?fit_ssm).
##'
##' @param x Argos data passed through prefilter()
##' @param time.step the regular time interval, in hours, to predict to.
##' Alternatively, a vector of prediction times, possibly not regular, must be
##' specified as a data.frame with id and POSIXt dates.
##' @param parameters a list of initial values for all model parameters and
##' unobserved states, default is to let sfilter specifiy these. Only play with
##' this if you know what you are doing...
##' @param fit.to.subset fit the SSM to the data subset determined by prefilter
##' (default is TRUE)
##' @param optim numerical optimizer to be used ("nlminb" or "optim")
##' @param verbose report progress during minimization
##' @param inner.control list of control settings for the inner optimization
##' (see ?TMB::MakeADFUN for additional details)
##'
##' @useDynLib crwHMM
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @importFrom stats approx cov sd predict nlminb optim na.omit
##' @importFrom dplyr mutate filter select full_join arrange lag bind_cols
##' @importFrom magrittr "%>%"
##' @importFrom tibble as_tibble
##' @importFrom sf st_crs st_coordinates st_geometry<- st_as_sf st_set_crs
##'
##' @examples
##' data(ellie)
##' pf <- prefilter(ellie, vmax=10, ang=c(15,25), min.dt=120)
##' out <- sfilter(pf, model="rw", time.step=24)
##'
##' @export

sfilter <-
  function(x,
           time.step = 6,
           parameters = NULL,
           fit.to.subset = TRUE,
           optim = c("nlminb", "optim"),
           verbose = FALSE,
           inner.control = NULL) {
    
    st <- proc.time()
    call <- match.call()
    optim <- match.arg(optim)
    
    if(is.null(time.step)) {
      print("\nNo time.step specified, using 6 h as a default time step")
      time.step <- 6
    } else if(length(time.step) > 1 & !is.data.frame(time.step)) {
      stop("\ntime.step must be a data.frame with id's when specifying multiple prediction times")
    } else if(length(time.step) > 1 & is.data.frame(time.step)) {
      if(sum(!names(time.step) %in% c("id","date")) > 0) stop("\n time.step names must be `id` and `date`")
    }
    
    ## drop any records flagged to be ignored, if fit.to.subset is TRUE
    ## add is.data flag (distinquish obs from reg states)
    if(fit.to.subset) xx <- x %>% filter(keep)
    else xx <- x
    
    prj <- st_crs(xx)
    loc <- as.data.frame(st_coordinates(xx))
    names(loc) <- c("x","y")
    st_geometry(xx) <- NULL
    d <- cbind(xx, loc) %>%
      mutate(isd = TRUE)
    
    if (length(time.step) == 1) {
      ## Interpolation times - assume on time.step-multiple of the hour
      tsp <- time.step * 3600
      tms <- (as.numeric(d$date) - as.numeric(d$date[1])) / tsp
      index <- floor(tms)
      ts <-
        data.frame(date = seq(
          trunc(d$date[1], "hour"),
          by = tsp,
          length.out = max(index) + 2
        ))
    } else {
      ts <- time.step %>%
        filter(id == unique(d$id)) %>%
        select(date)
    }
    
    ## merge data and interpolation times
    d.all <- full_join(d, ts, by = "date") %>%
      arrange(date) %>%
      mutate(isd = ifelse(is.na(isd), FALSE, isd)) %>%
      mutate(id = ifelse(is.na(id), na.omit(unique(id))[1], id))
    
    ## calc delta times in hours for observations & interpolation points (states)
    dt <- difftime(d.all$date, lag(d.all$date), units = "hours") %>%
      as.numeric() / 24
    
    #dt[1] <- 0.000001 # - 0 causes numerical issues in CRW model
    
    ## use approx & MA filter to obtain state initial values
    x.init <-
      approx(x = select(d, date, x),
             xout = d.all$date,
             rule = 2)$y
    x.init <-
      stats::filter(x.init, rep(1, 10) / 10) %>% as.numeric()
    x.init[1:4] <- x.init[5]
    x.init[which(is.na(x.init))] <-
      x.init[which(is.na(x.init))[1] - 1]
    
    y.init <-
      approx(x = select(d, date, y),
             xout = d.all$date,
             rule = 2)$y
    y.init <-
      stats::filter(y.init, rep(1, 10) / 10) %>% as.numeric()
    y.init[1:4] <- y.init[5]
    y.init[which(is.na(y.init))] <-
      y.init[which(is.na(y.init))[1] - 1]
    xs <- cbind(x.init, y.init)
    
    mu0 <- c(xs[1,1], xs[1,2])
    if (is.null(parameters)) {
      ## Estimate stochastic innovations
      es <- (xs[-1,] - xs[-nrow(xs),])/dt[-1]
      
      ## Estimate components of variance
      ar_es <- ar(es, order.max=1)
      sigma <- mean(sqrt(diag(var(es))))
      mu <- cbind(xs[,1], xs[,2])
      
      parameters <- list(
        log_beta= -0.75,
        log_sigma = log(pmax(1e-08, sigma)),
        mu=t(mu),
        l_psi = 0,
        l_tau = c(0, 0),
        l_rho_o = 0
      )
    }
    
    ## calculate prop'n of obs that are LS-derived
    d <- d %>% mutate(obs.type = factor(obs.type, levels = c("LS","KF"), labels = c("LS","KF")))
    pls <- table(d$obs.type)["LS"] / nrow(d)
    map <- {
      if (pls == 1) {
        list(
          l_psi = factor(NA),
          l_rho_o = factor(NA)
        )
      } else if (pls == 0) {
        list(
          l_tau = factor(c(NA, NA)),
          l_rho_o = factor(NA)
        )
      } else if (pls > 0 & pls < 1) {
        list(l_rho_o = factor(NA))
      }
    }
    
    ## TMB - data list
    data <- list(
      Y = rbind(d.all$x, d.all$y),
      mu0 = mu0,
      V0 = max(sqrt(diag(var(cbind(d.all$x, d.all$y),na.rm=T)))),
      dt = dt,
      isd = as.integer(d.all$isd),
      obs_mod = ifelse(d.all$obs.type == "LS", 0, 1),
      m = d.all$smin,
      M = d.all$smaj,
      c = d.all$eor,
      K = cbind(d.all$amf_x, d.all$amf_y)
    )
    
    ## TMB - create objective function
    if (is.null(inner.control)) {
      inner.control <- list(smartsearch = TRUE)
    }
    obj <-
      MakeADFun(
        data,
        parameters,
        map = map,
        random = c("mu"),
        DLL = "crwHMM",
        hessian = TRUE,
        silent = !verbose,
        inner.control = inner.control
      )
    #    newtonOption(obj, trace = verbose, smartsearch = TRUE)
    #    obj$env$inner.control$trace <- verbose
    #    obj$env$inner.control$smartsearch <- FALSE
    #    obj$env$inner.control$maxit <- 1
    obj$env$tracemgc <- verbose
    
    ## add par values to trace if verbose = TRUE
    myfn <- function(x) {
      print("pars:")
      print(x)
      obj$fn(x)
    }
    
    ## Minimize objective function
    message("Fitting model...")
    opt <-
      suppressWarnings(switch(optim,
                              nlminb = try(nlminb(obj$par, obj$fn, obj$gr))
                              , #myfn #obj$fn
                              optim = try(do.call(
                                optim,
                                args = list(
                                  par = obj$par,
                                  fn = obj$fn,
                                  gr = obj$gr,
                                  method = "L-BFGS-B"
                                )
                              ))))
    
    ## if error then exit with limited output to aid debugging
    rep <- suppressWarnings(try(sdreport(obj)))
    if (!inherits(opt, "try-error") & !inherits(rep, "try-error")) {
     # return(list(rep=rep, obj=obj, opt=opt) )
      ## Parameters, states and the fitted values
      message("Estmating locations...")
      fxd <- summary(rep, "report")
      tmp <- summary(rep, "random")
      alpha_idx <- rep(c(1:2), nrow(d.all))
      loc <- cbind(tmp[alpha_idx == 1,1], tmp[alpha_idx == 2,1],
                   tmp[alpha_idx == 1,2],tmp[alpha_idx == 2,2])
      colnames(loc) <- c("x", "y", "x.se", "y.se")
      loc <- as_tibble(loc)
      # vel <- cbind(tmp[alpha_idx == 2,1], tmp[alpha_idx == 4,1],
      #              tmp[alpha_idx == 2,2],tmp[alpha_idx == 4,2])
      # colnames(vel) <- c("u", "v", "u.se", "v.se")
      # vel <- as_tibble(vel)
      rdm <- loc %>% #bind_cols(loc, vel) %>%
        mutate(
          id = unique(d.all$id),
          date = d.all$date,
          isd = d.all$isd
        ) %>%
        # select(id, date, x, y, x.se, y.se, u, v, u.se, v.se, isd)
        select(id, date, x, y, x.se, y.se, isd)
      ## Fitted values (estimated locations at observation times)
      fd <- rdm %>%
        filter(isd) %>%
        select(-isd)
      ## Predicted values (estimated locations at regular time intervals, defined by `ts`)
      pd <- rdm %>%
        filter(!isd) %>%
        select(-isd)
      
      if (optim == "nlminb") {
        aic <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]]
      } else if (optim == "optim") {
        aic <- 2 * length(opt[["par"]]) + 2 * opt[["value"]]
      }
      
      out <- list(
        call = call,
        predicted = pd,
        fitted = fd,
        par = fxd,
        data = x,
        inits = parameters,
        pm = "crawl",
        ts = time.step,
        opt = opt,
        tmb = obj,
        rep = rep,
        aic = aic,
        time = proc.time() - st,
        
        report = obj$report()
      )
    } else {
      ## if optimiser fails
      out <- list(
        call = call,
        data = x,
        inits = parameters,
        pm = model,
        ts = time.step,
        tmb = obj,
        errmsg = opt
      )
    }
    
    class(out) <- append("crwHMM", class(out))
    return(out)
  }
