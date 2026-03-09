# Some code (for arguments and outputs) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet
#' @importFrom stats gaussian
grpnet.default <-
  function(x,
           y,
           group,
           family = c("gaussian", "binomial", "multinomial", "poisson",
                      "negative.binomial", "Gamma", "inverse.gaussian"),
           weights = NULL,
           beta = NULL, # aggiunto 12/03/24
           offset = NULL,
           alpha = 1,
           nlambda = 100,
           lambda.min.ratio = NULL,
           lambda = NULL,
           penalty.factor = NULL,
           penalty = c("LASSO", "MCP", "SCAD"),
           gamma = ifelse(penalty == "MCP", 3, 4),
           theta = 1,
           standardize = TRUE,
           intercept = TRUE,
           thresh = 1e-04,
           maxit = 1e05,
           ...){
    # group elastic net regularized regression (default)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-08-30


    ######***######   INITIAL CHECKS   ######***######

    ### get call
    grpnet.call <- match.call()

    ### check x and y
    x <- as.matrix(x) + 0.0
    nobs <- nrow(x)
    nvars <- ncol(x)
    if (is.null(lambda.min.ratio)) {
      lambda.min.ratio <- if (nobs < nvars) 0.05 else 0.0001
    }
    ny <- if(is.matrix(y)) nrow(y) else length(y)
    if(ny != nobs) stop("Inputs 'x' and 'y' must satisfy:\nnrow(x) == length(y)  or  nrow(x) == nrow(y)")

    ### x names
    xnames <- colnames(x)
    if(is.null(xnames)) xnames <- paste0("x", 1:nvars)

    ### check group
    if(missing(group)) {
      group <- ingroup <- 1:nvars
      gsize <- rep(1L, nvars)
      names(gsize) <- 1:nvars
      ngrps <- nvars
    } else {
      if(length(group) != nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) == length(group)")
      ingroup <- group
      group <- as.factor(group)
      ngrps <- nlevels(group)
      if(ngrps > nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) >= nlevels(group)")
      gsize <- table(group)
      group <- as.integer(group)
    }
    gnames <- names(gsize)

    ### check family
    family <- pmatch(as.character(family[1]), c("gaussian", "binomial", "multinomial", "poisson", "negative.binomial", "Gamma", "inverse.gaussian"))
    if(is.na(family)) stop("'family' not recognized")
    family <- c("gaussian", "binomial", "multinomial", "poisson", "negative.binomial", "Gamma", "inverse.gaussian")[family]
    if(family == "gaussian"){
      family <- gaussian()
    } else if(family == "binomial"){
      dr <- function(y, mu, wt){
        ylogy <- function(y, mu){
          res <- y * log(y / mu)
          res[is.nan(res)] <- 0
          res
        }
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        2 * wt * (ylogy(y, mu) + ylogy(1 - y, 1 - mu))
      }
      family <- list(family = "binomial",
                     linkinv = function(eta) {1 / (1 + exp(-eta))},
                     dev.resids = dr)
    } else if (family == "multinomial"){
      il <- function(eta){
        expeta <- exp(eta - apply(eta, 1, max))
        mu <- expeta / rowSums(expeta)
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        mu
      } # end il
      dr <- function(y, mu, wt) {
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        -2 * wt * rowSums(y * log(mu))
      }
      family <- list(family = "multinomial",
                     linkinv = il,
                     dev.resids = dr)
    } else if(family == "poisson"){
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        r <- mu * wt
        p <- which(y > 0)
        if(is.matrix(mu) && ncol(mu) > 1L){
          r[p,] <- (wt * (y * log(y/mu) - (y - mu)))[p,]
        } else {
          r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
        }
        2 * r
      }
      family <- list(family = "poisson",
                     linkinv = il,
                     dev.resids = dr)
    } else if(family == "negative.binomial"){
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        2 * wt * ( y * log(pmax(1, y) / mu) - (y + .Theta) * log((y + .Theta) / (mu + .Theta)) )
      }
      family <- list(family = "negative.binomial",
                     linkinv = il,
                     dev.resids = dr)
    } else if(family == "Gamma"){
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        -2 * wt * (log(y / mu) - (y - mu) / mu)
      }
      family <- list(family = "Gamma",
                     linkinv = il,
                     dev.resids = dr)
    } else if(family == "inverse.gaussian"){
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        wt * ( (y - mu)^2 / (y * mu^2) )
      }
      family <- list(family = "inverse.gaussian",
                     linkinv = il,
                     dev.resids = dr)
    } # end if(family == "gaussian")

    ### check weights
    if(is.null(weights)){
      noweights = TRUE
      weights <- wsqrt <- rep(1.0, nobs)
    } else {
      noweights = FALSE
      weights <- as.numeric(weights)
      if(length(weights) != nobs) stop("Inputs 'y' and 'weights' must satisfy:  length(y) == length(weights)")
      if(any(weights < 0)) stop("Input 'weights' must be non-negative")
      weights <- weights / mean(weights)
      wsqrt <- sqrt(weights)
    }


    ### check beta
    if (is.null(beta)){
      beta = numeric(nvars) # aggiunto 12/04/24
    }

    ### check y
    ylev <- NULL
    if(family$family == "gaussian"){
      y <- as.numeric(y)
    } else if(family$family == "binomial"){
      if(is.factor(y)){
        #if(nlevels(y) != 2L) stop("Input 'y' must be (coercible into) a factor with two levels when family = 'binomial'.")
        ylev <- levels(y)
        y <- ifelse(y == ylev[1], 0.0, 1.0)
      } else if(is.matrix(y)) {
        if(ncol(y) != 2L | any(y < 0)) stop("Input 'y' must be a matrix (# success, # failure) when family = 'binomial'")
        ytotal <- rowSums(y)
        weights <- weights * ytotal
        wsqrt <- sqrt(weights)
        y <- y[,1] / ytotal
      } else {
        y <- as.numeric(y)
        if(any(y < 0) | any(y > 1)) stop("Input 'y' must contain values between 0 and 1 when family = 'binomial'" )
      }
      if(is.null(ylev)) ylev <- c(0, 1)
    } else if(family$family == "multinomial"){
      yorig <- y
      if(is.factor(y)){
        yfac <- y
        yint <- as.integer(y)
        ylev <- levels(y)
        nlev <- nlevels(y)
        if(nlev < 3L) stop("Input 'y' must be (coercible into) a factor with 3 or more levels.")
        y <- matrix(0.0, nrow = nobs, ncol = nlev)
        colnames(y) <- ylev
        for(i in 1:nobs) y[i,yint[i]] <- 1.0
      } else if(is.matrix(y)){
        nlev <- ncol(y)
        ylev <- colnames(y)
        if(is.null(ylev)){
          ylev <- paste0("y", 1:nlev)
          colnames(y) <- ylev
        }
        y <- t(apply(y, 1, as.integer))
        if(min(y) < 0) stop("Input 'y' must be a matrix of counts (i.e, non-negative integers)")
        yrowsum <- rowSums(y)
        y <- y / yrowsum + 0.0
        weights <- weights * yrowsum
        wsqrt <- sqrt(weights)
      } else {
        stop("Invalid data input: 'y' must be a factor or matrix for when family = 'multinomial'.")
      }
    } else if(family$family == "poisson") {
      y <- as.integer(y)
      if(any(y < 0)) stop("Input 'y' must contain non-negative integers when family = 'poisson'.")
    } else if(family$family == "negative.binomial") {
      y <- as.integer(y)
      if(any(y < 0)) stop("Input 'y' must contain non-negative integers when family = 'negative.binomial'.")
    } else if(family$family == "Gamma") {
      y <- as.numeric(y)
      if(any(y <= 0)) stop("Input 'y' must contain non-negative numerics when family = 'Gamma'.")
    } else if(family$family == "inverse.gaussian") {
      y <- as.numeric(y)
      if(any(y <= 0)) stop("Input 'y' must contain non-negative numerics when family = 'inverse.gaussian'.")
    } else {
      warning("need to add checks here...")
    }

    ### check offset
    if(is.null(offset)){
      include.offset <- FALSE
      offset <- rep(0.0, nobs)
      if(family$family == "multinomial") offset <- matrix(offset, nrow = nobs, ncol = nlev)
    } else {
      include.offset <- TRUE
      if(family$family == "multinomial"){
        offset <- as.matrix(offset) + 0.0
        if(nrow(offset) != nobs) stop("Inputs 'x' and 'offset' must satisfy:  nrow(x) == nrow(offset)")
        if(ncol(offset) == 1L){
          offset <- matrix(offset, nrow = nobs, ncol = nlev)
        } else {
          if(ncol(offset) != nlev) stop("Input 'offset' must be a vector of length nobs\n or a matrix of dimension nobs x nlev")
        }
      } else {
        offset <- as.numeric(offset)
        if(length(offset) != nobs) stop("Inputs 'y' and 'offset' must satisfy:  length(y) == length(offset)")
      }
    }

    ### check alpha
    alpha <- as.numeric(alpha[1])
    if(alpha < 0 | alpha > 1) stop("Input 'alpha' must satisfy:  0 <= alpha <= 1")

    ### check lambda and nlambda
    if(is.null(lambda)){
      nolam <- TRUE
      nlambda <- as.integer(nlambda[1])
      if(nlambda <= 0) stop("Input 'nlambda' must be a positive integer")
      lambda <- rep(0.0, nlambda)
    } else {
      nolam <- FALSE
      lambda <- as.numeric(lambda)
      if(length(lambda) > 1L) lambda <- sort(lambda, decreasing = TRUE)
      if(any(lambda < 0)) stop("Input 'lambda' must contain non-negative values")
      nlambda <- length(lambda)
    } # end if(is.null(lambda))

    ### check lambda.min.ratio
    lambda.min.ratio <- as.numeric(lambda.min.ratio[1])
    if(lambda.min.ratio <= 0 | lambda.min.ratio >= 1) stop("Input 'lambda.min.ratio' must satisfy:  0 < lambda.min.ratio < 1")

    ### check penalty.factor
    if(is.null(penalty.factor)){
      penalty.factor <- sqrt(gsize)
    } else {
      penalty.factor <- as.numeric(penalty.factor)
      if(length(penalty.factor) != ngrps) stop("Inputs 'group' and 'penalty.factor' must satisfy:  nlevels(group) == length(penalty.factor)")
      if(any(penalty.factor < 0)) stop("Input 'penalty.factor' must be non-negative")
    }

    ### check penalty
    penalty <- penalty[1]
    if(is.integer(penalty)){
      if(penalty < 1L || penalty > 3L) stop("Invalid 'penalty': must be 1 (LASSO), 2 (MCP), or 3 (SCAD)")
    } else {
      penalty <- toupper(as.character(penalty))
      penalty <- pmatch(penalty, c("LASSO", "MCP", "SCAD"))
      if(is.na(penalty)) stop("Invalid 'penalty': must be 'LASSO', 'MCP', or 'SCAD'")
    }

    ### check gamma
    gamma <- as.numeric(gamma[1])
    if(penalty == "MCP"){
      if(gamma <= 1L) stop("Need gamma > 1 when penalty = 'MCP'")
    } else if(penalty == "SCAD"){
      if(gamma <= 2L) stop("Need gamma > 2 when penalty = 'SCAD'")
    }

    ### check theta
    if(family$family == "negative.binomial"){
      theta <- as.numeric(theta[1])
      if(theta <= 0) stop("Input 'theta' must be positive")
      .Theta <- theta
      env <- new.env(parent = .GlobalEnv)
      assign(".Theta", theta, envir = env)
    }

    ### standardize
    standardize <- as.logical(standardize[1])
    if(!any(standardize == c(TRUE, FALSE))) stop("Input 'standardize' must be TRUE or FALSE")

    ### check intercept
    intercept <- as.logical(intercept[1])
    if(!any(intercept == c(TRUE, FALSE))) stop("Input 'intercept' must be TRUE or FALSE")
    if(intercept) gnames <- c("(Intercept)", gnames)

    ### check thresh
    thresh <- as.numeric(thresh[1])
    if(thresh <= 0) stop("Input 'thresh' must be a positive numeric")

    ### check maxit
    maxit <- as.integer(maxit[1])
    if(maxit <= 0) stop("Input 'maxit' must be a positive integer")


    ######***######   WORK   ######***######

    ### check family
    if(family$family == "gaussian"){

      ## call fortran code
      res <- .Fortran("grpnet_gaussian",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = y,
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

    } else if(family$family == "binomial"){

      ## call fortran code
      res <- .Fortran("grpnet_binomial",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = y,
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

    } else if(family$family == "multinomial"){

      ## call fortran code
      res <- .Fortran("grpnet_multinom",
                      nobs = nobs,
                      nvars = nvars,
                      nresp = nlev,
                      x = x,
                      y = y,
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = matrix(0.0, nrow = nlev, ncol = nlambda),
                      betas = array(0.0, dim = c(nvars, nlev, nlambda)),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

      ## post-process coefs
      rownames(res$ibeta) <- ylev
      colnames(res$ibeta) <- paste0("s", 1:nlambda)
      betas <- vector("list", nlev)
      names(betas) <- ylev
      for(j in 1:nlev) {
        betas[[j]] <- res$betas[,j,]
        rownames(betas[[j]]) <- xnames
        colnames(betas[[j]]) <- paste0("s", 1:nlambda)
      }
      res$betas <- betas

    } else if(family$family == "poisson"){

      ## call fortran code
      res <- .Fortran("grpnet_poisson",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = as.numeric(y),
                      w = weights,
                      beta = beta,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      mus = matrix(0.0, nrow = nobs, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

    } else if(family$family == "negative.binomial"){

      ## call fortran code
      res <- .Fortran("grpnet_negbin",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = as.numeric(y),
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      theta = theta,
                      PACKAGE = "JZILGM")

    } else if(family$family == "Gamma"){

      ## call fortran code
      res <- .Fortran("grpnet_gamma",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = y,
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

    } else if(family$family == "inverse.gaussian"){

      ## call fortran code
      res <- .Fortran("grpnet_invgaus",
                      nobs = nobs,
                      nvars = nvars,
                      x = x,
                      y = y,
                      w = weights,
                      off = offset,
                      ngrps = ngrps,
                      gsize = gsize,
                      pw = penalty.factor,
                      alpha = alpha,
                      nlam = nlambda,
                      lambda = lambda,
                      lmr = lambda.min.ratio,
                      penid = penalty,
                      gamma = gamma,
                      eps = thresh,
                      maxit = maxit,
                      standardize = as.integer(standardize),
                      intercept = as.integer(intercept),
                      ibeta = rep(0.0, nlambda),
                      betas = matrix(0.0, nrow = nvars, ncol = nlambda),
                      iters = rep(0L, nlambda),
                      nzgrps = rep(0L, nlambda),
                      nzcoef = rep(0L, nlambda),
                      edfs = rep(0.0, nlambda),
                      devs = rep(0.0, nlambda),
                      nulldev = 0.0,
                      PACKAGE = "JZILGM")

    } # end if(family$family == "gaussian")


    ######***######   POST-PROCESSING   ######***######

    ## name coefficients
    if(family$family != "multinomial"){
      rownames(res$betas) <- xnames
      colnames(res$betas) <- paste0("s", 1:nlambda)
    }

    ## collect arguments
    args <- list(penalty.factor = penalty.factor,
                 penalty = penalty,
                 gamma = gamma,
                 theta = theta,
                 standardize = standardize,
                 intercept = intercept,
                 thresh = thresh,
                 maxit = maxit)

    ## intercept
    if(intercept){
      ingroup <- c(0, ingroup)
      ngrps <- ngrps + 1L
    }

    ## collect results
    res <- list(call = grpnet.call,
                a0 = res$ibeta,
                beta = res$betas,
                mu = res$mus,
                alpha = alpha,
                lambda = res$lambda,
                family = family,
                dev.ratio = 1 - res$devs / res$nulldev,
                nulldev = res$nulldev,
                df = res$edfs,
                nzgrp = res$nzgrps,
                nzcoef = res$nzcoef,
                ylev = ylev,
                nobs = nobs,
                group = ingroup,
                ngroups = ngrps,
                npasses = res$iters,
                offset = include.offset,
                args = args,
                term.labels = gnames)

    ### return results
    class(res) <- "grpnet"
    return(res)

  } # end grpnet.default
