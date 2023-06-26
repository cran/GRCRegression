ll4P <- function(pars, y, x1, scheme2, link.lambda, weights,
  only.value = FALSE){
  nn <- length(pars)
  NaNs <- list(f = NaN, gr = rep(NaN, nn),
    Hess = matrix(NaN, nrow = nn, ncol = nn))
  if(!all(is.finite(pars))) return(NaNs)

  etas <- c(x1 %*% pars)
  lambdas <- link.lambda$gInv(etas)
  qs <- PoisProb(lwrs = scheme2[y], uprs = scheme2[y + 1] - 1, lambdas)
  fValue <- sum(weights * log(qs)) # feature: NaN propagates to fValue!
  if(!is.finite(fValue)) return(NaNs)
  if(only.value) return(list(f = fValue))
  infs <- rep(Inf, length(lambdas))

  Dq.Dlambda <- dpois(scheme2[y] - 1, lambdas) -
    dpois(scheme2[y + 1] - 1, lambdas)
  Dlambda.Deta <- link.lambda$D.gInv(etas)
  tmp.Dl.Dbetas <- prot(Dq.Dlambda * Dlambda.Deta * weights / qs, infs)
  Dl.Dbetas <- c(tmp.Dl.Dbetas$d %*% x1)

  DDq.DDlambda <- dpois(scheme2[y] - 2, lambdas) -
    dpois(scheme2[y + 1] - 2, lambdas) - Dq.Dlambda
  HessDiag <- (DDq.DDlambda - Dq.Dlambda ^ 2 / qs) * Dlambda.Deta ^ 2 +
    Dq.Dlambda * link.lambda$DD.gInv(etas)
  HessDiag <- prot(HessDiag * weights / qs, infs)$d
  Hess <- makeX1D(x1, x1, HessDiag)
  return(list(f = fValue, gr = Dl.Dbetas, Hess = Hess,
    num.outliers = tmp.Dl.Dbetas$no))
}

mleInternal.P <- function(y, x1, scheme, link.lambda = link.log(),
  weights = rep(1, length(y)), xtol_rel = 1e-8, maxit = 100){
  stopifnot(nrow(x1) == length(y))
  stopifnot(min(diff(scheme)) > 0.5)
  scheme2 <- c(scheme, Inf)
  par.init <- double(ncol(x1))
  fun <- function(pars) ll4P(pars = pars, y = y, x1 = x1, scheme2 = scheme2,
    link.lambda = link.lambda, weights = weights, only.value = TRUE)$f
  GrHess <- function(pars) ll4P(pars = pars, y = y, x1 = x1, scheme2 = scheme2,
    link.lambda = link.lambda, weights = weights, only.value = FALSE)
  maxOut <- glmMaximizer(par.init = par.init, fun = fun, GrHess = GrHess,
    maxit = maxit, xtol_rel = xtol_rel)
  return(list(beta = maxOut$par, log.likelihood = maxOut$val, msg = maxOut$msg,
    num.iter = maxOut$num.iter, gr = maxOut$gr))
}
# tp <- genData.P(beta = c(-3, 2, 2, 4, -3), data.size= 9000, scheme=c(0, 2, 4, 10))
# output <- mleInternal.P(y = tp$y, x1 = tp$x, scheme = c(0, 2, 4, 10))

### GENERATING IMAT OBJECTS  #######################
make.Imat.4P <- function(beta, M, x1, link.lambda, weights){
  # force users to provide weights, to prevent mistake
  m <- nrow(x1)
  ptm <- proc.time()
  Imat <- list(M = M)
  Imat$key <- outer(0:M, c(0:(M + 1), Inf), paste, sep = ",")
  Imat$key[lower.tri(Imat$key, diag = TRUE)] <- NA
  Imat$parDim <- ncol(x1)
  etas <- c(x1 %*% beta)
  lambdas <- link.lambda$gInv(etas)
  infs <- rep(Inf, length.out = m)

  Dlambda.Deta <- link.lambda$D.gInv(etas)
  matColNms <- c(Imat$key)
  Imat$mat <- matrix(0, nrow = Imat$parDim ^ 2, ncol = 2 + M * (M + 5) / 2)
  colnames(Imat$mat) <- sort(matColNms[!is.na(matColNms)])
  noCollect <- 0
  for(rid in 0:M){
    for(cid in c((rid + 1):(M + 1), Inf)){
      # real block: {rid, rid + 1, ..., cid - 1}
      if((rid == 0) & (cid == Inf)) next
      qs <- PoisProb(lwrs = rid, uprs = cid - 1, lambdas = lambdas)
      Dq.Dlambda <- dpois(x = rid - 1, lambda = lambdas) -
        dpois(x = cid - 1, lambda = lambdas)
      tmp.diag <- prot((Dq.Dlambda * Dlambda.Deta) ^ 2 * weights / qs, infs)
      Imat$mat[, paste(rid, cid, sep = ",")] <- c(makeX1D(x1, x1, tmp.diag$d))
      noCollect <- noCollect + tmp.diag$no
    }
  }

  tmp.diagVanilla <- prot(Dlambda.Deta ^ 2 * weights / lambdas, infs)
  Imat$vanillaFisherInfoMat <- makeX1D(x1, x1, tmp.diagVanilla$d)

  temp <- matrix(dpois(x = rep((-1):M, each = m), lambda = lambdas), nrow = m)
  infos <- rowSums((temp[, -ncol(temp), drop = FALSE] -
    temp[, -1, drop = FALSE]) ^ 2 / temp[, -1, drop = FALSE])
  tmp.ResInfo <- prot((1 / lambdas - infos) * Dlambda.Deta ^ 2 * weights, infs)
  Imat$plainMat4lower <- c(makeX1D(x1, x1, tmp.ResInfo$d))
  Imat$num.outliers <- noCollect + tmp.diagVanilla$no + tmp.ResInfo$no
  Imat$time.consumed <- proc.time() - ptm
  return(Imat)
}
# set.seed(777)
# beta <- c(-3, 2, 2, 4, -3, -3, 1)
# tp <- genData.P(beta = beta, data.size= 9000, scheme=c(0, 2, 4, 10))
# link.lambda <- link.log()
# a <- make.Imat.4P(beta = beta, M = 4, x1 = tp$x, link.lambda = link.lambda, weights = rep(1, nrow(tp$x)))


###  GRCglm function  #######################
GRCglm.P <- function(y, x1, scheme, link.lambda = link.log,
  weights = rep(1, nrow(x1)), num.intercept = 1,
  xtol_rel = 1e-8, maxit = 100){
  if(typeof(link.lambda) == "closure") link.lambda <- link.lambda()
  # num.intercept: indicating whether x1 contains intercept. 1 if yes, 0 if no.
  ### PRE-PROCESSING OF DATA
  m <- nrow(x1)
  stopifnot(m == length(y))
  y <- c(y) # inherited from R lm, preventing y from having dimension 1.
  ## check singularity of x1
  svs <- svd(x = x1, nu = 0, nv = 0)$d
  if(max(svs) / min(svs) > 1e12){
    warning("coefficient matrix too singular, regression fails"); return(NULL)
  }
  rx <- ncol(x1)
  if(is.null(colnames(x1))) colnames(x1) <- paste0("x", 1:rx)
  ##############################################################################
  mleOutput <- mleInternal.P(y = y, x1 = x1, scheme = scheme,
    link.lambda = link.lambda, weights = weights,
    xtol_rel = xtol_rel, maxit = maxit)
  null.output <- mleInternal.P(y = y, x1 = matrix(1, nrow = m, ncol = 1),
    scheme = scheme, link.lambda = link.lambda, weights = weights,
    xtol_rel = xtol_rel, maxit = maxit)
  # START TO COMPUTE INFERENCE SCORES
  rpt <- list(coefficients = mleOutput$beta)
  rpt$mleOutput <- mleOutput
  rpt$beta <- mleOutput$beta
  names(rpt$coefficients) <- colnames(x1)
  names(rpt$beta) <- colnames(x1)
  rpt$fitting <- paste(deparse(match.call()), sep = "\n", collapse = "\n")
  rpt$log.likelihood <- mleOutput$log.likelihood
  l.saturated <- sum(weights * saturated.logL.P(scheme)[y])
  rpt$df.null <- m - num.intercept
  rpt$df.residual <- m - rx
  rpt$null.deviance <- 2 * (l.saturated - null.output$log.likelihood)
  rpt$deviance <- 2 * (l.saturated - mleOutput$log.likelihood)
  rpt$aic <- 2 * (rx - mleOutput$log.likelihood)
  rpt$bic <- rx * log(m) - 2 * mleOutput$log.likelihood
  rpt$McFaddenR2 <- 1 - mleOutput$log.likelihood / null.output$log.likelihood
  rpt$McFaddenAdjR2 <- 1 - (mleOutput$log.likelihood - rx) /
    null.output$log.likelihood
  ##############################################################################
  # feed parameters and postpone Fisher information matrix computing to summary.
  rpt$x1 <- x1
  rpt$y <- y
  rpt$weights <- weights
  rpt$scheme <- scheme
  rpt$link.lambda <- link.lambda
  class(rpt) <- "GRCglm.P"
  return(rpt)
}

summary.GRCglm.P <- function(object, level = 0.95, ...){
  stopifnot("GRCglm.P" %in% class(object))
  N <- length(object$scheme)
  object$lambdas <- object$link.lambda$gInv(c(object$x1 %*% object$beta))
  ##############################################################################
  # prepare Fisher information matrix
  # pmf: the probability matrix with m rows, and N := |G| columns.
  object$Imat <- make.Imat.4P(beta = object$beta, M = max(object$scheme),
    x1 = object$x1, link.lambda = object$link.lambda, weights = object$weights)
  object$FisherInfo <- rowSums(object$Imat$mat[, scheme2key(object$scheme),
    drop = FALSE])
  dim(object$FisherInfo) <- c(object$Imat$parDim, object$Imat$parDim)
  ##############################################################################
  object$var <- solvep(object$FisherInfo)
  object$stdErr <- sqrt(diag(object$var))
  names(object$stdErr) <- names(object$coefficients)
  interval <- c((1 - level) / 2, 1 - (1 - level) / 2)
  object$ci <- matrix(NA, nrow = length(object$coefficients), ncol = 2)
  for(i in 1:length(object$coefficients)){
    object$ci[i,] <- qnorm(interval, mean = object$coefficients[i],
      sd = object$stdErr[i])
  }
  colnames(object$ci) <- paste(round(100*interval, 1), "%")
  object$WOGRCStdErr <- sqrt(diag(solvep(object$Imat$vanillaFisherInfoMat)))
  object$zValue <- object$coefficients / object$stdErr
  object$sigLevel <- 2 * pnorm(- abs(object$zValue))
  class(object) <- c(class(object), "summary.GRCglm.P")
  return(object)
}

# library(MASS)
# source("GRCRegression_internal_functions.R")
# source("genData.R")
# set.seed(777)
# tp <- genData.P(beta = c(-3, 2, 2), data.size = 400, scheme = c(0:5, 7, 8), scope.lambda = c(1, 12))
# a <- summary(GRCglm.P(y = tp$y, x1 = tp$x, scheme = c(0:5, 7, 8)))
# b <- function(M)make.Imat.4P(beta = c(-3, 2, 2), M = M, x1 = tp$x, link.lambda = link.log(), weights = rep(1, nrow(tp$x)))
# d <- oracleEst(ImatMaker = b, scheme.list = list(c(0:4, 6, 7, 8), c(0:4, 7, 8, 9)))
