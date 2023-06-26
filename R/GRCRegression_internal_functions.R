link.log <- function(){
  l <- make.link("log")
  return(list(
    class = "0Inf",
    g = l$linkfun, gInv = l$linkinv,
    D.gInv = exp,
    DD.gInv = exp
  ))
}

solvep <- ginv # Maybe we need a finer one in the future. From Package "MASS"

pullUpEig <- function(A, threshold = 1e-7){
  # assume that A is symmetric. pull up its eigenvalues if they are too small
  temp <- eigen(A, symmetric = TRUE)
  if(all(temp$values > threshold)) return(A)
  temp$values <- pmax(temp$values, threshold * max(1, temp$values))
  return(temp$vectors %*% (temp$values * solve(temp$vectors)))
}
# pullUpEig(diag(1:7), threshold = 0.3)

PoisProb <- function(lwrs, uprs, lambdas){
  stopifnot(all(lwrs <= uprs))
  stopifnot(!any(is.infinite(lwrs) & is.infinite(uprs)))
  if(length(lwrs) == 1) lwrs <- rep(lwrs, length.out = length(lambdas))
  if(length(uprs) == 1) uprs <- rep(uprs, length.out = length(lambdas))
  useSum <- (uprs - lwrs < 50)
  feasible <- (lambdas < 1e10)

  probs <- rep(NaN, length.out = length(lambdas))
  if(sum(useSum & feasible) > 0.5)
    probs[useSum & feasible] <- unlist(lapply(which(useSum & feasible),
      function(i) sum(dpois(x = lwrs[i]:uprs[i], lambda = lambdas[i]))))

  if(sum((!useSum) & feasible) > 0.5){
    lwrs <- lwrs[(!useSum) & feasible] - 0.5
    uprs <- uprs[(!useSum) & feasible] + 0.5
    lambdas <- lambdas[(!useSum) & feasible]
    pt.lwr <- ppois(lwrs, lambdas, lower.tail = TRUE)
    pt.upr <- ppois(uprs, lambdas, lower.tail = TRUE)
    pf.lwr <- ppois(lwrs, lambdas, lower.tail = FALSE)
    pf.upr <- ppois(uprs, lambdas, lower.tail = FALSE)
    ind.f <- pf.lwr + pf.upr
    ind.t <- pt.lwr + pt.upr
    usef <- (ind.f < ind.t)
    temp <- pt.upr - pt.lwr
    temp[usef] <- pf.lwr[usef] - pf.upr[usef]
    probs[(!useSum) & feasible] <- temp
  }
  return(probs)
}
# PoisProb(lwrs = c(3,500,100, 1), uprs = c(5, 505, Inf, 100), lambdas = c(300, 100, 2, 300))

prot <- function(x, qs){ # force/protect NaN or Inf to zero
  inds <- ((abs(qs) < 1e-12) | (!is.finite(x)))
  x[inds] <- 0
  return(list(d = x, no = sum(inds)))
}
# prot(c(0:3, NA, NaN, Inf, -Inf, 3:0), 1:12)

makeX1D <- function(x1, u1, d){ # t(x1) %*% diag(d) %*% u1
  return(t(x1) %*% (d * u1))
}
# qq <- matrix(rnorm(49), nrow = 7)
# temp <- svd(qq <- qq %*% t(qq))
# qq - makeX1D(t(temp$u), t(temp$u), temp$d)

scheme2key <- function(scheme) return(paste0(scheme, ",", c(scheme[-1], Inf)))
# scheme2key(c(0,1,3,9,100))

saturated.logL.P <- function(scheme){
  # for i = 1, ..., |G|, find h[i], which is the maximum value of
  # log p(i | lambda), for group i.
  max.ll <- unlist(lapply(2:(length(scheme) - 1), function(i){
    optim.lambda <- exp(mean(log(scheme[i]:(scheme[i + 1] - 1))))
    return(log(PoisProb(lwrs = scheme[i], uprs = scheme[i + 1] - 1,
                        lambdas = optim.lambda)))
  }))
  return(c(0, max.ll, 0))
}

zLev.annotate <- function(x){
  temp <- c(" ", ".", "*", "**", "***")
  scr <- rowSums(outer(X = x, Y = c(0.001, 0.01, 0.05, 0.1), FUN = "<"))
  return(temp[scr + 1])
}

###PRINT GRCglm function
print.GRCglm.P <-
  function(x, digits = max(3L, getOption("digits") - 3L), ...)
  {
    cat("\nCall Fitting:  ", x$fitting, "\n\n", sep = "")
    if(length(x$coefficients) > 0.5) {
      cat("Coefficients:\n")
      print.default(format(x$coefficients, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    cat("Null Deviance:	   ",	format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nAIC:", format(signif(x$aic, digits)),
        "\tBIC:", format(signif(x$bic, digits)))
    cat("\n")
    cat("McFadden's R2:     ", format(signif(x$McFaddenR2, digits)), "\n")
    cat("McFadden's Adj R2: ", format(signif(x$McFaddenAdjR2, digits)))
    cat("\n")
    invisible(x)
  }

# tp <- genData.P(beta = c(-3, 2, 2, 4, -3, -3, 1), data.size = 9000, scheme = c(0, 2, 4, 10))
# output <- GRCglm.P(y = tp$y, x1 = tp$x, scheme = c(0, 2, 4, 10))

print.summary.GRCglm.P <-
  function(x, digits = max(3L, getOption("digits") - 3L), ...)
  {
    cat("\nCall Fitting:  ", x$fitting, "\n\n", sep = "")
    f <- function(u) format(u, digits = digits)

    if(length(x$coefficients) > 0.5) {
      cat("Coefficients:\n")
      temp <- cbind(f(x$coefficients), f(x$stdErr), f(x$zValue), f(x$sigLevel), f(x$ci))
      rownames(temp) <- names(x$coefficients)
      colnames(temp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", colnames(x$ci))
      temp[x$sigLevel < 2e-16, 4] <- "<2e-16"
      temp <- cbind(temp, zLev.annotate(x$sigLevel))
      print.default(temp, print.gap = 2, quote = FALSE)
      cat("---\nSignif. codes:  ", paste0(c(0, 0.001, 0.01, 0.05, 0.1), " ",
        intToUtf8(8216), c("***", "**", "*", ".", " "), intToUtf8(8217), " ",
        collapse = ""), "1", sep = "")
    } else cat("No coefficients\n\n")
    cat("\n\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    cat("Null Deviance:	   ",	format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nAIC:", format(signif(x$aic, digits)),
        "\tBIC:", format(signif(x$bic, digits)))
    cat("\n")
    cat("McFadden's R2:     ", format(signif(x$McFaddenR2, digits)), "\n")
    cat("McFadden's Adj R2: ", format(signif(x$McFaddenAdjR2, digits)))
    cat("\n")
    if(!is.null(x$optimScheme)){
      cat("\n\n======= ANALYSIS OF OPTIMAL GROUPING SCHEME =======")
      cat("\nOutput of the Searching Algorithm for Best Grouping Scheme\n")
      cat("-----------------------------------------------------\n")
      cat("CURRENT GROUPING SCHEME\n")
      cat("Scheme: ", x$scheme, "\n")
      cat("Score:  ", x$scoreFun(x$FisherInfo), "\n")
      cat("-----------------------------------------------------\n")
      cat("OPTIMAL GROUPING SCHEME\n")
      cat("Scheme: ", x$optimScheme$scheme, "\n")
      cat("Score:  ", x$optimScheme$score, "\n")
      cat("global optimal achieved: ", x$optimScheme$successful, "\n")
      cat("-----------------------------------------------------\n")
      cat("NO GROUPING\n")
      cat("Score:  ", x$scoreFun(x$Imat$vanillaFisherInfoMat), "\n")
      cat("-----------------------------------------------------\n")
      cat("Comparison of standard errors of each predictor:\n")
      tpCurStdErr <- f(x$stdErr)
      tpOptStdErr <- paste0(f(x$optimStdErr),
        " (", f(x$optimStdErr / x$stdErr * 100), "%)")
      tpWOGRCStdErr <- paste0(f(x$WOGRCStdErr),
        " (", f(x$WOGRCStdErr / x$stdErr * 100), "%)")
      char.opt <- 1 / diag(solve(x$optimScheme$fisherInfoMat))
      char.cur <- 1 / diag(solve(x$FisherInfo))
      char.inf <- 1 / diag(solve(x$Imat$vanillaFisherInfoMat))

      TotalErr    <- (char.inf - char.cur) / char.inf
      GroupsErr   <- (char.inf - char.opt) / char.inf
      GroupingErr <- (char.opt - char.cur) / char.inf
      FisherRatio <- char.cur / char.inf

      tempMat <- cbind(TotalErr, GroupsErr, GroupingErr, FisherRatio)
      rownames(tempMat) <- names(x$coefficients)
      x$GrpsGrpgErrs <- tempMat
      print.default(tempMat, print.gap = 2, quote = FALSE)
    }
    invisible(x)
  }
# tp <- genData.P(beta = c(-3, 2, 2, 4, -3, -3, 1), data.size = 9000, scheme = c(0, 1, 4, 6), scope.lambda = c(0.1, 12))
# a <- summary(GRCglm.P(y = tp$y, x1 = tp$x, scheme = c(0, 1, 4, 6)))

tryOp <- function(Op, fun, pars){ # Op: grad or hessian in package *pracma*
  for(k in c(5:9)){
    if(!all(is.finite(pars))) next
    ans <- try(Op(f = fun, x0 = pars, h = 10^(-k)), silent = TRUE)
    if("try-error" %in% class(ans)) next
    if(!all(is.finite(c(ans)))) next
    return(ans)
  }
  return(pars + NaN)
}
glmMaximizer <- function(par.init, fun, GrHess, maxit = 100, xtol_rel = 1e-8){
  par.new <- par.init
  stepSzs <- t(10 ^ (5:(-9)))
  nm <- function(x) sqrt(sum(c(x * x)))
  normalize <- function(x){
    the.norm <- nm(x)
    return(if(the.norm < 1e-6) c(x) else c(x / the.norm))
  }
  iter.count <- 1
  msg <- ""; newVal <- oldVal <- fun(par.init)
  while(iter.count < maxit){
    par.old <- par.new; oldVal <- newVal
    ana.GH <- GrHess(par.old)
    num.gr <- tryOp(grad, fun, par.old)
    if(!all(is.finite(num.gr))) {msg <- "gradient failed"; break}
    gradDir <- normalize(num.gr)
    hessCols <- c()
    if(all(is.finite(c(ana.GH$Hess)))){
      hessDir <- normalize(solve(pullUpEig(-ana.GH$Hess), ana.GH$gr))
      hessCols <- hessDir %*% stepSzs
      # Newton iteration requires Negative Hessian
    }
    newX <- cbind(gradDir %*% stepSzs, hessCols) + par.old
    ys <- apply(newX, 2, fun)

    if(sum(is.finite(ys)) < 0.5){msg <- "no feasible point ahead"; break}
    newX <- newX[, is.finite(ys), drop = FALSE]
    ys <- ys[is.finite(ys)]
    newVal <- max(ys)
    if(newVal <= oldVal){msg <- "seems peak achieved"; break}
    par.new <- newX[, which(ys == newVal)[1]]
    if(nm(par.new - par.old) < xtol_rel * nm(par.old)) break
    iter.count <- iter.count + 1
  }
  return(list(par = par.new, val = newVal,
    gr = grad(fun, par.new), #ok to fail
    msg = msg, num.iter = iter.count))
}
# par.init <- 1:5 + 3000
# fun <- function(x) sum(-x*x)
# GrHess <- function(x){
#   return(list(f = fun(x), gr = -2 * x, Hess = -2 * diag(5)))
# }
# aa <- glmMaximizer(par.init, fun, GrHess)
