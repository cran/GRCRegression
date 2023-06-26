### artificial data generators
# when generating data, always generate intercept.
mkblk <- function(ncol, data.size){ # gen random x1 block.
  x <- matrix(rnorm(data.size * (ncol - 1)), nrow = data.size)
  return(cbind(1, x))
}
# mkblk(ncol = 5, data.size = 7)

compScp <- function(vec, scope){ # see if vec falls inside scope
  return((vec >= min(scope)) & (vec <= max(scope)))
}
# compScp(vec = 1:10, scope = c(3, 5))

count2grc <- function(cnt, scheme){ # for making grc data from count data
  return(rowSums(outer(cnt, scheme - 0.5, FUN = ">")))
}
# count2grc(cnt = 0:15, scheme = c(0, 1, 2, 5))

genData.P <- function(beta, data.size, scheme,
  scope.lambda = c(0.01, 100), link.lambda = link.log, seed = list(no = 1),
  y.lowRatio = -1){
  if(typeof(link.lambda) == "closure") link.lambda <- link.lambda()
  if(is.null(seed$no)) stopifnot(!is.null(seed$x), !is.null(seed$y))
  num.success <- unused <- 0L
  dim.beta <- length(beta)
  x <- lambdas <- c()
  if(!is.null(seed$x)) set.seed(seed$x)
  while(num.success < data.size){
    x.temp <- mkblk(ncol = dim.beta, data.size = 2 * data.size)
    lambdas.temp <- link.lambda$gInv(c(x.temp %*% beta))
    inds <- compScp(vec = lambdas.temp, scope = scope.lambda)
    num.success <- num.success + sum(inds)
    unused <- unused + sum(!inds)
    x <- rbind(x, x.temp[inds, , drop = FALSE])
    lambdas <- c(lambdas, lambdas.temp[inds])
  }
  x <- x[1:data.size, , drop = FALSE]
  lambdas <- lambdas[1:data.size]
  if(!is.null(seed$y)) set.seed(seed$y)
  y <- rpois(n = data.size, lambda = lambdas)
  y <- count2grc(cnt = y, scheme = scheme)
  yRat.ckpt <- table(c(1:length(scheme), y)) - 1
  if(any(yRat.ckpt / length(y) < y.lowRatio))
    return(genData.P(beta = beta, data.size = data.size, scheme = scheme,
      scope.lambda = scope.lambda, link.lambda = link.lambda,
      y.lowRatio = y.lowRatio))
  return(list(x = x, y = y, unused = unused))
}
# genData.P(beta = c(2, -3, 1), data.size= 7, scheme=c(0, 2, 4, 10))
# genData.P(beta = c(2), data.size= 7, scheme=c(0, 2, 4, 6))
