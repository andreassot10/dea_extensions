dea_ar_type1_primal <- function (base, noutput, rts = 1, orientation = 1, 
  bounds = NULL, duals = TRUE) {
  
  require(lpSolveAPI)
  
  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)
  P <- as.data.frame(matrix(0, nrow = m, ncol = 2 * m - 2))
  Q <- as.data.frame(matrix(0, nrow = s, ncol = 2 * s - 2))
  if (s > 1 & m > 1) {
    re <- data.frame(matrix(0, nrow = n, 
      ncol = 1 + n + (2 * s - 2) + (2 * m - 2)))
    names(re) <- c('eff', paste('lambda', 1:n, sep = ''), 
      paste('tau.y', 1:(2 * s - 2), sep = ''), 
      paste('pi.x', 1:(2 * m - 2), sep = ''))
    P[1, which(1:ncol(P) %% 2 == 1)] <- bounds[1:(m - 1), 1]
    P[1, which(1:ncol(P) %% 2 != 1)] <- -bounds[1:(m - 1), 2]
    Q[1, which(1:ncol(Q) %% 2 == 1)] <- bounds[m:nrow(bounds), 1]
    Q[1, which(1:ncol(Q) %% 2 != 1)] <- -bounds[m:nrow(bounds), 2]
    k <- 1
    for (i in 2:nrow(P)) {
      P[i, k:(k + 1)] <- c(-1, 1)
      k <- k + 2
    }
    k <- 1
    for (i in 2:nrow(Q)) {
      Q[i, k:(k + 1)] <- c(-1, 1)
      k <- k + 2
    }
  } else if (s > 1 & m == 1) {
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + n + (2 * s - 2)))
    names(re) <- c('eff', paste('lambda', 1:n, sep = ''), 
      paste('tau.y', 1:(2 * s - 2), sep = ''))
    Q[1, which(1:ncol(Q) %% 2 == 1)] <- bounds[m:nrow(bounds), 1]
    Q[1, which(1:ncol(Q) %% 2 != 1)] <- -bounds[m:nrow(bounds), 2]
    k <- 1
    for (i in 2:nrow(Q)) {
      Q[i, k:(k + 1)] <- c(-1, 1)
      k <- k + 2
    }
  } else if (s == 1 & m > 1) {
    re <- data.frame(matrix(0, nrow = n, ncol = 1 + n + (2 * m - 2)))
    names(re) <- c('eff', paste('lambda', 1:n, sep = ''), 
      paste('pi.x', 1:(2 * m - 2), sep = ''))
    P[1, which(1:ncol(P) %% 2 == 1)] <- bounds[1:(m - 1), 1]
    P[1, which(1:ncol(P) %% 2 != 1)] <- -bounds[1:(m - 1), 2]
    k <- 1
    for (i in 2:nrow(P)) {
      P[i, k:(k + 1)] <- c(-1, 1)
      k <- k + 2
    }
  }
  P1 <- data.frame(P, matrix(0, nrow = nrow(P), ncol = ncol(Q)))
  Q <- data.frame(matrix(0, nrow = nrow(Q), ncol = ncol(P)), Q)
  P <- P1
  if (m > 1) {
    index <- 1 + which(rowSums(bounds)[1:(m - 1)] == 0)
    P[index, 1:(2 * m - 2)] <- 
      0 * P[index, 1:(2 * m - 2)]
  }
  if (s > 1) {
    index <- 1 + which(rowSums(bounds)[-(1:(m - 1))] == 0)
    Q[index, 1:(2 * s - 2)] <- 
      -abs(Q[index, 1:(2 * s - 2)])
  }
  names(P) <- NULL
  names(Q) <- NULL
  A.aux <- rbind(cbind(t(base[1:s]), Q), 
    cbind(-t(base[(s + 1):ncol(base)]), P))
  if (orientation == 1) {
    rhSide.aux <- c(rep(1, s), rep(0, m))
    thetaColumn <- c(rep(0, s), rep(1, m))
  } else {
    rhSide.aux <- c(rep(0, s), rep(1, m))
    thetaColumn <- c(rep(1, s), rep(0, m))
  }
  type <- c(rep('>=', s + m), '=')
  if (rts == 1) {
    A_finalRow <- rep(0, ncol(A.aux))
    rhSide.aux <- c(rhSide.aux, 0)
  } else {
    A_finalRow <- c(rep(1, n), rep(0, (2 * s - 2) + (2 * m - 2)))
    rhSide.aux <- c(rhSide.aux, 1)
  }
  A_finalRow <- c(0, A_finalRow)
  if (duals) {
    d.yx <- data.frame(matrix(0, nrow = nrow(base), ncol = s + m))
    names(d.yx) <- c(paste('dual.y', 1:s, sep = ''), 
      paste('dual.x', 1:m, sep = ''))
    re <- data.frame(re, d.yx)
  }
  lpmodel <- make.lp(
    nrow = m + s + 1, 
    ncol = 1 + n + (2 * s - 2) + (2 * m - 2)
  )
  set.objfn(
    lpmodel, 
    c(1, rep(0, n + (2 * s - 2) + (2 * m - 2)))
  )
  set.constr.type(
    lpmodel, 
    types = type, 
    constraints = 1:(m + s + 1)
  )
  for (i in 1:nrow(base)) {
    A <- data.frame(thetaColumn = unlist(base[i, ] * thetaColumn), A.aux)
    A <- rbind(A, A_finalRow)
    rhSide <- rhSide.aux
    rhSide[-length(rhSide.aux)] <- unlist(base[i, ] * 
      rhSide.aux[-length(rhSide.aux)])
    for (j in 1:ncol(A)) {
      set.column(lpmodel, column = j, x = A[, j])
    }
    set.rhs(lpmodel, b = rhSide)
    set.objfn(lpmodel, c(1, rep(0, n + (2 * s - 2) + (2 * m - 2))))
    solve(lpmodel)
    solution <- get.primal.solution(lpmodel)
    re[i, 1:ncol(A)] <- tail(solution, ncol(A))
    if (duals) {
      multipliers <- get.dual.solution(lpmodel)[2:(s + m + 1)]
      re[i, names(d.yx)] <- multipliers
    }
  }
  lambdas <- re[grep('lambda', names(re))]
  #if (m > 1 & s > 1) {
  slacks <- data.frame(-base[c(1:s)] + 
    as.matrix(lambdas) %*% as.matrix(base[c(1:s)]) + 
    t(as.matrix(Q[-c(1:(2 * m - 2))]) %*% t(re[grep('tau.y', names(re))])), 
    sweep(base[-c(1:s)], 1, re$eff, '*') - 
    as.matrix(lambdas) %*% as.matrix(base[-c(1:s)]) + 
    t(as.matrix(P[1:(2 * m - 2)]) %*% t(re[grep('pi.x', names(re))])))
  names(slacks) <- c(paste('slack.y', 1:s, sep = ''), 
    paste('slack.x', 1:m, sep = ''))
  #} else {
  # slacks <- 0
    #}
  re <- data.frame(re, slacks)
  return(re)
}