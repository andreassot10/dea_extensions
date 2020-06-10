#I've made corrections but still don't get same results with Despotis' original paper!
dea_gle <- function(base, noutput, eff) { 
  require(lpSolve)
	s <- noutput
	m <- ncol(base) - s
  n <- nrow(base)
	GLE <- list() 
	weights <- list()
	k <- 0
	aux <- apply(base, 2, max)
	baseNormalized <- sweep(base, 2, aux, '/') 
	for (t in (1 - seq(0, 1, 0.01))) {
		k <- k + 1
		if (t == 1) {
			f.dir <- c(rep("==", n), "==")
			f.rhs <- c(rep(0, n), 1)
			f.obj <- c(rep(1 / n, n), rep(0, s + m))
			f.con1 <- data.frame(cbind(diag(n), baseNormalized[, 1:s],
				-eff * baseNormalized[, (s + 1):(s + m)]))
			f.con2 <- c(rep(0, n), rep(1, s + m))
			names(f.con2) <- names(f.con1) ##
			f.con <- rbind(f.con1, f.con2)
			gle <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
			weights[[k]] <- gle[(n + 1):length(gle)]
			gle <- sweep(baseNormalized, 2, weights[[k]], '*') 
			re <- matrix(0, n, 1)
			re <- rowSums(data.frame(gle[, 1:s])) / 
			  rowSums(data.frame(gle[, (s + 1):ncol(gle)]))
			GLE[[k]] <- data.frame(re)
		} else if (t == 0) {
			zupper <- 1
			zlower <- 0
			f.dir <- c(rep("<=", n + n + n), rep(">=", s + m))
			f.rhs <- c(rep(0, n + n + n), rep(10 ^ -6, s + m))
			#f.dir <- c(rep("<=", n), rep(">=", s + m))
			#f.rhs <- c(rep(0, n), rep(10 ^ -6, s + m))
			f.obj <- rep(0, s + m) ##!!!Find out how to define dummy obj. function!!!
			f.con1 <- data.frame(cbind(baseNormalized[, 1:s], -eff * baseNormalized[, (s + 1):(s + m)]))
			f.con4 <- data.frame(cbind(diag(s), matrix(0, s, m)))
			f.con5 <- data.frame(cbind(matrix(0, m, s), diag(m)))
			names(f.con4) <- names(f.con1)
			names(f.con5) <- names(f.con1)
			while ((zupper - zlower) >= 10 ^ -6) {
				z <- mean(c(zupper, zlower))
				f.con2 <- data.frame(cbind(baseNormalized[, 1:s],
				  (-eff - z) * baseNormalized[, (s + 1):(s + m)]))
				f.con3 <- data.frame(cbind(-baseNormalized[, 1:s],
					(eff - z) * baseNormalized[, (s + 1):(s + m)]))
				names(f.con2) <- names(f.con1)
				names(f.con3) <- names(f.con1)
				f.con <- rbind(f.con1, f.con2, f.con3, f.con4, f.con5)
				#f.con <- rbind(f.con2, f.con3, f.con4)
				if (lp("min", f.obj, f.con, f.dir, f.rhs)$status != 0) {
					zlower <- z
				} else {
					zupper <- z
					gle <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
				}
			}
			weights[[k]] <- gle
			gle <- sweep(baseNormalized, 2, weights[[k]], '*') 
			re <- matrix(0, n, 1)
			re <- rowSums(data.frame(gle[, 1:s])) / 
			  rowSums(data.frame(gle[, (s + 1):ncol(gle)]))
			GLE[[k]] <- data.frame(re)
		} else {
			f.dir <- c(rep("==", n), rep("<=", n), "==")
			f.rhs <- c(rep(0, n + n), 1)
			f.obj <- c(rep(t * (1 / n), n), (1 - t), rep(0, s + m))
			f.con1 <- data.frame(cbind(diag(n), rep(0, n), baseNormalized [, 1:s],
				-eff * baseNormalized [, (s + 1):(s + m)]))
			f.con2 <- data.frame(cbind(diag(n), rep(-1, n), matrix(0, n, s + m)))
			names(f.con1) <- names(f.con2)
			f.con3 <- c(rep(0, n + 1), rep(1, s + m))
			names(f.con3) <- names(f.con2)
			f.con <- rbind(f.con1, f.con2, f.con3)
			gle <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
			weights[[k]] <- gle[(n + 2):length(gle)]
			gle <- sweep(baseNormalized, 2, weights[[k]], '*') 
			re <- matrix(0, n, 1)
			re <- rowSums(data.frame(gle[, 1:s])) / 
			  rowSums(data.frame(gle[, (s + 1):ncol(gle)]))
		}
		GLE[[k]] <- data.frame(re)
	}
	weights <- data.frame(matrix(unlist(weights), ncol = k, byrow = FALSE))
	#weights <- data.frame(unique(t(round(weights, 7))))
	weights <- data.frame((t(weights)))
	names(weights) <- c(paste("u", 1:s, sep = ""), paste("v", 1:m, sep = ""))
	GLE <- data.frame(matrix(unlist(GLE), ncol = k, byrow = F))
	#GLE <- t(unique(t(round(GLE, 3))))
	GLE <- t(t(GLE))
	#GLE[GLE >= 1] <- 1
	GLE <- data.frame(eff = eff, GLE, 
	  AVGLE = round(rowMeans(GLE, na.rm = TRUE), 3))
	GLE$q <- rowSums(ifelse(GLE[-c(1, ncol(GLE))] >= 1 | near(GLE[-c(1, ncol(GLE))], 1), 1, 0))
	GLE$RankingFactor <- GLE$AVGLE + GLE$q
	GLE <- list(GLE = GLE, weights = weights)
	return(GLE)
}