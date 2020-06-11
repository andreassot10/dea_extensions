dea_gle <- function(base, noutput, eff) { 
  require(lpSolve)
  require(dplyr)
  
	s <- noutput
	m <- ncol(base) - s
  n <- nrow(base)
	GLE <- list() 
	weights <- list()
	k <- 0
	baseNormalized <- sweep(base, 2, apply(base, 2, max), '/')
	
	for (t in (1 - seq(0, 1, 0.01))) {
		k <- k + 1
		if (t == 1) {
			f.dir <- c(rep("==", n), "==")
			f.rhs <- c(rep(0, n), 1)
			f.obj <- c(rep(1 / n, n), rep(0, s + m))
			f.con1 <- cbind(
			  diag(n), 
			  baseNormalized[, 1:s],
				-eff * baseNormalized[, (s + 1):(s + m)]
			  ) %>%
			  data.frame
			f.con2 <- c(rep(0, n), rep(1, s + m))
			names(f.con2) <- names(f.con1)
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
			f.obj <- rep(0, s + m)
			f.con1 <- cbind(
		    baseNormalized[, 1:s], 
			  -eff * baseNormalized[, (s + 1):(s + m)]
			  ) %>% 
			  data.frame
			f.con4 <- cbind(
			  diag(s), 
			  matrix(0, s, m)
			  ) %>%
			  data.frame
			f.con5 <- cbind(
			  matrix(0, m, s), 
			  diag(m)
			  ) %>%
			  data.frame
			names(f.con4) <- names(f.con1)
			names(f.con5) <- names(f.con1)
			while ((zupper - zlower) >= 10 ^ -6) {
				z <- mean(c(zupper, zlower))
				f.con2 <- cbind(
				  baseNormalized[, 1:s],
				  (-eff - z) * baseNormalized[, (s + 1):(s + m)]
				  ) %>%
				  data.frame
				f.con3 <- cbind(
				  -baseNormalized[, 1:s],
					(eff - z) * baseNormalized[, (s + 1):(s + m)]
					) %>%
				  data.frame
				names(f.con2) <- names(f.con1)
				names(f.con3) <- names(f.con1)
				f.con <- rbind(f.con1, f.con2, f.con3, f.con4, f.con5)
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
			f.con1 <- cbind(
			  diag(n), 
			  rep(0, n), 
			  baseNormalized [, 1:s],
				-eff * baseNormalized [, (s + 1):(s + m)]
				) %>%
			  data.frame
			f.con2 <- cbind(
			  diag(n), 
			  rep(-1, n), 
			  matrix(0, n, s + m)
			  ) %>%
			  data.frame
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
	weights <- weights %>%
	  unlist %>%
	  matrix(ncol = k, byrow = FALSE) %>%
	  data.frame
	#weights <- data.frame(unique(t(round(weights, 7))))
	weights <- weights %>%
	  t %>%
	  data.frame
	names(weights) <- c(paste("u", 1:s, sep = ""), paste("v", 1:m, sep = ""))
	GLE <- GLE %>%
	  unlist %>%
	  matrix(ncol = k, byrow = FALSE) %>%
	  data.frame
	names(GLE) <- paste0('t_', 1 - seq(0, 1, 0.01))
	rownames(weights) <- names(GLE)
	#GLE <- t(unique(t(round(GLE, 3))))
	GLE <- t(t(GLE))
	GLE[GLE > 1] <- 1
	GLE <- data.frame(
	  eff = eff, 
	  GLE, 
	  AVGLE = round(rowMeans(GLE, na.rm = TRUE), 3)
	)
	GLE$q <- rowSums(
	  ifelse(
	    GLE[-c(1, ncol(GLE))] >= 1 | 
	    near(GLE[-c(1, ncol(GLE))], 1), 
	   1, 0
	   )
	)
	GLE$RankingFactor <- GLE$AVGLE + GLE$q
	GLE <- list(
	  GLE = GLE, 
	  weights = weights
	)
	return(GLE)
}