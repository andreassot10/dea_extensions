dea_gle <- function(base, noutput, eff, t_values = seq(0, 1, 0.01)) { 
  
  require(lpSolve)
  require(dplyr)
  
	s <- noutput
	m <- ncol(base) - s
  n <- nrow(base)
	baseNormalized <- sweep(base, 2, apply(base, 2, max), '/')
	
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
	lin_prog <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
	weights_l1 <- lin_prog[(n + 1):length(lin_prog)]
	GLE_l1 <- baseNormalized %>% 
	  sweep(2, weights_l1, '*') %>%
	  mutate(
	    gle = 
	      rowSums(.[1:s]) / 
	      rowSums(.[(s + 1):ncol(.)])
	  ) %>%
	  select(gle)
	
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
	    weights_linf <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
	  }
	}
	GLE_linf <- baseNormalized %>% 
	  sweep(2, weights_linf, '*') %>%
	  mutate(
	    gle = 
	      rowSums(.[1:s]) / 
	      rowSums(.[(s + 1):ncol(.)])
	  ) %>%
	  select(gle)
	
	re <- t_values[c(-1, -length(t_values))] %>%
    mapply(
	    FUN = function(t) {
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
	      lin_prog <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
	      weights <- lin_prog[(n + 2):length(lin_prog)]
	      list(
	        baseNormalized %>% 
	          sweep(2, weights, '*') %>%
	          mutate(
	            gle = 
	              rowSums(.[1:s]) / 
	              rowSums(.[(s + 1):ncol(.)])
	          ) %>%
	          select(gle),
	        weights
	      )
	    }
    )
	GLE_rest <- re[1, ] %>%
	  bind_cols
	weights_rest <- do.call(rbind, re[2, ])
	
	weights <- rbind(
	  weights_linf, 
	  weights_rest, 
	  weights_l1
	) %>%
	as.data.frame
	names(weights) <- c(paste("u", 1:s, sep = ""), paste("v", 1:m, sep = ""))
	
	GLE <- cbind(
	  GLE_linf,
	  GLE_rest,
	  GLE_l1
	)
	names(GLE) <- paste0('t_', round(t_values, 2))
	rownames(weights) <- names(GLE)
	
	GLE[GLE > 1] <- 1
	
	GLE <- data.frame(
	  eff = eff, 
	  GLE, 
	  AVGLE = round(rowMeans(GLE, na.rm = TRUE), 3)
	) %>%
	  mutate(
	    q = rowSums(
	      ifelse(
	        .[-c(1, ncol(.))] >= 1 | 
	          near(.[-c(1, ncol(.))], 1), 
	        1, 0
	      )
	    ),
	    RankingFactor = AVGLE + q
	  )
	
	GLE <- list(
	  GLE = GLE, 
	  weights = weights
	)
	
	return(GLE)
}