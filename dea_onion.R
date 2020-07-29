dea_onion <- function(base, noutput, rts = 2, strong_eff = FALSE) {
  
  require(Benchmarking)
  require(dplyr)
  
  X = as.matrix(base[, (noutput + 1):ncol(base)])
  Y = as.matrix(base[, 1:noutput])
  
  base_row_names <- 1:nrow(X)
  tiers <- list()
  
  RTS <- ifelse(rts == 2, 'vrs', 'crs')
  
  RTS <- if (is.numeric(rts)) {
    case_when(
      rts == 2 ~ 'vrs', 
      rts == 1 ~ 'crs',
      TRUE ~ NA_character_
    )
  }
  
  if (is.na(RTS)) {
    stop(
      print(
        'You have supplied an erroneous numeric value for the returns-to-scale specification.\n
        Please set rts = 1 for CRS or rts = 2 for VRS.\n
        Otherwise, you must supply a string, which should be one of the following:\n
        "fdh", "vrs", "drs", "crs", "irs", "irs2", "add", "fdh+" or "vrs+".\n
        See Benchmarking::dea'
      )
    )
  }
  
  re <- dea(X, Y, RTS, SLACK = strong_eff)
  
  re <- near(re$objval, 1) & 
    ifelse(strong_eff, !re$slack, TRUE)
  
  tiers[[1]] <- base_row_names[re]
  
  k <- 1
  while (!all(is.na(tiers[[k]]))) {
    k <- k + 1
    X <- matrix(X[!(base_row_names %in% tiers[[k - 1]]), ], 
      ncol = ncol(base) - noutput)
    Y <- matrix(Y[!(base_row_names %in% tiers[[k - 1]]), ], ncol = noutput)
    
    base_row_names <- base_row_names[!(base_row_names %in% tiers[[k - 1]])]
    
    if (nrow(X) > 1) {
      re <- dea(X, Y, RTS, SLACK = strong_eff)
      
      re <- near(re$objval, 1) & 
        ifelse(strong_eff, !re$slack, TRUE)
      
      tiers[[k]] <- base_row_names[re]
    } else {
      tiers[[k]] <- NA
    }
    
  }
  
  tiers[[k]] <- base_row_names
  names(tiers) <- paste0('tier_', 1:k)
  tiers <- Filter(Negate(function(x) length(x) == 0), tiers)
  
  tiers <- do.call(rbind, lapply(tiers, as.data.frame)) %>%
    mutate(
      tier = gsub("\\..*", "", rownames(.)),
      tier = sub('tier_', '', tier)
    ) %>%
    rename_at(vars(-contains('tier')), ~ 'dmu')
  
  return(tiers)
}
