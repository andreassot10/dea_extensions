dea_onion <- function(base, noutput, rts = 2, strong_eff = FALSE,
  orientation = c('in', 'out'), direction = NULL,
  dea_model = c('ccr', 'bcc', 'ddf', 'mea', 'additive')) {
  
  require(Benchmarking)
  require(dplyr)
  
  base_row_names <- 1:nrow(base)
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
  
  if (dea_model %in% c('ccr', 'bcc')) {
    re <- dea(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation,
      SLACK = strong_eff
    )
    
    re <- near(re$objval, 1) & 
      ifelse(strong_eff, !re$slack, TRUE)
  } else if (dea_model == 'ddf') {
    re <- dea.direct(
      DIRECT = direction,
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation,
      SLACK = strong_eff
    )
    
    re <- near(re$objval, 0) & 
      ifelse(strong_eff, !re$slack, TRUE)
  } else if (dea_model == 'mea') {
    me <- mea(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation
    )
    
    re <- dea.direct(
      DIRECT = me$direct,
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation,
      SLACK = strong_eff
    )
    
    re <- !is.finite(re$objval) & 
      ifelse(strong_eff, !re$slack, TRUE)
  }  else if (dea_model == 'additive') {
    re <- !dea.add(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS
    )$slack
  }
  
  tiers[[1]] <- base_row_names[re]
  
  k <- 1
  while (!all(is.na(tiers[[k]]))) {
    k <- k + 1
    base <- filter(base, !(base_row_names %in% tiers[[k - 1]]))
    
    if (nrow(base) > 1) {
      if (dea_model %in% c('ccr', 'bcc')) {
        re <- dea(
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation,
          SLACK = strong_eff
        )
        
        re <- near(re$objval, 1) & 
          ifelse(strong_eff, !re$slack, TRUE)
      } else if (dea_model == 'ddf') {
        re <- dea.direct(
          DIRECT = direction,
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation,
          SLACK = strong_eff
        )
        
        re <- near(re$objval, 0) & 
          ifelse(strong_eff, !re$slack, TRUE)
      } else if (dea_model == 'mea') {
        me <- mea(
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation
        )
        
        re <- dea.direct(
          DIRECT = me$direct,
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation,
          SLACK = strong_eff
        )
        
        re <- !is.finite(re$objval) & 
          ifelse(strong_eff, !re$slack, TRUE)
      }  else if (dea_model == 'additive') {
        re <- !dea.add(
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS
        )$slack
      }
      
      tiers[[k]] <- base_row_names[re]
    } else {
      tiers[[k]] <- NA
    }
    
  }
  
  tiers[[k]] <- base_row_names[-unlist(tiers[-k])]
  names(tiers) <- paste0('tier_', 1:k)
  
  return(tiers)
}
