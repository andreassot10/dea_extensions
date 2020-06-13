dea_onion <- function(base, noutput, rts = 2, 
  orientation = c('in', 'out', 'non')) {
  
  require(Benchmarking)
  require(dplyr)
  
  base$row_names <- 1:nrow(base)
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
  
  try(
    if (is.na(RTS)) {
      stop(
        'You have supplied an erroneous numeric value for the returns-to-scale specification.\n
        Please set rts = 1 for CRS or rts = 2 for VRS.\n
        Otherwise, you must supply a string, which should be one of the following:\n
        "fdh", "vrs", "drs", "crs", "irs", "irs2", "add", "fdh+" or "vrs+".\n
        See ?Benchmarking::dea'
      )
    }
  )
  
  if (orientation == 'in') {
    re <- dea(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation,
      SLACK = TRUE
    )
    
    re <- near(re$eff, 1) %>% 
      rbind(!re$slack) %>%
      apply(2, all)
  } else if (orientation == 'out') {
    re <- dea(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS,
      ORIENTATION = orientation,
      SLACK = TRUE
    )
    
    re <- near(1 / re$eff, 1) %>% 
      rbind(!re$slack) %>%
      apply(2, all)
  } else if (orientation == 'non') {
    re <- !dea.add(
      X = as.matrix(base[, (noutput + 1):ncol(base)]), 
      Y = as.matrix(base[, 1:noutput]), 
      RTS = RTS
    )$slack
  } else {
    'wrong'
  }
  
  tiers[[1]] <- base_row_names[re]
  
  k <- 1
  while (!all(is.na(tiers[[k]]))) {
    k <- k + 1
    base <- filter(base, !(base_row_names %in% tiers[[k - 1]]))
    
    if (nrow(base) > 1) {
      if (orientation == 'in') {
        re <- dea(
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation,
          SLACK = TRUE
        )
        
        re <- near(re$eff, 1) %>% 
          rbind(!re$slack) %>%
          apply(2, all)
      } else if (orientation == 'out') {
        re <- dea(
          X = as.matrix(base[, (noutput + 1):ncol(base)]), 
          Y = as.matrix(base[, 1:noutput]), 
          RTS = RTS,
          ORIENTATION = orientation,
          SLACK = TRUE
        )
        
        re <- near(1 / re$eff, 1) %>% 
          rbind(!re$slack) %>%
          apply(2, all)
      } else if (orientation == 'non') {
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
