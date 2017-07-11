#' Rank and Replace method of adjustment for racial disparities 
#'
#' Idea: fit multivariate generalized linear model y ~ X'b, collect coeffs b    
#' The rank-and-replace method adjusts health status by ranking each sample by  
#' a summary index of health status and replacing the health status of each     
#' minority individual with that of the correspondingly ranked white,           
#' thus preserving the ranking of health status and its rank correlation with   
#' SES measures.                                                                
#' @param m a model.frame for the interior GLM function
#' @param x design matrix x, sans intercept
#' @param y target or dependent outcome
#' @param family generalized linear model family see help(glm), help(family), defaults to Gamma
#' @param link link function for generalized linear model
#' @keywords GLM, racial disparities, health disparities
#' @export
#' @examples
#' iomDisparity.glm()

iomDisparity.glm <- function(m, x, y, index, race,
                             family = Gamma, 
                             link = "log",
                             method = "glm.fit",
                             ...) {
  temp_id <- 1:dim(x)[1]
  x_id <- cbind(x, temp_id) # assign temporary id
  base_model <- glm(m, 
                    family = family(link = link), 
                    method = "glm.fit")
  base_coeffs <- base_model$coefficients
  # gather fitted values for just the index subset
  
  coeffs <- as.vector(base_coeffs[index + 1]) # + 1 account for intercept
  X_t <- t(data.matrix(x[index]))
  H_i <- apply(t(coeffs * X_t), 1, sum)
  H_i <- cbind(H_i, temp_id)
  
  # initiate transformation stage
  # add race
  H_i_plus_race <- as.data.frame(cbind(H_i, race = race))
  
  # split by race
  splitfunc1 <- function(x)  H_i_plus_race[H_i_plus_race$race == x, ]
  H_i_by_race <- lapply(unique(H_i_plus_race$race), splitfunc1)
  splitfunc2 <- function(i, y) {
    assign(paste0("race", i), y[[i]], envir = .GlobalEnv)
  }
  lapply(seq_along(H_i_by_race), splitfunc2, y = H_i_by_race)
  
  percRank <- function(x) trunc(rank(x)) / length(x) # calculates percentiles
  
  H_rank <- percRank(race1$H_i)
  race1 <- cbind(race1, H_rank)
  race1 <- race1[with(race1, order(H_rank)), ]
  
  H_rank <- percRank(race2$H_i)
  race2 <- cbind(race2, H_rank)
  race2 <- race2[with(race2, order(H_rank)), ]
  
  # TO DO: convert to Rcpp function
  LR <- 0 # last replaced
  for (i in 1:length(race2$H_rank)) {
    for (j in 1:(length(race1$H_rank) - LR)) {
      if (race2$H_rank[i] <= race1$H_rank[j + LR]) {
        race2$H_i[i] = race1$H_i[j + LR]
        LR = j + LR
        break
      }
    }
  }
  
  # H_i has been transformed for minorities, now to gather original coefficients
  sub1 <- race1[, c("H_i", "temp_id")]
  sub2 <- race2[, c("H_i", "temp_id")]
  transformed <- as.data.frame(rbind(sub1, sub2))
  
  # construct design matrix X
  new_x <- cbind(intercept = rep(1,dim(x)[1]), x_id[-index])
  sub_tmp <- subset(new_x, select = -c(temp_id))
  sub_mat <- t(data.matrix(sub_tmp))
  coeffs <- as.vector(base_coeffs[-(index + 1)])
  temp <- sub_mat * coeffs
  design_x <- as.data.frame(cbind(t(temp), temp_id))
  
  # collect X[-index] + T(=transformed)
  new_x <- merge(design_x, transformed, by = "temp_id")
  new_x$temp_id <- NULL
  
  new_y_i <- exp(apply(new_x, 1, sum))
  old_y_i <- exp(base_model$linear.predictors)
  
  new_expected_y <- mean(new_y_i, na.rm = TRUE)
  old_expected_y <- mean(old_y_i, na.rm = TRUE)
  
  return(list(summary(base_model), 
              new_expected_y = new_expected_y, 
              old_expected_y = old_expected_y))
}
