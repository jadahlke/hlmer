#' Center a vector of numeric values using cluster means.
#'
#' @param x Numeric vector.
#' @param cluster Vector of cluster identifiers.
#'
#' @return A cluster-mean centered vector.
#' @export
center_within <- function(x, cluster){
     x_center <- tapply(x, cluster, function(x) x - mean(x))
     cluster_lvls <- levels(factor(cluster))
     for(i in cluster_lvls) x[cluster == i] <- x_center[[i]]
     x
}

#' Obtain cluster means from a numeric vector.
#'
#' @param x Numeric vector.
#' @param cluster Vector of cluster identifiers.
#'
#' @return A vector of cluster means.
cluster_means <- function(x, cluster){
     x_center <- tapply(x, cluster, function(x) mean(x))
     cluster_lvls <- levels(factor(cluster))
     for(i in cluster_lvls) x[cluster == i] <- x_center[[i]]
     x
}

rel_lvl1 <- function(lmer_mod, data, predictors_lvl1 = NULL, cluster){
     se_mat <- se_lvl1(sigma = lmer_mod$sigma, data = data, predictors_lvl1 = predictors_lvl1, cluster = cluster)
     tao_vec <- attributes(lmer_mod$varcor[[1]])$stddev^2
     tao_mat <- matrix(tao_vec, nrow(se_mat), ncol(se_mat), TRUE)
     apply(tao_mat / (tao_mat + se_mat^2), 2, mean)
}

se_lvl1 <- function(sigma, data, predictors_lvl1 = NULL, cluster){
     if(!is.null(predictors_lvl1)){
          out <- t(simplify2array(by(data, data[,cluster], function(x){
               x1 <- cbind(`(Intercept)` = 1, as.matrix(x[,predictors_lvl1]))
               invmat <- solve(t(x1) %*% x1)
               (diag(sigma^2 * invmat))^.5
          })))
     }else{
          out <- cbind(`(Intercept)` = (sigma^2 / by(data, data[,cluster], nrow))^.5)
     }
     out
}



