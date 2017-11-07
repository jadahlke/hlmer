#' Create mixed-effects regression formula from multi-level formulas from up to 3 levels.
#'
#' @param eq_lvl1 Regression equation for level-1 effects.
#' @param eq_lvl2 List of regression equations for level-2 effects.
#' @param eq_lvl3 List of regression equations for level-3 effects.
#' @param cluster_lvl2 Name of grouping variable that defines level-2 clusters.
#' @param cluster_lvl3 Name of grouping variable that defines level-3 clusters.
#'
#' @return A list containing (1) a mixed-effects equation, (2) a full level-1 equation, and (3) a level-1 equation for random effects only.
#' @export
#'
#' @importFrom stats formula
#' @importFrom utils combn
#'
#' @examples
#' ## Level-1 equation:
#' eq_lvl1 <- "outcome ~ x1_lvl1 + x2_lvl1"
#'
#' ## Level-2 equations:
#' eq_lvl2 <- list("intercept ~ x1_lvl2 + x2_lvl2",
#'                 "x1_lvl1 ~ x1_lvl2 + x2_lvl2",
#'                 "x2_lvl1 ~ 1")
#'
#' ## Level-3 equations:
#' eq_lvl3 <- list("intercept ~ x1_lvl3",
#'                 "x1_lvl2 ~ x1_lvl3",
#'                 "x2_lvl2 ~ 1")
#'
#' ## Clustering variables:
#' cluster_lvl2 <- "cluster2"
#' cluster_lvl3 <- "cluster3"
#'
#' ## Two-level model with no level-2 explanatory variables:
#' create_regeq(eq_lvl1 = eq_lvl1, eq_lvl2 = NULL,
#'              cluster_lvl2 = cluster_lvl2)
#'
#' ## Two-level model with level-2 explanatory variables:
#' create_regeq(eq_lvl1 = eq_lvl1, eq_lvl2 = eq_lvl2,
#'              cluster_lvl2 = cluster_lvl2)
#'
#' ## Three-level model:
#' create_regeq(eq_lvl1 = eq_lvl1, eq_lvl2 = eq_lvl2, eq_lvl3 = eq_lvl3,
#'              cluster_lvl2 = cluster_lvl2, cluster_lvl3 = cluster_lvl3)
create_regeq <- function(eq_lvl1, eq_lvl2, eq_lvl3 = NULL, cluster_lvl2, cluster_lvl3 = NULL){
     if(is.null(cluster_lvl3)){
          lvl1 <- split_regeq(reg_eq = eq_lvl1)
          if(!is.null(eq_lvl2)){
               lvl2 <- clean_reglist(eq_lvl2)
          }else{
               lvl2 <- NULL
          }
          reg_eq <- paste(lvl1$y, paste(c(lvl1$x, lvl2$x), collapse = " + "), sep = " ~ ")
          if(is.null(lvl2)){
               reg_eq <- paste0(reg_eq, " + (1 | ", cluster_lvl2, ")")
          }else{
               reg_eq <- paste0(reg_eq, " + (", paste(lvl2$re, collapse = " + "), " | ", cluster_lvl2, ")")
          }
     }else{
          lvl1 <- split_regeq(reg_eq = eq_lvl1)
          if(!is.null(eq_lvl2)){
               lvl2 <- clean_reglist(eq_lvl2)
          }else{
               lvl2 <- NULL
          }
          if(!is.null(eq_lvl3)){
               lvl3 <- clean_reglist(reg_list = eq_lvl3)
          }else{
               lvl3 <- NULL
          }

          reg_eq <- paste(lvl1$y, paste(c(lvl1$x, lvl2$x, lvl3$x), collapse = " + "), sep = " ~ ")
          if(is.null(lvl3$re)){
               if(is.null(lvl2$re)){
                    reg_eq <- paste0(reg_eq, " + (1 | ", cluster_lvl2, ":", cluster_lvl3, ") + (1 | ", cluster_lvl3, ")")
               }else{
                    reg_eq <- paste0(reg_eq, " + (", paste(lvl2$re, collapse = " + "), " | ", cluster_lvl2, ":", cluster_lvl3, ") + (1 | ", cluster_lvl3, ")")
               }
          }else{
               if(is.null(lvl2$re)){
                    reg_eq <- paste0(reg_eq, " + (1 | ",
                                     cluster_lvl2, ":", cluster_lvl3, ") + (1 | ", cluster_lvl3, ")")
               }else{
                    reg_eq <- paste0(reg_eq, " + (", paste(lvl2$re, collapse = " + "), " | ",
                                     cluster_lvl2, ":", cluster_lvl3, ") + (",
                                     paste(lvl3$re, collapse = " + "), " | ", cluster_lvl3, ")")
               }
          }
     }
     if(is.null(lvl2$re)){
          lvl1_full <- paste(lvl1$y, paste(lvl1$x, collapse = " + "), sep = " ~ ")
          lvl1_re <- paste(lvl1$y, "1", sep = " ~ ")
     }else{
          lvl1_full <- paste(lvl1$y, paste(lvl1$x, collapse = " + "), sep = " ~ ")
          lvl1_re <- paste(lvl1$y, paste(lvl2$re, collapse = " + "), sep = " ~ ")
     }

     list(mixed = formula(reg_eq),
          lvl1_full = formula(lvl1_full),
          lvl1_re = formula(lvl1_re))
}


reg_terms <- function(x){
     x <- gsub(x = x, pattern = " ", replacement = "")
     x_split <- stringr::str_split(x, pattern = "\\*")

     term_list <- lapply(x_split, function(xi){
          xi <- sort(xi)

          combn_list <- list()
          for(i in 1:length(xi)){
               if(i == 1){
                    out <- as.list(1:length(xi))
               }else{
                    out <- data.frame(combn(1:length(xi), i))
                    out <- as.list(out)
               }
               combn_list <- append(combn_list, out)
          }
          as.character(unlist(lapply(combn_list, function(i) paste0(xi[i], collapse = ":"))))
     })

     term_vec <- unique(unlist(term_list))
     int_order <- stringr::str_count(term_vec, ":")
     term_dat <- dplyr::arrange(data.frame(order = int_order, term = term_vec), order, term)
     int_order <- term_dat$order
     term_vec <- as.character(term_dat$term)
     int_eq_vec <- gsub(x = term_vec, pattern = ":", replacement = "*")
     term_vec
}


split_regeq <- function(reg_eq){
     .split_regeq <- function(reg_eq){
          reg_eq <- gsub(x = reg_eq, pattern = " ", replacement = "")
          eq_split <- stringr::str_split(reg_eq, pattern = "~")[[1]]
          eq_lhs <- eq_split[1]
          eq_rhs <- eq_split[2]
          eq_rhs <- stringr::str_split(eq_rhs, pattern = "\\+")[[1]]
          list(y = gsub(x = eq_lhs, pattern = "\\*", replacement = ":"), x = reg_terms(eq_rhs))
     }

     if(!is.list(reg_eq)){
          out <- .split_regeq(reg_eq = reg_eq)
     }else{
          out <- lapply(reg_eq, .split_regeq)
          names(out) <- c(lapply(out, function(x) x[[1]]))
     }
     out
}


clean_reglist <- function(reg_list){
     mix_term <- NULL
     out <- split_regeq(reg_eq = reg_list)
     for(i in names(out)){
          if(length(out[[i]]$x) == 1){
               if(out[[i]]$x[1] == "0"){
                    out[[i]] <- NULL
               }else if(out[[i]]$x[1] == "1"){
                    out[[i]]$x <- NULL
               }
          }
          if(!is.null(out[[i]]$x))
               if(out[[i]]$y == "intercept"){
                    mix_term[[i]] <- out[[i]]$x
               }else{
                    if(length(out[[i]]$x) == 0){
                         mix_term[[i]] <- out[[i]]$y
                    }else{
                         mix_term[[i]] <- paste(out[[i]]$y, out[[i]]$x, sep = ":")
                    }
               }
     }
     y <- c(names(out))
     if(all(y != "intercept")) y <- c("0", y)
     if(any(y == "intercept")) y[which(y == "intercept")] <- "1"
     list(re = y, x = unique(as.character(unlist(mix_term))))
}


split_reglist <- function(reg_list){
     lapply(reg_list, function(x){
          y_out <- split_regeq(reg_eq = x)[[1]]
          x_out <- reg_terms(split_regeq(reg_eq = x)[[2]])
          out <- list()
          for(i in 1:length(y_out)){
               y_out[[i]] <- gsub(x = y_out[[i]], pattern = "\\*", replacement = ":")
               out[[y_out[[i]]]] <- list()
               out[[i]][["y"]] <- y_out[[i]]
               out[[i]][["x"]] <- x_out[[i]]
          }
          out
     })
}


