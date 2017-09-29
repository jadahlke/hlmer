#' lmer wrapper to produce HLM7-style output
#'
#'
#'
#' @param criterion_lvl1 Column label of \code{data} corresponding to criterion variable for level-1 observations.
#' @param cluster Column label of \code{data} corresponding to cluster/group identification variable.
#' @param predictors_lvl1 Column label of \code{data} corresponding to predictor variables for level-1 observations.
#' @param predictors_lvl2 Column label of \code{data} corresponding to predictor variables for level-2 observations.
#' @param criteria_lvl2 Optional list of level-1 statistics to be explained by level-2 predictors. Acceptable inputs are a scalar value "all" (indicates that all level 2 predictors should be used to predict all random effects)
#' or a list of vectors naming "Intercept" and/or \code{predictors_lvl1} random-effect variable names - vectors in this list must be named using \code{predictors_lvl2} variable names.
#' @param fixed_lvl1 Logical scalar or vector indicating which of the \code{predictors_lvl1} variables should be modeled as fixed effects (\code{TRUE}) or random effects (\code{FALSE}).
#' @param center_lvl1 Optional vector identifying the types of mean centering procedures to be applied to the \code{predictors_lvl1} variables. Options are "cluster" (within-cluster mean centering), "grand" (grand-mean centering), and "none" (no centering; default).
#' @param center_lvl2 Optional vector identifying the types of mean centering procedures to be applied to the \code{predictors_lvl2} variables. Options are "grand" (grand-mean centering) and "none" (no centering; default).
#' @param means_lvl1 Optional list of level-1 statistics to be explained by level-1 predictor mean values (only relevant when one or more level-1 predictors are cluster-mean centered). Acceptable formats for this argument are the same as for \code{criteria_lvl2},
#' with the exeption that lists supplied for this argument should be contain vectors named using \code{predictors_lvl1} variable names.
#' @param conf_level Confidence level to use in constructing confidence bounds around fixed effects.
#' @param model_type Numeric scalar indicating which of the following model types to run:
#' (0) the most complex model implied by the supplied arguments,
#' (1) unconditional model (i.e,. random-effects ANOVA),
#' (2) means-as-outcomes model (requrires that at least one level-2 predictor is avaialble, either from specification using the \code{predictors_lvl2} argument or from cluster-mean centering at least one level-1 predictor using the \code{center_lvl1} argument),
#' (3) random coefficients model (requires that at least one entry supplied for \code{fixed_lvl1} is \code{FALSE}),
#' (5) slopes and/or intercepts as outcomes model (requires at least one level-2 predictor or cluster-mean centered level-1 predictor and/or at least one cross-level interaction specified using \code{criteria_lvl2} and \code{means_lvl1}).
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#' @param ... Additional arugments to be passed to the \code{lme4::lmer()} function.
#'
#' @return Output from the \code{lme4::lmer()} function augmented with ICC statistics, random-effect reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
#' @export
#'
#' @import dplyr
#' @import lmerTest
#' @importFrom lme4 ranef
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @importFrom stats var
#' @importFrom stringr str_split
#'
#' @examples
#' ## Unconditional model (Raudenbush and Bryk Table 4.2):
#' hlmer(criterion_lvl1 = "MATHACH", cluster = "ID", model_type = 1, data = hsb)
#'
#' ## Means-as-outcomes model (Raudenbush and Bryk Table 4.3):
#' hlmer(criterion_lvl1 = "MATHACH", cluster = "ID", predictors_lvl2 = "MEANSES",
#' model_type = 2, data = hsb)
#'
#' ## Random-coefficients model (Raudenbush and Bryk Table 4.4):
#' hlmer(criterion_lvl1 = "MATHACH", cluster = "ID", predictors_lvl1 = "SES",
#' center_lvl1 = "cluster", model_type = 3, data = hsb)
#'
#' ## Random-coefficients model (Raudenbush and Bryk Table 4.5):
#' hlmer(criterion_lvl1 = "MATHACH", cluster = "ID", predictors_lvl1 = "SES",
#' predictors_lvl2 = "SECTOR", center_lvl1 = "cluster", criteria_lvl2 = "all",
#' means_lvl1 = "all", model_type = 4, data = hsb)
hlmer <- function(criterion_lvl1, cluster, predictors_lvl1 = NULL, predictors_lvl2 = NULL, criteria_lvl2 = "all",
                  fixed_lvl1 = FALSE, center_lvl1 = NULL, center_lvl2 = FALSE, means_lvl1 = NULL, conf_level = .95, model_type = 0, data, ...){
     call <- match.call()

     lmer_eq_lvl1 <- lmer_eq_lvl2 <- list()

     data <- data.frame(data)
     use_cols <- c(cluster, criterion_lvl1, predictors_lvl1, predictors_lvl2)
     data <- na.omit(data[,use_cols])

     eliminate_novariance <- function(x){
          if(any(zapsmall(apply(x[,predictors_lvl1], 2, var)) == 0)){
               x[0,]
          }else{
               x
          }
     }
     if(!is.null(predictors_lvl1)){
          orig_clusters <- unlist(unique(data[,cluster]))
          data <- data %>% group_by_(.dots = cluster) %>% do(eliminate_novariance(x = .))
          data <- data.frame(data)
          new_clusters <- unlist(unique(data[,cluster]))

          if(length(new_clusters) < length(orig_clusters))
               warning(length(orig_clusters) - length(new_clusters), " clusters had no variance on one or more variables and have been removed")
     }
     data[,cluster] <- factor(data[,cluster])

     if(!is.null(criteria_lvl2)){
          if(!is.list(criteria_lvl2) & length(criteria_lvl2) == 1){
               if(criteria_lvl2[1] == "all"){
                    criteria_lvl2 <- list()
                    for(i in predictors_lvl2) criteria_lvl2[[i]] <- c("Intercept", predictors_lvl1[!fixed_lvl1])
               }
          }else{
               criteria_lvl2 <- lapply(criteria_lvl2, function(x){
                    if(length(x) == 1  & x[1] == "all"){
                         c("Intercept", predictors_lvl1[!fixed_lvl1])
                    }else{
                         x
                    }
               })
               for(i in names(criteria_lvl2)){
                    if(length(criteria_lvl2[[i]]) == 1)
                         if(criteria_lvl2[[i]] == "none" | is.na(criteria_lvl2[[i]]))
                              criteria_lvl2[[i]] <- NULL
               }
               if(length(criteria_lvl2) == 0) criteria_lvl2 <- NULL
          }
     }


     if(!is.null(predictors_lvl1)){
          if(!is.null(center_lvl1)){
               lvl1_cwc <- lvl1_cgm <- lvl1_means <- NULL
               for(i in 1:length(center_lvl1)){
                    if(center_lvl1[i] == "cluster"){
                         lvl1_cwc <- cbind(lvl1_cwc, center_within(x = data[,predictors_lvl1[i]], cluster = data[,cluster]))
                         lvl1_means <- cbind(lvl1_means, cluster_means(x = data[,predictors_lvl1[i]], cluster = data[,cluster]))
                         lmer_eq_lvl1[[i]] <- c(orig = predictors_lvl1[i],
                                                lvl1 = paste0(predictors_lvl1[i], "_cwc"),
                                                means = paste0(predictors_lvl1[i], "_means"))
                    }
                    if(center_lvl1[i] == "grand"){
                         lvl1_cgm <- cbind(lvl1_cgm, data[,predictors_lvl1[i]] - mean(data[,predictors_lvl1[i]], na.rm = TRUE))
                         lmer_eq_lvl1[[i]] <- c(orig = predictors_lvl1[i],
                                                lvl1 = paste0(predictors_lvl1[i], "_cgm"),
                                                means = NA)
                    }
                    if(center_lvl1[i] == "none"){
                         lmer_eq_lvl1[[i]] <- c(orig = predictors_lvl1[i],
                                                lvl1 = predictors_lvl1[i],
                                                means = NA)
                    }
               }
               if(any(center_lvl1 == "cluster")){
                    colnames(lvl1_cwc) <- paste0(predictors_lvl1[center_lvl1 == "cluster"], "_cwc")
                    colnames(lvl1_means) <- paste0(predictors_lvl1[center_lvl1 == "cluster"], "_means")
                    data <- cbind(data, lvl1_cwc, lvl1_means)
               }
               if(any(center_lvl1 == "grand")){
                    colnames(lvl1_cgm) <- paste0(predictors_lvl1[center_lvl1 == "grand"], "_cgm")
                    data <- cbind(data, lvl1_cgm)
               }
               lmer_eq_lvl1 <- t(simplify2array(lmer_eq_lvl1))
          }else{
               lmer_eq_lvl1 <- cbind(orig = predictors_lvl1,
                                     lvl1 = predictors_lvl1,
                                     means = NA)
          }

          fixed_effects <- fe_lvl1 <- paste(lmer_eq_lvl1[,"lvl1"], collapse = " + ")
     }else{
          lmer_eq_lvl1 <- fixed_effects <- fe_lvl1 <- NULL
     }

     if(!is.null(predictors_lvl2)){
          if(any(center_lvl2)){
               lvl2_cgm <- NULL
               for(i in 1:length(center_lvl2)){
                    if(center_lvl2[i]){
                         lvl2_cgm <- cbind(lvl2_cgm, data[,predictors_lvl2[i]] - mean(data[,predictors_lvl2[i]], na.rm = TRUE))
                         lmer_eq_lvl2[[i]] <- c(orig = predictors_lvl2[i],
                                                lvl2 = paste0(predictors_lvl2[i], "_cgm"))
                    }else{
                         lmer_eq_lvl2[[i]] <- c(orig = predictors_lvl2[i],
                                                lvl2 = predictors_lvl2[i])
                    }
               }
               colnames(lvl2_cgm) <- paste0(predictors_lvl2[center_lvl2], "_cgm")
               data <- cbind(data, lvl2_cgm)
               lmer_eq_lvl2 <- t(simplify2array(lmer_eq_lvl2))
          }else{
               lmer_eq_lvl2 <- cbind(orig = predictors_lvl2,
                                     lvl2 = predictors_lvl2)
          }
     }else{
          lmer_eq_lvl2 <- NULL
     }

     if(!is.null(lmer_eq_lvl1) & !is.null(center_lvl1)){
          if(any(center_lvl1 == "cluster")){
               lmer_eq_lvl2 <- rbind(lmer_eq_lvl2,
                                     cbind(orig = lmer_eq_lvl1[center_lvl1 == "cluster", "means"],
                                           lvl2 = lmer_eq_lvl1[center_lvl1 == "cluster", "means"]))

          }
          if(!is.null(means_lvl1)){
               if(!is.list(means_lvl1) & length(means_lvl1) == 1){
                    if(means_lvl1[1] == "all"){
                         means_lvl1 <- list()
                         for(i in predictors_lvl1[center_lvl1 == "cluster"]) means_lvl1[[paste0(i, "_means")]] <- c("Intercept", predictors_lvl1[!fixed_lvl1])
                    }
               }else{
                    means_lvl1 <- lapply(means_lvl1, function(x){
                         if(length(x) == 1  & x[1] == "all"){
                              c("Intercept", predictors_lvl1[!fixed_lvl1])
                         }else{
                              x
                         }
                    })
                    for(i in names(means_lvl1)){
                         if(length(means_lvl1[[i]]) == 1)
                              if(means_lvl1[[i]] == "none" | is.na(means_lvl1[[i]]))
                                   means_lvl1[[i]] <- NULL
                    }
                    if(length(means_lvl1) == 0) means_lvl1 <- NULL
               }
               if(!is.null(means_lvl1)) criteria_lvl2 <- append(criteria_lvl2, means_lvl1)
          }
     }

     if(any(!fixed_lvl1) & !is.null(predictors_lvl1)){
          random_effects <- paste(lmer_eq_lvl1[!fixed_lvl1,"lvl1"], collapse = " + ")
          random_pred_lvl1 <- lmer_eq_lvl1[!fixed_lvl1,"lvl1"]
     }else{
          random_effects <- "1"
          random_pred_lvl1 <- NULL
     }

     if(!is.null(criteria_lvl2)){
          use_as_u0pred <- lapply(criteria_lvl2, function(x) any(x == "Intercept"))
          criteria_lvl2 <- lapply(criteria_lvl2, function(x) x[x != "Intercept"])

          if(!is.null(lmer_eq_lvl1) & !is.null(lmer_eq_lvl2)){
               if(!is.list(criteria_lvl2) & criteria_lvl2[1] == "all"){
                    xlvl_ints <- expand.grid(lmer_eq_lvl2[,2], lmer_eq_lvl1[,2])
                    xlvl_ints <- t(apply(xlvl_ints, 1, function(x) c(lvl2 = as.character(x[1]), lvl1 = as.character(x[2]))))
                    xlvl_ints <- apply(xlvl_ints, 1, function(x) paste(x, collapse = " * "))
                    xlvl_ints <- paste(xlvl_ints, collapse = " + ")
               }else{
                    .lvl2 <- names(criteria_lvl2)
                    rownames(lmer_eq_lvl2) <- lmer_eq_lvl2[,1]
                    names(criteria_lvl2) <- lmer_eq_lvl2[.lvl2,2]

                    rownames(lmer_eq_lvl1) <- lmer_eq_lvl1[,1]
                    criteria_lvl2 <- lapply(criteria_lvl2, function(x){lmer_eq_lvl1[x,2]})

                    xlvl_ints <- NULL
                    for(i in names(criteria_lvl2)) xlvl_ints <- c(xlvl_ints, paste(i, criteria_lvl2[[i]], sep = " * ", collapse = " + "))
                    xlvl_ints <- paste(xlvl_ints, collapse = " + ")
               }
          }else{
               xlvl_ints <- NULL
          }

          if(any(unlist(use_as_u0pred))){
               means_as_crit <- NULL
               for(i in names(criteria_lvl2)){
                    if(any(names(use_as_u0pred) == i)){
                         means_as_crit[i] <- use_as_u0pred[[i]]
                    }else{
                         means_as_crit[i] <- FALSE
                    }
                    fe_lvl2 <- lmer_eq_lvl2[,2][lmer_eq_lvl2[,"orig"] %in% names(means_as_crit)][means_as_crit]
               }
               fe_lvl2 <- paste(fe_lvl2, collapse = " + ")
               if(is.null(fixed_effects)){
                    fixed_effects <- fe_lvl2
               }else{
                    fixed_effects <- paste(fixed_effects, fe_lvl2, sep = " + ", collapse = " + ")
               }
          }else{
               fe_lvl2 <- NULL
          }
     }else{
          fe_lvl2 <- xlvl_ints <- NULL
     }

     if(is.null(fixed_effects)) fixed_effects <- "1"

     if(!is.null(lmer_eq_lvl1)) use_preds_lvl1 <- lmer_eq_lvl1[,"lvl1"]
     if(!is.null(lmer_eq_lvl2)) use_preds_lvl2 <- lmer_eq_lvl2[,"lvl2"]

     usable_mods <- c(null = TRUE,
                      meansout = !is.null(fe_lvl2),
                      randcoeff = !is.null(predictors_lvl1) & any(!fixed_lvl1),
                      intslopesout = !is.null(xlvl_ints))

     if(model_type == 0){
          maxmod <- max(which(usable_mods))
          usable_mods[-maxmod] <- FALSE
     }else{
          if(length(model_type) > 1) stop("Only one 'model_type' can be supplied", call. = FALSE)
          if(model_type > 4 | round(model_type) != model_type) stop("'model_type' must be an integer between 0 and 4", call. = FALSE)
          if(all(which(usable_mods) != model_type)) stop("The selected 'model_type' is not possible using the supplied data", call. = FALSE)
          usable_mods[-model_type] <- FALSE
     }

     if(usable_mods[1]){
          eq <- as.formula(paste0(criterion_lvl1, " ~ (1 | ", cluster, ")"))
     }else if(usable_mods[2]){
          eq <- as.formula(paste0(criterion_lvl1, " ~ ", fe_lvl2, " + (1 | ", cluster, ")"))
     }else if(usable_mods[3]){
          eq <- as.formula(paste0(criterion_lvl1, " ~ ", random_effects, " + (", random_effects," | ", cluster, ")"))
     }else if(usable_mods[4]){
          eq <- as.formula(paste0(criterion_lvl1, " ~ ", fixed_effects, " + ", xlvl_ints, " + (", random_effects, " | ", cluster, ")"))
     }

     mod <- lmerTest::lmer(eq, data = data, ...)
     if(usable_mods[1]){
          sum <- suppressMessages(summary(mod))
     }else{
          sum <- summary(mod)
     }
     tau <- attributes(sum$varcor[[1]])$stddev^2
     icc <- tau / (sum$sigma^2 + tau)
     ci <- cbind(`CI (Lower)` = sum$coefficients[,1] - qnorm((1 - conf_level) / 2, lower.tail = FALSE) * sum$coefficients[,2],
                 `CI (Upper)` = sum$coefficients[,1] + qnorm((1 - conf_level) / 2, lower.tail = FALSE) * sum$coefficients[,2])
     if(nrow(ci) == 1) rownames(ci) <- "(Intercept)"

     if(usable_mods[1]){
          chisq <- chisq_hlmer(model = mod, summary = sum, criterion_lvl1 = criterion_lvl1, cluster = cluster, data = data)
          rel <- rel_lvl1(lmer_mod = sum, cluster = cluster, data = data)
     }else if(usable_mods[2]){
          chisq <- chisq_hlmer(model = mod, summary = sum, criterion_lvl1 = criterion_lvl1,
                               predictors_lvl1 = NULL, predictors_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
          rel <- rel_lvl1(lmer_mod = sum, cluster = cluster, data = data)
     }else  if(usable_mods[3]){
          rel <- rel_lvl1(lmer_mod = sum, cluster = cluster, predictors_lvl1 = random_pred_lvl1, data = data)
          chisq <- chisq_hlmer(model = mod, summary = sum, criterion_lvl1 = criterion_lvl1,
                               predictors_lvl1 = use_preds_lvl1, predictors_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
     }else if(usable_mods[4]){
          if(any(!fixed_lvl1)){
               rel <- rel_lvl1(lmer_mod = sum, cluster = cluster, predictors_lvl1 = random_pred_lvl1, data = data)
          }else{
               rel <- rel_lvl1(lmer_mod = sum, cluster = cluster, data = data)
          }
          chisq <- chisq_hlmer(model = mod, summary = sum, criterion_lvl1 = criterion_lvl1,
                               predictors_lvl1 = use_preds_lvl1, predictors_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
     }

     out <- list(call = call,
                 model = mod,
                 summary = sum,
                 conf = ci,
                 reliability = rel,
                 icc = icc,
                 chisq_tau = chisq)

     class(out) <- "hlmerMod"
     out
}

print.hlmerMod <- function(x, ..., digits = 5){
     cat("Call: \n")
     print(x$call)

     cat("\n")
     print(x$summary)

     cat("\n")
     conf_level <- x$call$conf_level
     if(is.null(conf_level)) conf_level <- .95
     cat(paste0(round(conf_level * 100), "%"), "confidence intervals for fixed effects: \n")
     print(x$conf, digits = digits)

     cat("\n")
     cat("Reliability estimates for random effects: \n")
     print(x$reliability, digits = digits)

     cat("\n")
     cat("Intraclass correlation coefficients (ICCs) for random effects: \n")
     print(x$icc, digits = digits)

     cat("\n")
     cat("Approximate chi-square tests for random-effects variances: \n")
     print(x$chisq_tau, digits = digits)
}
