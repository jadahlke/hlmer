#' lmer wrapper to produce HLM7-style output
#'
#' This function is a wrapper for \code{lmer()} to formulate a mixed-effects linear model based on user-supplied variables.
#' In addition to the output from \code{lmer()}, \code{hlmer()} also provides estimates of ICC statistics, random-effect reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
#' A key feature of \code{hlmer()} is that it can automate the process of centering predictors and creating variables from the cluster means of level-1 predictors to be included in level 2 of the model.
#' The \code{model_type} argument can be used to specify one of four pre-formatted types of models (using arguments of 1, 2, 3, or 4) or to compute the complete model implied by other arguments (using a 0 argument).
#' Variables that are cluster-mean centered within this function will be labeled with a "cwc" suffix (abbreviation for centered within cluster) and
#' variables that are grand-mean centered within this function will be labeled with a "cgm" suffix (abbreviation for centered grand mean).
#'
#' @param y_lvl1 Column label of \code{data} corresponding to the criterion variable for level-1 observations.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#' @param x_lvl1 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-1 observations.
#' @param x_lvl2 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-2 observations.
#' @param y_lvl2 Optional list of level-1 statistics to be explained by level-2 predictors. Acceptable inputs are a scalar value of "all" (indicates that all level 2 predictors should be used to predict all random effects),
#' "intercepts" (indicates that all level 2 predictors should be used to predict random intercepts),
#' "slopes" (indicates that all level 2 predictors should be used to predict all random slopes),
#' or a list of vectors naming "Intercept" and/or \code{x_lvl1} random-effect variable names - vectors in this list must be named using \code{x_lvl2} variable names.
#' @param fixed_lvl1 Logical scalar or vector indicating which of the \code{x_lvl1} variables should be modeled as fixed effects (\code{TRUE}) or random effects (\code{FALSE}; default).
#' @param center_lvl1 Optional vector identifying the types of mean centering procedures to be applied to the \code{x_lvl1} variables. Options are "cluster" (within-cluster mean centering), "grand" (grand-mean centering), and "none" (no centering; default).
#' @param center_lvl2 Vector or scalar identifying the types of mean centering procedures to be applied to the \code{x_lvl2} variables. Options are \code{TRUE} (grand-mean centering) or \code{FALSE} (no centering; default).
#' @param y_lvl1means Optional list of level-1 statistics to be explained by level-1 predictor mean values (only relevant when one or more level-1 predictors are cluster-mean centered). Acceptable formats for this argument are the same as for \code{y_lvl2},
#' with the exeption that lists supplied for this argument should contain vectors named using \code{x_lvl1} variable names.
#' @param center_lvl1means Vector or scalar identifying the types of mean centering procedures to be applied to the means of \code{x_lvl1} variables. Options are \code{TRUE} (grand-mean centering) or \code{FALSE} (no centering; default).
#' If this argument is supplied as a vector, it should have as many elements as there are cluster-mean centered level-1 predictors.
#' @param conf_level Confidence level to use in constructing confidence bounds around fixed effects.
#' @param cred_level Credibility level to use in constructing credibility bounds (ranges of plausible values for random coefficients) around random effects.
#' @param model_type Numeric scalar indicating which of the following model types to run:
#' (0) the complete model implied by the supplied arguments (default),
#' (1) unconditional model (i.e,. random-effects ANOVA),
#' (2) means-as-outcomes model (requrires that at least one level-2 predictor is available, either from specification using the \code{x_lvl2} argument or from cluster-mean centering at least one level-1 predictor using the \code{center_lvl1} argument),
#' (3) random coefficients model (requires that at least one entry supplied for \code{fixed_lvl1} is \code{FALSE}),
#' (5) slopes and/or intercepts as outcomes model (requires at least one level-2 predictor or cluster-mean centered level-1 predictor and/or at least one cross-level interaction specified using \code{y_lvl2} and \code{y_lvl1means}).
#' @param remove_missing Logical scalar that determines whether cases with missing data should be omitted.
#' @param check_lvl1_variance Logical scalar that determines whether clusters with no variance in level-1 predictors should be screened out.
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#' @param ... Additional arugments to be passed to the \code{lme4::lmer()} function.
#'
#' @return Output from the \code{lmerTest::lmer()} function augmented with ICC statistics, random-effects reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
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
#' \dontrun{
#' ## Unconditional model (Raudenbush and Bryk Table 4.2):
#' hlmer(y_lvl1 = "MATHACH", cluster = "ID", model_type = 1, data = hsb)
#'
#' ## Means-as-outcomes model (Raudenbush and Bryk Table 4.3):
#' hlmer(y_lvl1 = "MATHACH", cluster = "ID", x_lvl2 = "MEANSES",
#' y_lvl2 = "intercepts", model_type = 2, data = hsb)
#'
#' ## Random-coefficients model (Raudenbush and Bryk Table 4.4):
#' hlmer(y_lvl1 = "MATHACH", cluster = "ID", x_lvl1 = "SES",
#' center_lvl1 = "cluster", model_type = 3, data = hsb)
#'
#' ## Slopes and intercepts as outcomes model (Raudenbush and Bryk Table 4.5):
#' hlmer(y_lvl1 = "MATHACH", cluster = "ID", x_lvl1 = "SES",
#' x_lvl2 = "SECTOR", center_lvl1 = "cluster", y_lvl2 = "all",
#' y_lvl1means = "all", model_type = 4, data = hsb)
#'
#' ## For a more complex model, we can specify which level-2 predictors
#' ## should be used to predict which level-1 random effects.
#' ## In this TIMSS example, the level-1 predictor is "self_efy" and
#' ## the level-2 predictors are "students" and "alg." If we want to
#' ## predict level-1 intercepts using both level-2 predictors, but we
#' ## only want to predict level-1 "self_efy" slopes using "students,"
#' ## we can specify that by passing a list to "y_lvl2" of the form:
#' ##
#' ## y_lvl2 = list(students = c("Intercept", "self_efy"),
#' ##               alg = "Intercept")
#' ##
#' ## This tells the program exactly which random effects to explain:
#' hlmer(y_lvl1 = "scores", cluster = "idteach",
#'       x_lvl1 = "self_efy", x_lvl2 = c("students", "alg"),
#'       center_lvl1 = "cluster", center_lvl2 = TRUE,
#'       y_lvl2 = list(students = c("Intercept", "self_efy"),
#'                     alg = "Intercept"),
#'       y_lvl1means = "intercepts", data = timss)
#' }
hlmer <- function(y_lvl1, cluster, x_lvl1 = NULL, x_lvl2 = NULL, y_lvl2 = NULL,
                  fixed_lvl1 = FALSE, center_lvl1 = NULL, center_lvl2 = FALSE,
                  y_lvl1means = NULL, center_lvl1means = FALSE, conf_level = .95, cred_level = .95,
                  model_type = 0, remove_missing = TRUE, check_lvl1_variance = TRUE, data, ...){
     .hlmer(y_lvl1 = y_lvl1, cluster = cluster, x_lvl1 = x_lvl1, x_lvl2 = x_lvl2, y_lvl2 = y_lvl2,
            fixed_lvl1 = fixed_lvl1, center_lvl1 = center_lvl1, center_lvl2 = center_lvl2,
            y_lvl1means = y_lvl1means, center_lvl1means = center_lvl1means, conf_level = conf_level, cred_level = cred_level,
            model_type = model_type, remove_missing = remove_missing, check_lvl1_variance = check_lvl1_variance, data, ...)
}


#' lmer wrapper to produce HLM7-style output from multi-level equations
#'
#' This function is a wrapper for \code{lmer()} to formulate a mixed-effects linear model based on user-supplied multi-level equations.
#' In addition to the output from \code{lmer()}, \code{hlmer2()} also provides estimates of ICC statistics, random-effect reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
#' To ensure that a simplified model (e.g., an unconditional model) run with this function will be based on the same cases as a more complex model, supply a \code{model_type} argument using the values outlined in the documentation for the \code{hlmer()} function.
#' IMPORTANT: Unlike the \code{hlmer()} function, this function requires that the user manually center predictors. To automate the centering process, see the \code{center_data()} function.
#'
#' @param conf_level Confidence level to use in constructing confidence bounds around fixed effects.
#' @param eq_lvl1 Regression equation (as a string) for level-1 models. Do not use parentheses in the model. When supplying an interaction between level-1 variables, use either "*" or ":" notation and manually supply the main effects. In the output, all within-level interactions will be indicated by a "_X_" symbol between the variables.
#' @param eq_lvl2 List of regression equations (as strings) for level-2 models, with outcomes labeled as "intercept" or as level-1 predictor names. If left NULL, random intercepts will be estimated, but all slopes will be treated as fixed effects.
#' To specify a fixed effect, simply omit the variable from the \code{eq_lvl2} list or provide a regression of the following form: "predictor ~ 0", where "predictor" is the name of a level-1 explanatory variable.
#' To specify a random effect with no level-2 explanatory variable, provide a regression of the following form: "predictor ~ 1".
#' When supplying an interaction between level-2 variables, use either "*" or ":" notation and manually supply the main effects. In the output, all within-level interactions will be indicated by a "_X_" symbol between the variables.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#' @param cred_level Credibility level to use in constructing credibility bounds (ranges of plausible values for random coefficients) around random effects.
#' @param remove_missing Logical scalar that determines whether cases with missing data should be omitted.
#' @param check_lvl1_variance Logical scalar that determines whether clusters with no variance in level-1 predictors should be screened out.
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#' @param ... Additional arugments to be passed to the \code{lme4::lmer()} function.
#'
#' @return Output from the \code{lmerTest::lmer()} function augmented with ICC statistics, random-effects reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
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
#' \dontrun{
#' ## Center the HSB database:
#' dat <- center_data(cluster = "ID", data = hsb)
#'
#' ## Unconditional model (Raudenbush and Bryk Table 4.2):
#' hlmer2(eq_lvl1 = "MATHACH ~ 1", cluster = "ID", data = dat)
#'
#' ## Means-as-outcomes model (Raudenbush and Bryk Table 4.3):
#' hlmer2(eq_lvl1 = "MATHACH ~ 1", eq_lvl2 = "intercept ~ MEANSES",
#' cluster = "ID", data = dat)
#'
#' ## Random-coefficients model (Raudenbush and Bryk Table 4.4):
#' hlmer2(eq_lvl1 = "MATHACH ~ SES_cwc",
#' eq_lvl2 = list("intercept ~ MEANSES", "SES_cwc ~ 1"),
#' cluster = "ID", data = dat)
#'
#' ## Slopes and intercepts as outcomes model (Raudenbush and Bryk Table 4.5):
#' hlmer2(eq_lvl1 = "MATHACH ~ SES_cwc",
#' eq_lvl2 = list("intercept ~ SECTOR + MEANSES", "SES_cwc ~ SECTOR + MEANSES"),
#' cluster = "ID", data = dat)
#'
#' ## For a more complex model, we can specify which level-2 predictors
#' ## should be used to predict which level-1 random effects.
#' ## In this TIMSS example, the level-1 predictor is "self_efy" and
#' ## the level-2 predictors are "students" and "alg." If we want to
#' ## predict level-1 intercepts using both level-2 predictors, but we
#' ## only want to predict level-1 "self_efy" slopes using "students,"
#' ## we can specify that with the following model:
#' ##
#' hlmer2(eq_lvl1 = "scores ~ self_efy",
#' eq_lvl2 = list("intercept ~ alg + students", "self_efy ~ students"),
#' cluster = "idteach", data = timss)
#' }
hlmer2 <- function(eq_lvl1, eq_lvl2 = NULL, cluster,
                   conf_level = 0.95, cred_level = 0.95,
                   remove_missing = TRUE, check_lvl1_variance = TRUE, data, ...){

     if(is.null(eq_lvl2)) eq_lvl2 <- "intercept ~ 1"

     eq_lvl1 <- gsub(x = eq_lvl1, pattern = " ", replacement = "")
     if(is.list(eq_lvl2)){
          eq_lvl2 <- lapply(eq_lvl2, function(x) gsub(x = x, pattern = " ", replacement = ""))
     }else{
          eq_lvl2 <- gsub(x = eq_lvl2, pattern = " ", replacement = "")
     }

     eq_lvl1_split <- stringr::str_split(eq_lvl1, pattern = "~")[[1]]
     eq_lvl2_split <- stringr::str_split(eq_lvl2, pattern = "~")

     eq_lvl1_lhs <- eq_lvl1_split[1]
     eq_lvl1_rhs <- eq_lvl1_split[2]
     eq_lvl1_rhs_split <- stringr::str_split(eq_lvl1_rhs, pattern = "\\+")[[1]]

     interaction_lvl1 <- grepl(pattern = "[*]", x = eq_lvl1_rhs_split) | grepl(pattern = "[:]", x = eq_lvl1_rhs_split)
     interaction_lvl1 <- eq_lvl1_rhs_split[interaction_lvl1]
     if(length(interaction_lvl1) > 0){
          interaction_lvl1 <- gsub(pattern = "[:]", x = interaction_lvl1, replacement = "*")
          eq_lvl1_rhs_split <- gsub(pattern = "[:]", x = eq_lvl1_rhs_split, replacement = "*")
          eq_lvl2_split <- lapply(eq_lvl2_split, function(x){c(gsub(pattern = "[:]", x = x[1], replacement = "*"), x[-1])})

          interaction_lvl1_label <- gsub(pattern = "[*]", x = interaction_lvl1, replacement = "_X_")
          eq_lvl1_rhs_split <- gsub(pattern = "[*]", x = eq_lvl1_rhs_split, replacement = "_X_")
          eq_lvl2_split <- lapply(eq_lvl2_split, function(x){c(gsub(pattern = "[*]", x = x[1], replacement = "_X_"), x[-1])})

          data_int <- list()
          for(i in 1:length(interaction_lvl1))
               data_int[[interaction_lvl1_label[i]]] <- eval(parse(text = interaction_lvl1[i]), data)

          data <- cbind(data, data.frame(data_int))
     }

     eq_lvl2_lhs <- lapply(eq_lvl2_split, function(x) x[1])
     eq_lvl2_rhs <- lapply(eq_lvl2_split, function(x) x[2])
     eq_lvl2_rhs_split <- lapply(eq_lvl2_rhs, function(x) stringr::str_split(x, pattern = "\\+")[[1]])
     names(eq_lvl2_rhs_split) <- unlist(eq_lvl2_lhs)

     y_lvl1 <- eq_lvl1_lhs
     x_lvl1 <- eq_lvl1_rhs_split

     x_lvl2 <- unique(unlist(eq_lvl2_rhs_split))
     x_lvl2 <- x_lvl2[x_lvl2 != "1"]
     x_lvl2 <- x_lvl2[x_lvl2 != "0"]
     for(i in names(eq_lvl2_rhs_split)) if(eq_lvl2_rhs_split[[i]][1] == "0") eq_lvl2_rhs_split[[i]] <- NULL

     interaction_lvl2 <- grepl(pattern = "[*]", x = x_lvl2)
     interaction_lvl2 <- x_lvl2[interaction_lvl2]
     if(length(interaction_lvl2) > 0){
          interaction_lvl2 <- gsub(pattern = "[:]", x = interaction_lvl2, replacement = "*")
          names(eq_lvl2_rhs_split) <- gsub(pattern = "[:]", x = names(eq_lvl2_rhs_split), replacement = "*")
          eq_lvl2_lhs <- lapply(eq_lvl2_lhs, function(x){gsub(pattern = "[:]", x = x, replacement = "*")})
          eq_lvl2_rhs_split <- lapply(eq_lvl2_rhs_split, function(x){gsub(pattern = "[:]", x = x, replacement = "*")})
          x_lvl2 <- gsub(pattern = "[:]", x = x_lvl2, replacement = "*")

          interaction_lvl2_label <- gsub(pattern = "[*]", x = interaction_lvl2, replacement = "_X_")
          names(eq_lvl2_rhs_split) <- gsub(pattern = "[*]", x = names(eq_lvl2_rhs_split), replacement = "_X_")
          eq_lvl2_lhs <- lapply(eq_lvl2_lhs, function(x){gsub(pattern = "[*]", x = x, replacement = "_X_")})
          eq_lvl2_rhs_split <- lapply(eq_lvl2_rhs_split, function(x){gsub(pattern = "[*]", x = x, replacement = "_X_")})
          x_lvl2 <- gsub(pattern = "[*]", x = x_lvl2, replacement = "_X_")

          data_int <- list()
          for(i in 1:length(interaction_lvl2))
               data_int[[interaction_lvl2_label[i]]] <- eval(parse(text = interaction_lvl2[i]), data)

          data <- cbind(data, data.frame(data_int))
     }

     y_lvl2 <- as.list(x_lvl2)
     names(y_lvl2) <- x_lvl2
     for(i in x_lvl2){
          y_lvl2[[i]] <- names(eq_lvl2_rhs_split)[unlist(lapply(eq_lvl2_rhs_split, function(x) any(x == i)))]
     }
     y_lvl2 <- lapply(y_lvl2, function(x){
          x_lo <- tolower(x)
          if(any(x_lo == "intercept")){
               x[x_lo == "intercept"] <- "Intercept"
          }
          x
     })

     fixed_lvl1 <- !(eq_lvl1_rhs_split %in% names(eq_lvl2_rhs_split))

     if(eq_lvl1_rhs_split[1] == "1"){
          x_lvl1 <- NULL
          fixed_lvl1 <- FALSE
     }
     if(length(x_lvl2) == 0) x_lvl2 <- NULL
     if(length(y_lvl2) == 0) y_lvl2 <- NULL

     out <- .hlmer(y_lvl1 = y_lvl1, cluster = cluster, x_lvl1 = x_lvl1, x_lvl2 = x_lvl2, y_lvl2 = y_lvl2,
                   fixed_lvl1 = fixed_lvl1, conf_level = conf_level, cred_level = cred_level,
                   center_lvl1 = NULL, center_lvl2 = FALSE, y_lvl1means = NULL, center_lvl1means = NULL,
                   remove_missing = remove_missing, check_lvl1_variance = check_lvl1_variance, data = data, ...)
     out$call <- match.call()
     out
}


.hlmer <- function(y_lvl1, cluster, x_lvl1 = NULL, x_lvl2 = NULL, y_lvl2 = NULL,
                  fixed_lvl1 = FALSE, center_lvl1 = NULL, center_lvl2 = FALSE,
                  y_lvl1means = NULL, center_lvl1means = FALSE, conf_level = .95, cred_level = .95,
                  model_type = 0, remove_missing = TRUE, check_lvl1_variance = TRUE, data, ...){
     call <- match.call()

     lmer_eq_lvl1 <- lmer_eq_lvl2 <- list()

     data <- data.frame(data)
     use_cols <- c(cluster, y_lvl1, x_lvl1, x_lvl2)
     if(remove_missing) data <- na.omit(data[,use_cols])

     if(check_lvl1_variance & !is.null(x_lvl1)){
          eliminate_novariance <- function(x){
               if(any(zapsmall(apply(x[,x_lvl1], 2, var)) == 0)){
                    x[0,]
               }else{
                    x
               }
          }

          orig_clusters <- unlist(unique(data[,cluster]))
          data <- data %>% group_by_(.dots = cluster) %>% do(eliminate_novariance(x = .))
          data <- data.frame(data)
          new_clusters <- unlist(unique(data[,cluster]))

          if(length(new_clusters) < length(orig_clusters))
               warning(length(orig_clusters) - length(new_clusters), " clusters had no variance on one or more variables and have been removed")
     }
     data[,cluster] <- factor(data[,cluster])

     if(!is.null(y_lvl2)){
          if(!is.list(y_lvl2) & length(y_lvl2) == 1){
               if(y_lvl2[1] == "all"){
                    y_lvl2 <- list()
                    for(i in x_lvl2) y_lvl2[[i]] <- c("Intercept", x_lvl1[!fixed_lvl1])
               }else if(y_lvl2[1] == "intercepts"){
                    y_lvl2 <- list()
                    for(i in x_lvl2) y_lvl2[[i]] <- c("Intercept")
               }else if(y_lvl2[1] == "slopes"){
                    y_lvl2 <- list()
                    for(i in x_lvl2) y_lvl2[[i]] <- c(x_lvl1[!fixed_lvl1])
               }
          }else{
               y_lvl2 <- lapply(y_lvl2, function(x){
                    if(length(x) == 1  & x[1] == "all"){
                         c("Intercept", x_lvl1[!fixed_lvl1])
                    }else{
                         x
                    }
               })
               for(i in names(y_lvl2)){
                    if(length(y_lvl2[[i]]) == 1)
                         if(y_lvl2[[i]] == "none" | is.na(y_lvl2[[i]]))
                              y_lvl2[[i]] <- NULL
               }
               if(length(y_lvl2) == 0) y_lvl2 <- NULL
          }
     }


     if(!is.null(x_lvl1)){
          if(!is.null(center_lvl1)){
               if(length(x_lvl1) > 1 & length(center_lvl1) == 1) center_lvl1 <- rep(center_lvl1, length(x_lvl1))
               lvl1_cwc <- lvl1_cgm <- lvl1_means <- NULL
               for(i in 1:length(center_lvl1)){
                    if(center_lvl1[i] == "cluster"){
                         lvl1_cwc <- cbind(lvl1_cwc, center_within(x = data[,x_lvl1[i]], cluster = data[,cluster]))
                         lvl1_means <- cbind(lvl1_means, cluster_means(x = data[,x_lvl1[i]], cluster = data[,cluster]))
                         lmer_eq_lvl1[[i]] <- c(orig = x_lvl1[i],
                                                lvl1 = paste0(x_lvl1[i], "_cwc"),
                                                means = paste0(x_lvl1[i], "_means"))
                    }
                    if(center_lvl1[i] == "grand"){
                         lvl1_cgm <- cbind(lvl1_cgm, data[,x_lvl1[i]] - mean(data[,x_lvl1[i]], na.rm = TRUE))
                         lmer_eq_lvl1[[i]] <- c(orig = x_lvl1[i],
                                                lvl1 = paste0(x_lvl1[i], "_cgm"),
                                                means = NA)

                    }
                    if(center_lvl1[i] == "none"){
                         lmer_eq_lvl1[[i]] <- c(orig = x_lvl1[i],
                                                lvl1 = x_lvl1[i],
                                                means = NA)
                    }
               }
               if(any(center_lvl1 == "cluster")){
                    colnames(lvl1_cwc) <- paste0(x_lvl1[center_lvl1 == "cluster"], "_cwc")
                    colnames(lvl1_means) <- paste0(x_lvl1[center_lvl1 == "cluster"], "_means")
                    data <- cbind(data, lvl1_cwc, lvl1_means)
               }
               if(any(center_lvl1 == "grand")){
                    colnames(lvl1_cgm) <- paste0(x_lvl1[center_lvl1 == "grand"], "_cgm")
                    data <- cbind(data, lvl1_cgm)
               }

               lmer_eq_lvl1 <- t(simplify2array(lmer_eq_lvl1))

               if(any(center_lvl1means)){
                    lvl1_means_cgm <- NULL
                    if(length(center_lvl1means) == 1) center_lvl1means <- rep(center_lvl1means, sum(center_lvl1 == "cluster"))
                    lvl1_means_orig <- lmer_eq_lvl1[center_lvl1 == "cluster","means"]
                    for(i in 1:length(center_lvl1means)){
                         if(center_lvl1means[i]){
                              lvl1_means_cgm <- cbind(lvl1_means_cgm, data[,lmer_eq_lvl1[center_lvl1 == "cluster","means"][i]] -
                                                           mean(data[,lmer_eq_lvl1[center_lvl1 == "cluster","means"][i]], na.rm = TRUE))
                              lmer_eq_lvl1[center_lvl1 == "cluster","means"][i] <- paste0(lmer_eq_lvl1[center_lvl1 == "cluster","means"][i], "_cgm")
                         }
                    }
                    colnames(lvl1_means_cgm) <- lmer_eq_lvl1[center_lvl1 == "cluster","means"]
                    data <- cbind(data, lvl1_means_cgm)
               }else{
                    lvl1_means_orig <- NULL
               }

          }else{
               lvl1_means_orig <- NULL
               lmer_eq_lvl1 <- cbind(orig = x_lvl1,
                                     lvl1 = x_lvl1,
                                     means = NA)
          }

          fixed_effects <- fe_lvl1 <- paste(lmer_eq_lvl1[,"lvl1"], collapse = " + ")
     }else{
          lvl1_means_orig <- lmer_eq_lvl1 <- fixed_effects <- fe_lvl1 <- NULL
     }

     if(!is.null(x_lvl2)){
          if(any(center_lvl2)){
               if(length(x_lvl2) > 1 & length(center_lvl2) == 1) center_lvl2 <- rep(center_lvl2, length(x_lvl2))
               lvl2_cgm <- NULL
               for(i in 1:length(center_lvl2)){
                    if(center_lvl2[i]){
                         lvl2_cgm <- cbind(lvl2_cgm, data[,x_lvl2[i]] - mean(data[,x_lvl2[i]], na.rm = TRUE))
                         lmer_eq_lvl2[[i]] <- c(orig = x_lvl2[i],
                                                lvl2 = paste0(x_lvl2[i], "_cgm"))
                    }else{
                         lmer_eq_lvl2[[i]] <- c(orig = x_lvl2[i],
                                                lvl2 = x_lvl2[i])
                    }
               }
               colnames(lvl2_cgm) <- paste0(x_lvl2[center_lvl2], "_cgm")
               data <- cbind(data, lvl2_cgm)
               lmer_eq_lvl2 <- t(simplify2array(lmer_eq_lvl2))
          }else{
               lmer_eq_lvl2 <- cbind(orig = x_lvl2,
                                     lvl2 = x_lvl2)
          }
     }else{
          lmer_eq_lvl2 <- NULL
     }

     if(!is.null(lmer_eq_lvl1) & !is.null(center_lvl1)){
          if(any(center_lvl1 == "cluster")){
               lmer_eq_meanslvl1 <- cbind(orig = lmer_eq_lvl1[center_lvl1 == "cluster", "means"],
                                          lvl2 = lmer_eq_lvl1[center_lvl1 == "cluster", "means"])
               if(!is.null(lvl1_means_orig)) lmer_eq_meanslvl1[,"orig"][center_lvl1 == "cluster"][center_lvl1means] <- lvl1_means_orig

               lmer_eq_lvl2 <- rbind(lmer_eq_lvl2, lmer_eq_meanslvl1)
               if(is.list(y_lvl1means)) names(y_lvl1means) <- paste0(names(y_lvl1means), "_means")
          }
          if(!is.null(y_lvl1means)){
               if(!is.list(y_lvl1means) & length(y_lvl1means) == 1){
                    if(y_lvl1means[1] == "all"){
                         y_lvl1means <- list()
                         for(i in x_lvl1[center_lvl1 == "cluster"]) y_lvl1means[[paste0(i, "_means")]] <- c("Intercept", x_lvl1[!fixed_lvl1])
                    }else if(y_lvl1means[1] == "intercepts"){
                         y_lvl1means <- list()
                         for(i in x_lvl1[center_lvl1 == "cluster"]) y_lvl1means[[paste0(i, "_means")]] <- c("Intercept")
                    }else if(y_lvl1means[1] == "slopes"){
                         y_lvl1means <- list()
                         for(i in x_lvl1[center_lvl1 == "cluster"]) y_lvl1means[[paste0(i, "_means")]] <- c(x_lvl1[!fixed_lvl1])
                    }
               }else{
                    y_lvl1means <- lapply(y_lvl1means, function(x){
                         if(length(x) == 1  & x[1] == "all"){
                              c("Intercept", x_lvl1[!fixed_lvl1])
                         }else{
                              x
                         }
                    })
                    for(i in names(y_lvl1means)){
                         if(length(y_lvl1means[[i]]) == 1)
                              if(y_lvl1means[[i]] == "none" | is.na(y_lvl1means[[i]]))
                                   y_lvl1means[[i]] <- NULL
                    }
                    if(length(y_lvl1means) == 0) y_lvl1means <- NULL
               }

               if(!is.null(y_lvl1means)) y_lvl2 <- append(y_lvl2, y_lvl1means)
          }
     }

     if(any(!fixed_lvl1) & !is.null(x_lvl1)){
          random_effects <- paste(lmer_eq_lvl1[!fixed_lvl1,"lvl1"], collapse = " + ")
          random_pred_lvl1 <- lmer_eq_lvl1[!fixed_lvl1,"lvl1"]
     }else{
          random_effects <- "1"
          random_pred_lvl1 <- NULL
     }

     if(!is.null(y_lvl2)){
          .lvl2 <- names(y_lvl2)
          rownames(lmer_eq_lvl2) <- lmer_eq_lvl2[,1]
          names(y_lvl2) <- lmer_eq_lvl2[.lvl2,2]

          use_as_u0pred <- lapply(y_lvl2, function(x) any(x == "Intercept"))
          y_lvl2 <- lapply(y_lvl2, function(x) x[x != "Intercept"])
          for(i in names(y_lvl2)) if(length(y_lvl2[[i]]) == 0) y_lvl2[[i]] <- NULL

          if(!is.null(lmer_eq_lvl1) & !is.null(lmer_eq_lvl2)){
               rownames(lmer_eq_lvl1) <- lmer_eq_lvl1[,1]
               y_lvl2 <- lapply(y_lvl2, function(x){lmer_eq_lvl1[x,2]})

               xlvl_ints <- NULL
               for(i in names(y_lvl2)) xlvl_ints <- c(xlvl_ints, paste(i, y_lvl2[[i]], sep = " * ", collapse = " + "))
               xlvl_ints <- paste(xlvl_ints, collapse = " + ")
          }else{
               xlvl_ints <- NULL
          }

          if(any(unlist(use_as_u0pred))){
               fe_lvl2 <- means_as_crit <- NULL
               for(i in names(use_as_u0pred)){
                    if(any(names(use_as_u0pred) == i)){
                         means_as_crit[i] <- use_as_u0pred[[i]]
                    }else{
                         means_as_crit[i] <- FALSE
                    }
               }
               fe_lvl2 <- lmer_eq_lvl2[,2][lmer_eq_lvl2[,"lvl2"] %in% names(means_as_crit)][means_as_crit]
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

     if(!is.null(lmer_eq_lvl1)){
          use_preds_lvl1 <- lmer_eq_lvl1[,"lvl1"]
     }else{
          use_preds_lvl1 <- NULL
     }
     if(!is.null(lmer_eq_lvl2)){
          use_preds_lvl2 <- lmer_eq_lvl2[,"lvl2"]
     }else{
          use_preds_lvl2 <- NULL
     }

     if(model_type != 0){
          usable_mods <- c(null = TRUE,
                           meansout = !is.null(fe_lvl2),
                           randcoeff = !is.null(x_lvl1) & any(!fixed_lvl1),
                           intslopesout = !is.null(xlvl_ints))
          if(length(model_type) > 1) stop("Only one 'model_type' can be supplied", call. = FALSE)
          if(model_type > 4 | round(model_type) != model_type) stop("'model_type' must be an integer between 0 and 4", call. = FALSE)
          if(all(which(usable_mods) != model_type)) stop("The selected 'model_type' is not possible using the supplied data", call. = FALSE)
          usable_mods[-model_type] <- FALSE
     }

     if(model_type == 0){
          eq <- paste0(y_lvl1, " ~ ")
          if(!is.null(fixed_effects)) eq <- paste0(eq, fixed_effects, " + ")
          if(!is.null(xlvl_ints)) eq <- paste0(eq, xlvl_ints, " + ")
          eq <- paste0(eq, "(")
          if(!is.null(random_effects)){
               eq <- paste0(eq, random_effects, " | ")
          }else{
               eq <- paste0(eq, "1 | ")
          }
          eq <- as.formula(paste0(eq, cluster, ")"))
     }else{
          if(usable_mods[1]){
               eq <- as.formula(paste0(y_lvl1, " ~ (1 | ", cluster, ")"))
          }else if(usable_mods[2]){
               eq <- as.formula(paste0(y_lvl1, " ~ ", fe_lvl2, " + (1 | ", cluster, ")"))
          }else if(usable_mods[3]){
               eq <- as.formula(paste0(y_lvl1, " ~ ", random_effects, " + (", random_effects," | ", cluster, ")"))
          }else if(usable_mods[4]){
               eq <- as.formula(paste0(y_lvl1, " ~ ", fixed_effects, " + ", xlvl_ints, " + (", random_effects, " | ", cluster, ")"))
          }
     }

     mod <- lmerTest::lmer(eq, data = data, ...)
     sum <- summary(mod)
     tau <- attributes(sum$varcor[[1]])$stddev^2
     icc <- as.numeric(tau[1] / (sum$sigma^2 + tau[1]))
     ci <- cbind(Estimate = sum$coefficients[,1],
                 `Std. Error` = sum$coefficients[,2],
                 `CI (Lower)` = sum$coefficients[,1] - qnorm((1 - conf_level) / 2, lower.tail = FALSE) * sum$coefficients[,2],
                 `CI (Upper)` = sum$coefficients[,1] + qnorm((1 - conf_level) / 2, lower.tail = FALSE) * sum$coefficients[,2])
     if(nrow(ci) == 1) rownames(ci) <- "(Intercept)"

     cv <- cbind(Estimate = sum$coefficients[rownames(sum$varcor[[1]]),1],
                 Variance = diag(sum$varcor[[1]]),
                 `CV (Lower)` = sum$coefficients[rownames(sum$varcor[[1]]),1] - qnorm((1 - cred_level) / 2, lower.tail = FALSE) * diag(sum$varcor[[1]])^.5,
                 `CV (Upper)` = sum$coefficients[rownames(sum$varcor[[1]]),1] + qnorm((1 - cred_level) / 2, lower.tail = FALSE) * diag(sum$varcor[[1]])^.5)
     if(nrow(cv) == 1) rownames(cv) <- "(Intercept)"

     if(model_type == 0){
          if(any(!fixed_lvl1)){
               rel <- rel_lvl1(summary = sum, cluster = cluster, x_lvl1 = random_pred_lvl1, data = data)
          }else{
               rel <- rel_lvl1(summary = sum, cluster = cluster, data = data)
          }
          chisq <- chisq_hlmer(model = mod, summary = sum, y_lvl1 = y_lvl1,
                               x_lvl1 = random_pred_lvl1, x_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
     }else{
          if(usable_mods[1]){
               rel <- rel_lvl1(summary = sum, cluster = cluster, data = data)
               chisq <- chisq_hlmer(model = mod, summary = sum, y_lvl1 = y_lvl1, cluster = cluster, data = data)
          }else if(usable_mods[2]){
               rel <- rel_lvl1(summary = sum, cluster = cluster, data = data)
               chisq <- chisq_hlmer(model = mod, summary = sum, y_lvl1 = y_lvl1,
                                    x_lvl1 = NULL, x_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
          }else  if(usable_mods[3]){
               rel <- rel_lvl1(summary = sum, cluster = cluster, x_lvl1 = random_pred_lvl1, data = data)
               chisq <- chisq_hlmer(model = mod, summary = sum, y_lvl1 = y_lvl1,
                                    x_lvl1 = random_pred_lvl1, x_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
          }else if(usable_mods[4]){
               if(any(!fixed_lvl1)){
                    rel <- rel_lvl1(summary = sum, cluster = cluster, x_lvl1 = random_pred_lvl1, data = data)
               }else{
                    rel <- rel_lvl1(summary = sum, cluster = cluster, data = data)
               }
               chisq <- chisq_hlmer(model = mod, summary = sum, y_lvl1 = y_lvl1,
                                    x_lvl1 = random_pred_lvl1, x_lvl2 = use_preds_lvl2, cluster = cluster, data = data)
          }
     }

     out <- list(call = call,
                 model = mod,
                 summary = sum,
                 conf = ci,
                 cred = cv,
                 reliability = rel,
                 icc = icc,
                 chisq_tau = chisq,
                 conf_level = conf_level, cred_level = cred_level,
                 data = data)

     class(out) <- c("hlmer", "hlmerMod")
     out
}

#' Print method for hlmer-class objects
#'
#' @param x hlmer-class object.
#' @param ... Additional arguments for \code{print()}.
#' @param digits Number of digits to which results should be printed.
#'
#' @return A formatted hlmer-class object.
#' @export
print.hlmer <- function(x, ..., digits = 5){
     print.hlmer.hlmerMod(x = x, ..., digits = digits)
}


#' Print method for hlmerMod-class objects
#'
#' @param x hlmerMod-class object.
#' @param ... Additional arguments for \code{print()}.
#' @param digits Number of digits to which results should be printed.
#'
#' @return A formatted hlmerMod-class object.
print.hlmer.hlmerMod <- function(x, ..., digits = 5){
     cat("Call: \n")
     print(x$call)

     cat("\n")
     print(x$summary)

     cat("\n")
     conf_level <- x$conf_level
     if(is.null(conf_level)) conf_level <- .95
     cat(paste0(round(conf_level * 100), "%"), "confidence intervals for fixed effects: \n")
     print(x$conf, digits = digits)

     cat("\n")
     cred_level <- x$cred_level
     if(is.null(cred_level)) cred_level <- .95
     cat(paste0(round(cred_level * 100), "%"), "credibility intervals for random effects: \n")
     print(x$cred, digits = digits)

     cat("\n")
     cat("Intraclass correlation coefficient (ICC): \n")
     print(x$icc, digits = digits)

     cat("\n")
     cat("Reliability estimates for random effects: \n")
     print(x$reliability, digits = digits)

     cat("\n")
     cat("Approximate chi-square tests for the variances of random effects: \n")
     print(x$chisq_tau, digits = digits)

     tau_dims <- dim(x$summary$varcor[[1]])
     n_random <- prod(tau_dims) - tau_dims[1] * (tau_dims[2] - 1) / 2 + 1
     n_fixed <- nrow(x$summary$coefficients)
     n_params <- n_random + n_fixed
     cat("\n")
     cat(paste0("Total number of parameters estimated (", n_random, " random + ", n_fixed, " fixed): ", n_params), "\n")
}


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
#' @export
cluster_means <- function(x, cluster){
     x_center <- tapply(x, cluster, function(x) mean(x))
     cluster_lvls <- levels(factor(cluster))
     for(i in cluster_lvls) x[cluster == i] <- x_center[[i]]
     x
}


#' Generate centered versions of all variables in a data set.
#'
#' Variables that are cluster-mean centered will be labeled with a "cwc" suffix (abbreviation for centered within cluster),
#' the cluster means of level-1 variables will be labeled with a "means" suffix, and
#' variables that are grand-mean centered will be labeled with a "cgm" suffix (abbreviation for centered grand mean).
#'
#' @param data Data frame, matrix, or tibble containing the data to center.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#'
#' @return An augmented database with variables indicating level-1 cluster means, cluster-mean centered level-1 variables, grand-mean centered level-1 variables, and grand-mean centered level-2 variables.
#' @export
#'
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' center_data(data = timss, cluster = "idteach")
#' }
center_data <- function(data, cluster){
     data <- data.frame(data)
     .data <- data
     numvars <- suppressWarnings(apply(data, 2, function(x) is.numeric(as.numeric(unlist(x)))))
     data <- data.frame(apply(data[,numvars], 2, function(x) as.numeric(unlist(x))))
     data[,cluster] <- NULL
     data <- cbind(.data[,cluster], data)
     colnames(data)[1] <- cluster

     cluster_vars <- t(simplify2array(by(data, data[,cluster], apply, 2, var, na.rm = TRUE)))
     cluster_vars[is.na(cluster_vars)] <- 0
     lvl1 <- apply(cluster_vars, 2, function(x) any(x > 0))

     dat_lvl1 <- data[,lvl1]
     dat_lvl2 <- data[,!lvl1]
     dat_lvl2[,cluster] <- NULL

     dat_lvl1_means <- apply(dat_lvl1, 2, function(x) cluster_means(unlist(x), cluster = unlist(data[,cluster])))
     dat_lvl1_cwc <- dat_lvl1 - dat_lvl1_means
     dat_lvl1_cgm <- apply(dat_lvl1, 2, function(x) x - mean(x, na.rm = TRUE))
     dat_lvl2_cgm <- apply(dat_lvl2, 2, function(x) x - mean(x, na.rm = TRUE))

     colnames(dat_lvl1_means) <- paste0(colnames(dat_lvl1_means), "_means")
     colnames(dat_lvl1_cwc) <- paste0(colnames(dat_lvl1_cwc), "_cwc")
     colnames(dat_lvl1_cgm) <- paste0(colnames(dat_lvl1_cgm), "_cgm")
     colnames(dat_lvl2_cgm) <- paste0(colnames(dat_lvl2_cgm), "_cgm")

     as_tibble(cbind(.data, dat_lvl1_means, dat_lvl1_cwc, dat_lvl1_cgm, dat_lvl2_cgm))
}



#' Estimate the reliability of random effects
#'
#' @param summary lmer summary model.
#' @param x_lvl1 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-1 observations.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#'
#' @return Vector of reliability estimates for random effects.
#' @export
rel_lvl1 <- function(summary, cluster, x_lvl1 = NULL, data){
     se_mat <- se_lvl1(sigma = summary$sigma, data = data, x_lvl1 = x_lvl1, cluster = cluster)
     tao_vec <- attributes(summary$varcor[[1]])$stddev^2
     tao_mat <- matrix(tao_vec, nrow(se_mat), ncol(se_mat), TRUE)
     out <- apply(tao_mat / (tao_mat + se_mat^2), 2, mean)
     names(out) <- rownames(summary$varcor[[1]])
     out
}

#' Estimate level-1 standard errors for random effects
#'
#' @param sigma Residual standard deviation.
#' @param x_lvl1 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-1 observations.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#'
#' @return Matrix of estimated level-1 standard errors for random effects.
#' @export
se_lvl1 <- function(sigma, cluster, x_lvl1 = NULL, data){
     if(!is.null(x_lvl1)){
          out <- t(simplify2array(by(data, data[,cluster], function(x){
               x1 <- cbind(`(Intercept)` = 1, as.matrix(x[,x_lvl1]))
               invmat <- solve(t(x1) %*% x1)
               (diag(sigma^2 * invmat))^.5
          })))
     }else{
          out <- cbind(`(Intercept)` = (sigma^2 / by(data, data[,cluster], nrow))^.5)
     }
     out
}


#' Chi-square estimates for random-effects variances
#'
#' @param model Raw lmer output model.
#' @param summary Optional lmer summary model.
#' @param y_lvl1 Column label of \code{data} corresponding to the criterion variable for level-1 observations.
#' @param cluster Column label of \code{data} corresponding to the cluster/group identification variable.
#' @param x_lvl1 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-1 observations.
#' @param x_lvl2 Column label(s) of \code{data} corresponding to the predictor variable(s) for level-2 observations.
#' @param data Data frame, matrix, or tibble containing the data to use in the linear model.
#'
#' @author Jeffrey Dahlke and Jonathan Brown
#'
#' @return Table of chi-square estimates for random-effects variances.
#' @export
chisq_hlmer <- function(model, summary = NULL, y_lvl1, cluster, x_lvl1 = NULL, x_lvl2 = NULL, data){
     mod <- model
     rm(model)
     if(is.null(summary)){
          smod <- summary(mod)
     }else{
          smod <- summary
     }
     rm(summary)

     ## Random effects residuals and coefficients
     re_u <- ranef(mod)[[cluster]]
     re_b <- coef(mod)[[1]]

     ## Data frame of sample sizes
     n <- data.frame(table(data[,cluster]))

     ## Mean criterion score by cluster
     ave <- aggregate(unlist(data[,y_lvl1]), by = list(Var1 = unlist(data[,cluster])), FUN = "mean")

     ## Matrix of level-1 standard errors
     se_mat <- se_lvl1(sigma = smod$sigma, data = data, x_lvl1 = x_lvl1, cluster = cluster)

     # Combine data into a usable data frame
     dat <- merge(n, ave, by = "Var1")
     dat <- cbind(dat, re_b, re_u, se_mat)
     fe_names <- colnames(re_b)
     re_names <- colnames(re_u)
     colnames(dat) <- c(cluster, "n", "mean", paste0("b_", fe_names), paste0("u_", re_names), paste0("se_", re_names))
     rownames(dat) <- NULL

     ## Add in columns for any necessary level-2 variables
     if(!is.null(x_lvl2)){
          lvl2 <- data[!duplicated(unlist(data[,cluster])),]
          dat <- suppressWarnings(full_join(dat, as_tibble(lvl2)[,c(cluster, x_lvl2)], by = cluster))
     }

     gammas <- smod$coefficients[,1]
     names(gammas) <- rownames(smod$coefficients)

     splitnames <- str_split(fe_names, pattern = ":")
     splitnames <- lapply(splitnames, function(x){
          if(length(x) == 1){
               c(x, NA)
          }else{
               x
          }
     })
     splitnames <- data.frame(t(simplify2array(splitnames)))
     splitnames <- splitnames[!is.na(splitnames[,2]),]
     if(nrow(splitnames) == 0) splitnames <- NULL

     if(!is.null(x_lvl1)){
          formula_lvl1 <- as.formula(paste(y_lvl1, "~", paste(x_lvl1, collapse = " + ")))
     }else{
          formula_lvl1 <- as.formula(paste(y_lvl1, "~ 1"))
     }
     lm_lvl1 <- by(data, data[,cluster], function(x){suppressWarnings(summary(lm(formula_lvl1, data = x))$coeff)})
     if(nrow(lm_lvl1[[1]]) == 1){
          lm_b_lvl1 <- unlist(lapply(lm_lvl1, function(x) x[,1]))
          lm_se_lvl1 <- unlist(lapply(lm_lvl1, function(x) x[,2]))
          lm_b_lvl1 <- data.frame(id = names(lm_b_lvl1), `(Intercept)` = lm_b_lvl1)
          lm_se_lvl1 <- data.frame(id = names(lm_se_lvl1), `(Intercept)` = lm_se_lvl1)
          colnames(lm_b_lvl1) <- colnames(lm_se_lvl1) <- c(cluster, "(Intercept)")
     }else{
          lm_b_lvl1 <- data.frame(t(simplify2array(lapply(lm_lvl1, function(x) x[,1]))))
          lm_se_lvl1 <- data.frame(t(simplify2array(lapply(lm_lvl1, function(x) x[,2]))))
          lm_b_lvl1 <- cbind(id = rownames(lm_b_lvl1), lm_b_lvl1)
          lm_se_lvl1 <- data.frame(id = rownames(lm_se_lvl1), lm_se_lvl1)
          colnames(lm_b_lvl1) <- colnames(lm_se_lvl1) <- c(cluster, rownames(lm_lvl1[[1]]))
     }

     chisq_out <- df <- NULL
     for(i in re_names){
          if(i == "(Intercept)"){
               xi <- fe_names[fe_names %in% x_lvl2]
               if(!is.null(x_lvl1)){
                    chisq_out[i] <- sum((lm_b_lvl1$`(Intercept)` - gammas[i] - gammas[xi] %*% t(dat[,xi]))^2 / dat[,paste0("se_", i)]^2)
               }else{
                    if(length(xi) > 0){
                         chisq_out[i] <- sum((dat$mean - gammas[i] - gammas[xi] %*% t(dat[,xi]))^2 / dat[,paste0("se_", i)]^2)
                    }else{
                         chisq_out[i] <- sum((dat$mean - gammas[i])^2 / dat[,paste0("se_", i)]^2)
                    }
               }
          }else{
               smod
               xi <- fe_names[fe_names %in% paste0(i, ":", x_lvl2)]
               xlvli <- gsub(x = xi, pattern = paste0(i, ":"), replacement = "")
               if(length(xi) > 0){
                    chisq_out[i] <- sum((lm_b_lvl1[,i] - gammas[i] - gammas[xi] %*% t(dat[,xlvli]))^2 / dat[,paste0("se_", i)]^2)
               }else{
                    chisq_out[i] <- sum((lm_b_lvl1[,i] - gammas[i])^2 / dat[,paste0("se_", i)]^2)
               }
          }
          df[i] <- smod$ngrps - length(xi) - 1
     }

     ## Combine results into a table of output
     cbind(tau = attributes(smod$varcor[[1]])$stddev^2, chisq = chisq_out, df = df, `Pr(>chisq)` = pchisq(q = chisq_out, df = df, lower.tail = FALSE))
}
