chisq_hlmer <- function(model, summary = NULL, criterion_lvl1, predictors_lvl1 = NULL, predictors_lvl2 = NULL, cluster, data){
     mod <- model
     rm(model)
     if(is.null(summary)){
          smod <- summary(mod)
     }else{
          smod <- summary
     }
     rm(summary)

     #Store EMPIRICAL BAYES (EB) adjusted random effects
     re_u <- ranef(mod)[[cluster]]
     re_b <- coef(mod)[[1]]

     ## Now create a data frame that combines all the information we want for the chi-square formula
     #Make a temp data frame of group IDs and sample sizes
     n = data.frame(table(data[,cluster]))

     #Obtain group means
     ave = aggregate(unlist(data[,criterion_lvl1]), by = list(Var1 = unlist(data[,cluster])), FUN = "mean")

     se_mat <- se_lvl1(sigma = smod$sigma, data = data, predictors_lvl1 = predictors_lvl1, cluster = cluster)


     #Formally combine group ID, group sample size, group mean, and group EB Random Effects, & change names for usage ease
     final <- merge(n, ave, by = "Var1")
     final <- cbind(final, re_b, re_u, se_mat)
     fe_names <- colnames(re_b)
     re_names <- colnames(re_u)
     colnames(final) <- c(cluster, "n", "mean", paste0("b_", fe_names), paste0("u_", re_names), paste0("se_", re_names))
     rownames(final) <- NULL

     #Also, add in a column for the level 2 predictor so that it can be easily called in the formula
     # Pull out level 2 data from the data frame and then match the cluster IDs to the "final2" data
     if(!is.null(predictors_lvl2)){
          lvl2 <- data[!duplicated(unlist(data[,cluster])),]
          final <- suppressWarnings(full_join(final, as_tibble(lvl2)[,c(cluster, predictors_lvl2)], by = cluster))
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

     if(!is.null(predictors_lvl1)){
          formula_lvl1 <- as.formula(paste(criterion_lvl1, "~", paste(predictors_lvl1, collapse = " + ")))
     }else{
          formula_lvl1 <- as.formula(paste(criterion_lvl1, "~ 1"))
     }
     lm_lvl1 <- by(data, data[,cluster], function(x){summary(lm(formula_lvl1, data = x))$coeff})
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
               xi <- fe_names[fe_names %in% predictors_lvl2]
               if(!is.null(predictors_lvl1)){
                    chisq_out[i] <- sum((lm_b_lvl1$`(Intercept)` - gammas[i] - gammas[xi] %*% t(final[,xi]))^2 / final[,paste0("se_", i)]^2)
               }else{
                    if(length(xi) > 0){
                         chisq_out[i] <- sum((final$mean - gammas[i] - gammas[xi] %*% t(final[,xi]))^2 / final[,paste0("se_", i)]^2)
                    }else{
                         chisq_out[i] <- sum((final$mean - gammas[i])^2 / final[,paste0("se_", i)]^2)
                    }
               }
          }else{
               smod
               xi <- fe_names[fe_names %in% paste0(i, ":", predictors_lvl2)]
               xlvli <- gsub(x = xi, pattern = paste0(i, ":"), replacement = "")
               if(length(xi) > 0){
                    chisq_out[i] <- sum((lm_b_lvl1[,i] - gammas[i] - gammas[xi] %*% t(final[,xlvli]))^2 / final[,paste0("se_", i)]^2)
               }else{
                    chisq_out[i] <- sum((lm_b_lvl1[,i] - gammas[i])^2 / final[,paste0("se_", i)]^2)
               }
          }
          df[i] <- smod$ngrps - length(xi) - 1
     }

     ## Combine results into a table of output
     cbind(tau = attributes(smod$varcor[[1]])$stddev^2, chisq = chisq_out, df = df, `Pr(>chisq)` = pchisq(q = chisq_out, df = df, lower.tail = FALSE))
}
