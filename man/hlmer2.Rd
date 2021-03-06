% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hlmer.R
\name{hlmer2}
\alias{hlmer2}
\title{lmer wrapper to produce HLM7-style output from multi-level equations}
\usage{
hlmer2(eq_lvl1, eq_lvl2 = NULL, cluster, conf_level = 0.95,
  cred_level = 0.95, remove_missing = TRUE,
  check_lvl1_variance = TRUE, data, ...)
}
\arguments{
\item{eq_lvl1}{Regression equation (as a string) for level-1 models. Do not use parentheses in the model. When supplying an interaction between level-1 variables, use either "*" or ":" notation and manually supply the main effects. In the output, all within-level interactions will be indicated by a "_X_" symbol between the variables.}

\item{eq_lvl2}{List of regression equations (as strings) for level-2 models, with outcomes labeled as "intercept" or as level-1 predictor names. If left NULL, random intercepts will be estimated, but all slopes will be treated as fixed effects.
To specify a fixed effect, simply omit the variable from the \code{eq_lvl2} list or provide a regression of the following form: "predictor ~ 0", where "predictor" is the name of a level-1 explanatory variable.
To specify a random effect with no level-2 explanatory variable, provide a regression of the following form: "predictor ~ 1".
When supplying an interaction between level-2 variables, use either "*" or ":" notation and manually supply the main effects. In the output, all within-level interactions will be indicated by a "_X_" symbol between the variables.}

\item{cluster}{Column label of \code{data} corresponding to the cluster/group identification variable.}

\item{conf_level}{Confidence level to use in constructing confidence bounds around fixed effects.}

\item{cred_level}{Credibility level to use in constructing credibility bounds (ranges of plausible values for random coefficients) around random effects.}

\item{remove_missing}{Logical scalar that determines whether cases with missing data should be omitted.}

\item{check_lvl1_variance}{Logical scalar that determines whether clusters with no variance in level-1 predictors should be screened out.}

\item{data}{Data frame, matrix, or tibble containing the data to use in the linear model.}

\item{...}{Additional arugments to be passed to the \code{lme4::lmer()} function.}
}
\value{
Output from the \code{lmerTest::lmer()} function augmented with ICC statistics, random-effects reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
}
\description{
This function is a wrapper for \code{lmer()} to formulate a mixed-effects linear model based on user-supplied multi-level equations.
In addition to the output from \code{lmer()}, \code{hlmer2()} also provides estimates of ICC statistics, random-effect reliability estimates, confidence intervals for fixed effects, and chi-square tests for random-effects variance.
To ensure that a simplified model (e.g., an unconditional model) run with this function will be based on the same cases as a more complex model, supply a \code{model_type} argument using the values outlined in the documentation for the \code{hlmer()} function.
IMPORTANT: Unlike the \code{hlmer()} function, this function requires that the user manually center predictors. To automate the centering process, see the \code{center_data()} function.
}
\examples{
\dontrun{
## Center the HSB database:
dat <- center_data(cluster = "ID", data = hsb)

## Unconditional model (Raudenbush and Bryk Table 4.2):
hlmer2(eq_lvl1 = "MATHACH ~ 1", cluster = "ID", data = dat)

## Means-as-outcomes model (Raudenbush and Bryk Table 4.3):
hlmer2(eq_lvl1 = "MATHACH ~ 1", eq_lvl2 = "intercept ~ MEANSES",
cluster = "ID", data = dat)

## Random-coefficients model (Raudenbush and Bryk Table 4.4):
hlmer2(eq_lvl1 = "MATHACH ~ SES_cwc",
eq_lvl2 = list("intercept ~ MEANSES", "SES_cwc ~ 1"),
cluster = "ID", data = dat)

## Slopes and intercepts as outcomes model (Raudenbush and Bryk Table 4.5):
hlmer2(eq_lvl1 = "MATHACH ~ SES_cwc",
eq_lvl2 = list("intercept ~ SECTOR + MEANSES", "SES_cwc ~ SECTOR + MEANSES"),
cluster = "ID", data = dat)

## For a more complex model, we can specify which level-2 predictors
## should be used to predict which level-1 random effects.
## In this TIMSS example, the level-1 predictor is "self_efy" and
## the level-2 predictors are "students" and "alg." If we want to
## predict level-1 intercepts using both level-2 predictors, but we
## only want to predict level-1 "self_efy" slopes using "students,"
## we can specify that with the following model:
##
hlmer2(eq_lvl1 = "scores ~ self_efy",
eq_lvl2 = list("intercept ~ alg + students", "self_efy ~ students"),
cluster = "idteach", data = timss)
}
}
