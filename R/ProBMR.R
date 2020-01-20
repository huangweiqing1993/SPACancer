#' Poisson Regression of Background Mutation Rate (BMR)
#'
#' \code{ProBMR} estimates BMRs using two mixed-effects Poisson regression models applied to link gene associated features and accommodate the overdispersion in mutation counts.
#'
#' @param m A mutation count matrix. Rows and columns correspond to genes and patients.
#' @param M A vector of gene coverage.
#' @param x A genomic covariates matrix. Rows and columns correspond to genes and features.
#'
#' @return A list with components:
#' \describe{
#' \item{bkgd1}{A list contains bmr (fitted values), fixef (fixed effects), ranef (random effects), VarCorr (Variance and Correlation Components) in the first regression.}
#' \item{bkgd2}{A list contains bmr (fitted values), ranef (random effects), VarCorr (Variance and Correlation Components) in the second regression.}
#' }
#'
#' @examples
#' # Load the example dataset:
#' data(mutation)
#'
#' fit <- ProBMR(m, M, x)
#'
#' @author Weiqing Huang
#'
#' @import lme4
#' @export ProBMR
ProBMR <- function(m, M, x) {
  ng <- nrow(m)
  np <- ncol(m)
  gp <- dimnames(m)
  data <- data.frame(m = c(m), M = rep(M, np), p = rep(1:np, each = ng))
  x <- kronecker(rep(1, np), x)
  fit <- glmer(m ~ x + (x | p), data, family = poisson, control = glmerControl(optimizer = "nloptwrap", calc.derivs = F), nAGQ = 0, offset = log(M))
  bmr <- matrix(fitted(fit), ng, dimnames = gp)
  bkgd1 <- list(bmr = bmr, fixef = fixef(fit), ranef = ranef(fit)$p, VarCorr = VarCorr(fit)$p)
  data <- data.frame(m = c(m), bmr = c(bmr), g = rep(1:ng, np))
  fit <- glmer(m ~ 0 + (1 | g), data, family = poisson, control = glmerControl(optimizer = "nloptwrap", calc.derivs = F), nAGQ = 0, offset = log(bmr))
  bmr <- matrix(fitted(fit), ng, dimnames = gp)
  bkgd2 <- list(bmr = bmr, ranef = ranef(fit)$g, VarCorr = VarCorr(fit)$g)
  list(bkgd1 = bkgd1, bkgd2 = bkgd2)
}
