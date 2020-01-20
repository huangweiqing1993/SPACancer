#' Identify driver genes
#'
#' \code{convolution_test} assigns a p value for each gene.
#'
#' @param bmr A Backround Mutation Rate (BMR) matrix. Rows and columns correspond to genes and patients.
#' @param n A mutation count matrix. Rows and columns correspond to genes and patients.
#'
#' @return A vector of p values
#'
#' @examples
#' # Load the example dataset:
#' data(mutation)
#'
#' convolution_test(bmr, n)
#'
#' @author Weiqing Huang
#'
#' @import data.table
#' @importFrom Rcpp evalCpp
#' @useDynLib SPACancer
#' @export convolution_test

convolution_test <- function(bmr, n) {
  ng <- nrow(bmr)
  gp <- dimnames(bmr)
  p0 <- matrix(dpois(0, bmr), ng, dimnames = gp)
  p1 <- 1 - p0
  s0 <- 0 * p0
  s1 <- round(-log(p1), 1) * 10
  n <- n != 0
  sob <- rowSums(s1 * n)
  p0 <- as.data.frame(t(p0))
  p1 <- as.data.frame(t(p1))
  s0 <- as.data.frame(t(s0))
  s1 <- as.data.frame(t(s1))
  mapply(function(p0, p1, s0, s1, sob) {
    l <- length(p0)
    data <- data.table(x = rep(1:l, 2), y = c(s0, s1), probability = c(p0, p1))
    res <- convolution(data, sob)$probability
    tail(res, 1)
  }, p0, p1, s0, s1, sob)
}

convolution <- function(dat, threshold=NA, verbose=F) {


  # make sure that the x’s have no ‘jumps’ (1:N) ----------------------------
  tmp_x <- data.table(x = dat[, unique(x)])
  setkey(tmp_x,x)
  tmp_x[, x_new := 1:.N]
  setkey(dat, x)
  dat <- tmp_x[dat][, .(x_new, probability, y)]
  setnames(dat, c("x", "probability", "y"))
  setkey(dat, x)


  # set threshold -----------------------------------------------------------
  max_possible_score <- dat[, max(y), by = x][, sum(V1)]
  if (is.na(threshold)) {
    thres <- max_possible_score
  } else {
    thres <- threshold
  }


  # set dummy vector of integers --------------------------------------------
  ints <- 0:(thres+1)


  # perform convolution -----------------------------------------------------
  out <- as.data.table(convolution_body(as.matrix(dat[, list(y, probability, x)]), threshold=(thres+1), integers=ints))
  setnames(out, c("y", "probability"))



  # return results ----------------------------------------------------------
  return(out)

}
