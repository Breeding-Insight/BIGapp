# function name: CalcBIC()

# Usage: CalcBIC(y, PC, K)

# Description:
#  Calculate BIC for GWAS model comparison.

# Arguments:
#  y: length N vector
#  PC: matrix of principal components with N rows and P columns
#  K: kinship matrix with N rows and N columns
#  [NOTE] Each row must represent the same genotype!

# Value:
# n.input: number of genotypes in the input data
# n.sample: sample size. When there are k missing values, n.sample = n.input - k
# BIC$withK: data frame of the BIC calculation results with kinship
# BIC$withoutK: data frame of the BIC calculation results without kinship
#  $n.PC: number of PCs
#  $LL: Log-likelihood
#  $n.par: number of parameters
#  $BIC: BIC (BIC = n.par * log(n.sample) - 2 * LL)

# Details:
#  This function calculate BIC with 1, 2, 3, ..., P number of PCs, with and
#  without kinship matrix. Result is returned as a list object. The mixed model
#  for GWAS was fitted by using mixed.solve() function. We can get estimated
#  parameters from this model fitting. Then, using the estimated parameters,
#  log-likelihood is calculated by using the equation (2) in Kang et al., 2008.
#  Finally, BIC is calculated by using the standard formula (BIC = K * log(N) - 2 * LL)


# Note for the users who were using GAPIT:
# The BIC in GAPIT uses an old/orginal formula (Schwarz 1978). In this case,
# BIC is "the larger, the better"
# After ~50 years of the work of Schwarz, now, many people uses a BIC formula
# which is "the smaller, the better". I also used this new, common definition.
# Thus, be careful with this difference.
# You can still cite (Schwarz 1978) for the reference of BIC.
# The old BIC and the new BIC have one-to-one correspondence. Thus, your result
# never change whichever definition you use.


#' function for BIC calculation
#'
#' @param y describe documentation
#' @param PC describe documentation
#' @param K describe documentation
#'
#' @import rrBLUP
#' @import MASS
CalcBIC <- function(y, PC, K) {
   # number of max. PC
   n.max.pc <- ncol(PC)

   # remove NA
   tf.na <- !is.na(y)
   y.common <- y[tf.na]
   K.common <- K[tf.na, tf.na]
   PC.common <- PC[tf.na, ]

   # ----- BIC with K ----- #
   # fit model by using mixed.solve function
   LL.vec <- rep(NA, n.max.pc + 1)
  # names(LL.vec) <- paste0("PC", 0:n.max.pc)
   names(LL.vec) <- paste0(0:n.max.pc)
   for ( k in 0:n.max.pc ) {
      if ( k == 0 ) { X.matrix <- NULL }
      if ( k >= 1 ) { X.matrix <- cbind(1, PC.common[, 1:k]) }
      ms <- mixed.solve(y = y.common, Z = NULL, K = K.common, X = X.matrix, method = "REML")
      n <- length(y.common)
      q <- length(ms$beta)
      H <- K.common + diag(nrow(K.common)) * (ms$Ve / ms$Vu)
      H.inv <- ginv(H)
      a <- (-1) * n * log(2 * pi * ms$Vu)
      b <- (-1) * as.numeric(determinant(H, logarithm = TRUE)$modulus)
      if ( k == 0 ) {
         X.matrix <- matrix(1, nrow = n, ncol = 1)
         c <- as.numeric((-1) * (1 / ms$Vu) * t(y.common - X.matrix %*% ms$beta) %*% H.inv %*% (y.common - X.matrix %*% ms$beta))
      } else {
         c <- as.numeric((-1) * (1 / ms$Vu) * t(y.common - X.matrix %*% ms$beta) %*% H.inv %*% (y.common - X.matrix %*% ms$beta))
      }
      LogLik.ML <- 0.5 * (a + b + c) # full likelihood based on REML
      LL.vec[k+1] <- LogLik.ML
   }
   n <- length(y.common)
   n.par.vec <- seq(from = 3, by = 1, length.out = n.max.pc + 1) # n.par starts from 3 & increases by 1
   BIC.vec <- -2 * LL.vec + n.par.vec * log(n)

   # ----- BIC without K ----- #
   LL.vec.woK <- rep(NA, n.max.pc + 1)
   #names(LL.vec.woK) <- paste0("PC", 0:n.max.pc)
   names(LL.vec.woK) <- paste0(0:n.max.pc)
   for ( k in 0:n.max.pc ) {
      if ( k == 0 ) {
         LL.vec.woK[k+1] <- as.numeric(logLik(lm(y.common ~ 1)))
      } else {
         df.lm <- data.frame("y" = y.common, PC.common[, 1:k])
         colnames(df.lm)[2:ncol(df.lm)] <- paste0("PC", 1:k)
         LL.vec.woK[k+1] <- as.numeric(logLik(lm(y ~., df.lm)))
      }
   }
   n.par.vec.woK <- seq(from = 2, by = 1, length.out = n.max.pc + 1) # n.par starts from 2
   BIC.vec.woK <- n.par.vec.woK * log(n) - 2 * LL.vec.woK

   # return
   ResList <- list("n.input" = length(y),
                   "n.sample" = length(y.common),
                   "BIC" = list("withK" = data.frame("n.PC" = names(LL.vec),
                                                     "LL" = LL.vec,
                                                     "n.par" = n.par.vec,
                                                     "BIC" = BIC.vec),
                                "withoutK" = data.frame("n.PC" = names(LL.vec.woK),
                                                        "LL" = LL.vec.woK,
                                                        "n.par" = n.par.vec.woK,
                                                        "BIC" = BIC.vec.woK)))
   return(ResList)
}
