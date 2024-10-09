
library(kpcalg)

kernelCItest <- function(x, y, S, suffStat, verbose = FALSE) {
  # Ensure x and y are numeric
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Extract necessary parameters from suffStat
  data <- suffStat$data
  method <- suffStat$ic.method
  sig <- suffStat$sig
  numPerms <- suffStat$p
  
  # Perform the CI test using Distance Covariance or HSIC
  if (length(S) == 0) { 
    # Marginal independence test
    if (method == 'hsic.gamma') {
      pval <- hsic.gamma(x = data[, x], y = data[, y], sig = sig)$p.value
    } else if (method == 'dcc.perm') {
      pval <- dcov.test(x = data[, x], y = data[, y], R = numPerms)$p.value
    }
  } else {
    # Conditional independence: Regress x and y on S and check residual independence
    residuals <- regrXonS(data[, c(x, y)], data[, S])
    resx <- residuals[, 1]
    resy <- residuals[, 2]
    
    if (method == 'hsic.gamma') {
      pval <- hsic.gamma(x = resx, y = resy, sig = sig)$p.value
    } else if (method == 'dcc.perm') {
      pval <- dcov.test(resx, resy, R = numPerms)$p.value
    }
  }
  return(pval)
}
