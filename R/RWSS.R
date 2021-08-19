RWSS <- function(x, lambda1, criterion1 = 0.01, b, lambda2 = lambda1, criterion2 = 0.05){
  m <- length(x)
  wi <- rep(1, m)
  ss <- smooth.spline(x, lambda = lambda1, all.knots = TRUE, w = wi)
  z <- ss$y
  r_org <- x - z
  sigma_MAV <- median(abs(r_org)) / 0.6745
  r <- r_org / sigma_MAV
  flag <- sum(abs(r[r<0]))
  niter = 0
  while ( 1 ){
    wi <- ifelse(r < 0, 1, ifelse(1 - (r / b) ^ 2 > 0, 1 - (r / b) ^ 2, 0))
    ss <- smooth.spline(x, lambda = lambda1, all.knots = TRUE, w = wi)
    z <- ss$y
    r_org <- x - z
    sigma_MAV <- median(abs(r_org[r_org < b * sigma_MAV])) / 0.6745
    r <- r_org / sigma_MAV
    if (abs(flag/sum(abs(r[r<0])) - 1) < criterion1){
      break
    }
    flag <- sum(abs(r[r<0]))
  }
  
  # WSS parts
  # pre-baseline 
  b <- z
  # initialize the weights
  weights <- rep(1, m)
  # redo the smoothing.spline fitting
  ss_re <- smooth.spline(b, lambda = lambda2, all.knots = TRUE, w = weights)
  z <- ss_re$y
  
  residuals <- log((b - z)^2)
  variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
  flag <- sum(abs(b - z))
  index = which((b-z) < 0)
  
  while(1){
    weights <- 1 / exp(variance$y)
    weights <- weights/sum(weights) # scale
    weights[index] = 1-weights[index]
    
    ss_re <- smooth.spline(b, lambda = lambda2, all.knots = TRUE, w = weights)
    res2 <- (b - ss_re$y)^2
    residuals <- log(ifelse(res2 == 0, min(res2[res2 != 0]), res2))
    variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
    if( abs(sum(abs(z - ss_re$y)) / flag - 1) < criterion2 ){break}
    flag <- sum(abs(z - ss_re$y))
    z <- ss_re$y
    
    index = which((b-z) < 0)
    
    if (length(index)<0.1*length(b)) {break}
  }
  return(z)
}
