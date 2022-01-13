RWSS_GCV <- function(x, lambda1, criterion1 = 0.01, b, lambda2 = lambda1, criterion2 = 0.05){
  m <- length(x)
  wi <- rep(1, m)
  ss <- smooth.spline(x, lambda = lambda1, all.knots = TRUE, w = wi)
  t.1 <- ss$y
  r_org <- x - t.1
  sigma_MAV <- median(abs(r_org)) / 0.6745
  r <- r_org / sigma_MAV
  flag <- sum(abs(r[r<0]))
  while ( 1 ){
    wi <- ifelse(r < 0, 1, ifelse(1 - (r / b) ^ 2 > 0, 1 - (r / b) ^ 2, 0))
    ss <- smooth.spline(t.1, all.knots = TRUE, w = wi, cv = FALSE)
    t.2 <- ss$y
    r_org <- t.1 - t.2
    sigma_MAV <- median(abs(r_org[r_org < b * sigma_MAV])) / 0.6745
    r <- r_org / sigma_MAV
    if (abs(flag/sum(abs(r[r<0])) - 1) < criterion1){
      break
    }
    flag <- sum(abs(r[r<0]))
    t.1 = t.2
  }
  
  # WSS parts
  # pre-baseline 
  b <- t.1
  # initialize the weights
  weights <- rep(1, m)
  # redo the smoothing.spline fitting
  ss_re <- smooth.spline(b, lambda = lambda2, all.knots = TRUE, w = weights)
  t_1 <- ss_re$y
  
  residuals <- log((b - t_1)^2)
  variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
  flag <- sum(abs(b - t_1))
  index = which((b-t_1) < 0)
  
  while(1){
    weights <- 1 / exp(variance$y)
    weights <- weights/sum(weights) # scale
    weights[index] = 1-weights[index]
    
    ss_re <- smooth.spline(t_1, lambda = lambda2, all.knots = TRUE, w = weights, cv=FALSE)
    res2 <- (t_1 - ss_re$y)^2
    residuals <- log(ifelse(res2 == 0, min(res2[res2 != 0]), res2))
    variance <- smooth.spline(residuals, cv = TRUE, all.knots = TRUE)
    
    if( abs(sum(abs(t_1 - ss_re$y)) / flag - 1) < criterion2 ){break}
    flag <- sum(abs(t_1 - ss_re$y))
    
    index = which((t_1-ss_re$y) < 0)
    t_1 <- ss_re$y
    
    
    if (length(index)<0.1*length(t_1)) {break}
  }
  return(t_1)
}
