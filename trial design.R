trial_design <- function(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  #phase 2 experimental treatment
  n_1 <- floor(n1*(1-a))
  #phase 2 control treatment
  n_00 <- floor((n1*2)*(1-a))
  #phase 3 experimental treatment 
  n_21 <- floor((n_1+n21)*(1-a))
  n_22 <- floor((n_1+n22)*(1-a))
  #phase 3 control treatment 
  n_01 <- floor((n_00+0.5*n21)*(1-a))
  n_02 <- floor((n_00+0.5*n22)*(1-a))
  # generate samples
  b <- matrix(0,nrow = 4, ncol = n_1)
  A <- sim_data(n, theta[1], mu[1], sigma, tau)
  b0 <- A[1:n_00,1]
  y01 <- A[1:n_01,2]
  y02 <- A[1:n_02,2]
  y1 <- matrix(0,nrow = 4, ncol = n_21)
  y2 <- matrix(0,nrow = 4, ncol = n_22)
  for (i in 1:4){
    d <- sim_data(n, theta[i+1], mu[i+1], sigma, tau)
    b[i,] <- d[1:n_1,1]
    y1[i,] <- d[1:n_21,2]
    y2[i,] <- d[1:n_22,2]
  }
  # interim analysis
  c <- rep(-1, 4)
  for (i in 1:4){
    if (chisq.test(rbind(c(sum(b[i,]), n_1-sum(b[i,])), c(sum(b0), n_00-sum(b0))))$p.value<p1){
      c[i] <- 1 
    }
  }
  e <- c
  #phase 3:
  d <- which(c==1, arr.ind=TRUE)
  #  1 success
  if (length(d)==1){
    e[d] <- t.test(y1[d,], y01, alternative = "two.sided")$p.value <  p2
  }
  # 2 and more sucess
  if (length(d)>1){
    for (j in d){
    e[j] <- t.test(y2[j,], y02, alternative = "two.sided")$p.value < p3
    }
  }
  return(e)
}



