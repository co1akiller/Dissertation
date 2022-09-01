
# At phase 2:
set.seed(NULL)
delta <- 0.17
delta2 <- 16
n <- 300
n1 <- 66
sigma <- 32
tau <- 10
#attrition
a <- 0.1
p1 <- 0.2
n21 <- 38
n22 <- 20
p2 <- 0.05
p3 <- 0.025
m <- 1000
M <- 100
theta_0 <- c(0.33,0.33,0.33,0.33,0.33)
mu_0 <- c(16,16,16,16,16)
theta_1 <- c(0.33,0.33+delta,0.33+delta,0.33+delta,0.33+delta)
mu_1 <- c(16,16+delta2,16+delta2,16+delta2,16+delta2)
theta <- c(0.33,0.33+delta,0.33,0.33,0.33)
mu <- c(16,16+delta2,16,16,16)
theta <- c(0.33,0.33+delta,0.33+delta,0.33+delta,0.33+delta)
mu <- c(16,16+delta2,16+delta2,16+delta2,16+delta2)

sim_data <- function(n, theta, mu, sigma, tau){
  b <- rbinom(n,1,theta)
  y <- rnorm(n,mean = mu + tau*b,sd=sigma)
  return(cbind(b,y))
}

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
  test_result <- rep(-1, 4)
  test_statistic <- rep(0,4)
  for (i in 1:4){
    test_statistic[i] <- chisq.test(rbind(c(sum(b[i,]), n_1-sum(b[i,])),
                                          c(sum(b0), n_00-sum(b0))))$p.value
    if (test_statistic[i]<p1){
      test_result[i] <- 1 
    }
  }
  if (length(which(test_result==1, arr.ind=TRUE))>2){
    test_result[which.min(test_statistic)] <- -1
    test_statistic[which.min(test_statistic)] <- 1
    if (length(which(test_result==1, arr.ind=TRUE))>2){
      test_result[which.min(test_statistic)]<- -1
    }
  }
  test_result_mu <- test_result
  #phase 3:
  test_result_reject <- which(test_result==1, arr.ind=TRUE)
  #  1 success
  if (length(test_result_reject)==1){
    test_result_mu[test_result_reject] <- t.test(y1[test_result_reject,], y01, alternative = "two.sided")$p.value <  p2
  }
  # 2 and more sucess
  if (length(test_result_reject)>1){
    for (j in test_result_reject){
      test_result_mu[j] <- t.test(y2[j,], y02, alternative = "two.sided")$p.value < p3
    }
  }
  return(test_result_mu)
}


trial_assess_mu <- function(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  f <- matrix(0,nrow=m,ncol=4)
  for (i in 1:m){
    f[i,] <- trial_design(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
  }
  k1 = k2 = 0
  True_null = False_null = numeric(0)
  for (j in 1:4){
    if (mu[j+1]==mu_0[j+1]){
      True_null <- append(True_null,j) 
    }
  }
  for (j in 1:4){
    if (mu[j+1]!=mu_0[j+1]){
      False_null <- append(False_null,j) 
    }
  }
  if (length(True_null)==0){
    fwer1 <- NA 
  }
  if (length(True_null)==1){
    f1 <- f[,True_null]
    for (i in 1:m){
      if (length(which(f1[i]==1, arr.ind=TRUE))>0){
        k1 <- k1+1
      }
    }  
  } 
  if (length(True_null)>1){
    f1 <- f[,True_null]
    for (i in 1:m){
      if (length(which(f1[i,]==1, arr.ind=TRUE))>0){
        k1 <- k1+1
      }
    }  
  } 
  if (length(False_null)==0){
    fwer2 <- NA 
  }
  if (length(False_null)==1){
    f1 <- f[,False_null]
    for (i in 1:m){
      if (length(which(f1[i]==1, arr.ind=TRUE))>0){
        k2 <- k2+1
      }
    }  
  } 
  if (length(False_null)>1){
    f1 <- f[,False_null]
    for (i in 1:m){
      if (length(which(f1[i,]==1, arr.ind=TRUE))>0){
        k2 <- k2+1
      }
    }  
  } 
  fwer1 <- k1/m
  fwer2 <- k2/m
  power <- 1-fwer2
  E1 <- numeric(0)
  P <- numeric(0)
  for (i in True_null){
    E1 <-append(E1,length(which(f[,i]==1, arr.ind=TRUE))/m)
  }
  for (i in False_null){
    P <-append(P,1- length(which(f[,i]!=1, arr.ind=TRUE))/m)
  }
  my_list <- list("fwer1" = fwer1,"power" = power, "error1"=E1, "power" = P)
  return(my_list)
}  
trial_assess_mu(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, tau, 
                n1, n21, n22, p1, p2, p3,a)

#sample size
sample_size <- function(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  f <- matrix(0,nrow=m,ncol=4)
  for (i in 1:m){
    f[i,] <- trial_design(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
  }
  N <- rep(0,m)
  for (i in 1:m){
    if (length(which(f[i,]==-1, arr.ind=TRUE))<3){
      N[i] <- 491
    }
    if (length(which(f[i,]==-1, arr.ind=TRUE))==3){
      N[i] <- 426
    }
    if (length(which(f[i,]==-1, arr.ind=TRUE))==4){
      N[i] <- 359
    }
  } 
  return(mean(N))
}
sample_size(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)



###################### Scenario 1
delta1_group <- seq(0.01,0.66,0.01) 
delta2_group <- 50*delta1_group
fwer1 <- rep(0,length(delta2_group))
power <- rep(0,length(delta2_group))
error1 <- matrix(0,nrow =3,ncol=length(delta2_group))
powerk <- rep(0,length(delta2_group))
N <- rep(0,length(delta2_group))
for (i in 1:length(delta2_group)){
  theta <- c(0.33,0.33+delta1_group[i],0.33,0.33,0.33)
  theta_0 <- c(0.33,0.33,0.33,0.33,0.33)
  theta_1 <- c(0.33,0.33+delta1_group[i],0.33+delta1_group[i]
               ,0.33+delta1_group[i],0.33+delta1_group[i])
  mu <- c(16,16+delta2_group[i],16,16,16)
  mu_0 <- c(16,16,16,16,16)
  mu_1 <- c(16,16+delta2_group[i],16+delta2_group[i]
            ,16+delta2_group[i],16+delta2_group[i])
  assess <- trial_assess_mu(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, 
                            tau, n1, n21, n22, p1, p2, p3,a)
  fwer1[i] <- assess$fwer1
  power[i] <- assess$power
  for (j in 1:3){
    error1[j,i] <- assess$error1[j]
  }
}
plot(delta2_group,fwer1, xlab = "Increase in mean value ",ylab = "FWER I")
plot(delta2_group,power ,xlab = "Increase in mean value ",ylab = "POWER")
plot(delta2_group,error1[1,], xlab = "Increase in mean value ",ylab = "PWER 2")
plot(delta2_group,error1[2,], xlab = "Increase in mean value ",ylab = "PWER 3")
plot(delta2_group,error1[3,], xlab = "Increase in mean value ",ylab = "PWER 4")

for (i in 1:length(delta2_group)){
  mu <- c(16,16+delta2_group[i],16,16,16)
  theta <- c(0.33,0.33+delta1_group[i],0.33,0.33,0.33)
  N[i] <- sample_size(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
}
plot(delta2_group, N, xlab = "Increase in mean value ",ylab = "Sample size")





###################### Scenario 2

# At phase 2:
theta_0 <- c(0.33,0.33,0.33,0.33,0.33,0.33)
mu_0 <- c(16,16,16,16,16,16)
theta_1 <- c(0.33,0.33+delta,0.33+delta,0.33+delta,0.33+delta,0.33+delta)
mu_1 <- c(16,16+delta2,16+delta2,16+delta2,16+delta2,16+delta2)
theta <- c(0.33,0.33+delta,0.33,0.33,0.33,0.33)
mu <- c(16,16+delta2,16,16,16,16)

sim_data <- function(n, theta, mu, sigma, tau){
  b <- rbinom(n,1,theta)
  y <- rnorm(n,mean = mu + tau*b,sd=sigma)
  return(cbind(b,y))
}


trial_design <- function(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  #phase 2 experimental treatment
  n_1 <- floor(n1*(1-a))
  #phase 2 control treatment
  n_00 <- floor((n1*sqrt(5))*(1-a))
  #phase 3 experimental treatment 
  n_21 <- floor((n_1+n21)*(1-a))
  n_22 <- floor((n_1+n22)*(1-a))
  #phase 3 control treatment 
  n_01 <- floor((n_00+0.5*n21)*(1-a))
  n_02 <- floor((n_00+0.5*n22)*(1-a))
  # generate samples
  b <- matrix(0,nrow = 5, ncol = n_1)
  A <- sim_data(n, theta[1], mu[1], sigma, tau)
  b0 <- A[1:n_00,1]
  y01 <- A[1:n_01,2]
  y02 <- A[1:n_02,2]
  y1 <- matrix(0,nrow = 5, ncol = n_21)
  y2 <- matrix(0,nrow = 5, ncol = n_22)
  for (i in 1:5){
    d <- sim_data(n, theta[i+1], mu[i+1], sigma, tau)
    b[i,] <- d[1:n_1,1]
    y1[i,] <- d[1:n_21,2]
    y2[i,] <- d[1:n_22,2]
  }
  # interim analysis
  test_result <- rep(-1, 5)
  test_statistic <- rep(0,5)
  for (i in 1:5){
    test_statistic[i] <- chisq.test(rbind(c(sum(b[i,]), n_1-sum(b[i,])),
                                          c(sum(b0), n_00-sum(b0))))$p.value
    if (test_statistic[i]<p1){
      test_result[i] <- 1 
    }
  }
  if (length(which(test_result==1, arr.ind=TRUE))>2){
    test_result[which.min(test_statistic)] <- -1
    test_statistic[which.min(test_statistic)] <- 1
    if (length(which(test_result==1, arr.ind=TRUE))>2){
      test_result[which.min(test_statistic)]<- -1
      test_statistic[which.min(test_statistic)] <- 1
      if (length(which(test_result==1, arr.ind=TRUE))>2){
        test_result[which.min(test_statistic)]<- -1
      }
    }
  }
  test_result_mu <- test_result
  #phase 3:
  test_result_reject <- which(test_result==1, arr.ind=TRUE)
  #  1 success
  if (length(test_result_reject)==1){
    test_result_mu[test_result_reject] <- t.test(y1[test_result_reject,], y01, alternative = "two.sided")$p.value <  p2
  }
  # 2 and more sucess
  if (length(test_result_reject)>1){
    for (j in test_result_reject){
      test_result_mu[j] <- t.test(y2[j,], y02, alternative = "two.sided")$p.value < p3
    }
  }
  return(test_result_mu)
}


trial_assess_mu <- function(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  f <- matrix(0,nrow=m,ncol=5)
  for (i in 1:m){
    f[i,] <- trial_design(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
  }
  k1 = k2 = 0
  True_null = False_null = numeric(0)
  for (j in 1:5){
    if (mu[j+1]==mu_0[j+1]){
      True_null <- append(True_null,j) 
    }
  }
  for (j in 1:5){
    if (mu[j+1]!=mu_0[j+1]){
      False_null <- append(False_null,j) 
    }
  }
  if (length(True_null)==0){
    fwer1 <- NA 
  }
  if (length(True_null)==1){
    f1 <- f[,True_null]
    for (i in 1:m){
      if (length(which(f1[i]==1, arr.ind=TRUE))>0){
        k1 <- k1+1
      }
    }  
  } 
  if (length(True_null)>1){
    f1 <- f[,True_null]
    for (i in 1:m){
      if (length(which(f1[i,]==1, arr.ind=TRUE))>0){
        k1 <- k1+1
      }
    }  
  } 
  if (length(False_null)==0){
    fwer2 <- NA 
  }
  if (length(False_null)==1){
    f1 <- f[,False_null]
    for (i in 1:m){
      if (length(which(f1[i]==1, arr.ind=TRUE))>0){
        k2 <- k2+1
      }
    }  
  } 
  if (length(False_null)>1){
    f1 <- f[,False_null]
    for (i in 1:m){
      if (length(which(f1[i,]==1, arr.ind=TRUE))>0){
        k2 <- k2+1
      }
    }  
  } 
  fwer1 <- k1/m
  fwer2 <- k2/m
  power <- 1-fwer2
  E1 <- numeric(0)
  P <- numeric(0)
  for (i in True_null){
    E1 <-append(E1,length(which(f[,i]==1, arr.ind=TRUE))/m)
  }
  for (i in False_null){
    P <-append(P,1- length(which(f[,i]!=1, arr.ind=TRUE))/m)
  }
  my_list <- list("fwer1" = fwer1,"power" = power, "error1"=E1, "power" = P)
  return(my_list)
}  


delta1_group <- seq(0.01,0.66,0.01) 
delta2_group <- 50*delta1_group
fwer1 <- rep(0,length(delta2_group))
power <- rep(0,length(delta2_group))
error1 <- matrix(0,nrow =4,ncol=length(delta2_group))
powerk <- rep(0,length(delta2_group))
N <- rep(0,length(delta2_group))
for (i in 1:length(delta2_group)){
  theta <- c(0.33,0.33+delta1_group[i],0.33,0.33,0.33,0.33)
  theta_0 <- c(0.33,0.33,0.33,0.33,0.33,0.33)
  theta_1 <- c(0.33,0.33+delta1_group[i],0.33+delta1_group[i]
               ,0.33+delta1_group[i],0.33+delta1_group[i],0.33+delta1_group[i])
  mu <- c(16,16+delta2_group[i],16,16,16,16)
  mu_0 <- c(16,16,16,16,16,16)
  mu_1 <- c(16,16+delta2_group[i],16+delta2_group[i]
            ,16+delta2_group[i],16+delta2_group[i],16+delta2_group[i])
  assess <- trial_assess_mu(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, 
                            tau, n1, n21, n22, p1, p2, p3,a)
  fwer1[i] <- assess$fwer1
  power[i] <- assess$power
  for (j in 1:4){
    error1[j,i] <- assess$error1[j]
  }
}
par(mfrow=c(1,1))
plot(delta2_group,fwer1, xlab = "Increase in mean value ",ylab = "FWER I")
plot(delta2_group,power ,xlab = "Increase in mean value ",ylab = "POWER")
plot(delta2_group,error1[1,], xlab = "Increase in mean value ",ylab = "PWER 2")
plot(delta2_group,error1[2,], xlab = "Increase in mean value ",ylab = "PWER 3")
plot(delta2_group,error1[3,], xlab = "Increase in mean value ",ylab = "PWER 4")
plot(delta2_group,error1[4,], xlab = "Increase in mean value ",ylab = "PWER 5")
sample_size <- function(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a){
  f <- matrix(0,nrow=m,ncol=5)
  for (i in 1:m){
    f[i,] <- trial_design(n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
  }
  N <- rep(0,m)
  for (i in 1:m){
    if (length(which(f[i,]==-1, arr.ind=TRUE))<4){
      N[i] <- 491+59
    }
    if (length(which(f[i,]==-1, arr.ind=TRUE))==4){
      N[i] <- 426+59
    }
    if (length(which(f[i,]==-1, arr.ind=TRUE))==5){
      N[i] <- 359+59
    }
  } 
  return(mean(N))
}

sample_size(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
for (i in 1:length(delta2_group)){
  mu <- c(16,16+delta2_group[i],16,16,16,16)
  theta <- c(0.33,0.33+delta1_group[i],0.33,0.33,0.33,0.33)
  N[i] <- sample_size(m,n, theta, mu, sigma, tau, n1, n21, n22, p1, p2, p3,a)
}
par(mfrow=c(1,1))
plot(delta2_group,N, xlab="Increase in mean value", ylab="Sample size" ) 


###################### Scenario 3

#rerun the codes before Scenario 1

h <- rep(0,100)
tau_group <- seq(1:100)
for (tau in tau_group){
  h[tau] <- sample_size(m,n, theta, mu, sigma,tau, n1, n21, n22, p1, p2, p3,a)
}
plot(tau_group,h, xlab="Increase in coefficient", ylab="Sample size" ) 

fwer1 <- rep(0,length(tau_group))
power <- rep(0,length(tau_group))
error1 <- matrix(0,nrow =3,ncol=length(tau_group))
powerk <- rep(0,length(tau_group))
for (tau in tau_group){
  theta_0 <- c(0.33,0.33,0.33,0.33,0.33)
  mu_0 <- c(16 + tau*theta_0[1],16+ tau*theta_0[2],16+ tau*theta_0[3],16+ tau*theta_0[4]
            ,16+ tau*theta_0[5])
  theta_1 <- c(0.33,0.33+delta,0.33+delta,0.33+delta,0.33+delta)
  mu_1 <- c(16+tau*theta_1[1],16+tau*theta_1[2],16+tau*theta_1[3],
            16+tau*theta_1[4],16+tau*theta_1[5])
  theta <- c(0.33,0.33+delta,0.33,0.33,0.33)
  mu <- c(16 + tau*theta[1],16+ tau*theta[2],16+ tau*theta[3],16+ tau*theta[4]
          ,16+ tau*theta[5])
  assess <- trial_assess_mu(m,n,theta_0, mu_0,theta_1,mu_1, theta, mu, sigma, 
                            tau, n1, n21, n22, p1, p2, p3,a)
  fwer1[tau] <- assess$fwer1
  power[tau] <- assess$power
  for (j in 1:3){
    error1[j,tau] <- assess$error1[j]
  }
}
plot(tau_group,fwer1, xlab = "Increase in coefficient ",ylab = "FWER I")
plot(tau_group,power ,xlab = "Increase in coefficient ",ylab = "POWER")
plot(tau_group,error1[1,], xlab = "Increase in coefficient ",ylab = "PWER 2")
plot(tau_group,error1[2,], xlab = "Increase in coefficient ",ylab = "PWER 3")
plot(tau_group,error1[3,], xlab = "Increase in coefficient ",ylab = "PWER 4")


