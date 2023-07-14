library(latex2exp)
library(MASS)
#social network
network <- function(phi, p1, J, t){
  big_boi_w <- vector(mode = "list", length = t)
  w <- matrix(0, J, J)
  for(s in 1:t){
    for(i in 1:J){
      k = i + 1
      while(k <= J){
        if(s >= 2){
          p12 = ((1-phi) * p1 * (1 - w[i,k])) + (1 - (1 - phi) * (1- p1)) * w[i,k]
          w[k,i] = w[i,k] = sample(c(1,0), 1, prob = c(p12, 1-p12))
          k = k + 1
        }
        else{
          w[k,i] = w[i,k] = sample(c(1,0), 1, prob = c(p1, 1-p1))
          k = k + 1
        }
      }
    }
    big_boi_w[[s]] = w
  }
  return(big_boi_w)
}


for(l in c(0.8, 0.9, 1.0)){
  for(k in c(0.8, 0.9, 1.0)){
    for(b in c(0.001, 0.005)){
      for(al in c(0.9, 0.95)){

#initialize variables
J <- 8
I <- diag(1, J)
w_bar <- matrix(NA, J, J)
sigma2 <- 1
time_obs <- 1000
beta <- 0.005
alpha <- 0.95
omega <- matrix(NA, J, J)
eta <- rep(0, J)
mean_eta <- rep(0, J)
miss <- 34
time_all <- miss + time_obs

mu <- matrix(NA, time_all, 2*J)
mu[1, seq(1, (2*J)-1, 2)] <- rnorm(J, 33.50589, .001)
mu[1, seq(2, 2*J, 2)] <- rnorm(J, -118.4686, .001)

# social network
super_w <- network(phi=1, p1=1, J, time_all)
w_matrix <- matrix(unlist(super_w),
                   nrow = J * time_all,
                   J, byrow = F)
w_stan <- array(w_matrix, dim = c(J, J, time_all))
w_stan <- aperm(w_stan, c(3, 1, 2))

for (t in 1:(time_all-1)){
  for (i in 1:J){
    for (j in 1:J){
      w_bar[i,j] <- w_stan[t,i,j] / sum(w_stan[t,i,])
      if (i == j){
        omega[i, j] <- sum(w_stan[t,i,])
      }
      else{
        omega[i, j] <- -alpha * w_stan[t,i,j]
      }
    }
  }
  
  # create Attraction matrix
  A <- (1 - beta)*I + beta*w_bar
  
  
  eta <- mvrnorm(1, mean_eta, sigma2 * solve(omega))
  
  # if(t <= 10){
  #   mu[t+1, seq(1, (2*J)-1, 2)] <- A%*%mu[t, seq(1,(2*J)-1,2)] + matrix(eta)
  #   mu[t+1, seq(2, 2*J, 2)] <- A%*%mu[t, seq(2,2*J,2)] + matrix(eta)
  # }
  mu[t+1, seq(1, (2*J)-1, 2)] <- A%*%mu[t, seq(1,(2*J)-1,2)] + matrix(eta)
  
  eta <- mvrnorm(1, mean_eta, sigma2 * solve(omega))
  mu[t+1, seq(2, 2*J, 2)] <- A%*%mu[t, seq(2,2*J,2)] + matrix(eta)
  
}



# delta priors
delta_alpha_mean <- c(1,1,1)
delta_beta_mean <- c(0,0,0)

delta_alpha_cov <- diag(1/2.2, 3)
delta_beta_cov <- diag(1/2.2, 3)




#create Z matrix
Z <- matrix(rep(diag(1, J), time_obs), nrow = J*time_obs, ncol = J, byrow = T)
Z <- rbind(Z[1:(J*15),], matrix(rep(0, J), miss, J), Z[-(1:(J*15)),])
#how many animals identified at time t
p <- c(rep(J-3, 10), rep(J, time_obs-10))
#p <- c(rep(J, time_obs))

#covariates
x1 <- c(rep(1, 330), rep(0, 670+miss))
x2 <- c(rep(0, 330), rep(1, 330), rep(0, 340+miss))
x3 <- c(rep(0, 660), rep(1, 340+miss))
x <- cbind(x1, x2, x3)

# fix mu
mu_lon <- matrix((t(mu[,seq(1,(2*J)-1,2)])), time_all*J, 1, byrow = T)
mu_lat <- matrix((t(mu[, seq(2, 2*J, 2)])), time_all*J, 1, byrow = T)

mu2 <- cbind(mu_lon, mu_lat)
mu2 <- mu2[-((J*15+1):((J*miss)+(J*15))),]

#partial missing observation
pm <- (10:1)*J
for(i in 1:10){
 mu2 <- mu2[-((pm[i]-2):pm[i]),]
 Z <- Z[-((pm[i]-2):pm[i]),]
}


# starting place
a <- cbind(rep(mean(mu2[1:20,1]), J), rep(mean(mu2[1:20,2]), J))
P <- diag(var(mu2[1:20,1]), J)


matplot(mu[,seq(1,(2*J)-1,2)], main="Long",
        lty = 1, type = "l")
matplot(mu[,seq(2,2*J,2)], main="Lat",
        lty = 1, type = "l")

library(rstan)


# super_w <- network(phi=l, p1=k, J, time_obs)
# w_matrix <- matrix(unlist(super_w),
#                    nrow = J * time_obs, 
#                    J, byrow = F)
# w_stan <- array(w_matrix, dim = c(J, J, time_obs))
# w_stan <- aperm(w_stan, c(3, 1, 2))
stan_data <- list(time_all = time_all, time_obs = time_obs, zrow = nrow(Z), 
                  murow = nrow(mu2), J = J, p = p,
                  Z = Z, a = a, P = P, mu = mu2, w = w_stan, 
                  x = x, delta_alpha_mean = delta_alpha_mean,
                  delta_beta_mean = delta_beta_mean, 
                  delta_alpha_cov = delta_alpha_cov,
                  delta_beta_cov = delta_beta_cov)
# stan_data <- list(time_obs = length(time_obs), zrow = nrow(Z_stan), 
#                       murow = nrow(mu), J = length(unique_id), p = p,
#                       Z = Z_stan, a = a, P = P, mu = cbind(dat$X, dat$Y), w = w_stan, x = x, delta_alpha_mean = delta_alpha_mean,
#                       delta_beta_mean = delta_beta_mean, delta_alpha_cov = delta_alpha_cov,
#                       delta_beta_cov = delta_beta_cov)

results <- stan(file = "~/Thesis/KF-Stan/Test/KF Stan7.stan",
                data = stan_data, chains = 4, cores = 4, iter = 1000, 
                algorithm = "NUTS")
results2 <- extract(results, permuted = FALSE, inc_warmup = FALSE)

save(results2, file = paste0("alpha", al, "beta", b, l, k, "simulation.rda"))

load(file = paste0("alpha0.95beta0.0050.80.8simulation.rda"))


start_sonar <- 330
end_sonar <- 660


# Used for creating Zaineb's plots
times <- 1:time_all
TM <- time_all
transp_level <- 0.02
expit <- function(x) 1/(1 + exp(-x))
setEPS()
postscript(file="simulation.eps", height = 4)
pdf(file="simulation.pdf", height = 4)
layout(matrix(1:2, 2, 1))
par(mar = c(1, 5.1, 3, 0.1), oma = c(1, 0, 0, 0))
for(param in c("alpha", "beta")){
  xlab <- ""; xaxt <- "n"; ylab <- TeX(paste0("$\\", param, "_t$"))
  chains <- expit(results2[, 1, paste0( "delta_", param, "[",1:3,"]")])
  main <- c("alignment", "attraction")[which(param == c("alpha", "beta"))]
  plot(range(times), 0:1, type = "n", ylab = ylab, 
       xlab = xlab, xaxt = xaxt, ylim = range(chains), cex.lab = 2, main = main)
  segments(x0 = times[1], y0 = chains[, 1], x1 = start_sonar, 
           col = scales::alpha("black", transp_level))
  segments(x0 = start_sonar, y0 = chains[, 2], x1 = end_sonar, 
           col = scales::alpha("dark blue", transp_level))
  segments(x0 = end_sonar, y0 = chains[, 3], x1 = times[TM], 
           col = scales::alpha("dark green", transp_level))
  abline(h=get(param), lty=2, lwd=3, col="black")
  # if(param == "pi"){
  #   lines(times, rowMeans(prop_edges) * diff(range(chains)) + min(chains), lwd = 2)
  #   # matplot(prop_edges * diff(range(chains)) + min(chains),
  #   #         type = "l", col = scales::alpha("gray", 0.01), lty = 1, add = T)
  # }
}
dev.off()
  }
}
}
}








