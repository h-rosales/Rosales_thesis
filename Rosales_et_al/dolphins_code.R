
library(sf)
library(latex2exp)

# creating data to input for stan

#read data
dat <- read.csv("~/Thesis/KF-Stan/Test/(New) CEE__2018_02_IDcor_alltracks.txt")


#convert character to POSIX
#dat$FILENAME <- gsub(".jpg", "", dat$FILENAME)
dat$FILENAME <- as.POSIXct(dat$FILENAME, format="%Y-%m-%d %H-%M-%S")

dat <- dat[dat$FILENAME >= "2018-06-15 16:28:54" & dat$FILENAME <= "2018-06-15 16:55:27", ]


#create time vector form beginning to end
time.1 <- as.POSIXct(min(dat$FILENAME), format="%Y-%m-%d %H:%M:%S"):as.POSIXct(max(dat$FILENAME), format="%Y-%m-%d %H:%M:%S")
time_all <- as.POSIXct(time.1, format="%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00")

#time when an animal was observed
time_obs <- sort(unique(dat$FILENAME))

#projection to meters
dolphins_sf <- st_as_sf(dat, coords = c("LONG", "LAT"), crs = "+proj=longlat")
dolphins_proj <- st_transform(dolphins_sf, crs = "+proj=utm +zone=10")
dat[, c("X", "Y")] <- st_coordinates(dolphins_proj)

#dat2 <- dat[order(dat$FILENAME, dat$ANIMAL_ID),]


#Compute Z transformation
unique_id <- sort(unique(dat$ANIMAL_ID))
Z <- vector(mode = "list", length = length(time_all))
names(Z) <- time_all
Z_t <- matrix(0, length(unique_id), length(unique_id))

for (i in 1:length(time_all)){
  if (time_all[i] %in% time_obs){
    Z_t <- as.integer(unique_id %in% dat$ANIMAL_ID[which(dat$FILENAME==time_all[i])])
    n <- matrix(0, nrow=sum(Z_t), ncol=length(unique_id))
    
    for (k in 1:sum(Z_t)){
      n[k, which(Z_t==1)[k]] <- 1
    }
    
    Z[[i]] <- n
  }
  else{
    Z[[i]] <- matrix(0, nrow = 1, ncol = length(unique_id))
  }
}

## adjust Z matrix for stan
#Z_stan <- do.call(rbind, Z)
Z_stan <- matrix(0, nrow = nrow(dat)+sum(!(time_all %in% time_obs)), ncol = length(unique_id))

k = 0
for(i in 1:length(time_all)){
  Z_stan[(k+1):(k+nrow(Z[[i]])),] <- Z[[i]]
  k = k + nrow(Z[[i]])
}


#observations
lat <- vector(mode = "list", length = length(time_obs))
long <- vector(mode = "list", length = length(time_obs))
names(lat) <- time_obs
names(long) <- time_obs

for (i in 1:length(time_obs)){
  y_pre <- dat[which(dat$FILENAME==time_obs[i]),c("ANIMAL_ID", "LAT", "LONG")]
  lat[[i]] <- y_pre[order(y_pre$ANIMAL_ID), "LAT"]
  long[[i]] <- y_pre[order(y_pre$ANIMAL_ID), "LONG"]
}

# fix observations lat and long to be in stan

total_lat <- unlist(lat, use.names = FALSE)

total_long <- unlist(long, use.names = FALSE)

#combining latitude and logitude values into one matrix
mu <- matrix(c(total_long, total_lat), nrow = length(total_lat), ncol = 2)


#same thing but for projects
#observations
lat_y <- vector(mode = "list", length = length(time_obs))
long_x <- vector(mode = "list", length = length(time_obs))
names(lat_y) <- time_obs
names(long_x) <- time_obs

for (i in 1:length(time_obs)){
  proj_pre <- dat[which(dat$FILENAME==time_obs[i]),c("ANIMAL_ID", "X", "Y")]
  lat_y[[i]] <- proj_pre[order(proj_pre$ANIMAL_ID), "Y"]
  long_x[[i]] <- proj_pre[order(proj_pre$ANIMAL_ID), "X"]
}

# fix observations lat and long to be in stan

total_lat_y <- unlist(lat_y, use.names = FALSE)

total_long_x <- unlist(long_x, use.names = FALSE)

#combining latitude and logitude values into one matrix
mu_proj <- matrix(c(total_long_x, total_lat_y), nrow = length(total_lat_y), ncol = 2)


#create p vector for indexing Z and mu
#this is number of animals identified at time t

p <- c(1:length(time_obs))

for (i in 1:length(time_obs)){
  p[i] <- length(lat[[i]])
}

# bottom code is just in case I need for all time observations top code is only
# for time that has been observed

# p <- c(1:length(time_all))
# 
# for (i in 1:length(time_all)){
#   p[i] <- sum(Z[[i]] == 1)
# }


#initialize values
##initial state values
#a <- matrix(c(921084, 3716583), nrow = length(unique_id), ncol = 2, byrow = T)
#P <- diag(10^6, length(unique_id))
P <- diag(10^5, length(unique_id))

a <- matrix(c(920964, 3716530), nrow = length(unique_id), ncol = 2, byrow = T)
# ##scalar parameters
# alpha <- 0.5
# beta <- 0.5
# 
# ##movement parameters
# sigma2 <- 1
# tau2 <- 1
# 
# ## identity matrix for easier calculation 
# I <- diag(1, length(unique_id))
```






```{r}
#social network
network <- function(phi, p1, J, t){
  big_boi_w <- vector(mode = "list", length = t)
  w <- matrix(0, J, J)
  for(s in 1:t){
    for(i in 1:J){
      k = i + 1
      while(k <= J){
        if(s >= 2){
          p12 = ((1-phi) * p1 * (1 - w[i,k])) + ((1 - (1 - phi) * (1 - p1)) * w[i,k])
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

# Sonar covariance
sonar_before <- sonar_during <- sonar_after <- rep(0, length(time_all))
start_sonar <- as.POSIXct("2018-06-15 16:42:50", origin = "1970-01-01 00:00:00")
end_sonar <- as.POSIXct("2018-06-15 16:52:50", origin = "1970-01-01 00:00:00")

#start_sonar <- as.POSIXct("2018-06-15 16:29:50", origin = "1970-01-01 00:00:00")
#end_sonar <- as.POSIXct("2018-06-15 16:30:45", origin = "1970-01-01 00:00:00")

sonar_before[1:sum(time_all < start_sonar)] <- rep(1, sum(time_all < start_sonar))

sonar_during[sum(time_all <= start_sonar):sum(time_all < end_sonar)] <- rep(1, sum(time_all >= start_sonar & time_all < end_sonar))

sonar_after[sum(time_all <= end_sonar):length(time_all)] <- rep(1, sum(time_all >= end_sonar))

x <- (cbind(sonar_before, sonar_during, sonar_after))



library(rstan)

# delta priors
delta_alpha_mean <- c(1,1,1)
delta_beta_mean <- c(0,0,0)

delta_alpha_cov <- diag(1/2.2, 3)
delta_beta_cov <- diag(1/2.2, 3)


#new stan6
set.seed(36)

#r data files
values1 <- c(0.7, 0.8, 0.9, 0.95, 0.99, 1.0)
values2 <- c(0.7, 0.8, 0.9, 0.95, 0.99, 1.0)


for(value_phi in values1){
  for(value_p1 in values2){
    
    super_w <- network(phi=value_phi, p1=value_p1, length(unique_id), length(time_all))

    w_matrix <- matrix(unlist(super_w),
                       nrow = length(unique_id) * length(time_all), 
                       length(unique_id), byrow = F)
    
    w_stan <- array(w_matrix, dim = c(length(unique_id), length(unique_id), length(time_all)))
    w_stan <- aperm(w_stan, c(3, 1, 2))
    
    stan_data <- list(time_all = length(time_all), time_obs = length(time_obs), zrow = nrow(Z_stan), 
                      murow = nrow(mu_proj), J = length(unique_id), p = p,
                      Z = Z_stan, a = a, P = P, mu = mu_proj, w = w_stan, x = x, delta_alpha_mean = delta_alpha_mean,
                      delta_beta_mean = delta_beta_mean, delta_alpha_cov = delta_alpha_cov,
                      delta_beta_cov = delta_beta_cov)
    
    results <- stan(file = "~/Thesis/KF-Stan/Test/KF Stan7.stan",
                    data = stan_data, chains = 4, cores = 4, iter = 1000, algorithm = "NUTS")
    #init=list(list(tau2=0.01, sigma2=8), list(tau2=0.1, sigma2=8), list(tau2=0.5, sigma2=8), list(tau2=1, sigma2=8)))
    results2 <- extract(results, permuted = FALSE, inc_warmup = FALSE)
    
    save(results2, file = paste0(value_phi, value_p1, "stan7.rda"))
    
  }
}

times <- time_all
TM <- length(times)
transp_level <- 0.02
expit <- function(x) 1/(1 + exp(-x))
setEPS()
postscript(file="1.01.0real.eps", height = 4)
pdf(file="1.01.0real.pdf", height = 4)
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
  # if(param == "pi"){
  #   lines(times, rowMeans(prop_edges) * diff(range(chains)) + min(chains), lwd = 2)
  #   # matplot(prop_edges * diff(range(chains)) + min(chains),
  #   #         type = "l", col = scales::alpha("gray", 0.01), lty = 1, add = T)
  # }
}
dev.off()



library(MASS)



```{r}
#posterior mean image
values <- c(0.7, 0.8, 0.9, 0.95, 0.99, 1.0)
sds_alpha32 <- sds_alpha21 <- sds_beta32 <- sds_beta21 <- means_beta32 <- means_beta21 <- means_alpha21 <- means_alpha32 <- matrix(NA, 6, 6)

dar <- 1
dar2 <- 1
expit <- function(x) 1/(1 + exp(-x))

#phi are rows and p1 are columns
for(phi in values){
  for(p1 in values){
    load(paste0("~/Thesis/KF-Stan/Test/post_dist/",phi, p1, "stan7.rda"))
    means_beta21[dar, dar2] <- mean(expit(results2[ , , 'delta_beta[2]']) - expit(results2[ , , 'delta_beta[1]']))
    means_beta32[dar, dar2] <- mean(expit(results2[ , , 'delta_beta[3]']) - expit(results2[, , 'delta_beta[2]']))
    means_alpha21[dar, dar2] <- mean(expit(results2[ , , 'delta_alpha[2]']) - expit(results2[ , , 'delta_alpha[1]']))
    means_alpha32[dar, dar2] <- mean(expit(results2[ , , 'delta_alpha[3]']) - expit(results2[ , , 'delta_alpha[2]']))
    
    sds_beta21[dar, dar2] <-  sd(expit(results2[ , , 'delta_beta[2]']) - expit(results2[ , , 'delta_beta[1]']))
    sds_beta32[dar, dar2] <-  sd(expit(results2[ , , 'delta_beta[3]']) - expit(results2[ , , 'delta_beta[2]']))
    sds_alpha21[dar, dar2] <-  sd(expit(results2[ , , 'delta_alpha[2]']) - expit(results2[ , , 'delta_alpha[1]']))
    sds_alpha32[dar, dar2] <-  sd(expit(results2[ , , 'delta_alpha[3]']) - expit(results2[ , , 'delta_alpha[2]']))
    #print(mean(results2[ , , 'delta_beta[2]']))
    dar2 <- dar2 + 1
  }
  dar <- dar + 1
  dar2 <- 1
}

#Mean
image(values, values, means_beta21, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "PRGn", rev = TRUE), zlim = c(-4.3, 4.3), main = "Means Diff Beta during and before")
image(values, values, means_beta32, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "PRGn", rev = TRUE), zlim = c(-4.3, 4.3), main = "Means Diff Beta after and during")

image(values, values, means_alpha21, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "YlOrRd", rev = TRUE), zlim = c(0, 7.42), main = "Means Diff alpha during and before")
image(values, values, means_alpha32, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "YlOrRd", rev = TRUE), zlim = c(0, 7.42), main = "Means Diff alpha after and during")

#SD
image(values, values, sds_beta21 / max(sds_beta21), xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Beta during and after")
image(values, values, sds_beta32 / max(sds_beta32), xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Beta after and during")


#plot(means_beta21[,4], ylim = c(0.1, 0.6))
#plot(means_beta21[,7], ylim = c(0.1, 0.6))

image(values, values, means_beta32, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "YlOrRd", rev = TRUE), zlim =  main = "Means Diff Beta after and during")

matplot(values, t(means_beta32), type = "l", xlab = "p1", col = hcl.colors(11, "PRGn", rev = TRUE))
matplot(values, t(means_alpha32), type = "l", xlab = "p1", col = hcl.colors(11))

t(means_beta32)
range(means_beta32)
```


```{r}
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


# setEPS()
# postscript(file="image_post_diff.eps", height = 4)
# pdf(file="image_post_diff.pdf", height = 4)

setEPS()
postscript(file=paste0("mean_beta_1_to_2.eps"), height = 4)
pdf(file="mean_beta_1_to_2.pdf", height = 4)
#mean beta 1 to 2
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, means_beta21, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "PRGn", rev = TRUE), zlim = range(means_beta21, means_beta32), main = "Means Diff Beta during and before")
# par(mar=c(2,3,1,3))
# image.scale(means_beta21, zlim = range(means_beta21, means_beta32), col = hcl.colors(16, "PRGn", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("mean_beta_2_to_3.eps"), height = 4)
pdf(file="mean_beta_2_to_3.pdf", height = 4)
#mean beta 2 to 3
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, means_beta32, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "PRGn", rev = TRUE), zlim = range(means_beta21, means_beta32), main = "Means Diff Beta after and during")
par(mar=c(2,3,1,3))
image.scale(means_beta32, zlim = range(means_beta21, means_beta32), col = hcl.colors(16, "PRGn", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("mean_alpha_1_to_2.eps"), height = 4)
pdf(file="mean_alpha_1_to_2.pdf", height = 4)

#mean alpha 1 to 2
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, means_alpha21, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "RdBu", rev = TRUE), zlim = range(means_alpha21, means_alpha32), main = "Means Diff Alpha during and before")
# par(mar=c(2,3,1,3))
# image.scale(means_alpha21, zlim = range(means_alpha21, means_alpha32), col = hcl.colors(16, "RdBu", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("mean_alpha_2_to_3.eps"), height = 4)
pdf(file="mean_alpha_2_to_3.pdf", height = 4)
#mean alpha 2 to 3
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, means_alpha32, xlab = "phi", ylab = "p1",
      col = hcl.colors(16, "RdBu", rev = TRUE), zlim = range(means_alpha21, means_alpha32), main = "Means Diff Alpha after and during")
par(mar=c(2,3,1,3))
image.scale(means_alpha32, zlim = range(means_alpha21, means_alpha32), col = hcl.colors(16, "RdBu", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("sd_beta_1_to_2.eps"), height = 4)
pdf(file="sd_beta_1_to_2.pdf", height = 4)
#SD beta 1 to 2
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, sds_beta21, xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Beta during and before")
# par(mar=c(2,3,1,3))
# image.scale(sds_beta21 / max(sds_beta21), zlim = c(0, 1), col = hcl.colors(16, "PuBu", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("sd_beta_2_to_3.eps"), height = 4)
pdf(file="sd_beta_2_to_3.pdf", height = 4)
#SD beta 3 to 2
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, sds_beta32, xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Beta after and during")
par(mar=c(2,3,1,3))
image.scale(sds_beta32 / max(sds_beta32), zlim = c(0, 1), col = hcl.colors(16, "PuBu", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("sd_alpha_1_to_2.eps"), height = 4)
pdf(file="sd_alpha_1_to_2.pdf", height = 4)
#SD alpha 1 to 2
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, sds_alpha21, xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Alpha during and before")
# par(mar=c(2,3,1,3))
# image.scale(sds_alpha21 / max(sds_alpha21), zlim = c(0, 1), col = hcl.colors(16, "PuBu", rev = TRUE), horiz = T)

dev.off()

setEPS()
postscript(file=paste0("sd_alpha_2_to_3.eps"), height = 4)
pdf(file="sd_alpha_2_to_3.pdf", height = 4)
#SD alpha 2 to 3
layout(matrix(c(1:2), nrow=2, ncol=1), widths=c(4,4,1), heights=c(4,1))
par(mar=c(5,5,1,2))
image(values, values, sds_alpha32, xlab = "phi", ylab = "p1", 
      col = hcl.colors(16, "PuBu", rev = TRUE), zlim = c(0, 1), main = "SD Diff Alpha after and during")
par(mar=c(2,3,1,3))
image.scale(sds_alpha32 / max(sds_alpha32), zlim = c(0, 1), col = hcl.colors(16, "PuBu", rev = TRUE), horiz = T)

dev.off()






