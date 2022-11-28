##--------------------------------------------------------------------------------------------------------
## SCRIPT : Killer Whales
##
## Authors : Matthieu Authier & Emmanuelle Cam
## Last update : 2022-09-16
##
## R version 4.0.5 (2021-03-31) -- "Shake and Throw"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

### get the missing packages
what_u_need <- c("tidyverse", "nimble", "MCMCvis", "glue")
cran_packages <- what_u_need[!(what_u_need %in% installed.packages())]
if(length(cran_packages) != 0) {
  lapply(cran_packages, install.packages, dependencies = TRUE)
}
lapply(what_u_need, library, character.only = TRUE)

### foreword
case_study <- "Crozet Killer Whale"
glue('This is the {case_study} case study\n',
     '\tIts primary aim is pedagogical and by no means prescriptive.\n',
     'The analysis provided is NOT necessarily correct given the assumptions we will make.'
     )
######### 

  # use your own path to locate the data file on your own computer
  # here the path locates the file in the 'data/killer whale' folder
  # perhaps your data file is another folder on your computer

######### Get the data file

orca <- readxl::read_excel(path = "killer whale/CMR-orcas-Matrix-Matriline-Size.xlsx",
                           sheet = 1
                           )
orca=CMR_orcas_Matrix_Matriline_Size
# Observations are called y

y <- orca[, 2:26]

View(y)

#Number of occasions

T <- ncol(y)

# Find the occasion when the orca was first contacted

first <- apply(y, 1, function(x) min(which(x !=0)))
first

# Let's get the number of individuals that have at least a detection.

## ----------------------------------------------------------------------------------------------
nobs <- sum(apply(y,1,sum) != 0)
nobs


# put 'NA's' before the individual was captured and marked (and released)

for (i in 1:nobs){
  if(first[i] > 1) y[i, 1:(first[i]-1)] <- NA
}

View(y)

# Then put the constants in a list. 
## ----------------------------------------------------------------------------------------------
my.constants <- list(N = nrow(y),      # nb of individuals
                     T = ncol(y),      # nb of occasions
                     first = first,    # occasions of first capture
                     PODSIZE = orca$SIZE,# pod size
                     HWI=as.numeric(orca$HWI) #HWI
                     )    
my.constants

# Put the data in a list. 
#Note that we add 1 to the data to have 1 for non-detections and 2 for detections. We do not keep 0's in the example, but only codes >0
#You may use the coding you prefer of course, you will just need to adjust the $\Omega$ and $\Gamma$ matrices in the model above.  
## ----------------------------------------------------------------------------------------------
my.data <- list(y = y + 1)

# latent states must be initialized
# Building initial values for latent state
# this can be particularly painful
# individuals should not be immortal, but death should only be considered after the last resighting
# Constraints have to be imposed on initial latent states that are....
#....really painful to specify.....

zinits <- function(y) {
  T <- length(y)
  first <- min(which(y != 0)) # first resighting
  last <- max(which(y == 1)) # last resighting
  z <- rep(NA, T) # create an object of length T that contains NA's
  if(last == T) {
    z[first:T] <- 1 # if resighted at the last occasion: the individual remains alive from first to last
  } else {
    if(last == (T - 1)) { # if resighted for the last time at T-1
      z[first:(T - 1)] <- 1 # the individual remains alive from first to T-1
      z[T] <- sample(1:2, size = 1) # after that, the individual either survives or dies at T with prob =0.5
    } else{ # other situations
      death <- sample((last + 1):T, size = 1) # select the occurrence of latent state 'dead' at random between last1 and T
      z[first:(death - 1)] <- 1 # before death the individual remains alive
      z[death:T] <- 2 # the individual remains dead up to T
    }
  }
  return(z)
}

############################################################################
#POLYNOMIAL MODEL


initial.values <- function() list(muphi = rnorm(1),
                                  betaphi1 = rnorm(1),
                                  betaphi2 = rnorm(1),
                                  mup = rnorm(1),
                                  eps = runif(1),
                                  z = t(apply(y, 1, zinits))
)
# model with logit(phi) = intercept + betaphi1*PODSIZE + betaphi2*PODSIZE^2
hmm.phiLINptRE <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # survival process
  for(i in 1:N) {
    # polynomial relationship with podsize
    logit(phi[i]) <- muphi + (betaphi1) * (PODSIZE[i] - 1) + (betaphi2) * (PODSIZE[i]-1)^2
    gamma[1,1,i] <- phi[i]      
    gamma[1,2,i] <- 1 - phi[i]  
    gamma[2,1,i] <- 0           
    gamma[2,2,i] <- 1           
  }
  # detection process
  for (t in 1:(T-1)){
    logit(p[t]) <- eps[t] # eps is the random effect (epsilon)
    eps[t] ~ dnorm(mup, sd = sdeps) 
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  # priors
  muphi ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  betaphi1 ~ dnorm(0, sd= 1.0)
  betaphi2 ~dpois(0) 
  mup ~ dnorm(0, sd= 1.5)
  sdeps ~ dunif(0,10) #prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (t in (first[i]+1):T){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i])
      y[i, t] ~ dcat(omega[z[i, t], 1:2, t-1])
    }
  }
})

parameters.to.save <- c("phi", "muphi", "betaphi1","betaphi2", 
                        "mup", "p", "eps", "sdeps"
)

n.iter <- 1000
n.burnin <- 500
nc <- 2
nt <- 2

orca_mcmc2 <- nimbleMCMC(code = hmm.phiLINptRE, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = nc,
                        thin = nt,
                        summary = TRUE, WAIC = TRUE,
)

###
MCMCsummary(object = orca_mcmc2$samples, 
            round = 2
)

MCMCplot(orca_mcmc2$samples)

MCMCtrace(orca_mcmc2$samples, 
          pdf=FALSE, ind=TRUE,
          params = c('muphi', 'betaphi1',"betaphi2", 'mup', 'sdeps')
)

### posterior analysis: extract joint posterior distribution
n_final <- (n.iter - n.burnin) / nt

posterior <- orca_mcmc2$samples %>% 
  do.call("rbind", .) %>% 
  as.data.frame() %>%
  mutate(chain_id = rep(1:nc, each = n_final),
         iteration = rep(1:n_final, nc)
  )

### prediction for pod size between 1 and 11
podsize <- 0:10
X <- matrix(1, nrow = length(podsize), ncol = 3)
X[, 3] <- podsize

### 
phi_pred <- posterior %>%
  select("muphi", "betaphi1","betaphi2") %>%
  as.matrix() %*% t(X) %>% # matrix computation
  plogis() %>% # back to the probability scale
  as.data.frame() 

names(phi_pred) <- paste("phi", podsize, sep = "")

phi_pred %>%
  as.data.frame() %>%
  pivot_longer(col = starts_with("phi"),
               values_to = "value",
               names_to = "param"
  ) %>%
  mutate(pod = str_sub(param, start = 4),
         pod = as.numeric(pod)
  ) %>%
  group_by(param, pod) %>%
  summarize(posterior_mean = mean(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)
  ) %>%
  ggplot(aes(x = pod, y = posterior_mean, 
             ymin = lower, ymax = upper)
  ) +
  geom_ribbon(fill = "midnightblue", alpha = 0.2) +
  geom_line(color = "midnightblue") +
  scale_x_continuous(name = "PodSize", breaks = 0:20,
                     labels = 1:21
  ) +
  scale_y_continuous(name = quote(phi), 
                     breaks = seq(0, 1, 0.1)
  ) +
  coord_cartesian(ylim = c(0.0, 1.0)) +
  theme_bw()

#####SHIT MODEL


initial.values <- function() list(muphi = rnorm(1),
                                  betaphi1 = rnorm(1),
                                  betaphi2 = rnorm(1),
                                  mup = rnorm(1),
                                  eps = runif(1),
                                  z = t(apply(y, 1, zinits))
)
# model with logit(phi) = intercept + betaphi1*PODSIZE + betaphi2*PODSIZE^2
hmm.phiLINptRE <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # survival process
  for(i in 1:N) {
    # polynomial relationship with podsize
    logit(phi[i]) <- muphi + (betaphi1) * (PODSIZE[i] - 1) + (betaphi2) * (PODSIZE[i] - 1)^2
    gamma[1,1,i] <- phi[i]      
    gamma[1,2,i] <- 1 - phi[i]  
    gamma[2,1,i] <- 0           
    gamma[2,2,i] <- 1           
  }
  # detection process
  for (t in 1:(T-1)){
    logit(p[t]) <- eps[t] # eps is the random effect (epsilon)
    eps[t] ~ dnorm(mup, sd = sdeps) 
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  # priors
  muphi ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  betaphi1 ~ dnorm(0, sd= 1.0)
  betaphi2 ~dnorm(0, sd= 5.0) 
  mup ~ dnorm(0, sd= 1.5)
  sdeps ~ dunif(0,10) #prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (t in (first[i]+1):T){
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:2, i])
      y[i, t] ~ dcat(omega[z[i, t], 1:2, t-1])
    }
  }
})

parameters.to.save <- c("phi", "muphi", "betaphi1","betaphi2", 
                        "mup", "p", "eps", "sdeps"
)

n.iter <- 10000
n.burnin <- 5000
nc <- 2
nt <- 4

orca_mcmc2 <- nimbleMCMC(code = hmm.phiLINptRE, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = nc,
                         thin = nt,
                         summary = TRUE, WAIC = TRUE,
)

###
MCMCsummary(object = orca_mcmc2$samples, 
            round = 2
)

MCMCplot(orca_mcmc2$samples)

MCMCtrace(orca_mcmc2$samples, 
          pdf=FALSE, ind=TRUE,
          params = c('muphi', 'betaphi1',"betaphi2", 'mup', 'sdeps')
)

### posterior analysis: extract joint posterior distribution
n_final <- (n.iter - n.burnin) / nt

posterior <- orca_mcmc2$samples %>% 
  do.call("rbind", .) %>% 
  as.data.frame() %>%
  mutate(chain_id = rep(1:nc, each = n_final),
         iteration = rep(1:n_final, nc)
  )

### prediction for pod size between 1 and 11
podsize <- 0:10
X <- matrix(1, nrow = length(podsize), ncol = 3)
X[, 2] <- podsize

### 
phi_pred <- posterior %>%
  select("muphi", "betaphi1","betaphi2") %>%
  as.matrix() %*% t(X) %>% # matrix computation
  plogis() %>% # back to the probability scale
  as.data.frame() 

names(phi_pred) <- paste("phi", podsize, sep = "")

phi_pred %>%
  as.data.frame() %>%
  pivot_longer(col = starts_with("phi"),
               values_to = "value",
               names_to = "param"
  ) %>%
  mutate(pod = str_sub(param, start = 4),
         pod = as.numeric(pod)
  ) %>%
  group_by(param, pod) %>%
  summarize(posterior_mean = mean(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975)
  ) %>%
  ggplot(aes(x = pod, y = posterior_mean, 
             ymin = lower, ymax = upper)
  ) +
  geom_ribbon(fill = "midnightblue", alpha = 0.2) +
  geom_line(color = "midnightblue") +
  scale_x_continuous(name = "PodSize", breaks = 0:20,
                     labels = 1:21
  ) +
  scale_y_continuous(name = quote(phi), 
                     breaks = seq(0, 1, 0.1)
  ) +
  coord_cartesian(ylim = c(0.0, 1.0)) +
  theme_bw()


