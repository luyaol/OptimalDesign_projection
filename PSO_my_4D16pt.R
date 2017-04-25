require(lhs)

####################################################################
## Define function for minimum distance
min_d = function(D){
  D = as.matrix(D)
  size <- nrow(D)
  dimsn <- ncol(D)
  distM <- dist(D)
  min_d <- min(distM)#/sqrt(dimsn)
  return(min_d)
}

fctrn_roshan_paper = function(Design){
  try_d = as.matrix(Design)
  p = ncol(try_d)
  n = nrow(try_d)
  ctrn_tmp = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      ctrn_tmp = ctrn_tmp +  exp(-2*sum(log(abs((try_d[i,]-try_d[j,])))))
    }
  }
  ctrn = ctrn_tmp
  return(-ctrn)
}


ctrn_my = function(Design,p_act=.25){
  nsample = nrow(Design)
  ndim = ncol(Design)
  try_d = Design
  ctrn_new = 1/0.5855625*p_act^3*(1-p_act)*(min_d(try_d[,-3]) + min_d(try_d[,-1]) + min_d(try_d[,-2]) + min_d(try_d[,-4])) + 
    15 * p_act*(1-p_act)^3*(min_d(try_d[,1]) + min_d(try_d[,2])+min_d(try_d[,3])+min_d(try_d[,4])) +
    3 * p_act^2*(1-p_act)^2*(min_d(try_d[,c(1,2)]) + min_d(try_d[,c(2,3)])+min_d(try_d[,c(3,4)]) +
                               min_d(try_d[,c(2,4)]) +min_d(try_d[,c(1,4)]) + min_d(try_d[,c(1,3)])) +
    p_act^4*min_d(try_d)
  return(ctrn_new)
}


######################################################################
## Start Particle Swarm algorithm

## Design settings:
D = 4
nsample = 16


## Settings:
theta <- .5
alpha <-2
beta <-2

## Set seed for PSO
N_seed <- 5000

## Set loop for PSO
nloop <- 500


###################################################################
## Function to find updated best design

find.c <- function(Design, num_seed, n_sample){
  c <- rep(0,num_seed)
  for(k in 1:num_seed){
    Ds <- matrix(Design[,k],nrow=n_sample,byrow=FALSE)
    c[k] <- ctrn_my(Ds)
  }
  ind <- which.max(c)
  out_put <- list(c=c,ind=ind)
  return(out_put)
}

check.bound <- function(x,bd = 0.25){
  x <- (x > bd) * bd + ( x < -bd ) * ( -bd ) + (x>=-bd)*(x<=bd) * x
  return(x)
}

check.dist <- function(D1,n_s=nsample){
  Design <- matrix(D1,nrow=n_s,byrow=F)
  return(ctrn_my(Design))
}


PSO_find_Mm = function(Seed){
  Design_seed <- matrix(runif(D*nsample*N_seed),ncol=N_seed)
  # Design_seed <- matrix(0,ncol=N_seed,nrow=D*nsample)
  # set.seed(Seed)
  # for(i in 1:N_seed){
  #   LHS_D <- as.vector(randomLHS(nsample,D))
  #   Design_seed[,i] <- LHS_D
  # }
  
  
  ##################################################################
  ##################################################################
  ## Start PSO
  
  Z0 <- Design_seed
  
  ###################
  ## Initial values:
  
  c0 <- rep(0,N_seed)
  
  for(j in 1:N_seed){
    D0 <- matrix(Z0[,j],nrow=nsample,byrow=FALSE)
    c0[j] <- ctrn_my(D0)
  }
  
  
  init_opt_ind <- which.max(c0)
  
  init_c <- find.c(Design = Design_seed , num_seed=N_seed,n_sample = nsample)
  
  g0 <- Z0[,init_opt_ind]
  p0 <- Z0
  expand_pop = N_seed / 2
  expand_size = 0.2
  update_count=0
  for(i in 1:nloop){
    
    if(i > 1){
      
      optm_c_p <- find.c(p,num_seed=N_seed,n_sample=nsample)
      V <- check.bound(theta * V + alpha * runif(1) * (g - Z) + beta * runif(1) * (p - Z))
      Z_ori = Z
      Z <- V + Z
      Z <- (Z<0)*0 + (Z>1)*1 + (Z>=0)*(Z<=1)*Z
      V = Z - Z_ori
      optm_c <- find.c(Z,num_seed=N_seed,n_sample=nsample)
      p <- (optm_c$c > optm_c_p$c) * Z + (optm_c$c <= optm_c_p$c) * p
      update_ans = optm_c$c[optm_c$ind] > check.dist(g)
      g <- (optm_c$c[optm_c$ind] > check.dist(g)) * Z[,optm_c$ind] + (optm_c$c[optm_c$ind] <= check.dist(g)) * g
             }
    
    if(i==1){
      V <- check.bound(alpha * runif(1) * (g0 - Z0) + beta * runif(1) * (p0 - Z0))
      Z <- V + Z0
      Z <- (Z<0)*0 + (Z>1)*1 + (Z>=0)*(Z<=1)*Z
      V <- Z - Z0
      optm_c <- find.c(Z,num_seed=N_seed,n_sample=nsample)
      p <- (optm_c$c > init_c$c) * Z + (optm_c$c <= init_c$c) * Z0
      g <- (optm_c$c[optm_c$ind] > check.dist(g0)) * Z[,optm_c$ind] + (optm_c$c[optm_c$ind] <= check.dist(g0)) * g0
    }
    if(i>1){
      #(i%%2==0){
      print(paste("optimal design criteria is ",check.dist(g)))
      print(paste("best criteria in the pop is ", optm_c$c[optm_c$ind]))
      print(paste("update? ", update_ans))
      print("loop ");print(i)
      print(Sys.time())
      #DG <- matrix(g,ncol=2,byrow=F)
      #plot(DG)
      #DG <- matrix(g,ncol=3,byrow=F)
      #plot(DG[,1:2]);plot(DG[,2:3])
    }
    
  }
  DG = matrix(g,ncol=D,byrow=F)
  return(DG)
}

best_Design_my = PSO_find_Mm(Seed=rnorm(1)*10000)


#save.image(file="PSO1.RData")


