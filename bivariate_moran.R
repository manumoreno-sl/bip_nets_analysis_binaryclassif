# Written by RogerioJB
# source https://stackoverflow.com/questions/45177590/map-of-bivariate-spatial-correlation-in-r-bivariate-lisa

# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  yp <- (y - mean(y, na.rm=T))/sd(y, na.rm=T)
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}


# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 1000){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}

# Adjacency Matrix (Queen)
moran_calc= function(A,z1,autocorrelation=TRUE,z2,alpha=0.05,degree=2){
  require(Matrix)
  # A : Igraph's network as adjacency matrix, of type sparse
  # z1 : Variable for each node
  # autocorrelation 
  # z2 : Variable to make bivariate correlation when autocorrelation=FALSE
  
  if(degree==2){
  A =A %*% A
  }
  
  W  <- as.matrix(A/Matrix::rowSums(A))
  W[which(is.na(W))] <- 0
  

  if(autocorrelation==TRUE){
  x=z1
  y=z1
  }else{
  x=z1
  y=z2
  }
  alpha=alpha/length(z1)
  m   <- moran_I(x, y, W) # m[[1]]  global value  # m[[2]]  local values
  m_i <- m[[2]]  # local values
  local_sims <- simula_moran(x, y, W)$local_sims
  intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=c(alpha/2, 1-alpha/2),na.rm=T)))
  sig        <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
  names=names(sig)
  return(list(glob=m[[1]][1,1],loc=m_i,sig=unname(sig),v_name=names))
}
# Identifying the significant values 




