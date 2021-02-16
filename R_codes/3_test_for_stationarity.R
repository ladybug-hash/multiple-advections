
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))

TT <- 5

RESID_MAT <- matrix(, ncol = TT * 550 * 2, nrow = 40)
RESID_OUT_MAT <- matrix(, ncol = 550 * 2, nrow = 40)

for(yr in 1980:2019){

	dat_temp <- read.table(paste(root, 'Data/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat2_temp <- read.table(paste(root, 'Data/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	dat_temp[which(dat_temp < -25)] <- -25
	dat2_temp[which(dat2_temp < -25)] <- -25

	DAT1 <- dat_temp
	DAT2 <- dat2_temp

	locs <- read.table(paste(root, 'Data/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	load(paste(root, 'Data/covariates_' , yr, '.RData', sep = ''))

	locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
	locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / 100
	std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

	sim_grid_locations <- locs_demean

	Z1 <- DAT1
	Z2 <- DAT2

	n <- nrow(sim_grid_locations)
	nn <- nrow(locs_demean)

	start_index <- 139

	X_list <- list()

	for(VAR in 1:2){

		X1 <- X2 <- X3 <- X4 <- NULL
		for(ll in 1:TT){
			X1 <- c(X1, data_array[start_index + ll, , 1, VAR])
			X2 <- c(X2, data_array[start_index + ll, , 2, VAR])
		}
		X_list[[VAR]] <- X <- cbind(rep(1, length(X1)), (X1 - mean(X1)) / sd(X1), (X2 - mean(X2)) / sd(X2), rep(std_locs[, 1], TT), rep(std_locs[, 2], TT))

	}

	X_full <- rbind(cbind(X_list[[1]], matrix(0, ncol = ncol(X_list[[1]]), nrow = nrow(X_list[[1]]))), cbind(matrix(0, ncol = ncol(X_list[[1]]), nrow = nrow(X_list[[1]])), X_list[[2]]))


	R1 <- R2 <- NULL
	for(start_index in 139:143){

		R1 <- c(R1, Z1[start_index + 1, ])
		R2 <- c(R2, Z2[start_index + 1, ])
	}

	Z_rand_sample <- matrix( c(R1, R2), nrow = 1)

	beta_GLS <- solve(t(X_full) %*% X_full) %*% t(X_full) %*% t(Z_rand_sample)
	resid <- t(Z_rand_sample) - X_full %*% beta_GLS

	RESID_MAT[yr - 1979, ] <- resid

	X0_list <- list()

	for(VAR in 1:2){

		X1 <- X2 <- X3 <- X4 <- NULL
		for(ll in (TT + 1):(TT + 1)){
			X1 <- c(X1, data_array[start_index + ll, , 1, VAR])
			X2 <- c(X2, data_array[start_index + ll, , 2, VAR])
		}
		X0_list[[VAR]] <- X <- cbind(rep(1, length(X1)), (X1 - mean(X1)) / sd(X1), (X2 - mean(X2)) / sd(X2), std_locs)

	}

	X0_full <- rbind(cbind(X0_list[[1]], matrix(0, ncol = ncol(X0_list[[1]]), nrow = nrow(X0_list[[1]]))), cbind(matrix(0, ncol = ncol(X0_list[[1]]), nrow = nrow(X0_list[[1]])), X0_list[[2]]))


	R1 <- R2 <- NULL
	for(start_index in 144:144){

		R1 <- c(R1, Z1[start_index + 1, ])
		R2 <- c(R2, Z2[start_index + 1, ])
	}

	Z_out <- matrix( c(R1, R2), nrow = 1)
	resid_out <- t(Z_out) - X0_full %*% beta_GLS
	RESID_OUT_MAT[yr - 1979, ] <- resid_out
}

n <- nn <- nrow(locs)

compute_splag <- function(LOCS){

        N <- nrow(LOCS)
        xlags <- ylags <- matrix(, ncol = N, nrow = N)

        for(i in 1:N){
                xlags[i, ] <- LOCS[i, 1] - LOCS[, 1]
                ylags[i, ] <- LOCS[i, 2] - LOCS[, 2]
        }

        splag <- list()

        for(ll in 1:nrow(lag_targ)){
                locs_sub <- NULL

                for(i in 1:N){
                        locs_temp <- which(xlags[i, ] == lag_targ[ll, 1] & ylags[i, ] == lag_targ[ll, 2])
                        if(length(locs_temp) != 0) locs_sub <- rbind(locs_sub, cbind(rep(i, length(locs_temp)), locs_temp))
                }
                splag[[ll]] <- locs_sub
        }

        return(splag)

}

compute_splag_binning <- function(LOCS){

        N <- nrow(LOCS)
        xlags <- ylags <- matrix(, ncol = N, nrow = N)

        for(i in 1:N){
                xlags[i, ] <- LOCS[i, 1] - LOCS[, 1]
                ylags[i, ] <- LOCS[i, 2] - LOCS[, 2]
        }

        splag <- list()

        for(ll in 1:nrow(lag_targ)){
                locs_sub <- NULL

                for(i in 1:N){
                        locs_temp <- which(round(xlags[i, ], 0) == lag_targ[ll, 1] & round(ylags[i, ], 0) == lag_targ[ll, 2])
                        if(length(locs_temp) != 0) locs_sub <- rbind(locs_sub, cbind(rep(i, length(locs_temp)), locs_temp))
                }
                splag[[ll]] <- locs_sub
        }

        return(splag)

}


get_A<-function(Z,lag,tlag,index)
{ 
  nTime = dim(Z)[1]
  A1 = A2 = 0
  Sn = dim(lag)[1]
  for(s in 1:Sn)
  { 
    for(t in index)
    { 
      A1 = A1 + Z[t,lag[s,1]]*Z[t+tlag,lag[s,2]]
      A2 = A2 + Z[t+nTime/2,lag[s,1]]*Z[t+tlag+nTime/2,lag[s,2]]
    }
  }
  A1 = A1/Sn/length(index)
  A2 = A2/Sn/length(index)
  return(list(A1,A2))
}


get_G<-function(Z,splag,ulag,index)
{
  m = length(splag)
  G = integer(2*m)

  for(slag in 1:length(splag))
  {
    tmp = get_A(Z,splag[[slag]],ulag[slag],index);
    G[slag] = tmp[[1]]
    G[slag + m] = tmp[[2]]
  }

  return(G)
}


get_A_space<-function(Z1, Z2,lag1, lag2, index)
{
  nTime = dim(Z1)[1]
  A1 = A2 = 0
  Sn1 = dim(lag1)[1]
  Sn2 = dim(lag2)[1]
  for(s in 1:Sn1)
  {
    for(t in index)
    {
      A1 = A1 + Z1[t, lag1[s,1]] * Z1[t, lag1[s,2]]
    }
  }
  for(s in 1:Sn2)
  {
    for(t in index)
    {
      A2 = A2 + Z2[t, lag2[s,1]] * Z2[t, lag2[s,2]]
    }
  }
  A1 = A1/Sn1/length(index)
  A2 = A2/Sn2/length(index)
  return(list(A1,A2))
}

get_G_space<-function(Z1, Z2,splag1, splag2, index)
{
  m = length(splag1)
  G = integer(2*m)

  for(slag in 1:length(splag1))
  {
    tmp = get_A_space(Z1, Z2,splag1[[slag]],splag2[[slag]],index);
    G[slag] = tmp[[1]]
    G[slag + m] = tmp[[2]]
  }

  return(G)
}

get_A_space_cross<-function(Z1_1, Z1_2, Z2_1, Z2_2, lag1, lag2, index)
{
  nTime = dim(Z1_1)[1]
  A1 = A2 = 0
  Sn1 = dim(lag1)[1]
  Sn2 = dim(lag2)[1]
  for(s in 1:Sn1)
  {
    for(t in index)
    {
      A1 = A1 + Z1_1[t, lag1[s,1]] * Z2_1[t, lag1[s,2]]
    }
  }
  for(s in 1:Sn2)
  {
    for(t in index)
    {
      A2 = A2 + Z1_2[t, lag2[s,1]] * Z2_2[t, lag2[s,2]]
    }
  }
  A1 = A1/Sn1/length(index)
  A2 = A2/Sn2/length(index)
  return(list(A1,A2))
}

get_G_space_cross <- function(Z1_1, Z1_2, Z2_1, Z2_2, splag1, splag2, index)
{
  m = length(splag1)
  G = integer(2*m)

  for(slag in 1:length(splag1))
  {
    tmp = get_A_space_cross(Z1_1, Z1_2, Z2_1, Z2_2, splag1[[slag]],splag2[[slag]],index);
    G[slag] = tmp[[1]]
    G[slag + m] = tmp[[2]]
  }

  return(G)
}


stationary.test_space<-function(Z1, Z2, splag1, splag2)
{
  m = length(splag1)

  nTime = dim(Z1)[1]
  index = 1:nTime

  G = get_G_space(Z1, Z2, splag1, splag2, index) * sqrt(length(index))

  G_ = NULL

  sample_size <- 5

  #for(i in 1:(floor((nTime/2-max(ulag))/10)))
  for(i in 1:(floor(nTime/sample_size)))
  {
    index = ((i-1)*sample_size+1):(i*sample_size)

    G_ = cbind(G_,get_G_space(Z1, Z2, splag1, splag2,index)*sqrt(length(index)))

    #index = index + max(ulag)

    #G_ = cbind(G_,get_G_space(Z1, Z2, splag1, splag2,index)*sqrt(length(index)))
  }

  Sigma = cov(t(G_))

  X=cbind(diag(m),-diag(m))

  ts = t(X%*%G)%*%solve(X%*%Sigma%*%t(X),X%*%G)
  return(1-pchisq(ts,df=m))

}

stationary.test_space_cross<-function(Z1_1, Z1_2, Z2_1, Z2_2, splag1, splag2)
{
  m = length(splag1)

  nTime = dim(Z1_1)[1]
  index = 1:nTime

  G = get_G_space_cross(Z1_1, Z1_2, Z2_1, Z2_2, splag1, splag2, index) * sqrt(length(index))

  G_ = NULL

  sample_size <- 5

  #for(i in 1:(floor((nTime/2-max(ulag))/10)))
  for(i in 1:(floor(nTime/sample_size)))
  {
    index = ((i-1)*sample_size+1):(i*sample_size)

    G_ = cbind(G_,get_G_space_cross(Z1_1, Z1_2, Z2_1, Z2_2, splag1, splag2,index)*sqrt(length(index)))

    #index = index + max(ulag)

    #G_ = cbind(G_,get_G_space(Z1, Z2, splag1, splag2,index)*sqrt(length(index)))
  }

  Sigma = cov(t(G_))

  X=cbind(diag(m),-diag(m))

  ts = t(X%*%G)%*%solve(X%*%Sigma%*%t(X),X%*%G)
  return(1-pchisq(ts,df=m))

}

stationary.test<-function(Z,splag,ulag)
{
  m = length(splag)

  nTime = dim(Z)[1]
  index = 1:(nTime/2 - max(ulag))

  G = get_G(Z,splag,ulag,index) * sqrt(length(index))

  G_ = NULL
  for(i in 1:(floor((nTime/2-max(ulag))/10)))
  {
    index = ((i-1)*10+1):(i*10)

    G_ = cbind(G_,get_G(Z,splag,ulag,index)*sqrt(length(index)))

    index = index + max(ulag)

    G_ = cbind(G_,get_G(Z,splag,ulag,index)*sqrt(length(index)))
  }

  Sigma = cov(t(G_))

  X=cbind(diag(m),-diag(m))

  ts = t(X%*%G)%*%solve(X%*%Sigma%*%t(X),X%*%G)
  return(1-pchisq(ts,df=m))

}

#################################################################################

lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1)), ncol = 2, byrow = T)
lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1), c(0, 2), c(2, 0), c(2, 2)), ncol = 2, byrow = T)
lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1), c(0, 2), c(2, 0), c(2, 2), c(0, 3), c(3, 0), c(3, 3)), ncol = 2, byrow = T)

#ind1 <- which(locs[, 1] <= 44)

ind1 <- which(locs[, 2] <= 24.5)

locs1 <- locs_demean[ind1, ]
locs2 <- locs_demean[-ind1, ]

Z <- Z2 <- NULL

for(tt in 1:40){
	Z_temp <- matrix(RESID_MAT[tt, 1:(n * TT)], ncol = n , byrow = T)
	Z2_temp <- matrix(RESID_MAT[tt, n * TT + 1:(n * TT)], ncol = n , byrow = T)
	Z <- rbind(Z, Z_temp)	
	Z2 <- rbind(Z2, Z2_temp)	
}

ZA <- Z[, ind1]
ZB <- Z[, -ind1]
Z2A <- Z2[, ind1]
Z2B <- Z2[, -ind1]
for(zz in 1:nrow(Z)){
        ZA[zz, ] <- ZA[zz, ] - mean(ZA[zz, ])
        ZB[zz, ] <- ZB[zz, ] - mean(ZB[zz, ])
        Z2A[zz, ] <- Z2A[zz, ] - mean(Z2A[zz, ])
        Z2B[zz, ] <- Z2B[zz, ] - mean(Z2B[zz, ])
}

splag1 <- compute_splag_binning(locs1)
splag2 <- compute_splag_binning(locs2)

stationary.test_space(ZA, ZB, splag1, splag2)
stationary.test_space(Z2A, Z2B, splag1, splag2)
stationary.test_space_cross(ZA, ZB, Z2A, Z2B, splag1, splag2)

############################################################ ROTATED ############################################################

#unique(sort(abs(round(xlags, 0))))

p_space_old <- read.table(paste(root, 'Results/multivariate_stationary_real_M3_A_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
R <-  p_space_old[6:8]

ROTATED <- cbind(R[1] * cos(R[3]) * locs_demean[, 1] + R[1] * sin(R[3]) * locs_demean[, 2], -R[2] * sin(R[3]) * locs_demean[, 1] + R[2] * cos(R[3]) * locs_demean[, 2])

lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1)), ncol = 2, byrow = T)
lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1), c(0, 2), c(2, 0), c(2, 2)), ncol = 2, byrow = T)
lag_targ <- 2 * matrix(c(c(0, 1), c(1, 0), c(1, 1), c(0, 2), c(2, 0), c(2, 2), c(0, 3), c(3, 0), c(3, 3)), ncol = 2, byrow = T)

#ind1 <- which(locs[, 1] <= 44)

ind1 <- which(locs[, 2] <= 24.5)

locs1 <- ROTATED[ind1, ]
locs2 <- ROTATED[-ind1, ]

Z <- Z2 <- NULL

for(tt in 1:40){
	Z_temp <- matrix(RESID_MAT[tt, 1:(n * TT)], ncol = n , byrow = T)
	Z2_temp <- matrix(RESID_MAT[tt, n * TT + 1:(n * TT)], ncol = n , byrow = T)
	Z <- rbind(Z, Z_temp)	
	Z2 <- rbind(Z2, Z2_temp)	
}

ZA <- Z[, ind1]
ZB <- Z[, -ind1]
Z2A <- Z2[, ind1]
Z2B <- Z2[, -ind1]
for(zz in 1:nrow(Z)){
        ZA[zz, ] <- ZA[zz, ] - mean(ZA[zz, ])
        ZB[zz, ] <- ZB[zz, ] - mean(ZB[zz, ])
        Z2A[zz, ] <- Z2A[zz, ] - mean(Z2A[zz, ])
        Z2B[zz, ] <- Z2B[zz, ] - mean(Z2B[zz, ])
}

splag1 <- compute_splag_binning(locs1)
splag2 <- compute_splag_binning(locs2)

stationary.test_space(ZA, ZB, splag1, splag2)
stationary.test_space(Z2A, Z2B, splag1, splag2)
stationary.test_space_cross(ZA, ZB, Z2A, Z2B, splag1, splag2)


