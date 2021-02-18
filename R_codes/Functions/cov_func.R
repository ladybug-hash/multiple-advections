#---------STATIONARY------------#

#-------UNIVARIATE----------#

nonfrozen_matern_uni <- function(theta, wind_mu, wind_var, max_time_lag = 0, LOCS){

	locs1 <- cbind(LOCS, rep(0, nrow(LOCS)))

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

	mu_V <- wind_mu
	sigma_V <- wind_var

	output11 <- list()

	for(tt in 0:max_time_lag){

		locs2 <- cbind(LOCS, rep(tt, nrow(LOCS)))

		Sigma_temp1 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2, beta, nu), v_mean = mu_V, v_var = sigma_V)
		output11[[tt + 1]] <- Sigma_temp1

	}

	S <- toeplitz_mat(output11)

  	return(S)
}

matern_multi_step_one <- function(theta, DIST){

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

	S <-  sigma2 * Matern(DIST, range = beta, nu = nu)
	
  	return(S)
}

#-------MULTIVARIATE----------#

purely_spatial_parsimonious_matern <- function(theta, q = 2, DIST){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	S <- matrix(NA,  q * nrow(DIST), q * nrow(DIST))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- (i - 1) * nrow(DIST) + 1:nrow(DIST)
	      		temp1 <- (j - 1) * nrow(DIST) + 1:nrow(DIST)
	      
	      		if(i == j){
		
				temp2 <- ifelse(DIST != 0, sigma2[i] * (DIST / beta)^nu[i] * besselK(DIST / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (DIST / beta)^nu3 * besselK(DIST / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}

purely_spatial_parsimonious_matern_variable_asymmetry <- function(theta, q = 2, DIST, DIST12){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	S <- matrix(NA,  q * nrow(DIST), q * nrow(DIST))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- (i - 1) * nrow(DIST) + 1:nrow(DIST)
	      		temp1 <- (j - 1) * nrow(DIST) + 1:nrow(DIST)
	      
	      		if(i == j){
		
				temp2 <- ifelse(DIST != 0, sigma2[i] * (DIST / beta)^nu[i] * besselK(DIST / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (DIST12 / beta)^nu3 * besselK(DIST12 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}


nonfrozen_matern_cov_multi_advec_added_dimensions_small_scale <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, wind_var12, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	LOCS1 <- LOCS1_FULL <- cbind(LOCS, rep(1, n), rep(1, n))
	LOCS2 <- LOCS2_FULL <- cbind(LOCS, rep(0, n), rep(1, n))

	for(tt in 2:TT){
		LOCS1_FULL <- rbind(LOCS1_FULL, cbind(LOCS, rep(1, n), rep(tt, n)))
		LOCS2_FULL <- rbind(LOCS2_FULL, cbind(LOCS, rep(0, n), rep(tt, n)))
	}

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	output11 <- output22 <- list()

	for(tt in 1:TT){

		locs2 <- cbind(LOCS, rep(1, n), rep(tt, n))

		Sigma_temp1 <- nonfrozen_rcpp_added_dimensions(Loc1 = LOCS1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = wind_mu1, v_var = wind_var1)

		locs2 <- cbind(LOCS, rep(0, n), rep(tt, n))

		Sigma_temp2 <- nonfrozen_rcpp_added_dimensions(Loc1 = LOCS1, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = wind_mu2, v_var = wind_var2)

		output11[[tt]] <- Sigma_temp1
		output22[[tt]] <- Sigma_temp2
	}

	Sigma11 <- toeplitz_mat(output11)
	Sigma22 <- toeplitz_mat(output22)

	Sigma12 <- nonfrozen_rcpp_multi_cross_added_dimensions(Loc1 = LOCS1_FULL, Loc2 = LOCS2_FULL, param = c(sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = c(wind_mu1, wind_mu2), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var2)))

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	return(S)
}

nonfrozen_matern_cov_multi_small_scale <- function(theta, wind_mu, wind_var, max_time_lag = 0, q = 2, LOCS, k = c(0, 0)){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	locs1 <- cbind(LOCS, rep(0, nrow(LOCS)))

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	mu_V <- wind_mu
	sigma_V <- wind_var

	output11 <- output22 <- output12 <- list()

	for(tt in 0:max_time_lag){

		locs2 <- cbind(LOCS, rep(tt, nrow(LOCS)))

		Sigma_temp1 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = mu_V, v_var = sigma_V)
		Sigma_temp2 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = mu_V, v_var = sigma_V)
		locs2 <- cbind(LOCS[, 1] - k[1], LOCS[, 2] - k[2], rep(tt, nrow(LOCS)))
		Sigma_temp12 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = mu_V, v_var = sigma_V)

		output11[[tt + 1]] <- Sigma_temp1
		output22[[tt + 1]] <- Sigma_temp2
		output12[[tt + 1]] <- Sigma_temp12

	}

	Sigma11 <- toeplitz_mat(output11)
	Sigma22 <- toeplitz_mat(output22)
	Sigma12 <- toeplitz_mat(output12)

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	return(S)
}

nonfrozen_matern_cov_multi_advec_small_scale <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, wind_var12, max_time_lag = 0, q = 2, LOCS){

	start_time <- theta[7]

	if(is.na(start_time))	start_time <- 0

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	locs1 <- locs <- cbind(LOCS, rep(0 + start_time, nrow(LOCS)))

	if(max_time_lag > 0){
		for(tt in 1:max_time_lag){
			locs <- rbind(locs, cbind(LOCS, rep(0 + start_time + tt, nrow(LOCS))))
		}
	}

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	output11 <- output22 <- output12 <- list()

	for(tt in 0:max_time_lag){

		locs2 <- cbind(LOCS, rep(0 + start_time + tt, nrow(LOCS)))

		Sigma_temp1 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = wind_mu1, v_var = wind_var1)
		Sigma_temp2 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = wind_mu2, v_var = wind_var2)

		output11[[tt + 1]] <- Sigma_temp1
		output22[[tt + 1]] <- Sigma_temp2

	}

	Sigma11 <- toeplitz_mat(output11)
	Sigma22 <- toeplitz_mat(output22)
	Sigma12 <- nonfrozen_rcpp_multi_cross(Loc1 = locs, Loc2 = locs, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = c(wind_mu1, wind_mu2), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var2)))

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	return(S)
}

nonfrozen_matern_cov_multi <- function(theta, wind_mu, wind_var, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	locs1 <- cbind(LOCS, rep(0, n))

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	mu_V <- wind_mu
	sigma_V <- wind_var

	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	clusterExport(cl, c("root", "LOCS", "TT", "n", "nu", "beta", "sigma2", "rho", "mu_V", "sigma_V", "locs1"), envir = environment())
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = '')))

	output <- foreach(tt=0:(TT - 1), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep(tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = mu_V, v_var = sigma_V)

		return(Sigma_temp)

	}

	Sigma11 <- toeplitz_mat(output)

	output <- foreach(tt=0:(TT - 1), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep(tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = mu_V, v_var = sigma_V)

		return(Sigma_temp)

	}

	Sigma22 <- toeplitz_mat(output)


	output <- foreach(tt=0:(TT - 1), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep(tt, n))

			
		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = mu_V, v_var = sigma_V)

		return(Sigma_temp)

	}

	Sigma12 <- toeplitz_mat(output)

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	stopCluster(cl)

	return(S)
}

nonfrozen_matern_cov_multi_advec <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, wind_var12, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	locs1 <- locs <- cbind(LOCS, rep(1, n))

	for(tt in 2:TT){
		locs <- rbind(locs, cbind(LOCS, rep( tt, n)))
	}

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	clusterExport(cl, c("root", "LOCS", "TT", "n", "nu", "beta", "sigma2", "rho", "wind_mu1", "wind_mu2", "wind_var1", "wind_var2", "wind_var12", "locs1"), envir = environment())
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = '')))

	output <- foreach(tt = 1:TT, .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep( tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = wind_mu1, v_var = wind_var1)

		return(Sigma_temp)

	}

	Sigma11 <- toeplitz_mat(output)

	output <- foreach(tt = 1:TT, .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep( tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = wind_mu2, v_var = wind_var2)

		return(Sigma_temp)

	}

	Sigma22 <- toeplitz_mat(output)

	Sigma12 <- nonfrozen_rcpp_multi_cross(Loc1 = locs, Loc2 = locs, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = c(wind_mu1, wind_mu2), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var2)))

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	stopCluster(cl)

	return(S)
}

nonfrozen_matern_cov_multi_advec_added_dimensions <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, wind_var12, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	LOCS1 <- LOCS1_FULL <- cbind(LOCS, rep(.2, n), rep(1, n))
	LOCS2 <- LOCS2_FULL <- cbind(LOCS, rep(0, n), rep(1, n))

	for(tt in 2:TT){
		LOCS1_FULL <- rbind(LOCS1_FULL, cbind(LOCS, rep(0.2, n), rep(tt, n)))
		LOCS2_FULL <- rbind(LOCS2_FULL, cbind(LOCS, rep(0, n), rep(tt, n)))
	}

	rho_temp <- theta[6]
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2

	rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	clusterExport(cl, c("root", "LOCS", "TT", "n", "nu", "beta", "sigma2", "rho", "wind_mu1", "wind_mu2", "wind_var1", "wind_var2", "wind_var12", "LOCS1", "LOCS2", "LOCS1_FULL", "LOCS2_FULL"), envir = environment())
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = '')))

	output <- foreach(tt = 1:TT, .packages = "MASS", .noexport = "nonfrozen_rcpp_added_dimensions") %dopar% {

		locs2 <- cbind(LOCS, rep(0.2, n), rep(tt, n))

		Sigma_temp <- nonfrozen_rcpp_added_dimensions(Loc1 = LOCS1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = wind_mu1, v_var = wind_var1)

		return(Sigma_temp)

	}

	Sigma11 <- toeplitz_mat(output)

	output <- foreach(tt = 1:TT, .packages = "MASS", .noexport = "nonfrozen_rcpp_added_dimensions") %dopar% {

		locs2 <- cbind(LOCS, rep(0, n), rep(tt, n))
		Sigma_temp <- nonfrozen_rcpp_added_dimensions(Loc1 = LOCS2, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = wind_mu2, v_var = wind_var2)

		return(Sigma_temp)

	}

	Sigma22 <- toeplitz_mat(output)

	Sigma12 <- nonfrozen_rcpp_multi_cross_added_dimensions(Loc1 = LOCS1_FULL, Loc2 = LOCS2_FULL, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = c(wind_mu1, wind_mu2), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var2)))

	S <- rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

	stopCluster(cl)

	return(S)
}


nonfrozen_lmc_cov_small_scale <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, coef_mat1, coef_mat2, max_time_lag, q = 2, LOCS){
  
	nu <- theta[5:6]
	beta <- theta[3:4]
	var <- theta[1:2]
	  
	n <- nrow(LOCS)
	
	locs1 <- cbind(LOCS, rep(0, n))

	output11 <- output22 <- output12 <- list()

	for(tt in 0:max_time_lag){

		locs2 <- cbind(LOCS, rep(tt, n))

		Sigma_temp1 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(var[1], beta[1], nu[1]), v_mean = wind_mu1, v_var = wind_var1)
		Sigma_temp2 <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(var[2], beta[2], nu[2]), v_mean = wind_mu2, v_var = wind_var2)
		output11[[tt + 1]] <- Sigma_temp1
		output22[[tt + 1]] <- Sigma_temp2

	}

	Sigma11 <- toeplitz_mat(output11)
	Sigma22 <- toeplitz_mat(output22)

 	S <- rbind(cbind(coef_mat1[1, 1] * Sigma11 + coef_mat2[1, 1] * Sigma22, coef_mat1[1, 2] * Sigma11 + coef_mat2[1, 2] * Sigma22), cbind(coef_mat1[2, 1] * Sigma11 + coef_mat2[2, 1] * Sigma22, coef_mat1[2, 2] * Sigma11 + coef_mat2[2, 2] * Sigma22))
 
  	return(S)
}

nonfrozen_lmc_cov <- function(theta, wind_mu1, wind_mu2, wind_var1, wind_var2, max_time_lag, q = 2, LOCS){
  
	nu <- theta[5:6]
	beta <- theta[3:4]
	var <- theta[1:2]
	  
	alpha <- matrix(c(theta[7], theta[8], theta[9], theta[10]), ncol=2, byrow=T)
  
	n <- nrow(LOCS)
	
	locs1 <- cbind(LOCS, rep(0, n))

	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	clusterExport(cl, c("root", "LOCS", "TT", "n", "nu", "beta", "var", "wind_mu1", "wind_mu2", "wind_var1", "wind_var2", "locs1"), envir = environment())
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
	clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = '')))

	output <- foreach(tt=0:(TT - 1), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep(tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(var[1], beta[1], nu[1]), v_mean = wind_mu1, v_var = wind_var1)

		return(Sigma_temp)

	}

	Sigma11 <- toeplitz_mat(output)

	output <- foreach(tt=0:(TT - 1), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {

		locs2 <- cbind(LOCS, rep(tt, n))

		Sigma_temp <- nonfrozen_rcpp(Loc1 = locs1, Loc2 = locs2, param = c(var[2], beta[2], nu[2]), v_mean = wind_mu2, v_var = wind_var2)

		return(Sigma_temp)

	}

	Sigma22 <- toeplitz_mat(output)


  	S <- list()

	S[[1]] <- Sigma11
	S[[2]] <- Sigma22
	
  	S1 <- rbind(cbind(alpha[1,1]^2 * S[[1]] + alpha[1,2]^2 * S[[2]], alpha[1,1] * alpha[2,1] * S[[1]] + alpha[1,2] * alpha[2,2] * S[[2]]), cbind(alpha[1,1] * alpha[2,1] * S[[1]] + alpha[1,2] * alpha[2,2] * S[[2]], alpha[2,1]^2 * S[[1]] + alpha[2,2]^2 * S[[2]]))
  
  	return(S1)
}

