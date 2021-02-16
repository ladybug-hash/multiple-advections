
workstation = F

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'multiple-advections/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

ANISOTROPY <- T

M1 <- F
M1_A <- F

M2 <- F
M2_A <- F
M2_A_VAR_ASYM <- F

M3 <- F
M3_A <- T

M5 <- F
M5_A <- F

M6 <- F
M6_A <- F
M6_A_VAR_ASYM <- F

FULL <- T

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

##############################################     NOW DO MLE USING THE OBTAINED VALUES AS STARTING POINT    ###################################################


n <- nn <- 550

args <- commandArgs(trailingOnly = TRUE)
set <- as.numeric(args[1])

cat("  SET: ", set, '\n')

if(FULL){
	subs <- 1:40
}else{
	set.seed(set)
	subs <- sample(1:40, 30)
}

REPLICATES <- length(subs)

index_in <- c(1:(n * TT), n * (TT + 1) + 1:(n * TT))

index_out_temp <- seq(1, nn * (TT + 1) * 2, 1) %in% index_in
index_out <- seq(1, nn * (TT + 1) * 2, 1)[!index_out_temp]

index_in_uni <- 1:(n * TT)
index_out_temp <- seq(1, nn * (TT + 1), 1) %in% index_in
index_out_uni <- seq(1, nn * (TT + 1), 1)[!index_out_temp] 

Z_rand_sample <- RESID_MAT[subs, ]

Z_out <- RESID_OUT_MAT[subs, ]

NEGLOGLIK_SPACE <- function(p){

	theta <- exp(p[1:5])

	cat(p, '\n')

	out <- 0

	for(ii in 1:2){

		Sigma <- matern_multi_step_one(theta[c(ii, 3, ii + 3)], DIST = dist0)
		cholmat <- t(chol(Sigma))
		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample0[[ii]]))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample0[[ii]])
			out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
			
		}
	}
	return(out)
}

NEGLOGLIK_SPACE_A <- function(p){

	theta <- exp(p[1:5])

	cat(p, '\n')

	ROTATED <- cbind(p[6] * cos(p[8]) * sim_grid_locations[, 1] + p[6] * sin(p[8]) * sim_grid_locations[, 2], -p[7] * sin(p[8]) * sim_grid_locations[, 1] + p[7] * cos(p[8]) * sim_grid_locations[, 2])

	dist0 <- as.matrix(dist(ROTATED, diag = TRUE, upper = TRUE))	

	out <- 0

	for(ii in 1:2){

		Sigma <- matern_multi_step_one(theta[c(ii, 3, ii + 3)], DIST = dist0)
		cholmat <- t(chol(Sigma))
		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample0[[ii]]))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample0[[ii]])
			out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
			
		}
	}
	return(out)
}

dist0 <- as.matrix(dist(sim_grid_locations, diag = TRUE, upper = TRUE))	

Z_rand_sample0 <- list()

for(ii in 1:2){

	R <- Z_rand_sample[, (ii - 1) * n * TT + 1:(n * TT)]

	Z_rand_sample0_temp <- matrix(, ncol = n, nrow = TT * REPLICATES)

	for(rr in 1:REPLICATES){
		for(tt in 1:TT){
			Z_rand_sample0_temp[ (rr - 1) * TT + tt, ] <- R[rr, (tt - 1) * n + 1:n] 
		}
	}

	Z_rand_sample0[[ii]] <- Z_rand_sample0_temp
}

if(!M5_A & !M5){
	if(ANISOTROPY){
		init <- c(-2.163, -1.567, -0.4808, 0.2985, 0.349, 0.474, 1.108, -0.7664)
		fit1 <- optim(par = init, fn = NEGLOGLIK_SPACE_A, control = list(trace = 5, maxit = 2000)) #

		p_space <- fit1$par

		ROTATED <- cbind(p_space[6] * cos(p_space[8]) * sim_grid_locations[, 1] + p_space[6] * sin(p_space[8]) * sim_grid_locations[, 2], -p_space[7] * sin(p_space[8]) * sim_grid_locations[, 1] + p_space[7] * cos(p_space[8]) * sim_grid_locations[, 2])
	}else{
		init <- c(-2.163, -1.567, -0.4808, 0.2985, 0.349)
		fit1 <- optim(par = init, fn = NEGLOGLIK_SPACE, control = list(trace = 5, maxit = 2000)) #
		p_space <- fit1$par
	}
}

if(M6_A_VAR_ASYM){

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M6  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	locs_s <- locs_t <- NULL

	for(tt in 1:TT){

		locs_s <- rbind(locs_s, ROTATED)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	NEGLOGLIK_ST <- function(p){
	
		cat(p, '\n')

		if(abs(p[1]) > 1 | p[3] < 0 | p[3] > 1 | p[4] < 0 | p[4] > 1){
			return(Inf)
		}else{

			locs_s <- NULL

			for(tt in 1:TT){
				locs_s <- rbind(locs_s, cbind(ROTATED[, 1] - p[5], ROTATED[, 2] - p[6]))
			}

			dist_s12 <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	

			theta <- c(exp(p_space[1:5]), p[1])

			dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
			DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
			DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

			DIST12 <- dist_s12 / (dist_temp1)^(p[4] / 2)

			Sigma <- purely_spatial_parsimonious_matern_variable_asymmetry(theta, DIST = DIST0, DIST12 = DIST12)  / DIST_TEMP

			cholmat <- tryCatch(t(chol(Sigma)) , error = function(a) numeric(0) )

			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				return(out)
			}
		}
	}

	start_time = Sys.time()

	init <- c(0.0142, -0.248, 0.217, 0.3674, 0, 0)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	cat(fit_st$par, '\n')

	p <- fit_st$par

	locs_s <- locs_t <- NULL

	for(tt in 1:(TT + 1)){

		locs_s <- rbind(locs_s, ROTATED)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	locs_s <- NULL

	for(tt in 1:(TT + 1)){
		locs_s <- rbind(locs_s, cbind(ROTATED[, 1] - p[5], ROTATED[, 2] - p[6]))
	}

	dist_s12 <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	

	theta <- c(exp(p_space[1:5]), p[1])

	dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
	DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
	DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)
	DIST12 <- dist_s12 / (dist_temp1)^(p[4] / 2)

	Sigma <- purely_spatial_parsimonious_matern_variable_asymmetry(theta, DIST = DIST0, DIST12 = DIST12)  / DIST_TEMP

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	end_time = Sys.time()

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_A_VAR_ASYM_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M6_A_VAR_ASYM_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_A_VAR_ASYM_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_A_VAR_ASYM_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_A_VAR_ASYM_set_', set, sep = ''), '\n')
	}

}


if(M6_A){

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M6  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	locs_s <- locs_t <- NULL

	for(tt in 1:TT){

		locs_s <- rbind(locs_s, ROTATED)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	NEGLOGLIK_ST <- function(p){
	
		cat(p, '\n')

		if(abs(p[1]) > 1 | p[3] < 0 | p[3] > 1 | p[4] < 0 | p[4] > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])

			dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
			DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
			DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

			Sigma <- purely_spatial_parsimonious_matern(theta, DIST = DIST0)  / DIST_TEMP

			cholmat <- tryCatch(t(chol(Sigma)) , error = function(a) numeric(0) )

			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				return(out)
			}
		}
	}

	start_time = Sys.time()

	init <- c(0.0142, -0.248, 0.217, 0.3674)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	cat(fit_st$par, '\n')

	p <- fit_st$par

	locs_s <- locs_t <- NULL

	for(tt in 1:(TT + 1)){

		locs_s <- rbind(locs_s, ROTATED)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	theta <- c(exp(p_space[1:5]), p[1])

	dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
	DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
	DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

	Sigma <- purely_spatial_parsimonious_matern(theta, DIST = DIST0)  / DIST_TEMP

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	end_time = Sys.time()

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_A_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M6_A_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_A_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_A_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_A_set_', set, sep = ''), '\n')
	}

}

if(M6){

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M6  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	locs_s <- locs_t <- NULL

	for(tt in 1:TT){

		locs_s <- rbind(locs_s, sim_grid_locations)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	NEGLOGLIK_ST <- function(p){
	
		cat(p, '\n')

		if(abs(p[1]) > 1 | p[3] < 0 | p[3] > 1 | p[4] < 0 | p[4] > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])

			dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
			DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
			DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

			Sigma <- purely_spatial_parsimonious_matern(theta, DIST = DIST0)  / DIST_TEMP

			cholmat <- tryCatch(t(chol(Sigma)) , error = function(a) numeric(0) )

			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				return(out)
			}
		}
	}

	start_time = Sys.time()

	init <- c(0.0142, -0.248, 0.217, 0.3674)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	cat(fit_st$par, '\n')

	p <- fit_st$par

	locs_s <- locs_t <- NULL

	for(tt in 1:(TT + 1)){

		locs_s <- rbind(locs_s, sim_grid_locations)
		locs_t <- c(locs_t, rep(tt, n))

	}

	dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
	dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

	theta <- c(exp(p_space[1:5]), p[1])

	dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * p[3]) + 1)
	DIST_TEMP <- rbind(cbind(dist_temp1, dist_temp1), cbind(dist_temp1, dist_temp1))
	DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

	Sigma <- purely_spatial_parsimonious_matern(theta, DIST = DIST0)  / DIST_TEMP

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	end_time = Sys.time()

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M6_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M6_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M6_set_', set, sep = ''), '\n')
	}
}


if(M5_A){

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M5  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_ST <- function(p){

		cat(p, '\n')

		theta <- exp(p[1:4])
		mu1 <- p[5:6]
		mu2 <- p[7:8]
		wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
		wind_var1 <- t(wind_var_chol) %*% wind_var_chol
		wind_var1 <- 9 * wind_var1 / (0.077 * 10000)

		wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
		wind_var2 <- t(wind_var_chol) %*% wind_var_chol
		wind_var2 <- 9 * wind_var2 / (0.077 * 10000)

		coef_chol <- matrix(c(p[15], p[16], 0, p[17]), ncol = 2, byrow = T)
		A1 <- t(coef_chol) %*% coef_chol
		coef_chol <- matrix(c(p[18], p[19], 0, p[20]), ncol = 2, byrow = T)
		A2 <- t(coef_chol) %*% coef_chol

		ROTATED <- cbind(p[21] * cos(p[23]) * sim_grid_locations[, 1] + p[21] * sin(p[23]) * sim_grid_locations[, 2], -p[22] * sin(p[23]) * sim_grid_locations[, 1] + p[22] * cos(p[23]) * sim_grid_locations[, 2])

		Sigma <- tryCatch(nonfrozen_lmc_cov_small_scale(theta = c(1, 1, theta), wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var1, wind_var2 = wind_var2, coef_mat1 = A1, coef_mat2 = A2, max_time_lag = TT - 1, LOCS = ROTATED) , error = function(a) numeric(0) )

		if(length(Sigma) == 0){
			return(Inf)
		}else{
			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				
				return(out)
			}
		}
	}

	start_time = Sys.time()

	#init <- c(log(0.20549312), log(0.20549312), 0.032623811, 0.321388974, -0.073211462, -0.005373963, 0.025145157, 0.015633703, 4.685808, 0.02526665, 4.766496, 4.685808, 0.02526665, 4.766496, 1, 0, 1, 1, 0, 1, 1.609055, -0.04168374, 1.53141749)
	init <- c(-0.3308062, -0.3308062, 0.32623811, 0.321388974, 0.000001, 0.000001, 0.000001, 0.000001, 4.685808, 0, 4.766496, 4.685808, 0, 4.766496, 1, 0, 1, 1, 0, 1, 0.7375330, 1.687539, -0.7756992)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	end_time = Sys.time()

	cat(fit_st$par, '\n')

	p <- fit_st$par

	theta <- exp(p[1:4])
	mu1 <- p[5:6]
	mu2 <- p[7:8]
	wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
	wind_var1 <- t(wind_var_chol) %*% wind_var_chol
	wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
	wind_var2 <- t(wind_var_chol) %*% wind_var_chol

	coef_chol <- matrix(c(p[15], p[16], 0, p[17]), ncol = 2, byrow = T)
	A1 <- t(coef_chol) %*% coef_chol
	coef_chol <- matrix(c(p[18], p[19], 0, p[20]), ncol = 2, byrow = T)
	A2 <- t(coef_chol) %*% coef_chol

	ROTATED <- cbind(p[21] * cos(p[23]) * sim_grid_locations[, 1] + p[21] * sin(p[23]) * sim_grid_locations[, 2], -p[22] * sin(p[23]) * sim_grid_locations[, 1] + p[22] * cos(p[23]) * sim_grid_locations[, 2])

	Sigma <- nonfrozen_lmc_cov_small_scale(theta = c(1, 1, theta), wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var1, wind_var2 = wind_var2, coef_mat1 = A1, coef_mat2 = A2, max_time_lag = TT, LOCS = ROTATED)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)


	if(FULL){
		write.table(c(fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M5_A_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M5_A_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M5_A_FULL', sep = ''), '\n')
	}else{
		write.table(c(fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M5_A_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M5_A_set_', set, sep = ''), '\n')
	}
}

if(M5){

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M5  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_ST <- function(p){

		cat(p, '\n')

		theta <- exp(p[1:4])
		mu1 <- p[5:6]
		mu2 <- p[7:8]
		wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
		wind_var1 <- t(wind_var_chol) %*% wind_var_chol
		wind_var1 <- 9 * wind_var1 / (0.077 * 10000)

		wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
		wind_var2 <- t(wind_var_chol) %*% wind_var_chol
		wind_var2 <- 9 * wind_var2 / (0.077 * 10000)

		coef_chol <- matrix(c(p[15], p[16], 0, p[17]), ncol = 2, byrow = T)
		A1 <- t(coef_chol) %*% coef_chol
		coef_chol <- matrix(c(p[18], p[19], 0, p[20]), ncol = 2, byrow = T)
		A2 <- t(coef_chol) %*% coef_chol

		Sigma <- tryCatch(nonfrozen_lmc_cov_small_scale(theta = c(1, 1, theta), wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var1, wind_var2 = wind_var2, coef_mat1 = A1, coef_mat2 = A2, max_time_lag = TT - 1, LOCS = sim_grid_locations) , error = function(a) numeric(0) )

		if(length(Sigma) == 0){
			return(Inf)
		}else{
			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				
				return(out)
			}
		}
	}

	start_time = Sys.time()

	init <- c(log(0.20549312), log(0.20549312), 0.032623811, 0.321388974, -0.073211462, -0.005373963, 0.025145157, 0.015633703, 4.685808, 0.02526665, 4.766496, 4.685808, 0.02526665, 4.766496, 1, 0, 1, 1, 0, 1)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	cat(fit_st$par, '\n')

	p <- fit_st$par

	theta <- exp(p[1:4])
	mu1 <- p[5:6]
	mu2 <- p[7:8]
	wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
	wind_var1 <- t(wind_var_chol) %*% wind_var_chol
	wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
	wind_var2 <- t(wind_var_chol) %*% wind_var_chol

	coef_chol <- matrix(c(p[15], p[16], 0, p[17]), ncol = 2, byrow = T)
	A1 <- t(coef_chol) %*% coef_chol
	coef_chol <- matrix(c(p[18], p[19], 0, p[20]), ncol = 2, byrow = T)
	A2 <- t(coef_chol) %*% coef_chol

	Sigma <- nonfrozen_lmc_cov_small_scale(theta = c(1, 1, theta), wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var1, wind_var2 = wind_var2, coef_mat1 = A1, coef_mat2 = A2, max_time_lag = TT, LOCS = sim_grid_locations)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	end_time = Sys.time()

	if(FULL){
		write.table(c(fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M5_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M5_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M5_FULL', sep = ''), '\n')
	}else{
		write.table(c(fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M5_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M5_set_', set, sep = ''), '\n')
	}
}

if(M1_A){
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_ST <- function(p){

		cat(p, '\n')

		out <- 0

		for(ii in 1:2){

			theta <- exp(c(p_space[ii], p_space[3], p_space[ii + 3]))
			mu <- p[(ii - 1) * 2 + 1:2]
			wind_var_chol <- matrix(c(p[(ii - 1) * 3 + 5], p[(ii - 1) * 3 + 6], 0, p[(ii - 1) * 3 + 7]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_uni(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = ROTATED)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample_temp[[ii]]))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample_temp[[ii]])
				out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
				
			}
		}
		return(out)
	}

	start_time = Sys.time()

	Z_rand_sample_temp <- Z_out_temp <- list()

	for(ii in 1:2){

		Z_rand_sample_temp[[ii]] <- Z_rand_sample[, (ii - 1) * n * TT + 1:(n * TT)]
		Z_out_temp[[ii]] <- Z_out[, (ii - 1) * n + 1:n]
	}

	init <- c(0, 0, 0, 0, 5.45967, 0, 5.45967, 5.45967, 0, 5.45967)
	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	end_time = Sys.time()

	p <- fit_st$par

	MSE_temp <- NULL

	for(ii in 1:2){
		theta <- exp(c(p_space[ii], p_space[3], p_space[ii + 3]))
		mu <- p[(ii - 1) * 2 + 1:2]
		wind_var_chol <- matrix(c(p[(ii - 1) * 3 + 5], p[(ii - 1) * 3 + 6], 0, p[(ii - 1) * 3 + 7]), ncol = 2, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol

		wind_var <- 9 * wind_var / (0.077 * 10000)

		Sigma <- nonfrozen_matern_uni(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = ROTATED)

		C11_m <- Sigma[index_in_uni, index_in_uni]
		C22_m <- Sigma[index_out_uni, index_out_uni]
		C12_m <- Sigma[index_in_uni, index_out_uni]

		pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample_temp[[ii]])
		
		MSE_temp <- c(MSE_temp, (pred - t(Z_out_temp[[ii]]))^2)

	}

	MSE <- mean(MSE_temp)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M1_A_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M1_A_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M1_A_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M1_A_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M1_A_set_', set, sep = ''), '\n')
	}
}

if(M1){
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_ST <- function(p){

		cat(p, '\n')

		out <- 0

		for(ii in 1:2){

			theta <- exp(c(p_space[ii], p_space[3], p_space[ii + 3]))
			mu <- p[(ii - 1) * 2 + 1:2]
			wind_var_chol <- matrix(c(p[(ii - 1) * 3 + 5], p[(ii - 1) * 3 + 6], 0, p[(ii - 1) * 3 + 7]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_uni(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = sim_grid_locations)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample_temp[[ii]]))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample_temp[[ii]])
				out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
				
			}
		}
		return(out)
	}

	start_time = Sys.time()

	Z_rand_sample_temp <- Z_out_temp <- list()

	for(ii in 1:2){

		Z_rand_sample_temp[[ii]] <- Z_rand_sample[, (ii - 1) * n * TT + 1:(n * TT)]
		Z_out_temp[[ii]] <- Z_out[, (ii - 1) * n + 1:n]
	}

	init <- c(0, 0, 0, 0, 7.45967, 0, 7.45967, 7.45967, 0, 7.45967)
	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	end_time = Sys.time()

	p <- fit_st$par

	MSE_temp <- NULL

	for(ii in 1:2){
		theta <- exp(c(p_space[ii], p_space[3], p_space[ii + 3]))
		mu <- p[(ii - 1) * 2 + 1:2]
		wind_var_chol <- matrix(c(p[(ii - 1) * 3 + 5], p[(ii - 1) * 3 + 6], 0, p[(ii - 1) * 3 + 7]), ncol = 2, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol

		wind_var <- 9 * wind_var / (0.077 * 10000)

		Sigma <- nonfrozen_matern_uni(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = sim_grid_locations)

		C11_m <- Sigma[index_in_uni, index_in_uni]
		C22_m <- Sigma[index_out_uni, index_out_uni]
		C12_m <- Sigma[index_in_uni, index_out_uni]

		pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample_temp[[ii]])
		
		MSE_temp <- c(MSE_temp, (pred - t(Z_out_temp[[ii]]))^2)

	}

	MSE <- mean(MSE_temp)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M1_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M1_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M1_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M1_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M1_set_', set, sep = ''), '\n')
	}
}

if(M3_A){
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M3  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_ST <- function(p){
		cat(p, '\n')

		if(abs(p[1]) > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1], 0)
			mu1 <- p[2:3]
			mu2 <- p[4:5]
			wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT - 1, LOCS = ROTATED)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				
			}

			return(out)
		}

	}

	start_time = Sys.time()

	init <- c(0.7, 0, 0, 0, 0, 5.748917, 0.6336984,  5.935040,  0.9178897, 6.3098361, -0.449795,  6.3433873, 2.210878, -0.3663625, 2.3823858)
	#init <- c(0.7872109, -0.07675721, -0.04057635, -0.05899823, -0.01698744, 6.800345, 0.6538839, 6.715755, 0.09911932, 7.370236, -0.7061986, 6.868115, 0.9774688, -0.3456212, 1.121637)
	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	end_time = Sys.time()

	p <- fit_st$par

	theta <- c(exp(p_space[1:5]), p[1], 0)
	mu1 <- p[2:3]
	mu2 <- p[4:5]
	wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- 9 * wind_var / (0.077 * 10000)

	Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT, LOCS = ROTATED)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M3_A_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M3_A_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M3_A_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M3_A_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M3_A_set_', set, sep = ''), '\n')
	}

}

if(M3){
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M3  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	if(FULL){
		p1 <- read.table(paste(root, 'Results/multivariate_stationary_real_M1_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	}else{
		p1 <- read.table(paste(root, 'Results/multivariate_stationary_real_M1_set_', set, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	}

	NEGLOGLIK_ST <- function(p){
		cat(p, '\n')

		if(abs(p[1]) > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1], 0)
			mu1 <- p1[9:10]
			mu2 <- p1[11:12]
			wind_var_chol <- matrix(c(p1[13], p1[14], p[2], p[3], 0, p1[15], p[4], p[5], 0, 0, p1[16], p1[17], 0, 0, 0, p1[18]), ncol = 4, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)
				
			}

			return(out)
		}

	}

	start_time = Sys.time()

	init <- rep(0, 5)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

	end_time = Sys.time()

	p <- fit_st$par

	theta <- c(exp(p_space[1:5]), p[1], 0)
	mu1 <- p1[9:10]
	mu2 <- p1[11:12]
	wind_var_chol <- matrix(c(p1[13], p1[14], p[2], p[3], 0, p1[15], p[4], p[5], 0, 0, p1[16], p1[17], 0, 0, 0, p1[18]), ncol = 4, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- 9 * wind_var / (0.077 * 10000)

	Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT, LOCS = sim_grid_locations)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	if(FULL){
		write.table(c(p1, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M3_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M3_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M3_FULL', sep = ''), '\n')
	}else{
		write.table(c(p1, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M3_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M3_set_', set, sep = ''), '\n')
	}

}

if(M2_A_VAR_ASYM){

	NEGLOGLIK_ST <- function(p){
		cat(p, '\n')

		if(abs(p[1]) > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])
			mu <- p[2:3]
			wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = ROTATED, k = p[7:8])

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)

				return(out)
			}
		}

	}

	start_time = Sys.time()

	init <- c(0, 0, 0, 7.45967, 0, 7.45967, 0, 0)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	end_time = Sys.time()

	cat(fit_st$par, '\n')

	p <- fit_st$par

	theta <- c(exp(p_space[1:5]), p[1])
	mu <- p[2:3]
	wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- 9 * wind_var / (0.077 * 10000)

	Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = ROTATED)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_A_VAR_ASYM_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M2_A_VAR_ASYM_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_A_VAR_ASYM_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_A_VAR_ASYM_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_A_VAR_ASYM_set_', set, sep = ''), '\n')
	}

}

if(M2_A){

	NEGLOGLIK_ST <- function(p){
		cat(p, '\n')

		if(abs(p[1]) > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])
			mu <- p[2:3]
			wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = ROTATED)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)

				return(out)
			}
		}

	}

	start_time = Sys.time()

	init <- c(0, 0, 0, 5.45967, 0, 5.45967)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 1000)) #

	end_time = Sys.time()

	cat(fit_st$par, '\n')

	p <- fit_st$par

	theta <- c(exp(p_space[1:5]), p[1])
	mu <- p[2:3]
	wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- 9 * wind_var / (0.077 * 10000)

	Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = ROTATED)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_A_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M2_A_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_A_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_A_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_A_set_', set, sep = ''), '\n')
	}

}

if(M2){

	NEGLOGLIK_ST <- function(p){
		cat(p, '\n')

		if(abs(p[1]) > 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])
			mu <- p[2:3]
			wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			wind_var <- 9 * wind_var / (0.077 * 10000)

			Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = sim_grid_locations)

			cholmat <- tryCatch(t(chol(Sigma)), error = function(a) numeric(0) )
			if( length(cholmat) == 0){
				return(Inf)
			}else{
				z <- forwardsolve(cholmat, t(Z_rand_sample))
				logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
				out  <- 1/2 * logsig + 1/2 * sum(z^2)

				return(out)
			}
		}

	}

	start_time = Sys.time()

	init <- c(0, 0, 0, 7.45967, 0, 7.45967)

	fit_st <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

	end_time = Sys.time()

	cat(fit_st$par, '\n')

	p <- fit_st$par

	theta <- c(exp(p_space[1:5]), p[1])
	mu <- p[2:3]
	wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- 9 * wind_var / (0.077 * 10000)

	Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = sim_grid_locations)

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	if(FULL){
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_FULL', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		write.table(t(pred), file = paste(root, 'Results/multivariate_stationary_real_M2_FULL_PRED', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_FULL', sep = ''), '\n')
	}else{
		write.table(c(p_space, fit_st$par, -fit_st$value, MSE, as.numeric(end_time - start_time, units = "secs")), file = paste(root, 'Results/multivariate_stationary_real_M2_set_', set, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_real_M2_set_', set, sep = ''), '\n')
	}

}


