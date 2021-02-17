
workstation = F

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'multiple-advections/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

M3 = T
M2 = F

N <- 23
nn <- n <- N^2
TT <- 5
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
locs <- expand.grid(grid_x, grid_y) %>% as.matrix()


args <- commandArgs(trailingOnly = TRUE)
set <- as.numeric(args[1])

rho_val = as.numeric(args[3])

config = as.numeric(args[2])

k = as.numeric(args[4])

cat(config, '\n')

if(config == 1){
	VAR_MAT_MARGIN <- k * matrix(c(1, 0, 0.9, 0, 0, 1, 0, 0.9, 0.9, 0, 1, 0, 0, 0.9, 0, 1), ncol = 4, nrow = 4, byrow = T)
}else if(config == 2){
	VAR_MAT_MARGIN <- k * matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), ncol = 4, nrow = 4, byrow = T)
}else if(config == 3){
	VAR_MAT_MARGIN <- k * matrix(c(1, 0, -0.9, 0, 0, 1, 0, -0.9, -0.9, 0, 1, 0, 0, -0.9, 0, 1), ncol = 4, nrow = 4, byrow = T)
}

THETA_SPACE <- c(1, 1, 0.23, 0.5, 1, rho_val)
MU1 <- c(0.1, 0.1)
MU2 <- c(0.1, 0.1)

cov1 <- nonfrozen_matern_cov_multi_advec_small_scale(theta = THETA_SPACE, wind_mu1 = MU1, wind_mu2 = MU2, wind_var1 = VAR_MAT_MARGIN[1:2, 1:2], wind_var2 = VAR_MAT_MARGIN[3:4, 3:4], wind_var12 = VAR_MAT_MARGIN[1:2, 3:4], max_time_lag = TT, LOCS = locs)	

set.seed(set)
r1 <- rmvn(30, rep(0, ncol(cov1)), cov1, ncores = 25)

REPLICATES <- nrow(r1)

sim_grid_locations <- locs

n <- nrow(sim_grid_locations)

index_in <- c(1:(n * TT), n * (TT + 1) + 1:(n * TT))

index_out_temp <- seq(1, nn * (TT + 1) * 2, 1) %in% index_in
index_out <- seq(1, nn * (TT + 1) * 2, 1)[!index_out_temp]

R1 <- r1[, index_in]

Z_rand_sample <- r1[, index_in]

Z_out <- r1[, index_out]

if(M2){

	start_time = Sys.time()

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_SPACE <- function(p){

		theta <- exp(p[1:5])

		cat(theta, '\n')

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
		if(ii == 1) R <- R1[, 1:(n * TT)] 	else	R <- R1[, n * TT + 1:(n * TT)]
		Z_rand_sample0_temp <- matrix(, ncol = n, nrow = TT * REPLICATES)

		for(rr in 1:REPLICATES){
			for(tt in 1:TT){
				Z_rand_sample0_temp[ (rr - 1) * TT + tt, ] <- R[rr, (tt - 1) * n + 1:n] 
			}
		}

		Z_rand_sample0[[ii]] <- Z_rand_sample0_temp
	}

	init <- c(0, 0, log(0.05), 0, 0)
	#init <- log(THETA_SPACE[1:5])

	fit1 <- optim(par = init, fn = NEGLOGLIK_SPACE, control = list(trace = 5, maxit = 2000)) #

	p_space <- fit1$par
	
	NEGLOGLIK_ST <- function(p){

		if(p[1] <= -1 | p[1] >= 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1])
			mu <- p[2:3]
			wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			cat(p[1:3], wind_var[1, 1], wind_var[2, 2], wind_var[1, 2], '\n')

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

	init <- c(0, 0.0001, 0.0001, 0.1, 0, 0.1)
	#init <- c(rho_val, (MU1[1] + MU2[1]) / 2, (MU1[2] + MU2[2]) / 2, t(chol(VAR_MAT_MARGIN[1:2, 1:2]))[lower.tri(VAR_MAT_MARGIN[1:2, 1:2], diag = T)])

	fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	p <- fit1$par

	theta <- c(exp(p_space[1:5]), p[1])

	mu <- p[2:3]
	wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT, LOCS = locs)

	C11_m_M3 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M3 <- C12_m <- Sigma[index_in, index_out]

	M3_temp = t(C12_m_M3) %*% solve(C11_m_M3)

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	MSE <- mean((pred - t(Z_out))^2)

	pred_var <- sum(diag(t(C12_m) %*% solve(C11_m) %*% C12_m))

	end_time = Sys.time()

	total_time <- as.numeric(end_time - start_time, units = "secs")

	write.table(matrix(c("M2", set, config, THETA_SPACE, k, p_space, p, fit1$value, MSE, pred_var, total_time), nrow = 1), file = paste(root, 'Results/multivariate_stationary_simulation_study', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

	cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_simulation_study', sep = ''), '\n')

}


if(M3){
	start_time = Sys.time()

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M3  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	NEGLOGLIK_SPACE <- function(p){

		theta <- exp(p[1:5])

		cat(theta, '\n')

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
		if(ii == 1) R <- R1[, 1:(n * TT)] 	else	R <- R1[, n * TT + 1:(n * TT)]
		Z_rand_sample0_temp <- matrix(, ncol = n, nrow = TT * REPLICATES)

		for(rr in 1:REPLICATES){
			for(tt in 1:TT){
				Z_rand_sample0_temp[ (rr - 1) * TT + tt, ] <- R[rr, (tt - 1) * n + 1:n] 
			}
		}

		Z_rand_sample0[[ii]] <- Z_rand_sample0_temp
	}

	init <- c(0, 0, log(0.05), 0, 0)
	#init <- log(THETA_SPACE[1:5])

	fit1 <- optim(par = init, fn = NEGLOGLIK_SPACE, control = list(trace = 5, maxit = 2000)) #

	p_space <- fit1$par

	NEGLOGLIK_ST <- function(p){
		cat(p[1:5], '\n')

		if(p[1] <= -1 | p[1] >= 1){
			return(Inf)
		}else{

			theta <- c(exp(p_space[1:5]), p[1], 0)
			mu1 <- p[2:3]
			mu2 <- p[4:5]
			wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol

			Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)	

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

	#init <- c(rho_val, MU1, MU2, t(chol(VAR_MAT_MARGIN))[lower.tri(VAR_MAT_MARGIN, diag = T)])
	init <- c(rho_val, 0.1, 0.1, 0.1, 0.1, 1, 0, 0.9, 0, 1, 0, 0.9, 0.4358899, 0, 0.4358899)

	fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 2000)) #

	p <- fit1$par

	theta <- c(exp(p_space), p[1], 0)
	mu1 <- p[2:3]
	mu2 <- p[4:5]
	wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol

	wind_var <- VAR_MAT_MARGIN

	Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT, LOCS = locs)	

	C11_m_M2 <- C11_m <- Sigma[index_in, index_in]
	C22_m <- Sigma[index_out, index_out]
	C12_m_M2 <- C12_m <- Sigma[index_in, index_out]

	pred <- t(C12_m) %*% solve(C11_m) %*% t(Z_rand_sample)

	M2_temp = t(C12_m_M2) %*% solve(C11_m_M2)

	MSE <- mean((pred - t(Z_out))^2)

	pred_var <- sum(diag(t(C12_m) %*% solve(C11_m) %*% C12_m))

	end_time = Sys.time()

	total_time <- as.numeric(end_time - start_time, units = "secs")

	write.table(matrix(c("M3", set, config, THETA_SPACE, k, p_space, p, fit1$value, MSE, pred_var, total_time), nrow = 1), file = paste(root, 'Results/multivariate_stationary_simulation_study', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

	cat('DONE. Saved in ', paste(root, 'Results/multivariate_stationary_simulation_study', sep = ''), '\n')
}

