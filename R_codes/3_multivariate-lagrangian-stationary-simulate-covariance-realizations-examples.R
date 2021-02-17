#frozen

workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
sourceCpp(file=paste(root, "R_codes/Functions/distR.cpp",sep=''))

N <- 50
n <- N^2
TT <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

ref_loc <- n / 2 + floor(N / 2)

nonfrozen = T

example = 1

for(example in 1:3){

	if(example == 1){

		Sigma_V <- list()
		Sigma_V[[1]] <- matrix(c(1, 0, 0, 1), ncol = 2, nrow = 2)
		Sigma_V[[2]] <- 0.1 * matrix(c(1, 0, 0, 1), ncol = 2, nrow = 2)
		Sigma_V[[3]] <- 0.1 * matrix(c(1, 0.9, 0.9, 1), ncol = 2, nrow = 2)

		for(EXAMPLE in 1:3){

			cat('COMPUTING COVARIANCE MATRIX ... ', '\n')

			cov1 <- nonfrozen_matern_cov_multi(c(1, 1, 0.23, 0.5, 1, 0.5), wind_mu = c(0.1001, 0.1001), wind_var = Sigma_V[[EXAMPLE]], max_time_lag = TT - 1, LOCS = sim_grid_locations)

			write.table(cov1[c(ref_loc, n * TT + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_covariance_single_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

			cat('GENERATING REALIZATIONS ... ', '\n')

			set.seed(12345)
			r1 <- rmvn(10, rep(0, ncol(cov1)), cov1, ncores = 25)

			write.table(r1, file = paste(root, 'Data/multivariate_stationary_realizations_single_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		}

	}else if(example == 2){

		if(nonfrozen){	

			for(EXAMPLE in 1:3){

				if(EXAMPLE == 1){
					A11 <- 1
					A12 <- sqrt(1 - A11^2)
					A22 <- sqrt(0.5)
					A21 <- sqrt(1 - A22^2)
					WIND_MU1 <- WIND_MU2 <- c(0.1001, 0.1001)
				}else if(EXAMPLE == 2){
					A11 <- 1
					A12 <- sqrt(1 - A11^2)
					A22 <- sqrt(0.5)
					A21 <- sqrt(1 - A22^2)
					WIND_MU1 <- c(0.1001, 0.1001)
					WIND_MU2 <- c(-0.1001, -0.1001)
				}else if(EXAMPLE == 3){
					A11 <- sqrt(0.6)
					A12 <- 
					A22 <- sqrt(0.1)
					A21 <- sqrt(1 - A22^2)
					WIND_MU1 <- c(0.1001, 0.1001)
					WIND_MU2 <- c(-0.1001, -0.1001)
				}else{
					A11 <- 1
					A12 <- A13 <- A14 <- sqrt(1 - A11^2)
					A22 <- A23 <- A24 <- sqrt(0.25)
					A21 <- sqrt(1 - A22^2 - A23^2 - A24^2)
					WIND_MU1 <- c(0.1001, 0.1001)
					WIND_MU2 <- c(-0.1001, -0.1001)
					WIND_MU3 <- c(0.1001, -0.1001)
					WIND_MU4 <- c(-0.1001, 0.1001)
				}

				cat('COMPUTING COVARIANCE MATRIX ... ', '\n')
				if(EXAMPLE <= 3){	
					cov1 <- nonfrozen_lmc_cov(theta = c(1, 1, 0.23, 0.23, 1, 1, A11, A12, A21, A22), wind_mu1 = WIND_MU1, wind_mu2 = WIND_MU2, wind_var1 = 0.001 * diag(2), wind_var2 = 0.001 * diag(2), max_time_lag = TT - 1, LOCS= sim_grid_locations)
				}else{
					cov1 <- nonfrozen_lmc_cov_R3(theta = c(1, 1, 0.23, 0.23, 1, 1, A11, A12, A21, A22), wind_mu1 = WIND_MU1, wind_mu2 = WIND_MU2, wind_mu3 = WIND_MU3, wind_mu4 = WIND_MU4, wind_var1 = 0.001 * diag(2), wind_var2 = 0.001 * diag(2), wind_var3 = 0.001 * diag(2), wind_var4 = 0.001 * diag(2), max_time_lag = TT - 1, LOCS= sim_grid_locations)
				}
				write.table(cov1[c(ref_loc, n * TT + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_covariance_lmc_multiple_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

				cat('GENERATING REALIZATIONS ... ', '\n')

				set.seed(12345)
				r1 <- rmvn(10, rep(0, ncol(cov1)), cov1, ncores = 25)

				write.table(r1, file = paste(root, 'Data/multivariate_stationary_realizations_lmc_multiple_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
			}

		}else{
			cov1 <- lmc_cov(theta = c(1, 1, 0.1, 0.23, 0.5, 1, A11, A12, A21, A22), wind1 = c(0.1001, 0.1001), wind2 = c(-0.1001, -0.1001), max_time_lag = TT - 1, LOCS = sim_grid_locations)
		}
	}else if(example == 'matern_multi_advec'){
		if(nonfrozen){	

			for(EXAMPLE in 1:3){
				cat('EXAMPLE: ', EXAMPLE, '\n')

				if(EXAMPLE == 1){
					VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, 0.9, 0, 0, 1, 0, 0.9, 0.9, 0, 1, 0, 0, 0.9, 0, 1), ncol = 4, nrow = 4, byrow = T)
				}else if(EXAMPLE == 2){
                                        VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), ncol = 4, nrow = 4, byrow = T)
                                }else{
					VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, -0.9, 0, 0, 1, 0, -0.9, -0.9, 0, 1, 0, 0, -0.9, 0, 1), ncol = 4, nrow = 4, byrow = T)
				}

				cat('COMPUTING COVARIANCE ... ', '\n')

				cov1 <- nonfrozen_matern_cov_multi_advec(theta = c(1, 1, 0.23, 0.5, 1, 0.5), wind_mu1 = c(0.1001, 0.1001), wind_mu2 = c(-0.1001, 0.1001), wind_var1 = VAR_MAT_MARGIN[1:2, 1:2], wind_var2 = VAR_MAT_MARGIN[3:4, 3:4], wind_var12 = VAR_MAT_MARGIN[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)	
				#cov1 <- nonfrozen_matern_cov_multi_advec_small_scale_NEW(theta = c(1, 1, 0.5, 0.5, 1, 0.5), wind_mu1 = c(0.1001, 0.1001), wind_mu2 = c(-0.1001, 0.1001), wind_var1 = VAR_MAT_MARGIN[1:2, 1:2], wind_var2 = VAR_MAT_MARGIN[3:4, 3:4], wind_var12 = VAR_MAT_MARGIN[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)	
				write.table(cov1[c(ref_loc, n * TT + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_multiple_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
				write.table(cov1[c(ref_loc, n + ref_loc, n * 2 + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

				cat('GENERATING REALIZATIONS ... ', '\n')

				set.seed(12345)
				r1 <- rmvn(10, rep(0, ncol(cov1)), cov1, ncores = 25)

				write.table(r1, file = paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
			}
		}else{
			cov1 <- frozen_matern_cov_multi_advec(c(1, 1, 0.23, 1, 1, 0.5), wind1 = c(0.1001, 0.1001), wind2 = c(-0.1001, -0.1001), max_time_lag = TT - 1, LOCS = sim_grid_locations)
		}
	}else if(example == 'matern_multi_added_dimensions'){
		if(nonfrozen){	

			for(EXAMPLE in 1:3){
				cat('EXAMPLE: ', EXAMPLE, '\n')
				if(EXAMPLE == 1){
					VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, 0, 0.9999, 0, 0, 0, 1, 0, 0, 0.9999, 0, 0, 0, 1, 0, 0, 0.9, 0.9999, 0, 0, 1, 0, 0, 0, 0.9999, 0, 0, 1, 0, 0, 0, 0.9, 0, 0, 1), ncol = 6, nrow = 6, byrow = T)
				}else if(EXAMPLE == 2){
					VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, 0, 0.9999, 0, 0, 0, 1, 0, 0, 0.9999, 0, 0, 0, 1, 0, 0, 0, 0.9999, 0, 0, 1, 0, 0, 0, 0.9999, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), ncol = 6, nrow = 6, byrow = T)
                                }else{
					VAR_MAT_MARGIN <- 0.1 * matrix(c(1, 0, 0, 0.99, 0, 0, 0, 1, 0, 0, 0.99, 0, 0, 0, 1, 0, 0, -0.9, 0.99, 0, 0, 1, 0, 0, 0, 0.99, 0, 0, 1, 0, 0, 0, -0.9, 0, 0, 1), ncol = 6, nrow = 6, byrow = T)
				}
				cat('COMPUTING COVARIANCE ... ', '\n')
				cov1 <- nonfrozen_matern_cov_multi_advec_added_dimensions(theta = c(1, 1, 0.23, 0.5, 1, 0.5), wind_mu1 = c(0.1001, 0.1001, 0), wind_mu2 = c(0.1001, 0.1001, 0), wind_var1 = VAR_MAT_MARGIN[1:3, 1:3], wind_var2 = VAR_MAT_MARGIN[4:6, 4:6], wind_var12 = VAR_MAT_MARGIN[1:3, 4:6], max_time_lag = TT - 1, LOCS = sim_grid_locations)	
				write.table(cov1[c(ref_loc, n * TT + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
				write.table(cov1[c(ref_loc, n + ref_loc, n * 2 + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
				set.seed(12345)
				cat('GENERATING REALIZATIONS ... ', '\n')
				r1 <- rmvn(10, rep(0, ncol(cov1)), cov1, ncores = 25)
				write.table(r1, file = paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
			}
		}
	}

}

jpeg(file = paste(root, 'thesis-defense/Figures/covariance.jpg', sep = ''), width = 1000, height = 1000)

par(mfrow = c(3, 3))

for(i in 1:2){
	for(tt in 1:3){

		image.plot(matrix(cov1[(i - 1) * n * TT + ref_loc, (i - 1) * n * TT + (tt - 1) * n + 1:n], N, N), zlim = c(0, 1))

	}
}

for(tt in 1:3){

	image.plot(matrix(cov1[ref_loc + n * TT, (tt - 1) * n + 1:n], N, N), zlim = c(0, 1))

}

dev.off()



jpeg(file = paste(root, 'thesis-defense/Figures/realizations.jpg', sep = ''), width = 1200, height = 1000)

par(mfrow = c(2, 3))

for(tt in 1:3){

	image.plot(matrix(r1[3, (tt - 1) * n + 1:n], N, N), zlim = range(r1))

}

for(tt in 1:3){

	image.plot(matrix(r1[3, n * TT + (tt - 1) * n + 1:n], N, N), zlim = range(r1))

}

dev.off()

jpeg(file = paste(root, 'thesis-defense/Figures/covariance.jpg', sep = ''), width = 1000, height = 1000)

COV <- cov1 - cov2

image.plot(COV^2)

dev.off()

