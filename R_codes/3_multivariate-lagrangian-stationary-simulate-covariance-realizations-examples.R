
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'multiple-advections/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

N <- 50
n <- N^2
TT <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

ref_loc <- n / 2 + floor(N / 2)

nonfrozen = T

M3 <- F
M4 <- T

if(M3){
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
}

if(M4){
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
			cov1 <- nonfrozen_matern_cov_multi_advec_added_dimensions(theta = c(1, 1, 0.23, 0.5, 1, 0.5), wind_mu1 = c(0.1001, 0.1001, 0.0501), wind_mu2 = c(0.1001, 0.1001, -0.0501), wind_var1 = VAR_MAT_MARGIN[1:3, 1:3], wind_var2 = VAR_MAT_MARGIN[4:6, 4:6], wind_var12 = VAR_MAT_MARGIN[1:3, 4:6], max_time_lag = TT - 1, LOCS = sim_grid_locations)	
			write.table(cov1[c(ref_loc, n * TT + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
			write.table(cov1[c(ref_loc, n + ref_loc, n * 2 + ref_loc), ], file = paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
			set.seed(12345)
			cat('GENERATING REALIZATIONS ... ', '\n')
			r1 <- rmvn(10, rep(0, ncol(cov1)), cov1, ncores = 25)
			write.table(r1, file = paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example', EXAMPLE, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
		}
	}
}

