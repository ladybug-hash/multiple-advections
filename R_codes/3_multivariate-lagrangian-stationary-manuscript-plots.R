#plots

workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

load_data = F
fig1 = F
fig1_new = F
fig2_new = F
fig1_newnew = T
fig2_newnew = T
fig2 = F
plot_realizations = F

plot_multi_advec = F
plot_nonfrozen_multi_advec = F

fig3 = F
fig3_new = F
fig4 = F
fig5 = F

beamer = F

if(load_data){

	cat('READING DATA ... ', '\n')

	if(!beamer){
	}else{

	}
}

N <- 50
n <- N^2
TT <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
dat3 <- sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

cat('PLOTTING ... ', '\n')

if(fig1){

	DAT1 <- read.table(paste(root, 'Data/multivariate_stationary_covariance_single_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT2 <- read.table(paste(root, 'Data/multivariate_stationary_covariance_single_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT3 <- read.table(paste(root, 'Data/multivariate_stationary_covariance_single_example3', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	r1 <- read.table(paste(root, 'Data/multivariate_stationary_realizations_single_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r2 <- read.table(paste(root, 'Data/multivariate_stationary_realizations_single_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r3 <- read.table(paste(root, 'Data/multivariate_stationary_realizations_single_example3', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	
	DAT <- list()
	DAT[[1]] <- DAT1
	DAT[[2]] <- DAT2
	DAT[[3]] <- DAT3

	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), length(DAT)))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2
	REALIZATIONS_LIST[,, 3] <- r3

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])
	
	mod_labels <- c('(a)', '(b)', '(c)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	 #jpeg(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_covariance.jpg', sep = ''), width = 700, height = 600, res = 5, units = 'in')
	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_covariance_single.pdf', sep = ''), width = 9, height = 8)

	split.screen( rbind(c(0.13, 0.38, 0.67, 0.97), c(0.53, 0.93, 0.67, 0.97),
			      c(0.13, 0.38, 0.36, 0.66), c(0.53, 0.93, 0.36, 0.66),
			      c(0.13, 0.38, 0.05, 0.35), c(0.53, 0.93, 0.05, 0.35), c(0.38,0.45,0.05,0.95), c(0.93,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 2, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 2, 3 ), screen = 4 )
	split.screen( figs = c( 3, 3 ), screen = 5 )
	split.screen( figs = c( 2, 3 ), screen = 6 )
	 	
 	ll <- 1 
	
	for(cov_mod in 1:length(DAT)){

		for(tt in 1:TT){

			screen((cov_mod - 1) * 15 + 8 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		matrix_vals <- matrix(DAT[[cov_mod]][ll, (tt - 1) * n + 1:n], N, N)
	      
	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[11]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
	      		}

			if(cov_mod == 1){
				mtext(paste('u = ', tt - 1, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
			}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 11 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		matrix_vals <- matrix(DAT[[cov_mod]][ll + 1, n * TT + (tt - 1) * n + 1:n], N, N)
	      
	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[22]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
				mtext(mod_labels[cov_mod], side = 2, line = 4.5, adj = 0.5, cex = 2, font = 2, col = 4)
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
	      		}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 14 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		matrix_vals <- matrix(DAT[[cov_mod]][ll , n * TT + (tt - 1) * n + 1:n], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[12]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
		
	      		}

			if(cov_mod == length(DAT)){

	      			mtext(expression(h[x]), side = 1, line = 1.2, adj = 0.5, cex = 0.75, font=2)
	      			axis(1, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			}
	    	}
	    
	    	for(tt in 1:TT){
	     		screen((cov_mod - 1) * 15 + 17 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
	      
	      		if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(Z[1]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
	      		}
			if(cov_mod == 1){
                        	mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
                	}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 20 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(1.5,0.5,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, n * TT + (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
	      
	      		if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(Z[2]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
	      		}
			if(cov_mod == length(DAT)){

	     			mtext(expression(s[x]), side = 1, line = 1.2, adj = 0.5, cex = 0.75, font=2)
	      			axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
			}
	    	}
	}
	  
	screen(7)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.78,0.78,0.99,0.99)
	
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1, length.out = 3))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.41,0.41,0.62,0.62)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.04,0.04,0.25,0.25)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	screen(8)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.78,0.78,0.99,0.99)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.41,0.41,0.62,0.62)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.04,0.04,0.25,0.25)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	close.screen( all=TRUE)
	dev.off()
	  
}

if(plot_nonfrozen_multi_advec){
		DAT1 <- read.table(paste(root, 'Data/multivariate_stationary_multi_advec_covariance_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
		DAT2 <- read.table(paste(root, 'Data/multivariate_stationary_multi_advec_covariance_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
		r1 <- read.table(paste(root, 'Data/multivariate_stationary_multi_advec_realizations_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
		r2 <- read.table(paste(root, 'Data/multivariate_stationary_multi_advec_realizations_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	DAT <- list()
	DAT[[1]] <- DAT1
	DAT[[2]] <- DAT2

	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), length(DAT)))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	  
	  new_grid_x <- matrix(dat3[, 1], N, N)
	  new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	 #jpeg(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_covariance.jpg', sep = ''), width = 700, height = 600, res = 5, units = 'in')
	 pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_multi_advec_covariance.pdf', sep = ''), width = 9, height = 8)
	  split.screen( rbind(c(0.13, 0.52, 0.67, 0.97), c(0.54, 0.93, 0.67, 0.97),
			      c(0.13, 0.52, 0.36, 0.66), c(0.54, 0.93, 0.36, 0.66),
			      c(0.13, 0.52, 0.05, 0.35), c(0.54, 0.93, 0.05, 0.35), c(0.38,0.45,0.05,0.95), c(0.93,0.99,0.05,0.95)))
	  
	  split.screen( figs = c( 3, 3 ), screen = 1 )
	  split.screen( figs = c( 3, 3 ), screen = 2 )
	  split.screen( figs = c( 3, 3 ), screen = 3 )
	  split.screen( figs = c( 3, 3 ), screen = 4 )
	  split.screen( figs = c( 3, 3 ), screen = 5 )
	  split.screen( figs = c( 3, 3 ), screen = 6 )

	for(cov_mod in 1:length(DAT)){
	cov_part <- 1	 	
		for(ll in 1:3){
		    for(tt in 1:TT){

		      screen((cov_part - 1) * 18 + (ll - 1) * 3 + 8 + tt + (cov_mod - 1) * 9 )
		      par(pty="s") 
		      par(mai=c(0.02,0.02,0.02,0.02))
		     par(mgp=c(2,0.5,0))

		      matrix_vals <- matrix(DAT[[cov_mod]][ll, (tt - 1) * n + 1:n], N, N)
		      
		      poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
		      
		      if(tt == 1){
			mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
			mtext(expression(C[11]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
			
			axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			
		      }
			if(ll == 1){
				mtext(paste('u = ', tt - 1, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
			}
		    }
		}
	   cov_part <- 2 
		for(ll in 1:3){
		    for(tt in 1:TT){

		      screen((cov_part - 1) * 18 + (ll - 1) * 3 + 8 + tt + (cov_mod - 1) * 9 )
		      par(pty="s") 
		      par(mai=c(0.02,0.02,0.02,0.02))
		     par(mgp=c(2,0.5,0))

		      matrix_vals <- matrix(DAT[[cov_mod]][ll + 3, n * TT + (tt - 1) * n + 1:n], N, N)
		      
		      poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
		      
		      if(tt == 1){
			mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
			mtext(expression(C[11]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
			
			axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			
		      }
		    }
		}

		cov_part <- 3

		for(ll in 1:3){
		    for(tt in 1:TT){

		      screen((cov_part - 1) * 18 + (ll - 1) * 3 + 8 + tt + (cov_mod - 1) * 9 )
		      par(pty="s") 
		      par(mai=c(0.02,0.02,0.02,0.02))
		     par(mgp=c(2,0.5,0))

		      matrix_vals <- matrix(DAT[[cov_mod]][ll, n * TT + (tt - 1) * n + 1:n], N, N)
		      
		      poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
		      
		      if(tt == 1){
			mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
			mtext(expression(C[11]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
			
			axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			
		      }
		    }
		}
	   } 

	  screen(7)
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.78,0.78,0.99,0.99)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1, length.out = 3))
	  
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.41,0.41,0.62,0.62)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.04,0.04,0.25,0.25)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	  screen(8)
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.78,0.78,0.99,0.99)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.41,0.41,0.62,0.62)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	  x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	  y1 <- c(0.04,0.04,0.25,0.25)
	  
	  legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	  close.screen( all=TRUE)
	  dev.off()
	  
}

if(plot_multi_advec){

	DAT3_cross <- read.table(paste(root, 'Data/multivariate_stationary_crosscovariance_example3', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_crosscovariance.pdf', sep = ''), width = 9, height = 8)

	par(mfrow = c(3, 3))

	for(tt in 1:TT){
		image.plot(matrix(DAT3_cross[1, n * TT + n * (tt - 1) + 1:n], N, N))
	}
	for(tt in 1:TT){
		image.plot(matrix(DAT3_cross[2, n * TT + n * (tt - 1) + 1:n], N, N))
	}
	for(tt in 1:TT){
		image.plot(matrix(DAT3_cross[3, n * TT + n * (tt - 1) + 1:n], N, N))
	}

	dev.off()
}

if(fig2){

	DAT1 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_multiple_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT2 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_multiple_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT1_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT2_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT3 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_added_dim_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT4 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_covariance_matern_added_dim_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT3_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT4_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	r1 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r2 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r3 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r4 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
		
	DAT <- list()
	DAT[[1]] <- DAT1
	DAT[[2]] <- DAT2
	DAT[[3]] <- DAT3
	DAT[[4]] <- DAT4
	
	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), 4))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2
	REALIZATIONS_LIST[,, 3] <- r3
	REALIZATIONS_LIST[,, 4] <- r4

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_multiple.pdf', sep = ''), width = 9, height = 11)
	split.screen( rbind(c(0.13, 0.38, 0.74, 0.96), c(0.53, 0.93, 0.74, 0.96),
			      c(0.13, 0.38, 0.51, 0.73), c(0.53, 0.93, 0.51, 0.73),
			      c(0.13, 0.38, 0.28, 0.5), c(0.53, 0.93, 0.28, 0.5),
			      c(0.13, 0.38, 0.05, 0.27), c(0.53, 0.93, 0.05, 0.27), c(0.38,0.45,0.05,0.95), c(0.93,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 2, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 2, 3 ), screen = 4 )
	split.screen( figs = c( 3, 3 ), screen = 5 )
	split.screen( figs = c( 2, 3 ), screen = 6 )
	split.screen( figs = c( 3, 3 ), screen = 7 )
	split.screen( figs = c( 2, 3 ), screen = 8 )
	 	
 	ll <- 1 

	for(cov_mod in 1:4){

		for(tt in 1:TT){

	      		screen((cov_mod - 1) * 15 + 10 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		matrix_vals <- matrix(DAT[[cov_mod]][ll, (tt - 1) * n + 1:n], N, N)
	      
	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[11]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
		
	      		}
			if(cov_mod == 1){
                        	mtext(paste('u = ', tt - 1, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
                	}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 13 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		matrix_vals <- matrix(DAT[[cov_mod]][ll + 1, n * TT + (tt - 1) * n + 1:n], N, N)
	      
	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[22]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
				mtext(mod_labels[cov_mod], side = 2, line = 4.5, adj = 0.5, cex = 2, font = 2, col = 4)
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
		
	      		}
	      
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 16 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      		#matrix_vals <- matrix(DAT[[cov_mod]][ll , n * TT + (tt - 1) * n + 1:n], N, N)
			if(cov_mod == 1){
	      			matrix_vals <- matrix(DAT1_cross[tt, n * TT + (tt - 1) * n + 1:n], N, N)
			}else if(cov_mod == 2){
	      			matrix_vals <- matrix(DAT2_cross[tt, n * TT + (tt - 1) * n + 1:n], N, N)
			}else if(cov_mod == 3){
	      			matrix_vals <- matrix(DAT3_cross[tt, n * TT + (tt - 1) * n + 1:n], N, N)
			}else{
	      			matrix_vals <- matrix(DAT4_cross[tt, n * TT + (tt - 1) * n + 1:n], N, N)
			}
	      
	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0,1))
			abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
	      
	      		if(tt == 1){
				mtext(expression(h[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(C[12]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
	      		}

			if(cov_mod == length(DAT)){

	      			mtext(expression(h[x]), side = 1, line = 1.5, adj = 0.5, cex = 0.75, font=2)
	      
	      			axis(1, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 19 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
	      
	      		if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(Z[1]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
		
				axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
	      		}
			if(cov_mod == 1){
                        	mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
                	}
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 15 + 22 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.3,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, n * TT + (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
	      
	      		if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				mtext(expression(Z[2]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
				axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
	      		}
			if(cov_mod == length(DAT)){

	     			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 0.75, font=2)
	      
	      			axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
			}
	    	}
	}
		  
	screen(9)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.85,0.85,0.99,0.99)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1)
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.57,0.57,0.71,0.71)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.3,0.3,0.44,0.44)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.02,0.02,0.16,0.16)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))
	  
	screen(10)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.85,0.85,0.99,0.99)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.57,0.57,0.71,0.71)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.3,0.3,0.44,0.44)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	x1 <- c(0.025,0.09,0.09,0.025) + 0.05
	y1 <- c(0.02,0.02,0.16,0.16)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))
	  
	close.screen( all=TRUE)
	dev.off()
	  
}


if(fig1_new){

	r1 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r2 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r3 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r4 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
		
	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), 4))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2
	REALIZATIONS_LIST[,, 3] <- r3
	REALIZATIONS_LIST[,, 4] <- r4

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_multiple_NEW.pdf', sep = ''), width = 10, height = 7)

	split.screen( rbind(c(0.08, 0.49, 0.72, 0.92), c(0.52, 0.93, 0.72, 0.92),
			      c(0.08, 0.49, 0.51, 0.71), c(0.52, 0.93, 0.51, 0.71),
			      c(0.08, 0.49, 0.28, 0.48), c(0.52, 0.93, 0.28, 0.48),
			      c(0.08, 0.49, 0.07, 0.27), c(0.52, 0.93, 0.07, 0.27), c(0.93,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 1, 3 ), screen = 1 )
	split.screen( figs = c( 1, 3 ), screen = 2 )
	split.screen( figs = c( 1, 3 ), screen = 3 )
	split.screen( figs = c( 1, 3 ), screen = 4 )
	split.screen( figs = c( 1, 3 ), screen = 5 )
	split.screen( figs = c( 1, 3 ), screen = 6 )
	split.screen( figs = c( 1, 3 ), screen = 7 )
	split.screen( figs = c( 1, 3 ), screen = 8 )
	 	
	for(cov_mod in 1:4){

	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 6 + 9 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.5,0))

	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
	      
	      		if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1, adj = 0.5, cex = 0.75, font=2)
				text(-0.6, 0.475, mod_labels[cov_mod], col = 'blue', xpd = NA, cex = 1.5)
				axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
	      		}
	      		if(tt == 2 & cov_mod == 1){
				mtext(expression(Z[1]), side = 3, line = 1.5, adj = 0.5, cex = 1.5, font=2, col = 'blue')
	      		}
			if(cov_mod == 1){
                        	mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
                	}
			if(cov_mod == dim(REALIZATIONS_LIST)[3]){

	     			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 0.75, font=2)
	      
	      			axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
			}
			
	    	}
	    
	    	for(tt in 1:TT){
	      		screen((cov_mod - 1) * 6 + 12 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.3,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, n * TT + (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
			if(cov_mod == 1){
                        	mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1)
                	}
	      		if(tt == 2 & cov_mod == 1){
				mtext(expression(Z[2]), side = 3, line = 1.5, adj = 0.5, cex = 1.5, font=2, col = 'blue')
	      		}
	      
			if(cov_mod == dim(REALIZATIONS_LIST)[3]){

	     			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 0.75, font=2)
	      
	      			axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
			}
	    	}
	}

	screen(9)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.07
	y1 <- c(0.35,0.35,0.65,0.65)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5), CEX = 1)
		  
	close.screen( all=TRUE)
	dev.off()
	  
}

if(fig2_new){

	DAT1_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT2_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT3_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT4_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	DAT <- list()
	DAT[[1]] <- DAT1_cross
	DAT[[2]] <- DAT2_cross
	DAT[[3]] <- DAT3_cross
	DAT[[4]] <- DAT4_cross

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_covariance_multiple_NEW.pdf', sep = ''), width = 20, height = 5.5)

	split.screen( rbind(c(0.04, 0.265, 0.1, 0.88), c(0.275, 0.485, 0.1, 0.88), c(0.515, 0.725, 0.1, 0.88), c(0.735, 0.945, 0.1, 0.88),
			      c(0.96,0.99,0.05,0.93)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 3, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 3, 3 ), screen = 4 )
	 	
	for(cov_mod in 1:4){

		for(ll in 1:3){
			for(tt in 1:TT){
				screen((cov_mod - 1) * 9 + (ll - 1) * 3  + 5 + tt)
				par(pty="s") 
				par(mai=c(0.02,0.02,0.02,0.02))
				par(mgp=c(2,0.5,0))

				matrix_vals <- matrix(DAT[[cov_mod]][ll, n * TT + (tt - 1) * n + 1:n], N, N)

				poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0, 0.5))
				abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
		      
				if(tt == 1 & cov_mod == 1){
					mtext(expression(h[y]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
					axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 1, tcl=-0.3)
				}
				if(ll == 3){
					mtext(expression(h[x]), side = 1, line = 2, adj = 0.5, cex = 1.5, font=2)
					axis(1, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 1, tcl=-0.3)
				}
				if(tt == 2 & ll == 1){
					mtext(mod_labels[cov_mod], side = 3, line = 2, adj = 0.5, cex = 2, col = 'blue')
				}
				if(ll == 1){
					label1 <- bquote(paste(t[1], "=", .(tt), sep=""))
					mtext(label1, side = 3, line = 0, cex = 1.5)
				}
				if(tt == 3  & cov_mod == 4){
					label1 <- bquote(paste(t[2], "=", .(ll), sep=""))
					text(1.2, 0.475, label1, xpd = NA, cex = 1.5, srt = 270)
				}
			}
		}
	}

	screen(5)
	x1 <- c(0.025,0.09,0.09,0.025) + 0.15
	y1 <- c(0.3,0.3,0.7,0.7)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1.5)
		  
	close.screen( all=TRUE)
	dev.off()
	  
}

if(fig1_newnew){

	r1 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r2 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r3 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_multiple_example3', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r4 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r5 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	r6 <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_realizations_matern_added_dim_example3', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
		
	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), 6))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2
	REALIZATIONS_LIST[,, 3] <- r3
	REALIZATIONS_LIST[,, 4] <- r4
	REALIZATIONS_LIST[,, 5] <- r5
	REALIZATIONS_LIST[,, 6] <- r6

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

	mod_labels <- c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_example_multiple.pdf', sep = ''), width = 14, height = 9)

	split.screen( rbind(c(0.08, 0.47, 0.73, 0.93), c(0.54, 0.93, 0.73, 0.93),
			      c(0.08, 0.47, 0.48, 0.68), c(0.54, 0.93, 0.48, 0.68),
			      c(0.08, 0.47, 0.27, 0.47), c(0.54, 0.93, 0.27, 0.47),
			      c(0.08, 0.47, 0.06, 0.26), c(0.54, 0.93, 0.06, 0.26), c(0.93,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 1, 3 ), screen = 1 )
	split.screen( figs = c( 1, 3 ), screen = 2 )
	split.screen( figs = c( 1, 3 ), screen = 3 )
	split.screen( figs = c( 1, 3 ), screen = 4 )
	split.screen( figs = c( 1, 3 ), screen = 5 )
	split.screen( figs = c( 1, 3 ), screen = 6 )
	split.screen( figs = c( 1, 3 ), screen = 7 )
	split.screen( figs = c( 1, 3 ), screen = 8 )

	cov_mod = 1

	for(tt in 1:TT){
		screen(9 + tt)
		par(pty="s") 
		par(mai=c(0.02,0.02,0.02,0.02))
		par(mgp=c(2,0.5,0))

      
		matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
      
		if(tt == 1){
			mtext(expression(s[y]), side = 2, line = 1.5, adj = 0.5, cex = 1.2, font=2)
			axis(2, at = seq(0, 1, length.out = 3), labels = seq(0, 1, length.out = 3), cex.axis = 1, tcl=-0.3)
		}
		if(tt == 2 & cov_mod == 1){
			mtext(expression(Z[1]), side = 3, line = 1, adj = 0.5, cex = 1.75, font=2, col = 'blue')
		}
		mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1.3)
		if(cov_mod == dim(REALIZATIONS_LIST)[3]){

			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 1.2, font=2)
      
			axis(1, at = seq(0, 1, length.out = 3), labels = seq(0, 1, length.out = 3), cex.axis = 1, tcl=-0.3)
		}
		
	}
	
	for(cov_mod in 1:3){

	    	for(tt in 1:TT){
			screen((cov_mod - 1) * 6 + 15 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.3,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, n * TT + (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
			if(tt == 1){
				mtext(expression(s[y]), side = 2, line = 1.5, adj = 0.5, cex = 1.2, font=2)
				text(-0.6, 0.475, mod_labels[cov_mod], col = 'blue', xpd = NA, cex = 2)
				axis(2, at = seq(0, 1, length.out = 3), labels = seq(0, 1, length.out = 3), cex.axis = 1, tcl=-0.3)
			}
	      		if(tt == 2 & cov_mod == 1){
				mtext(expression(Z[2]), side = 3, line = 0, adj = 0.5, cex = 1.75, font=2, col = 'blue')
	      		}
	      
			if(cov_mod == 3){
	     			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 1.2, font=2)
	      
	      			axis(1, at = seq(0, 1, length.out = 3), labels = seq(0, 1, length.out = 3), cex.axis = 1, tcl=-0.3)

			}
	    	}
	}

	cov_mod = 4

	for(tt in 1:TT){
		screen(12 + tt)
		par(pty="s") 
		par(mai=c(0.02,0.02,0.02,0.02))
		par(mgp=c(2,0.5,0))

      
		matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
      
		if(tt == 2 & cov_mod == 4){
			mtext(expression(Z[1]), side = 3, line = 1, adj = 0.5, cex = 1.75, font=2, col = 'blue')
		}
		mtext(paste('t = ', tt, sep = ''), side = 3, line = 0, adj = 0.5, cex = 1.3)
		if(cov_mod == dim(REALIZATIONS_LIST)[3]){

			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 1.2, font=2)
      
			axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 0.6, tcl=-0.3)
		}
		
	}
	
	for(cov_mod in 4:6){

	    	for(tt in 1:TT){
			screen((cov_mod - 4) * 6 + 18 + tt)
	      		par(pty="s") 
	      		par(mai=c(0.02,0.02,0.02,0.02))
	     		par(mgp=c(2,0.3,0))
	      
            		matrix_vals <- matrix(REALIZATIONS_LIST[8, n * TT + (tt - 1) * n + 1:n, cov_mod], N, N)

	      		poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = zlim_range1)
			if(tt == 1){
				text(-0.3, 0.475, mod_labels[cov_mod], col = 'blue', xpd = NA, cex = 2)
			}
	      		if(tt == 2 & cov_mod == 4){
				mtext(expression(Z[2]), side = 3, line = 0, adj = 0.5, cex = 1.75, font=2, col = 'blue')
	      		}
	      
			if(cov_mod == dim(REALIZATIONS_LIST)[3]){

	     			mtext(expression(s[x]), side = 1, line = 1.5, adj = 0.5, cex = 1.2, font=2)
	      
	      			axis(1, at = seq(0, 1, length.out = 3), labels = seq(0, 1, length.out = 3), cex.axis = 1, tcl=-0.3)
			}
	    	}
	}

	screen(9)
	x1 <- c(0.425,0.49,0.49,0.425)
	y1 <- c(0.35,0.35,0.65,0.65)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5), CEX = 1)
		  
	close.screen( all=TRUE)
	dev.off()
	  
}

if(fig2_newnew){

	DAT1_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT2_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT3_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_multiple_example3', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT4_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example1', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT5_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT6_cross <- read.table(paste(root, 'Data/multivariate_stationary_nonfrozen_crosscovariance_matern_added_dim_example3', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	DAT <- list()
	DAT[[1]] <- DAT1_cross
	DAT[[2]] <- DAT2_cross
	DAT[[3]] <- DAT3_cross
	DAT[[4]] <- DAT4_cross
	DAT[[5]] <- DAT5_cross
	DAT[[6]] <- DAT6_cross

	mod_labels <- c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')
	  
	new_grid_x <- matrix(dat3[, 1], N, N)
	new_grid_y <- matrix(dat3[, 2], N, N, byrow = F)

	pdf(file = paste(root, 'thesis-defense/Figures/multivariate_stationary_sec2_covariance_multiple.pdf', sep = ''), width = 20, height = 15)

	split.screen( rbind(c(0.05, 0.33, 0.54, 0.94), c(0.35, 0.63, 0.54, 0.94), c(0.65, 0.93, 0.54, 0.94), c(0.05, 0.33, 0.09, 0.48), c(0.35, 0.63, 0.09, 0.48), c(0.65, 0.93, 0.09, 0.48), c(0.96,0.99,0.05,0.93)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 3, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 3, 3 ), screen = 4 )
	split.screen( figs = c( 3, 3 ), screen = 5 )
	split.screen( figs = c( 3, 3 ), screen = 6 )
	 	
	for(cov_mod in 1:6){

		for(ll in 1:3){
			for(tt in 1:TT){
				screen((cov_mod - 1) * 9 + (ll - 1) * 3  + 7 + tt)
				par(pty="s") 
				par(mai=c(0.02,0.02,0.02,0.02))
				par(mgp=c(2,0.5,0))

				matrix_vals <- matrix(DAT[[cov_mod]][ll, n * TT + (tt - 1) * n + 1:n], N, N)

				poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab="", yaxt = 'n', xaxt = 'n', zlim = c(0, 0.5))
				abline(h = 0.5, v = 0.5, lty = 3, col = '#FFFFFF', lwd = 2)
		      
				if(tt == 1 & (cov_mod == 1 | cov_mod == 4)){
					mtext(expression(h[y]), side = 2, line = 2, adj = 0.5, cex = 1.5, font=2)
					axis(2, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 1.5, tcl=-0.3)
				}
				if(ll == 3 & cov_mod > 3){
					mtext(expression(h[x]), side = 1, line = 2.5, adj = 0.5, cex = 1.5, font=2)
					axis(1, at = seq(0.1, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 1.5, tcl=-0.3)
				}
				if(tt == 2 & ll == 1 & cov_mod < 4){
					mtext(mod_labels[cov_mod], side = 3, line = 2, adj = 0.5, cex = 2, col = 'blue')
				}
				if(tt == 2 & ll == 1 & cov_mod > 3){
					mtext(mod_labels[cov_mod], side = 3, line = 0.5, adj = 0.5, cex = 2, col = 'blue')
				}
				if(ll == 1 & cov_mod < 4){
					label1 <- bquote(paste(t[1], "=", .(tt), sep=""))
					mtext(label1, side = 3, line = 0, cex = 1.75)
				}
				if(tt == 3  & (cov_mod == 3 || cov_mod == 6)){
					label1 <- bquote(paste(t[2], "=", .(ll), sep=""))
					text(1.2, 0.475, label1, xpd = NA, cex = 1.75, srt = 270)
				}
			}
		}
	}

	screen(7)
	x1 <- c(0.01,0.15,0.15,0.01) + 0.09
	y1 <- c(0.3,0.3,0.7,0.7)
	  
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1.5)
		  
	close.screen( all=TRUE)
	dev.off()
	  
}


if(fig3){

	TT <- 3

	radian <- seq(0, 2 * pi, length = 50)
	radius <- 1
	hx = radius *  cos(radian)
	hy = radius *  sin(radian)

	sim_grid_locations <- cbind(c(0, hx), c(0, hy))
	n <- nrow(sim_grid_locations)

	DAT <- list()	

	sim_vals_temp <- read.table(paste(root, 'Results/simulation-crosscovariance-values', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT[[1]] <- sim_vals_temp[1:(nrow(sim_vals_temp) / 2), c(1, 3, 5, 11, 13, 15, 21, 23, 25)]
	DAT[[2]] <- sim_vals_temp[(nrow(sim_vals_temp) / 2) + 1:(nrow(sim_vals_temp) / 2), c(1, 3, 5, 11, 13, 15, 21, 23, 25)]
	DAT[[4]] <- DAT[[3]] <- DAT[[1]]

	lims <- max(abs(DAT[[1]]))

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	
	pdf(file = paste(root, 'thesis-defense/Figures/simulation1.pdf', sep = ''), width = 20, height = 20)

	split.screen( rbind(c(0.06, 0.50, 0.54, 0.93), c(0.53, 0.97, 0.54, 0.93),
			      c(0.06, 0.50, 0.1, 0.49), c(0.53, 0.97, 0.1, 0.49), c(0.99,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 3, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 3, 3 ), screen = 4 )
	 	
	for(cov_mod in 1:4){

		for(ll in 1:9){

			screen(5 + (cov_mod - 1) * 9 + ll)
			par(pty="s") 
			par(mai=c(0.02,0.02,0.02,0.02))
			plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), xaxt = 'n', yaxt = 'n§')

			abline(h = 0, v = 0, lty = 3)

			lines(x = DAT[[cov_mod]][, ll] * cos(radian), y = DAT[[cov_mod]][, ll] * sin(radian), col = 1)
	
			if(ll %in% c(1, 4, 7) & (cov_mod == 1 | cov_mod == 3)){

				mtext(expression(h[x]), side = 1, line = 1.5, adj = 0.5, cex = 0.75, font=2)
				#axis(1, at = seq(0, 0.9, length.out = 3), labels = seq(-0.5, 0.5, length.out = 3), cex.axis = 0.65, tcl=-0.3)
			}
		}
	}

	close.screen( all=TRUE)
	dev.off()
	  
}

if(fig3_new){

	TT <- 3

	radian <- seq(0, 2 * pi, length = 50)
	radius <- 1
	hx = radius * cos(radian)
	hy = radius * sin(radian)

	sim_grid_locations <- cbind(c(0, hx), c(0, hy))
	n <- nrow(sim_grid_locations)

	DAT <- list()	

	sim_vals_temp <- read.table(paste(root, 'Results/simulation-crosscovariance-values', sep = ''), sep = " ", header = FALSE) %>% as.matrix()
	DAT[[1]] <- sim_vals_temp[1:(nrow(sim_vals_temp) / 2), c(1, 3, 5, 11, 13, 15, 21, 23, 25)]
	DAT[[2]] <- sim_vals_temp[(nrow(sim_vals_temp) / 2) + 1:(nrow(sim_vals_temp) / 2), c(1, 3, 5, 11, 13, 15, 21, 23, 25)]

	lims <- max(abs(DAT[[1]]))

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')
	
	pdf(file = paste(root, 'thesis-defense/Figures/simulation1.pdf', sep = ''), width = 11, height = 6)

	split.screen( rbind(c(0.08, 0.52, 0.1, 0.9), c(0.55, 0.99, 0.1, 0.9), c(0.99,0.99,0.05,0.95)))
	  
	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 3, 3 ), screen = 2 )
	 	
	for(cov_mod in 1:2){
		for(ll in 1:9){

			screen(3 + (cov_mod - 1) * 9 + ll)
			par(pty="s") 
			par(mai=c(0.02,0.02,0.02,0.02))

			if(ll %in% c(1, 4) & cov_mod == 1){
				plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), xaxt = 'n', yaxt = 'n')
				mtext(expression(h[y]), side = 2, line = 2, adj = 0.5, cex = 0.9, font=2)
				axis(2, at = seq(-lims + 0.03, lims - 0.03, length.out = 3), labels = seq(-round(lims - 0.03, 2), round(lims - 0.03, 2), length.out = 3), cex.axis = 0.9)
			}else if(ll == 7 & cov_mod == 1){
				plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), xaxt = 'n', yaxt = 'n')
				axis(2, at = seq(-lims + 0.03, lims - 0.03, length.out = 3), labels = seq(-round(lims - 0.03, 2), round(lims - 0.03, 2), length.out = 3), cex.axis = 0.9)
				axis(1, at = seq(-lims + 0.03, lims - 0.03, length.out = 3), labels = seq(-round(lims - 0.03, 2), round(lims - 0.03, 2), length.out = 3), cex.axis = 0.9)
				mtext(expression(h[x]), side = 1, line = 2, adj = 0.5, cex = 0.9, font=2)
				mtext(expression(h[y]), side = 2, line = 2, adj = 0.5, cex = 0.9, font=2)
				labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= 0.9")
				mtext(labs, side = 2, line = 3.3, adj = 0.5, cex = 1, font=1)
			}else if(ll == 7 & cov_mod == 2){
                                plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), yaxt = 'n', xaxt = 'n')
				axis(1, at = seq(-lims + 0.03, lims - 0.03, length.out = 3), labels = seq(-round(lims - 0.03, 2), round(lims - 0.03, 2), length.out = 3), cex.axis = 0.9)
				mtext(expression(h[x]), side = 1, line = 2, adj = 0.5, cex = 0.9, font=2)
                        }else if(ll %in% c(8, 9)){
                                plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), yaxt = 'n', xaxt = 'n')
				axis(1, at = seq(-lims + 0.03, lims - 0.03, length.out = 3), labels = seq(-round(lims - 0.03, 2), round(lims - 0.03, 2), length.out = 3), cex.axis = 0.9)
				mtext(expression(h[x]), side = 1, line = 2, adj = 0.5, cex = 0.9, font=2)
                        }else{
                                plot(0, 0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), xaxt = 'n', yaxt = 'n')
			}
			if(ll == 1){
				labs <- bquote(rho[bold(V[22]) *"|"* bold(V[11])] ~ "= -0.9")
				mtext(labs, side = 3, line = 0, adj = 0.5, cex = 1, font=1)
			}
			if(ll == 2){
				labs <- bquote(rho[bold(V[22]) *"|"* bold(V[11])] ~ "= 0")
				mtext(labs, side = 3, line = 0, adj = 0.5, cex = 1, font=1)
			}
			if(ll == 3){
				labs <- bquote(rho[bold(V[22]) *"|"* bold(V[11])] ~ "= 0.9")
				mtext(labs, side = 3, line = 0, adj = 0.5, cex = 1, font=1)
			}
			if(ll == 4 & cov_mod == 1){
				labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= 0")
				mtext(labs, side = 2, line = 3.3, adj = 0.5, cex = 1, font=1)
			}
			if(ll == 1 & cov_mod == 1){
				labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= -0.9")
				mtext(labs, side = 2, line = 3.3, adj = 0.5, cex = 1, font=1)
			}
			if(ll == 2 & cov_mod == 1){
				label1 <- bquote(paste(t[1], " = ", 1, ', ', t[2], " = ", 1, sep=""))
				mtext(label1, side = 3, line = 1.5, adj = 0.5, cex = 1.5, font=1, col = 4)
			}
			if(ll == 2 & cov_mod == 2){
				label1 <- bquote(paste(t[1], " = ", 2, ', ', t[2], " = ", 2, sep=""))
				mtext(label1, side = 3, line = 1.5, adj = 0.5, cex = 1.5, font=1, col = 4)
			}
			abline(h = 0, v = 0, lty = 3)
			lines(x = DAT[[cov_mod]][, ll] * cos(radian), y = DAT[[cov_mod]][, ll] * sin(radian), col = 1)
		}
	}

	close.screen( all=TRUE)
	dev.off()
}



if(fig4){

	pdf(file = paste(root, 'thesis-defense/Figures/simulation2.pdf', sep = ''), width = 16, height = 4)
	
	split.screen( rbind(c(0.03,0.5,0.1,0.9), c(0.52,0.99,0.1,0.9), c(0.99,0.99,0.1,0.95)))
	split.screen( figs = c( 1, 3 ), screen = 1 ) 
	split.screen( figs = c( 1, 3 ), screen = 2 ) 

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config1_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	#sigs_est <- matrix(, ncol = 3, nrow = nrow(p))
	#for(samps in 1:nrow(p)){
	#	wind_var_chol <- matrix(c(p[samps, 9], p[samps, 10], 0, p[samps, 11]), ncol = 2, byrow = T)
	#	wind_var <- t(wind_var_chol) %*% wind_var_chol
	#	sigs_est[samps,] <- c(wind_var[1, 1], wind_var[1, 2], wind_var[2, 2])
	#}

	screen(4)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	#mtext("Parameter Estimate", side = 2, line = 3, cex = 1.5)
	#labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= -0.9" * " and " * rho[bold(V[22]) *"|"* bold(V[11])] ~ "= -0.9")
	mtext("k = 1", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config3_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	screen(5)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	#labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= -0.9" * " and " * rho[bold(V[22]) *"|"* bold(V[11])] ~ "= 0")
	mtext("k = 2", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config5_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	screen(6)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	#labs <- bquote(rho[bold(V[11]) *"|"* bold(V[22])] ~ "= -0.9" * " and " * rho[bold(V[22]) *"|"* bold(V[11])] ~ "= 0.9")
	mtext("k = 3", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MLOE', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(7)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n', ylim = c(0, max(DAT1)))
	mtext("MLOE", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MMOM', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(8)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n', ylim = c(min(DAT1), 0))
	mtext("MMOM", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MSE', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(9)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n')
	mtext("MSE", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	close.screen( all=TRUE)
	dev.off()
}


if(fig5){

	pdf(file = paste(root, 'thesis-defense/Figures/simulation3.pdf', sep = ''), width = 16, height = 4)
	
	split.screen( rbind(c(0.03,0.5,0.1,0.9), c(0.52,0.99,0.1,0.9), c(0.99,0.99,0.1,0.95)))
	split.screen( figs = c( 1, 3 ), screen = 1 ) 
	split.screen( figs = c( 1, 3 ), screen = 2 ) 

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config1_mu_set2_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	screen(4)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	mtext("k = 1", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config3_mu_set2_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	screen(5)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	mtext("k = 2", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_parameter_estimation_multi_config5_mu_set2_NEW', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	p <- DAT1
	theta <- cbind(exp(p[, 1:5]), p[, 6])
	mu <- p[, 7:8]
	std_theta <- theta - matrix(c(1, 1, 0.23, 0.5, 1, 0.5), nrow = nrow(theta), ncol = ncol(theta), byrow = T)
	std_theta <- std_theta / apply(theta, 2, sd)

	screen(6)
	par(mai=c(0.05,0.05,0.05,0.05))
	boxplot(cbind(std_theta), ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-12, 65))
	axis(1, at = 1:6, labels = c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]), expression(theta[6])))
	mtext(expression(bold(Theta)), side = 1, line = 3, cex = 2)
	mtext("k = 3", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MLOE_mu_set2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(7)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n', ylim = c(0, max(DAT1)))
	mtext("MLOE", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MMOM_mu_set2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(8)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n', ylim = c(min(DAT1), 0))
	mtext("MMOM", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	DAT1 <- read.table(paste(root, 'Results/multivariate_stationary_simulation_MSE_mu_set2', sep = ''), sep = " ", header = FALSE) %>% as.matrix()

	screen(9)
	par(mai=c(0.05, 0.4, 0.05, 0.4))
	boxplot(DAT1, ylab = '', xlab = '', xaxt = 'n')
	mtext("MSE", side = 3, line = 1, adj = 0.5, cex = 1.5, font=1)
	axis(1, at = 1:3, labels = c("k = 1", "k = 2", "k = 3"))

	close.screen( all=TRUE)
	dev.off()
}


