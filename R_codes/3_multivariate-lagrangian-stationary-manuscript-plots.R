
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'multiple-advections/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

fig1_newnew = T
fig2_newnew = T

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

	pdf(file = paste(root, 'Figures/multivariate_stationary_sec2_example_multiple.pdf', sep = ''), width = 14, height = 9)

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

	pdf(file = paste(root, 'Figures/multivariate_stationary_sec2_covariance_multiple.pdf', sep = ''), width = 20, height = 15)

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

