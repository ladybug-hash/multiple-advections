
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'multiple-advections/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R",sep=''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R",sep=''))

yr <- 2019

dat <- read.table(paste(root, 'Data/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat2 <- read.table(paste(root, 'Data/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat3 <- read.table(paste(root, 'Data/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

dat[which(dat < -25)] <- -25
dat2[which(dat2 < -25)] <- -25

pdf(file = paste(root, 'Figures/spacetime-maps-residuals-manuscript-NEW.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in 139:139){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	zlim_range1 <- zlim_range2 <- range(c(dat[hr_index, ], dat2[hr_index, ]))

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	hr_label <- c('0:00', '3:00', '6:00', '9:00', '12:00', '15:00', '18:00', '21:00')
	for(hr in hr_index){
		
		hr_count <- hr_count + 1
		
		for(variable in 1:2){
			
			screen((variable - 1) * 5 + 2 + hr_count)
			par(pty = 's')
			par(mai=c(0.2,0.2,0.2,0.2))
			
			if(hr_count == 1 & variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], nx = 25, ny = 25, zlim = zlim_range2, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('925 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(hr_count == 1 & variable == 1){
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], nx = 25, ny = 25, zlim = zlim_range1, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('880 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], nx = 25, ny = 25, zlim = zlim_range2, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}else{
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], nx = 25, ny = 25, zlim = zlim_range1, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
			points(cbind(dat3[190, 1], dat3[190, 2]), col = 'black', pch = 4, cex = 4, lwd = 6)
			points(cbind(dat3[422, 1], dat3[422, 2]), col = 'black', pch = 4, cex = 4, lwd = 6)
			
			if(hr_count == 1){
				mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			}

			if(variable == 1){
				if(mod(hr, 8) == 0){
					mtext(hr_label[1], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
				}else{
					mtext(hr_label[mod(hr, 8) + 1], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
				}
			}else{
				mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
			}
		}				
	}
	screen(2)
	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.15,0.15,0.8,0.8)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range2[1], zlim_range2[2], length.out = 5), 0), CEX = 2)

	close.screen( all=TRUE)

}

dev.off()


