
workstation = F

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

############################################     SIMULATION STUDY 1     ############################################

RHO_VAL <- c(-0.9, -0.6, -0.3, 0.3, 0.6, 0.9)

MSE <- array(, dim = c(50, 12, 3))

for(aa in 1:3){
	for(set in 1:50){
		for(rr in 1:length(RHO_VAL)){
			p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_M1_config', aa,'_set_', set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			MSE[set, (rr - 1) * 2 + 1, aa] <- p[length(p)]

			p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_M2_config', aa,'_set_', set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			MSE[set, (rr - 1) * 2 + 2, aa] <- p[length(p)]
		}
	}
}

pdf(file = paste(root, 'thesis-defense/Figures/simulation-parameter-estimates1.pdf', sep = ''), width = 15, height = 6)

split.screen( rbind(c(0.05,0.98,0.05,0.97), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 )

screen(3)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[,,1], ylim = range(MSE), xaxt = 'n', xlab = '', col = rep(c("#d9534f", "#3fb1e2"), 6))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(a)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)
mtext('MSE', side = 2, line = 3, cex = 1.5)

screen(4)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[,,2], ylim = range(MSE), xaxt = 'n', yaxt = 'n', xlab = '', col = rep(c("#d9534f", "#3fb1e2"), 6))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(b)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)

screen(5)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[,,3], ylim = range(MSE), xaxt = 'n', yaxt = 'n', xlab = '', col = rep(c("#d9534f", "#3fb1e2"), 6))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(c)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)

close.screen( all=TRUE)
dev.off()


############################################     SIMULATION STUDY 2     ############################################

RHO_VAL <- 0.9

TRUE_PARAM <- rep(c(1, 1, 0.05, 1, 1, RHO_VAL, 0.1, 0.1, 0.1, 0.1, 0), each  = 3)

PARAM <- PARAM_TEMP <- matrix(, nrow = 50, ncol = 11 * 3)

for(aa in 1:3){

	for(set in 1:50){
		p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_2_M3_config', aa,'_set_', 150 + set, '_rho_', RHO_VAL, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
		wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol

		PARAM_TEMP[set, aa - 1 + seq(1, 33, by = 3)] <- c(exp(p[1:5]), p[6:8], wind_var[1, 1], wind_var[2, 2], wind_var[1, 2])
	}
}


for(aa in 1:33){
	PARAM[, aa] <- (PARAM_TEMP[, aa] - TRUE_PARAM[aa])  / sd(PARAM_TEMP[, aa])
}


LABELS <- c(bquote(underline(hat(sigma)[11]^2)), bquote(underline(hat(sigma)[22]^2)), bquote(underline(hat(a))), bquote(underline(hat(nu)[11])), bquote(underline(hat(nu)[22])), bquote(underline(hat(rho))), bquote(underline(hat(mu)[1])), bquote(underline(hat(mu)[2])), bquote(underline(hat(Sigma)[paste(1, ",",1, sep = "")])), bquote(underline(hat(Sigma)[paste(2, ",",2, sep = "")])), bquote(underline(hat(Sigma)[paste(1, ",", 2, sep = "")])), expression(hat(Sigma)[paste(1, ",", 2, sep = "")]))

pdf(file = paste(root, 'thesis-defense/Figures/simulation-parameter-estimates2.pdf', sep = ''), width = 15, height = 6)

boxplot(PARAM, medcol = rep(c("#00FFFF", "#f97432", "#916fcd"), 3), col = rep(c("#00FFFF", "#f97432", "#916fcd"), 3), xaxt = 'n')

axis(1, at = seq(2, 33, by = 3), labels = LABELS[1:11], cex.axis = 1.5, tick = F, line=0.5)
mtext('Standardized Estimates', side = 2, line = 3, cex = 1.5)

dev.off()

pdf(file = paste(root, 'thesis-defense/Figures/parameter-estimates.pdf', sep = ''), width = 25, height = 10)

split.screen( rbind(c(0.05,0.99,0.55,0.95), c(0.05,0.42,0.05,0.45), c(0.4,0.99,0.05,0.45)))
split.screen( figs = c( 1, 4 ), screen = 1 )
#split.screen( rbind(c(0.02,0.32,0.5,0.95), c(0.34,0.64,0.5,0.95), c(0.66,0.99,0.5,0.95), c(0.02,0.32,0.05,0.5), c(0.34,0.64,0.05,0.5), c(0.66,0.99,0.05,0.5)))
#split.screen( rbind(c(0.02,0.17,0.1,0.95), c(0.17,0.32,0.1,0.95), c(0.27,0.42,0.1,0.95), c(0.41, 0.55,0.1,0.95), c(0.5,0.69,0.1,0.95), c(0.71,0.99,0.1,0.95)))

screen(4)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
#boxplot(p_vals[, 1:2], outwex = 0.5, boxwex = 0.2)
boxplot(p_vals[, 1:2], boxwex = 0.1, at = c(0.1, 0.25), xlim = c(0, 0.4), xaxt = 'n', frame = F, cex.axis = 2)
axis(1, at = c(0.1, 0.25), labels = c(expression(hat(sigma)[11]^2), expression(hat(sigma)[22]^2)), cex.axis = 2, tick = F, col.axis = 4)

screen(5)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(p_vals[, 3], boxwex = 0.05, at = 0.04, xlim = c(0, 0.1), xaxt = 'n', frame = F, cex.axis = 2)
axis(1, at = 0.04, labels = expression(hat(a)), cex.axis = 2, tick = F, col.axis = 4)

screen(6)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(p_vals[, 4:5], boxwex = 0.1, at = c(0.1, 0.25), xlim = c(0, 0.4), xaxt = 'n', frame = F, cex.axis = 2)
axis(1, at = c(0.1, 0.25), labels = c(expression(hat(nu)[11]), expression(hat(nu)[22])), cex.axis = 2, tick = F, col.axis = 4)

screen(7)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(p_vals[, 6], boxwex = 0.05, at = 0.04, xlim = c(0, 0.1), xaxt = 'n', frame = F, cex.axis = 1.75)
axis(1, at = 0.04, labels = expression(hat(rho)), cex.axis = 2, tick = F, col.axis = 4)

screen(2)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(p_vals[, 7:10], boxwex = 0.65, at = seq(0.5, 3.5, length.out = 4), xlim = c(0, 5), xaxt = 'n', frame = F, cex.axis = 1.5)
axis(1, at = seq(0.5, 3.5, length.out = 4), labels = c(expression(hat(mu)[1]), expression(hat(mu)[2]), expression(hat(mu)[3]), expression(hat(mu)[4])), cex.axis = 2, tick = F, col.axis = 4)

screen(3)
#par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(p_vals[, 11:20], ylim = c(-0.045, 0.26), boxwex = 4, at = seq(1, 50, length.out = 10), xlim = c(0, 52), xaxt = 'n', frame = F, cex.axis = 2)
axis(1, at = seq(1, 50, length.out = 10), labels = c(expression(hat(Sigma)[paste(1, ",",1, sep = "")]), expression(hat(Sigma)[paste(2, ",",2, sep = "")]), expression(hat(Sigma)[paste(3, ",",3, sep = "")]), expression(hat(Sigma)[paste(4, ",",4, sep = "")]), expression(hat(Sigma)[paste(1, ",",2, sep = "")]), expression(hat(Sigma)[paste(1, ",",3, sep = "")]), expression(hat(Sigma)[paste(1, ",",4, sep = "")]), expression(hat(Sigma)[paste(2, ",",3, sep = "")]), expression(hat(Sigma)[paste(2, ",",4, sep = "")]), expression(hat(Sigma)[paste(3, ",",4, sep = "")])), cex.axis = 2, tick = F, col.axis = 4)

close.screen( all=TRUE)
dev.off()


############################################     SIMULATION STUDY 2 BOXPLOT 2     ############################################

MSE <- array(, dim = c(50, length(RHO_VAL) * 2, 3))

for(aa in 1:3){

	for(set in 1:50){
		for(rr in 1:length(RHO_VAL)){

			if(rr < 4){
				p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_2_M3_config', aa,'_set_', set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			}else{
				p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_2_M3_config', aa,'_set_', 100 + set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			}
			MSE[set, (rr  - 1) * 2 + 1, aa] <- p[12]

			if(rr < 4){
				p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_2_M2_config', aa,'_set_', set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			}else{
				p <- read.table(paste(root, 'Results/multivariate_stationary_simulated_MSE_2_M2_config', aa,'_set_', 100 + set, '_rho_', RHO_VAL[rr], sep = ''), header = FALSE, sep = " ") %>% as.matrix()
			}
			MSE[set, (rr  - 1) * 2 + 2, aa] <- p[21]
		}
	}
}

pdf(file = paste(root, 'thesis-defense/Figures/simulation-parameter-estimates3.pdf', sep = ''), width = 15, height = 6)

split.screen( rbind(c(0.05,0.98,0.05,0.97), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 )

screen(3)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[, , 1], ylim = range(MSE), xaxt = 'n', xlab = '', col = rep(c("#93c44b", "#3fb1e2"), 3))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(d)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)
mtext('MSE', side = 2, line = 3, cex = 1.5)

screen(4)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[, , 2], ylim = range(MSE), xaxt = 'n', yaxt = 'n', xlab = '', col = rep(c("#93c44b", "#3fb1e2"), 3))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(e)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)

screen(5)
par(pty="s") 
par(mai=c(0.2,0.2,0.2,0.2))
boxplot(MSE[, , 3], ylim = range(MSE), xaxt = 'n', yaxt = 'n', xlab = '', col = rep(c("#93c44b", "#3fb1e2"), 3))
axis(1, at = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5), labels = RHO_VAL, cex.axis = 1.2, tick = F, line=0.5)
mtext('(f)', side = 3, line = 1, cex = 1.5)
mtext(expression(rho), side = 1, line = 3, cex = 1.5)

close.screen( all=TRUE)
dev.off()
