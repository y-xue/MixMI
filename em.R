library(parallel)
library(GPfit)
source("adaptive_gd.R")
source("util.R")
source('mixtureMITemporalConfig.R')
source('GP.R')

L <- function(pi1, pi2, pi3, w1, w2, w3, s, m_reg1, sig_reg1, m_reg2, sig_reg2, m_gp, s2_gp, epsilon=1e-8) {
	p1 = dnorm(s,m_reg1,sig_reg1)
	loglik_1 = 0
	if (w1 == 0) {
		loglik_1 = log_with_limits(pi1 * p1, epsilon)
	}
	else {
		loglik_1 = w1 * log_with_limits(pi1 * p1 / w1, epsilon)
	}

	p2 = dnorm(s,m_reg2,sig_reg2)
	loglik_2 = 0
	if (w2 == 0) {
		loglik_2 = log_with_limits(pi2 * p2, epsilon)
	}
	else {
		loglik_2 = w2 * log_with_limits(pi2 * p2 / w2, epsilon)
	}

	if (is.na(s2_gp)) {
		# to be changed
		p3 = 0
	} else {
		p3 = dnorm(s,m_gp,sqrt(s2_gp))
	}
	loglik_3 = 0
	if (w3 == 0) {
		loglik_3 = log_with_limits(pi3 * p3, epsilon)
	}
	else {
		loglik_3 = w3 * log_with_limits(pi3 * p3 / w3, epsilon)
	}

	res = loglik_1 + loglik_2 + loglik_3

	return(res)
}

L_rr <- function(pi1, pi2, w1, w2, s, m1, sig1, m2, sig2, epsilon=1e-8) {

	p1 = dnorm(s,m1,sig1)
	loglik_1 = 0
	if (w1 == 0) {
		loglik_1 = log_with_limits(pi1 * p1, epsilon)
	}
	else {
		loglik_1 = w1 * log_with_limits(pi1 * p1 / w1, epsilon)
	}

	p2 = dnorm(s,m2,sig2)
	loglik_2 = 0
	if (w2 == 0) {
		loglik_2 = log_with_limits(pi2 * p2, epsilon)
	}
	else {
		loglik_2 = w2 * log_with_limits(pi2 * p2 / w2, epsilon)
	}

	res = loglik_1 + loglik_2

	return(res)
}

em_rrg_obs_only <- function(S,Z,Yreg,Ygp,xte_vec_tr,xtr_vec_tr,t,r_v_tr,mix_model_num,w1,w2,w3,pi1,pi2,pi3,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn)
{
	print("em_rrg_obs_only")

	epsilon = 1e-8

	N <- length(S)

	param <- list(lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,l,pi1,pi2,pi3,w1,w2,w3,-Inf,Inf,mix_model_num)
	names(param) <- c('lr_beta1','lr_sigma1','lr_beta2','lr_sigma2','ll','pi1','pi2','pi3','w1','w2','w3','loglik','abs_error','mix_model_num')

	print(sprintf("prev pi1: %s",pi1))
	print(sprintf("prev pi2: %s",pi2))
	print(sprintf("prev pi3: %s",pi3))

	if (pi3 == 0) {
		pi3 = 0.005
		pi1 = pi1/(pi1+pi2) * (1-pi3)
		pi2 = 1-pi1-pi3
		w3 = rep(pi3,N)
		w1 = w1/(w1+w2) * (1-pi3)
		w2 = 1-w1
	}

	M <- rep(0,N)
	K <- rep(0,N)
	sig2vec <- rep(0,N)

	Rinv_lst = mclapply(1:N, function(i) Rinverse(l,xtr_vec_tr[i,][r_v_tr[i,]]), mc.cores=num_cores)
	M = unlist(mclapply(1:N,function (i) yhat(l,xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
	sig2vec = unlist(mclapply(1:N, function(i) sig2(l,Ygp[i,-t][r_v_tr[i,]],xtr_vec_tr[i,][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
	K = unlist(mclapply(1:N, function(i) s2(l,sig2vec[i],xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))

	loglik = sum(unlist(mclapply(1:N, function(i) L(pi1,pi2,pi3,w1[i],w2[i],w3[i],S[i],Z[i,]%*%lr_beta1,lr_sigma1,Yreg[i,]%*%lr_beta2,lr_sigma2,M[i],K[i]), mc.cores=num_cores)))
	prev_loglik <- loglik + tolerance + 1
	param$loglik = loglik

	ll = l
	j <- 1

	print(sprintf("loglik: %s",loglik))
	print(sprintf("num predictions: %s",length(S)))
	print(sprintf("abs error: %s",sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Yreg%*%lr_beta2)+pi3*M)))))
	GP_pred_error = sum(abs(S - M))
	print(sprintf("GP abs error: %s",GP_pred_error))

	while ((abs(loglik-prev_loglik) > tolerance && j < em_max_iter)) {
		print(sprintf("iter: %d",j))
		# E-step
		# 0.4s
		# ptm <- proc.time()
		for (i in 1:N) {
			ry = Ygp[i,-t][r_v_tr[i,]]
			
			if (length(ry) == 0 || length(unique(ry)) == 1) {
				if (length(ry) == 0 || unique(ry) != S[i]) {
					preg1 <- pi1 * dnorm(S[i],Z[i,]%*%lr_beta1,lr_sigma1)
					preg2 <- pi2 * dnorm(S[i],Yreg[i,]%*%lr_beta2,lr_sigma2)
					pGP = 1e-8
				} else {
					preg1 = 1e-8
					preg2 = 1e-8
					pGP = 1-preg1-preg2
				}
			} else {
				preg1 <- pi1 * dnorm(S[i],Z[i,]%*%lr_beta1,lr_sigma1)
				preg2 <- pi2 * dnorm(S[i],Yreg[i,]%*%lr_beta2,lr_sigma2)
				pGP <- pi3 * dnorm(S[i],M[i],sqrt(K[i]))
			}
			
			if (round(preg1 + preg2 + pGP, 8) == 0) {
				preg1 = pi1
				preg2 = pi2
				pGP = pi3
			}

			w1[i] <-  preg1 / (preg1 + preg2 + pGP)
			w2[i] <- preg2 / (preg1 + preg2 + pGP)
			w3[i] <- pGP / (preg1 + preg2 + pGP)

			if (w1[i] < epsilon) {
				if (w2[i] < w3[i]) {
					w3[i] = w3[i] - (epsilon-w1[i])
				} else {
					w2[i] = w2[i] - (epsilon-w1[i])
				}
				w1[i] = epsilon
			}

			if (w2[i] < epsilon) {
				if (w1[i] < w3[i]) {
					w3[i] = w3[i] - (epsilon-w2[i])
				} else {
					w1[i] = w1[i] - (epsilon-w2[i])
				}
				w2[i] = epsilon
			}

			if (w3[i] < epsilon) {
				if (w1[i] < w2[i]) {
					w2[i] = w2[i] - (epsilon-w3[i])
				} else {
					w1[i] = w1[i] - (epsilon-w3[i])
				}
				w3[i] = epsilon
			}
		}

		# M-step
		pi1 <- (sum(w1)) / N
		pi2 <- (sum(w2)) / N
		pi3 <- (sum(w3)) / N

		print(sprintf("pi1: %s",pi1))
		print(sprintf("pi2: %s",pi2))
		print(sprintf("pi3: %s",pi3))

		sw = sqrt(w1)
		swZ = sw * Z
		swS = sw * S
		ztwz <- t(swZ) %*% swZ
		ridge <- 0.00001
		pen <- ridge * diag(ztwz)
		if (length(pen)==1) pen <- matrix(pen)
		v <- solve(ztwz + diag(pen))
		lr_beta1 <- v %*% t(swZ) %*% swS
		residuals <- swS - swZ %*% lr_beta1
		lr_sigma1 <- sqrt((sum(residuals^2))/sum(w1))
		# print(proc.time() - ptm)

		sw2 = sqrt(w2)
		swY = sw2 * Yreg
		swS = sw2 * S
		ztwz <- t(swY) %*% swY
		ridge <- 0.00001
		pen <- ridge * diag(ztwz)
		if (length(pen)==1) pen <- matrix(pen)
		v <- solve(ztwz + diag(pen))
		lr_beta2 <- v %*% t(swY) %*% swS
		residuals <- swS - swY %*% lr_beta2
		lr_sigma2 <- sqrt((sum(residuals^2))/sum(w2))

		ll <- Adam_one_obs_only(w3,pi3,Ygp,S,xte_vec_tr,xtr_vec_tr,Rinv_lst,t,r_v_tr,ll,step,gd_precision,gd_miter,j)

		Rinv_lst = mclapply(1:N, function(i) Rinverse(ll,xtr_vec_tr[i,][r_v_tr[i,]]), mc.cores=num_cores)
		M = unlist(mclapply(1:N,function (i) yhat(ll,xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
		sig2vec = unlist(mclapply(1:N, function(i) sig2(ll,Ygp[i,-t][r_v_tr[i,]],xtr_vec_tr[i,][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
		K = unlist(mclapply(1:N, function(i) s2(ll,sig2vec[i],xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))

		prev_loglik <- loglik
		loglik = sum(unlist(mclapply(1:N, function(i) L(pi1,pi2,pi3,w1[i],w2[i],w3[i],S[i],Z[i,]%*%lr_beta1,lr_sigma1,Yreg[i,]%*%lr_beta2,lr_sigma2,M[i],sqrt(K[i])), mc.cores=num_cores)))
		print(sprintf("loglik: %s",loglik))
		
		pred_error = sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Yreg%*%lr_beta2)+pi3*M)))
		print(sprintf("abs error: %s",pred_error))

		GP_pred_error = sum(abs(S - M))
		print(sprintf("GP abs error: %s",GP_pred_error))

		if (GP_pred_error < pred_error) {
			param$mix_model_num = 0
			pred_error = GP_pred_error
		} else {
			param$mix_model_num = 1
		}

		if (pred_error < param$abs_error) {
			param$abs_error = pred_error
		} else {
			break
		}

		if (loglik > param$loglik) {
			param$loglik <- loglik
			param$w1 <- w1
			param$w2 <- w2
			param$w3 <- w3
			param$pi1 <- pi1
			param$pi2 <- pi2
			param$pi3 <- pi3
			param$lr_beta1 <- lr_beta1
			param$lr_sigma1 <- lr_sigma1
			param$lr_beta2 <- lr_beta2
			param$lr_sigma2 <- lr_sigma2
			param$ll <- ll
		}

		j = j + 1
	}

	return(param)
}

em_double_reg <- function(S,Z,Y,T,t,w1,w2,pi1,pi2,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,em_max_iter,tolerance) {
	print("em_double_reg")

	epsilon = 1e-8

	N = length(S)
	# initialization
	param = list(lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,pi1,pi2,w1,w2,-Inf,Inf,2)
	names(param) = c('lr_beta1','lr_sigma1','lr_beta2','lr_sigma2','pi1','pi2','w1','w2','loglik','abs_error','mix_model_num')

	print(sprintf("prev pi1: %s",pi1))
	print(sprintf("prev pi2: %s",pi2))

	# 1s
	loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],S[i],Z[i,]%*%lr_beta1,lr_sigma1,Y[i,]%*%lr_beta2,lr_sigma2), mc.cores=num_cores)))
	prev_loglik = loglik + tolerance + 1
	param$loglik = loglik

	j = 1

	print(sprintf("loglik: %s",loglik))
	print(sprintf("num predictions: %s",length(S)))
	print(sprintf("abs error: %s",sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Y%*%lr_beta2))))))

	while ((abs(loglik-prev_loglik) > tolerance && j < em_max_iter)) {
		print(sprintf("iter: %d",j))
		# E-step
		# 0.4s
		# ptm <- proc.time()
		for (i in 1:N) {

			p1 = pi1 * dnorm(S[i],Z[i,]%*%lr_beta1,lr_sigma1)
			p2 = pi2 * dnorm(S[i],Y[i,]%*%lr_beta2,lr_sigma2)

			if (round(p1+p2,8) == 0) {
				p1 = pi1
				p2 = pi2
			}

			w1[i] = p1 / (p1 + p2)
			w2[i] = p2 / (p1 + p2)

			if (w1[i] < epsilon) {
				w1[i] = epsilon
				w2[i] = 1 - epsilon
			}

			if (w2[i] < epsilon) {
				w2[i] = epsilon
				w1[i] = 1 - epsilon
			}
		}

		# M-step
		pi1 <- (sum(w1)) / N
		pi2 <- (sum(w2)) / N

		print(sprintf("pi1: %s",pi1))
		print(sprintf("pi2: %s",pi2))

		if (pi1 == 0 || round(pi1,8) == 0) {
			return(param)
		}
		if (pi2 == 0 || round(pi2,8) == 0) {
			return(param)
		}

		# 0.1
		# ptm <- proc.time()
		sw1 = sqrt(w1)
		swZ = sw1 * Z
		swS = sw1 * S
		ztwz <- t(swZ) %*% swZ
		ridge <- 0.00001
		pen <- ridge * diag(ztwz)
		if (length(pen)==1) pen <- matrix(pen)
		v <- solve(ztwz + diag(pen))
		lr_beta1 <- v %*% t(swZ) %*% swS
		residuals <- swS - swZ %*% lr_beta1
		lr_sigma1 <- sqrt((sum(residuals^2))/sum(w1))
		# print(proc.time() - ptm)

		# ptm <- proc.time()
		sw2 = sqrt(w2)
		swY = sw2 * Y
		swS = sw2 * S
		ztwz <- t(swY) %*% swY
		ridge <- 0.00001
		pen <- ridge * diag(ztwz)
		if (length(pen)==1) pen <- matrix(pen)
		v <- solve(ztwz + diag(pen))
		lr_beta2 <- v %*% t(swY) %*% swS
		residuals <- swS - swY %*% lr_beta2
		lr_sigma2 <- sqrt((sum(residuals^2))/sum(w2))
		

		prev_loglik <- loglik
		# print(prev_loglik)
		# loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],S[i],Z[i,],Y[i,],lr_beta1,lr_sigma1,lr_beta2,lr_sigma2), mc.cores=num_cores)))
		loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],S[i],Z[i,]%*%lr_beta1,lr_sigma1,Y[i,]%*%lr_beta2,lr_sigma2), mc.cores=num_cores)))
		print(sprintf("loglik: %s",loglik))
		pred_error = sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Y%*%lr_beta2))))
		print(sprintf("abs error: %s",pred_error))

		# if (pred_error < param$abs_error) {
		# 	param$abs_error = pred_error

		# 	if (loglik > param$loglik) {
		# 		param$loglik <- loglik
		# 		param$w1 <- w1
		# 		param$w2 <- w2
		# 		param$pi1 <- pi1
		# 		param$pi2 <- pi2
		# 		param$lr_beta1 <- lr_beta1
		# 		param$lr_sigma1 <- lr_sigma1
		# 		param$lr_beta2 <- lr_beta2
		# 		param$lr_sigma2 <- lr_sigma2
		# 	}
		# }

		if (pred_error < param$abs_error) {
			param$abs_error = pred_error
		} else {
			break
		}

		if (loglik > param$loglik) {
			param$loglik <- loglik
			param$w1 <- w1
			param$w2 <- w2
			param$pi1 <- pi1
			param$pi2 <- pi2
			param$lr_beta1 <- lr_beta1
			param$lr_sigma1 <- lr_sigma1
			param$lr_beta2 <- lr_beta2
			param$lr_sigma2 <- lr_sigma2
		}

		

		j = j + 1
	}

	return(param)
}
