library(parallel)
library(GPfit)
source("adaptive_gd.R")
source("util.R")
source('mixtureMITemporalConfig.R')
source('GP.R')

L <- function(pi1, pi2, pi3, w1, w2, w3, U1, U2, U3, S1, S2, S3, s, Z, beta_reg1, sig_reg1, Y, beta_reg2, sig_reg2, Ygp, m_gp, s2_gp, epsilon=1e-8) {
	p1 = dnorm(s,Z%*%beta_reg1,sig_reg1) * dmvnorm(Z,U1,S1)
	loglik_1 = 0
	if (w1 == 0) {
		loglik_1 = log_with_limits(pi1 * p1, epsilon)
	}
	else {
		loglik_1 = w1 * log_with_limits(pi1 * p1 / w1, epsilon)
	}

	p2 = dnorm(s,Y%*%beta_reg2,sig_reg2) * dmvnorm(Y,U2,S2)
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
		p3 = dnorm(s,m_gp,sqrt(s2_gp)) * dmvnorm(Ygp,U3,S3)
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

L_rr <- function(pi1, pi2, w1, w2, U1, U2, S1, S2, s, Z, beta1, sig1, Y, beta2, sig2, epsilon=1e-8) {

	p1 = dnorm(s,Z%*%beta1,sig1) * dmvnorm(Z,U1,S1)
	loglik_1 = 0
	if (w1 == 0) {
		loglik_1 = log_with_limits(pi1 * p1, epsilon)
	}
	else {
		loglik_1 = w1 * log_with_limits(pi1 * p1 / w1, epsilon)
	}

	p2 = dnorm(s,Y%*%beta2,sig2) * dmvnorm(Y,U2,S2)
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

get_w <- function(N,t,S,Ygp,r_v_tr,Z,Yreg,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,M,K,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3,epsilon=1e-8) {
	w1 = rep(0,N); w2 = rep(0,N); w3 = rep(0,N)
	for (i in 1:N) {
		ry = Ygp[i,-t][r_v_tr[i,]]
		
		if (length(ry) == 0 || length(unique(ry)) == 1) {
			if (length(ry) == 0 || unique(ry) != S[i]) {
				preg1 <- pi1[i] * dnorm(S[i],Z[i,]%*%lr_beta1,lr_sigma1) * dmvnorm(Z[i,],U1,S1)
				preg2 <- pi2[i] * dnorm(S[i],Yreg[i,]%*%lr_beta2,lr_sigma2) * dmvnorm(Yreg[i,],U2,S2)
				pGP = 1e-8
			} else {
				# preg1 = 1e-8
				# preg2 = 1e-8
				# pGP = 1-preg1-preg2
				preg1 = pi1
				preg2 = pi2
				pGP = pi3
			}
		} else {
			preg1 <- pi1[i] * dnorm(S[i],Z[i,]%*%lr_beta1,lr_sigma1) * dmvnorm(Z[i,],U1,S1)
			preg2 <- pi2[i] * dnorm(S[i],Yreg[i,]%*%lr_beta2,lr_sigma2) * dmvnorm(Yreg[i,],U2,S2)
			pGP <- pi3[i] * dnorm(S[i],M[i],sqrt(K[i])) * dmvnorm(Ygp[i,-t],U3,S3)
		}

		if (round(preg1 + preg2 + pGP, 8) == 0) {
			preg1 = pi1[i]
			preg2 = pi2[i]
			pGP = pi3[i]
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
	res = list(w1,w2,w3)
	names(res) = c("w1","w2","w3")

	return(res)
}

get_w_rr <- function(N,S,Z,Y,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,pi1,pi2,U1,U2,S1,S2,epsilon=1e-8) {
	w1 = rep(0,N); w2 = rep(0,N)
	for (p in 1:N) {

        p1 = pi1 * dnorm(S[p],Z[p,]%*%lr_beta1,lr_sigma1) * dmvnorm(Z[i,],U1,S1)
        p2 = pi2 * dnorm(S[p],Y[p,]%*%lr_beta2,lr_sigma2) * dmvnorm(Y[i,],U2,S2)

        if (round(p1+p2,8) == 0) {
            p1 = pi_1
            p2 = pi_2
        }

        w1[p] = p1 / (p1 + p2)
        w2[p] = p2 / (p1 + p2)

        if (w1[p] < epsilon) {
            w1[p] = epsilon
            w2[p] = 1 - epsilon
        }

        if (w2[p] < epsilon) {
            w2[p] = epsilon
            w1[p] = 1 - epsilon
        }
    }
    res = list(w1,w2)
	names(res) = c("w1","w2")

	return(res)
}

get_ww <- function(N,t,Ygp,Z,Yreg,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3,epsilon=1e-8) {
	w1 = rep(0,N); w2 = rep(0,N); w3 = rep(0,N)
	for (i in 1:N) {
		
		preg1 <- pi1[i] * dmvnorm(Z[i,],U1,S1)
		preg2 <- pi2[i] * dmvnorm(Yreg[i,],U2,S2)
		pGP <- pi3[i] * dmvnorm(Ygp[i,-t],U3,S3)

		if (round(preg1 + preg2 + pGP, 8) == 0) {
			preg1 = pi1[i]
			preg2 = pi2[i]
			pGP = pi3[i]
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
	res = list(w1,w2,w3)
	names(res) = c("w1","w2","w3")

	return(res)
}

get_ww_rr <- function(N,Z,Y,pi1,pi2,U1,U2,S1,S2,epsilon=1e-8) {
	w1 = rep(0,N); w2 = rep(0,N)
	for (p in 1:N) {

        p1 = pi1 * dmvnorm(Z[i,],U1,S1)
        p2 = pi2 * dmvnorm(Y[i,],U2,S2)

        if (round(p1+p2,8) == 0) {
            p1 = pi_1
            p2 = pi_2
        }

        w1[p] = p1 / (p1 + p2)
        w2[p] = p2 / (p1 + p2)

        if (w1[p] < epsilon) {
            w1[p] = epsilon
            w2[p] = 1 - epsilon
        }

        if (w2[p] < epsilon) {
            w2[p] = epsilon
            w1[p] = 1 - epsilon
        }
    }
    res = list(w1,w2)
	names(res) = c("w1","w2")

	return(res)
}

em_rrg_obs_only <- function(S,Z,Yreg,Ygp,xte_vec_tr,xtr_vec_tr,t,r_v_tr,mix_model_num,w1,w2,w3,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn)
{
	print("em_rrg_obs_only")

	epsilon = 1e-8

	N <- length(S)

	param <- list(lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,l,pi1,pi2,pi3,w1,w2,w3,-Inf,Inf,mix_model_num)
	names(param) <- c('lr_beta1','lr_sigma1','lr_beta2','lr_sigma2','ll','pi1','pi2','pi3','w1','w2','w3','loglik','abs_error','mix_model_num')

	print(sprintf("prev pi1: %s",pi1))
	print(sprintf("prev pi2: %s",pi2))
	print(sprintf("prev pi3: %s",pi3))

	M <- rep(0,N)
	K <- rep(0,N)
	sig2vec <- rep(0,N)

	Rinv_lst = mclapply(1:N, function(i) Rinverse(l,xtr_vec_tr[i,][r_v_tr[i,]]), mc.cores=num_cores)
	M = unlist(mclapply(1:N,function (i) yhat(l,xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
	sig2vec = unlist(mclapply(1:N, function(i) sig2(l,Ygp[i,-t][r_v_tr[i,]],xtr_vec_tr[i,][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
	K = unlist(mclapply(1:N, function(i) s2(l,sig2vec[i],xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))

	loglik = sum(unlist(mclapply(1:N, function(i) L(pi1,pi2,pi3,w1[i],w2[i],w3[i],S[i],U1,U2,U3,S1,S2,S3,Z[i,],lr_beta1,lr_sigma1,Yreg[i,],lr_beta2,lr_sigma2,Ygp[i,-t],M[i],K[i]), mc.cores=num_cores)))
	prev_loglik <- loglik + tolerance + 1
	param$loglik = loglik

	ll = l
	j <- 1

	wws = get_ww(N,t,Ygp,Z,Yreg,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
	ww1 = wws$w1; ww2 = wws$w2; ww3 = wws$w3

	print(sprintf("loglik: %s",loglik))
	print(sprintf("num predictions: %s",length(S)))
	print(sprintf("abs error: %s",sum(abs(S - (ww1*(Z%*%lr_beta1)+ww2*(Yreg%*%lr_beta2)+ww3*M)))))
	GP_pred_error = sum(abs(S - M))
	print(sprintf("GP abs error: %s",GP_pred_error))

	while ((abs(loglik-prev_loglik) > tolerance && j < em_max_iter)) {
		print(sprintf("iter: %d",j))
		# E-step
		# 0.4s
		# ptm <- proc.time()
		ws = get_w(N,t,S,Ygp,r_v_tr,Z,Yreg,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,M,K,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
		w1 = ws$w1; w2 = ws$w2; w3 = ws$w3

		# M-step
		# estimate pi
		pi1 <- (sum(w1)) / N
		pi2 <- (sum(w2)) / N
		pi3 <- (sum(w3)) / N

		print(sprintf("pi1: %s",pi1))
		print(sprintf("pi2: %s",pi2))
		print(sprintf("pi3: %s",pi3))

		if (pi1 == 0 || round(pi1,8) == 0) {
			return(param)
		}
		if (pi2 == 0 || round(pi2,8) == 0) {
			return(param)
		}
		if (pi3 == 0 || round(pi3,8) == 0) {
			return(param)
		}

		# estimate U, S
		U1 = apply(w1*Z,2,sum)/sum(w1)
		U2 = apply(w2*Yreg,2,sum)/sum(w2)
		U3 = apply(w3*Ygp[,-t],2,sum)/sum(w3)

		S1 = Reduce('+',lapply(1:N,function(ri) {w1[ri]*(Z[ri,]-U1)%*%t(Z[ri,]-U1)})) / sum(w1)
		S2 = Reduce('+',lapply(1:N,function(ri) {w2[ri]*(Yreg[ri,]-U2)%*%t(Yreg[ri,]-U2)})) / sum(w2)
		S3 = Reduce('+',lapply(1:N,function(ri) {w3[ri]*(Ygp[ri,-t]-U3)%*%t(Ygp[ri,-t]-U3)})) / sum(w3)

		# estimate regression parameters
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

		ll <- Adam_one_obs_only(w3,pi3,U3,S3,Ygp,S,xte_vec_tr,xtr_vec_tr,Rinv_lst,t,r_v_tr,ll,step,gd_precision,gd_miter,j)

		Rinv_lst = mclapply(1:N, function(i) Rinverse(ll,xtr_vec_tr[i,][r_v_tr[i,]]), mc.cores=num_cores)
		M = unlist(mclapply(1:N,function (i) yhat(ll,xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
		sig2vec = unlist(mclapply(1:N, function(i) sig2(ll,Ygp[i,-t][r_v_tr[i,]],xtr_vec_tr[i,][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))
		K = unlist(mclapply(1:N, function(i) s2(ll,sig2vec[i],xte_vec_tr[i],xtr_vec_tr[i,][r_v_tr[i,]],Ygp[i,-t][r_v_tr[i,]],Rinv=Rinv_lst[[i]]), mc.cores=num_cores))

		wws = get_ww(N,t,Ygp,Z,Yreg,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
		ww1 = wws$w1; ww2 = wws$w2; ww3 = wws$w3

		prev_loglik <- loglik
		loglik = sum(unlist(mclapply(1:N, function(i) L(pi1,pi2,pi3,w1[i],w2[i],w3[i],S[i],U1,U2,U3,S1,S2,S3,Z[i,],lr_beta1,lr_sigma1,Yreg[i,],lr_beta2,lr_sigma2,Ygp[i,-t],M[i],K[i]), mc.cores=num_cores)))
		print(sprintf("loglik: %s",loglik))
		
		pred_error = sum(abs(S - (ww1*(Z%*%lr_beta1)+ww2*(Yreg%*%lr_beta2)+ww3*M)))
		print(sprintf("abs error: %s",pred_error))

		pi_pred_error = sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Yreg%*%lr_beta2)+pi3*M)))
		print(sprintf("abs error using pi: %s",pi_pred_error))

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
			param$U1 <- U1
			param$U2 <- U2
			param$U3 <- U3
			param$S1 <- S1
			param$S2 <- S2
			param$S3 <- S3
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

em_double_reg <- function(S,Z,Y,T,t,w1,w2,pi1,pi2,U1,U2,S1,S2,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,em_max_iter,tolerance) {
	print("em_double_reg")

	epsilon = 1e-8

	N = length(S)
	# initialization
	param = list(lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,pi1,pi2,w1,w2,-Inf,Inf,2)
	names(param) = c('lr_beta1','lr_sigma1','lr_beta2','lr_sigma2','pi1','pi2','w1','w2','loglik','abs_error','mix_model_num')

	print(sprintf("prev pi1: %s",pi1))
	print(sprintf("prev pi2: %s",pi2))

	# 1s
	loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],U1,U2,S1,S2,S[i],Z[i,],lr_beta1,lr_sigma1,Y[i,],lr_beta2,lr_sigma2), mc.cores=num_cores)))
	prev_loglik = loglik + tolerance + 1
	param$loglik = loglik

	j = 1

	get_ww_rr(N,Z,Y,pi1,pi2,U1,U2,S1,S2)
	ww1 = wws$w1; ww2 = wws$w2

	print(sprintf("loglik: %s",loglik))
	print(sprintf("num predictions: %s",length(S)))
	print(sprintf("abs error: %s",sum(abs(S - (ww1*(Z%*%lr_beta1)+ww2*(Y%*%lr_beta2))))))

	while ((abs(loglik-prev_loglik) > tolerance && j < em_max_iter)) {
		print(sprintf("iter: %d",j))
		# E-step
		# 0.4s
		# ptm <- proc.time()
		ws = get_w_rr(N,S,Z,Y,lr_beta1,lr_sigma1,lr_beta2,lr_sigma2,pi1,pi2,U1,U2,S1,S2)
		w1 = ws$w1; w2 = ws$w2

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

		# estimate U, S
		U1 = apply(w1*Z,2,sum)/sum(w1)
		U2 = apply(w2*Y,2,sum)/sum(w2)

		S1 = Reduce('+',lapply(1:N,function(ri) {w1[ri]*(Z[ri,]-U1)%*%t(Z[ri,]-U1)})) / sum(w1)
		S2 = Reduce('+',lapply(1:N,function(ri) {w2[ri]*(Y[ri,]-U2)%*%t(Y[ri,]-U2)})) / sum(w2)
		
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
		
		wws = get_ww_rr(N,Z,Y,pi1,pi2,U1,U2,S1,S2)
        ww1 = wws$w1; ww2 = wws$w2

		prev_loglik <- loglik
		# print(prev_loglik)
		# loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],S[i],Z[i,],Y[i,],lr_beta1,lr_sigma1,lr_beta2,lr_sigma2), mc.cores=num_cores)))
		loglik = sum(unlist(mclapply(1:N, function(i) L_rr(pi1,pi2,w1[i],w2[i],U1,U2,S1,S2,S[i],Z[i,],lr_beta1,lr_sigma1,Y[i,],lr_beta2,lr_sigma2), mc.cores=num_cores)))
		print(sprintf("loglik: %s",loglik))
		pred_error = sum(abs(S - (ww1*(Z%*%lr_beta1)+ww2*(Y%*%lr_beta2))))
		print(sprintf("abs error: %s",pred_error))
		pi_pred_error = sum(abs(S - (pi1*(Z%*%lr_beta1)+pi2*(Y%*%lr_beta2))))
		print(sprintf("abs error using pi: %s",pi_pred_error))

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
