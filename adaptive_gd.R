source('GP.R')
source('util.R')
source('mixtureMITemporalConfig.R')
library(parallel)

L_dl <- function(wgp,l,yi,s,xte,xtr,Rinv,nug_thres=20) {
	n = length(xtr)

	if (n == 0 || length(yi) == 0 || length(unique(yi)) == 1) {
		return(0)
	}
	else {

		m1n = matrix(rep(1,n),n,1)

		ITRI = (t(m1n)%*%Rinv%*%m1n)[1,]

		R_dl = cov_func_dl(l,xtr, xtr)

		Rinv_dl = -1 * Rinv %*% R_dl %*% Rinv

		ITRI_dl = t(m1n) %*% Rinv_dl %*% m1n

		rT = cov_func(l,xte,xtr)

		r = t(rT)

		rT_dl = cov_func_dl(l,xte,xtr)

		r_dl = cov_func_dl(l,xtr,xte)

		rT_Rinv_I = rT%*%Rinv%*%m1n

		rT_Rinv_I_dl = (rT_dl %*% Rinv + rT %*% Rinv_dl) %*% m1n

		fc = 1 - rT_Rinv_I

		gc = ITRI

		fc_dl = -1 * rT_Rinv_I_dl

		gc_dl = ITRI_dl

		part1 = (fc / gc) %*% t(m1n) + rT

		part1_dl = ((fc_dl * gc - gc_dl * fc) / gc^2) %*% t(m1n) + rT_dl

		yhat_dl = (part1_dl %*% Rinv + part1 %*% Rinv_dl) %*% yi

		rT_Rinv_r = rT%*%Rinv%*%r
		rT_Rinv_r_dl = rT_dl %*% Rinv %*% r + rT %*% Rinv_dl %*% r + rT %*% Rinv %*% r_dl

		fc1 = (1-t(m1n)%*%Rinv%*%r)^2

		fc1_dl = 2 * (1-t(m1n)%*%Rinv%*%r) * (-1 * t(m1n) %*% (Rinv_dl %*% r + Rinv %*% r_dl))

		gc1 = ITRI

		gc1_dl = ITRI_dl

		part2 = fc1 / gc1

		part2_dl = (fc1_dl * gc1 - gc1_dl * fc1) / gc1^2

		part3 = 1 - rT_Rinv_r + part2

		part3_dl = -1 * rT_Rinv_r_dl + part2_dl

		ITRY = (t(m1n)%*%Rinv%*%yi)

		ITRY_dl = t(m1n) %*% Rinv_dl %*% yi

		part4 = yi - m1n %*% (ITRY / ITRI)

		part4_dl = -1 * m1n %*% ((ITRY_dl * ITRI - ITRI_dl * ITRY) / ITRI^2)

		sigma2 = 1/n * (t(part4) %*% Rinv %*% part4)

		sigma2_dl = 1/n * (t(part4_dl) %*% Rinv %*% part4 + t(part4) %*% Rinv_dl %*% part4 + t(part4) %*% Rinv %*% part4_dl)

		s2 = sigma2 * part3

		s2_dl = sigma2 * part3_dl + sigma2_dl * part3

		f.T <- s2

		f_dl.T <- s2_dl

		g.T <- simple_yhat(l, xtr, yi, xte)

		g_dl.T <- yhat_dl

		Ldl = (wgp*(-1/(2*f.T) * f_dl.T + (s-g.T)*g_dl.T/f.T + f_dl.T*(s-g.T)^2/(2*f.T^2)))

	}

	return(Ldl)
}

L_GP <- function(wgp,pi_gp,U3,S3,X,l,Yobs,Y,s,xte,xtr,Rinv,epsilon=1e-8) {		
	if (length(Yobs) == 0 || length(unique(Yobs)) == 1) {
		L2 = pi_gp
		loglik = wgp * log_with_limits(L2,epsilon)
	}
	else {
		# gp_sigma2 <- sig2(l,Yobs,xtr,Rinv=Rinv)
		# m <- yhat(l,xte,xtr,Yobs,Rinv=Rinv)
		# k <- s2(l,gp_sigma2,xte,xtr,Yobs,Rinv=Rinv)
		gp_pred = simple_GP_pred(l,xtr,Yobs,xte)
		m = gp_pred$yhat
		k = gp_pred$mse
		# loglik <- loglik + wgp[i]*log(pi_gp*dnorm(S[i],m,sqrt(k))/wgp[i])
		p2 = dnorm(s,m,sqrt(k)) * dmvnorm(X,U3,S3)
		# p2 = round(p2,digits=8)
		if (wgp == 0) {
			loglik = log_with_limits(pi_gp*p2,epsilon)
		}
		else {
			L2 = pi_gp * p2 / wgp
			loglik = wgp * log_with_limits(L2,epsilon)
		}
	}
	return(loglik)
}

simple_L <- function(wgp,pi_gp,l,xtr,ytr,xte,epsilon=1e-8) {
	gp_pred = simple_GP_pred(l,xtr,ytr,xte)
	m = gp_pred$yhat
	k = gp_pred$mse
	# print(m)
	# print(k)
	p2 = dnorm(S[x],m,sqrt(k)) * dmvnorm(X[x,],U3,S3)
	# print(dmvnorm(X[x,],U3,S3))
	# print(dnorm(S[x],m,sqrt(k)))
	# print(p2)
	L2 = pi_gp * p2 / wgp
	# print(L2)
	# loglik = wgp * log_with_limits(L2,epsilon)
	loglik = wgp * log(L2)
	return(loglik)	
}

dff <- function(wgp,pi_gp,l,s,xtr,ytr,xte,epsilon=1e-8) {
	Rinv = solve(cov_func(l,xtr,xtr))
	gx = L_dl(wgp,l,ytr,s,xte,xtr,Rinv)

	gx0 = simple_L(wgp,pi_gp,l+1e-4,xtr,ytr,xte)
	gx1 = simple_L(wgp,pi_gp,l-1e-4,xtr,ytr,xte)

	rgx = (gx0 - gx1)/(2*1e-4)

	return(gx-rgx)
}

Adam_one_obs_only <- function(wgp,pi_gp,U3,S3,X,Y,S,xte_vec,xtr_vec,Rinv_lst,xstar,r_v,l,tao,e,miter,emiter) {
	N = length(S)
	beta1 = 0.9
	beta2 = 0.999
	m = 0
	v = 0
	prev_m = m
	prev_v = v
	epsilon = 1e-8

	loglik_window = rep(0,5)

	res = l
	max_loglik = -Inf

	for (i in 1:miter) {
		print(i)
		prev_l <- l

		g <- sum(unlist(mclapply(1:N, function(x) L_dl(wgp[x],prev_l,Y[x,-xstar][r_v[x,]],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv_lst[[x]]), mc.cores=num_cores)))
		
		lst1 = unlist(mclapply(1:N, function(x) {
			y = Y[x,-xstar][r_v[x,]]
			
			Rinv = solve(cov_func(prev_l,xtr_vec[x,][r_v[x,]],xtr_vec[x,][r_v[x,]]))
			gx = L_dl(wgp[x],prev_l,Y[x,-xstar][r_v[x,]],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv)
			
			gx0 = simple_L(wgp[x],pi_gp,prev_l+1e-4,xtr_vec[x,][r_v[x,]],y,xte_vec[x])
			gx1 = simple_L(wgp[x],pi_gp,prev_l-1e-4,xtr_vec[x,][r_v[x,]],y,xte_vec[x])

			rgx = (gx0 - gx1)/(2*1e-4)
			
			diff = gx-rgx
			diff
			
		},mc.cores=20))


		lst2 = sapply(1:N, function(x) {
			Rinv = solve(cov_func(prev_l,xtr_vec[x,][r_v[x,]],xtr_vec[x,][r_v[x,]]))
			gx = L_dl(wgp[x],prev_l,Y[x,-xstar][r_v[x,]],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv)
			gx0 = simple_L(wgp[x],pi_gp,prev_l+1e-4,xtr_vec[x,][r_v[x,]],Y[x,-xstar][r_v[x,]],xte_vec[x])
			gx1 = simple_L(wgp[x],pi_gp,prev_l-1e-4,xtr_vec[x,][r_v[x,]],Y[x,-xstar][r_v[x,]],xte_vec[x])

			rgx = (gx0 - gx1)/(2*1e-4)
			
			diff = gx-rgx
			diff
		})

		lst3 = sapply(1:N, function(x) dff(wgp[x],pi_gp,prev_l,S[x],xtr_vec[x,][r_v[x,]],Y[x,-xstar][r_v[x,]],xte_vec[x]))

		lst4 = rep(0,N)
		for (x in 1:N) {
			lst4[x] = dff(wgp[x],pi_gp,prev_l,S[x],xtr_vec[x,][r_v[x,]],Y[x,-xstar][r_v[x,]],xte_vec[x])
		}

		lst5 = rep(0,N)
		for (x in 1:5) {
			
			y = Y[x,-xstar][r_v[x,]]
			
			if (length(y) >= 3 && length(unique(y)) != 1) {
				Rinv = solve(cov_func(prev_l,xtr_vec[x,][r_v[x,]],xtr_vec[x,][r_v[x,]]))
				gx = L_dl(wgp[x],prev_l,Y[x,-xstar][r_v[x,]],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv)
				
				gx0 = simple_L(wgp[x],pi_gp,prev_l+1e-4,xtr_vec[x,][r_v[x,]],y,xte_vec[x])
				gx1 = simple_L(wgp[x],pi_gp,prev_l-1e-4,xtr_vec[x,][r_v[x,]],y,xte_vec[x])

				rgx = (gx0 - gx1)/(2*1e-4)
			
			}
			diff = gx-rgx
			lst5[x] = diff
		}

		print(lst5[1:5])
		print(lst4[1:5])
		print(lst3[1:5])
		print(lst2[1:5])
		print(lst1[1:5])

		print(sum(lst5-lst4))
		print(sum(lst5-lst3))
		print(sum(lst5-lst2))
		print(sum(lst5-lst1))

		# for (x in 1:N) {
		# 	L_dl(wgp[x],prev_l,Y[x,-xstar][r_v[x,]],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv_lst[[x]])
		# }

		m = beta1 * prev_m + (1-beta1) * g
		v = beta2 * prev_v + (1-beta2) * g^2

		m_hat = m / (1-beta1^i)
		v_hat = v / (1-beta2^i)

		l <- prev_l + m_hat * tao / (sqrt(v_hat) + epsilon)

		# ptm <- proc.time()
		Rinv_lst = mclapply(1:N, function(x) Rinverse(l,xtr_vec[x,][r_v[x,]]), mc.cores=num_cores)

		loglik = sum(unlist(mclapply(1:N, function(x) L_GP(wgp[x],pi_gp,U3,S3,X[x,],l,Y[x,-xstar][r_v[x,]],Y[x,-xstar],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv_lst[[x]]), mc.cores=num_cores)))
		# loglik = 0
		# for (x in 1:N) {
		# 	loglik = loglik + L_GP(wgp[x],pi_gp,U3,S3,X[x,],l,Y[x,-xstar][r_v[x,]],Y[x,-xstar],S[x],xte_vec[x],xtr_vec[x,][r_v[x,]],Rinv_lst[[x]])
		# }

		print(round(loglik,1))

		if (loglik > max_loglik) {
			res = l
			max_loglik = loglik
		}

		if (i <= 5) {
			loglik_window[i] = loglik
		}
		else {
			loglik_window = c(loglik_window,loglik)[-1]
		}

		if (i >= 5 && (abs(max(loglik_window) - min(loglik_window)) < e
			 || abs(sum(loglik_window[-1]-loglik_window[-5])) < e)) {
			break
		}
	}

	return(res)
}