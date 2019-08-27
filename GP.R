library(GPfit)
library(parallel)

cov_func <- function(l, p.vec, q.vec, power=1.95) {
	lenp <- length(p.vec)
	lenq <- length(q.vec)
	M <- matrix(rep(0,lenp*lenq),lenp,lenq)
	for (i in 1:lenp) {
		for (j in 1:lenq) {
			M[i,j] <- exp(-10^l*(abs(p.vec[i]-q.vec[j]))^power)
		}
	}
	return(M)
}

cov_func_dl <- function(l, p.vec, q.vec, power=1.95) {
	lenp <- length(p.vec)
	lenq <- length(q.vec)
	M <- matrix(rep(0,lenp*lenq),lenp,lenq)
	for (i in 1:lenp) {
		for (j in 1:lenq) {
			M[i,j] <- -1 * exp(-10^l*(abs(p.vec[i]-q.vec[j]))^power) * (abs(p.vec[i]-q.vec[j]))^power * log(10) * 10^l
		}
	}
	return(M)
}

simple_yhat <- function(beta, X, Y, xnew, nug_thres=20, power=1.95, M=1) {
	if (is.matrix(X) == FALSE){
		X = as.matrix(X)
	}
	if (is.matrix(xnew) == FALSE){
		xnew = as.matrix(xnew)
	}

	n = nrow(X)

	One = rep(1,n);
	R = cov_func(beta,X,X);
	temp = eigen(R,symmetric = TRUE, only.values = TRUE);
	eig_val = temp$values;
	condnum = kappa(R,triangular = TRUE,exact=TRUE);
	max_eigval = eig_val[1];
	delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*(exp(nug_thres)-1))));

	LO = diag(n);
	Sig = R + delta*LO;

	L = chol(Sig);

	Sig_invOne = solve(L,solve(t(L),One));	
	Sig_invY = solve(L,solve(t(L),Y));

	xn = matrix(xnew,nrow=1);
	r = exp(-(abs(X-as.matrix(rep(1,n))%*%(xn))^power)%*%(10^beta));
	yhat = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invY;

	return(yhat)
}

simple_GP_pred <- function(beta, X, Y, xnew, nug_thres=20, power=1.95, M=1) {
	# Adapted code from GPfit[1].
	# [1] MacDonald, B., Ranjan, P., & Chipman, H. (2013).
	# GPfit: an R package for Gaussian process model fitting 
	# using a new optimization algorithm. arXiv preprint arXiv:1305.0759.

	res = list(0,0)
	names(res) = c("yhat","mse")
	if (length(Y) == 0 || length(unique(Y)) == 1) {
		res$yhat = unique(Y)
		res$mse = NA
	} else {

		if (is.matrix(X) == FALSE){
			X = as.matrix(X)
		}
		if (is.matrix(xnew) == FALSE){
			xnew = as.matrix(xnew)
		}

		n = nrow(X)

		One = rep(1,n);
		R = cov_func(beta,X,X);
		delta = 0

		LO = diag(n);
		Sig = R + delta*LO;

		L = chol(Sig);

		Sig_invOne = solve(L,solve(t(L),One));	
		Sig_invY = solve(L,solve(t(L),Y));	

		mu_hat = solve(t(One)%*%Sig_invOne,t(One)%*%Sig_invY);
		Sig_invb = solve(L,solve(t(L),(Y-One%*%mu_hat)));
		sig2 = t(Y-One%*%mu_hat)%*%Sig_invb/n;

		if (delta == 0){
			Sig_invLp = solve(L,solve(t(L),LO));
		} else {
			s_Onei = One;
			s_Yi = Y;
			s_Li = LO;
			Sig_invOne = matrix(0, ncol = 1, nrow = n);
			Sig_invY = matrix(0, ncol = 1, nrow = n);
			Sig_invLp = matrix(0, ncol = n, nrow = n)
			for (it in 1:M)
			{
				s_Onei = solve(L,solve(t(L),delta*s_Onei));
				Sig_invOne = Sig_invOne + s_Onei/delta;

				s_Yi = solve(L,solve(t(L),delta*s_Yi));
				Sig_invY = Sig_invY + s_Yi/delta;

				s_Li = solve(L,solve(t(L),delta*s_Li));
				Sig_invLp = Sig_invLp + s_Li/delta;
			}
		}

		xn = matrix(xnew,nrow=1);
		r = exp(-(abs(X-as.matrix(rep(1,n))%*%(xn))^power)%*%(10^beta));
		res$yhat = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invY;


		if (delta == 0) {
			Sig_invr = solve(L,solve(t(L),r));
		} else {
			s_ri = r;
			Sig_invr = matrix(0, ncol = 1, nrow = n);
			for (it in 1:M)
			{
				s_ri = solve(L,solve(t(L),delta*s_ri));
				Sig_invr = Sig_invr + s_ri/delta;
			}
		}
		cp_delta_r = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invr;

		cp_delta_Lp = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invLp;
		mse = sig2*(1-2*cp_delta_r+cp_delta_Lp%*%R%*%t(cp_delta_Lp));

		res$mse = mse*(mse>0)

		if (res$mse <= 0) {
			print("mse == 0. set to 1e-8")
			res$mse = 1e-8
			
		}
	}
	return(res)
}

Rinverse <- function(l, gpx, nug_thres=20) {
	n = length(gpx)
	if (n == 0) {
		return(NULL)
	}

	R = cov_func(l,gpx,gpx)

	temp = eigen(R,symmetric = TRUE, only.values = TRUE);
	eig_val = temp$values;

	condnum = kappa(R,triangular = TRUE,exact=TRUE);
	max_eigval = eig_val[1]
	delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*
			(exp(nug_thres)-1))))

	LO = diag(n);
	Sig = R + delta*LO;

	L = chol(Sig)

	Rinv <- solve(L, solve(t(L)))

	return(Rinv)
}

fit_gp <- function(xtr,ytr) {
	if (length(unique(ytr)) == 1) {
		ytr[1] = ytr[1] + 0.000001
		print("unique")
    }
	
	if (length(ytr) <= 1) {
		gpmod = list(NA)
		names(gpmod) = c("beta")
	}
	else {
		gpmod = GP_fit(xtr,ytr)
	}
	
	return(gpmod)
}

gp_predict_one_rt <- function(l, xtr, ytr, xte, Rinv=NULL) {
	res = list(0,0,FALSE)
	names(res) = c("pred","mse","idx")
	
	if (length(unique(ytr)) == 1) {
		res$pred = ytr[1]
		res$mse = 0
		res$idx = TRUE
	}
	else {
		gp_pred = simple_GP_pred(l,xtr,ytr,xte)
		res$pred = gp_pred$yhat
		res$mse = gp_pred$mse
		res$idx = FALSE
	}
	return(res)
}
