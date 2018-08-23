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
		# temp = eigen(R,symmetric = TRUE, only.values = TRUE);
		# eig_val = temp$values;
		# condnum = kappa(R,triangular = TRUE,exact=TRUE);
		# max_eigval = eig_val[1];
		# delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*(exp(nug_thres)-1))));

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
		## Adding in the itterative approach section
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
			## if delta != 0, start itterations
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
			gpmod = GP_fit(X,Y)
			gpmod$beta=beta
			gpmod$sig2=sig2
			mse = predict.GP(gpmod, xnew)$MSE

			print(sprinf("mse: %s",res$mse))
			print(sprinf("mse from GP_fit: %s",mse))
			# print("mse <= 0. set to 1e-8")
			# res$mse = 1e-8
			res$mse = mse
		}
	}
	return(res)
}

Rinverse <- function(l, gpx, nug_thres=20) {
	n = length(gpx)
	if (n == 0) {
		return(NULL)
	}

	# print(gpx)

	R = cov_func(l,gpx,gpx)

	# print(R)
	
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

# sig2 <- function(l, Y, gpx, Rinv=NULL, nug_thres=20) {
# 	if (length(Y) == 0 || length(unique(Y)) == 1) {
# 		return(NA)
# 	}

# 	n = length(gpx)

# 	if (is.null(Rinv)) {
# 		Rinv = Rinverse(l,gpx)
# 	}

# 	m1n = matrix(rep(1,length(gpx)),length(gpx),1)
# 	ITRI = (t(m1n)%*%Rinv%*%m1n)[1,]
# 	ITRY = (t(m1n)%*%Rinv%*%Y)[1,]
	
# 	part4 = Y - m1n %*% (ITRY / ITRI)
# 	sigma2 = 1/n * (t(part4) %*% Rinv %*% part4)

# 	return(sigma2)
# }

# s2 <- function(l,gpsig2,i,i.bar,Y,Rinv=NULL,nug_thres=20) {
# 	if (is.na(gpsig2)) {
# 		return(NA)
# 	}

# 	ki <- cov_func(l,i,i.bar)
	
# 	if (is.null(Rinv)) {
# 		Rinv = Rinverse(l,i.bar)
# 	}
# 	kiT <- t(ki)
# 	m1n <- matrix(rep(1,length(i.bar)),length(i.bar),1)

# 	mse = gpsig2 * (1 - ki%*%Rinv%*%kiT + (1 - t(m1n)%*%Rinv%*%kiT)^2/(t(m1n)%*%Rinv%*%m1n))

# 	# to change
# 	if (mse <= 0) {
# 		print(mse)
# 		print(l)
# 		print(i.bar)
# 		print(Y)
# 		print(i)
# 		print(gpsig2)
# 		print(Rinv)
# 		gpmod = GP_fit(i.bar,Y)
# 		gpmod$beta=l
# 		gpmod$sig2=gpsig2
# 		mse = predict.GP(gpmod, i)$MSE
# 	}

# 	return(mse)
# }

# yhat <- function(l,i,i.bar,Yi_1,Rinv=NULL) {
# 	if (length(Yi_1) == 0) {
# 		return(Inf)
# 	}
# 	if (length(unique(Yi_1)) == 1) {
# 		return(unique(Yi_1))
# 	}

# 	if (is.null(Rinv)) {
# 		Rinv = Rinverse(l,i.bar)
# 	}
	
# 	ki <- cov_func(l,i,i.bar)

# 	m1n <- matrix(rep(1,length(i.bar)),length(i.bar),1)

# 	return((((1-ki%*%Rinv%*%m1n) / (t(m1n)%*%Rinv%*%m1n)) %*% t(m1n) + ki) %*% Rinv %*% Yi_1)
# }

fit_gp <- function(xtr,ytr) {
	if (length(unique(ytr)) == 1) {
		ytr[1] = ytr[1] + 0.000001
		print("unique")
    }
	
	if (length(ytr) <= 1) {
		# gpmod = list(-3)
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
	# else if (length(ytr) == 2) {
 #        res$pred = approx(xtr,ytr,xte,rule=2)$y
 #        res$mse = 0
 #        res$idx = TRUE
 #    }
	else {
		# res$pred = yhat(l,xte,xtr,ytr,Rinv=Rinv)
		# res$mse = s2(l,sig2(l, ytr, xtr, Rinv=Rinv),xte,xtr,ytr,Rinv=Rinv)
		gp_pred = simple_GP_pred(l,xtr,ytr,xte)
		res$pred = gp_pred$yhat
		res$mse = gp_pred$mse
		res$idx = FALSE
	}
	return(res)
}
