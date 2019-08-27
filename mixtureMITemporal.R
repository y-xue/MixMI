source('em.R')
source('regression.R')
source('GP.R')
library(parallel)

initial_impute_sample <- function(y, ry, x=NULL, ...)
{   
    return(sample(y[ry], size=sum(!ry), replace=TRUE))
}

initialize_imp <- function(data, m, visit_col_sequence, r, prev_data=NA)
{
  imp <- vector("list", ncol(data))
  nmis <- apply(is.na(data),2,sum)

  if(m > 0) {
    
    for(j in visit_col_sequence) {
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(!r[,j]), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[r[,j] == FALSE], 1:m)
      y <- data[,j]
      ry <- r[,j]
      for(i in 1:m) {
        if (nmis[j]<nrow(data)) {
          if (is.data.frame(prev_data)) {
            imp[[j]][,i] <- prev_data[[j]][!ry]
          }
          else {
            imp[[j]][,i] <- initial_impute_sample(y, ry, x)
          }
        }
        else imp[[j]][,i] <- rnorm(nrow(data))
      }
    }
  }

  return(imp)
}

initialize_imp_tensor <- function(pv_tensor,m) {
  imp_tensor <- vector('list', length(pv_tensor))
  for (i in 1:length(pv_tensor)) {
    visit_col_sequence <- (1:ncol(pv_tensor[[i]]))[apply(is.na(pv_tensor[[i]]),2,any)]
    r <- (!is.na(pv_tensor[[i]]))
    imp_tensor[[i]] <- initialize_imp(pv_tensor[[i]],m,visit_col_sequence,r)
  }
  return(imp_tensor)
}

sampler <- function(pv_tensor, prt_m, artificial_prt_tensor, ori_tensor, model_type, out_cdn, gpmodel_dir, m, maxit, obs_only, imp_tensor, r_list, r_vlist, predictor_matrix_list, visit_col_sequence_list, em_max_iter, tolerance, step, gd_miter, gd_precision, ridge, observation, printFlag, ...)
{
    print("sampler")
    print(out_cdn)
    print(model_type)
    print(sprintf("obs_only: %s", obs_only))

    num_time_point <- length(pv_tensor)
    num_pt <- dim(pv_tensor[[1]])[1]
    num_var <- dim(pv_tensor[[1]])[2]

    if (maxit < 1) {
        return(imp_tensor)
    }

    # 0 for GP in RRG, 1 for RRG, 2 for RR
    mix_model_num = 1

    pi3 = 0.1
    pi1 = (1-pi3)/2
    pi2 = (1-pi3)/2
    
    pi_1_m = array(pi1,dim=c(14,num_time_point,m))
    pi_2_m = array(pi2,dim=c(14,num_time_point,m))
    pi_3_m = array(pi3,dim=c(14,num_time_point,m))

    w1_lst = list()
    w2_lst = list()
    w3_lst = list()
    for (i in 1:m) {
        w1_lst[[i]] = list()
        w2_lst[[i]] = list()
        w3_lst[[i]] = list()
        for (v in 2:num_var) {
            # assume the first variable (time) is excluded
            for (t in 1:num_time_point) {
                w1_lst[[i]][[(v-1)*num_time_point+t]] = 0
                w2_lst[[i]][[(v-1)*num_time_point+t]] = 0
                w3_lst[[i]][[(v-1)*num_time_point+t]] = 0
            }
        }
    }

    mix_model_1_param = list(pi_1_m,pi_2_m,pi_3_m,w1_lst,w2_lst,w3_lst)
    names(mix_model_1_param) = c("pi_1_m","pi_2_m","pi_3_m","w1_lst","w2_lst","w3_lst")

    pi_1_m = array(0.5,dim=c(14,num_time_point,m))
    pi_2_m = array(0.5,dim=c(14,num_time_point,m))
    mix_model_2_param = list(pi_1_m,pi_2_m,w1_lst,w2_lst)
    names(mix_model_2_param) = c("pi_1_m","pi_2_m","w1_lst","w2_lst")

    rt = TRUE
    if (is.null(prt_m)) {
        rt = FALSE
        print("time index")
        prt_m = matrix(NA, num_pt, num_time_point)
        for (p in 1:num_pt) {
            prt_m[p,] = c(1:num_time_point)
        }
    }
    prt_m_norm = t(apply(prt_m,1,function(row) (row-min(row))/(max(row)-min(row))))

    for (k in 1:maxit) {
        print(sprintf("iter: %s",k))
        for (i in 1:m) {
            print(sprintf("imp: %s",i))

            # fill the data with the last set of imputations
            for (t in 1:num_time_point) {
                r <- r_list[[t]]
                visit_col_sequence <- visit_col_sequence_list[[t]]

                for (v in visit_col_sequence) {
                    pv_tensor[[t]][!r[,v],v] <- imp_tensor[[t]][[v]][,i]
                }
            }

            for (v in 2:num_var) {
                # assume the first variable (time) is excluded
                pt_list <- vector('list',num_time_point)
                for (t in 1:num_time_point) {
                    pt_list[[t]] <- pv_tensor[[t]][,v]
                }
                pt_df <- do.call(cbind,pt_list)

                for (t in 1:num_time_point) {
                    print(sprintf("v %s, t %s",v,t))

                    data <- pv_tensor[[t]]
                    r <- r_list[[t]]
                    visit_col_sequence <- visit_col_sequence_list[[t]]
                    predictor_matrix <- predictor_matrix_list[[t]]

                    x1 <- data[, predictor_matrix[v, ] == 1, drop = FALSE]
                    y <- data[, v]
                    ry <- r[, v]

                    ori_y = NULL
                    if (!is.null(ori_tensor)) {
                        ori_y = ori_tensor[[t]][,v]
                    }
                    
                    if (is.null(artificial_prt_tensor)) {
                        print("prt_m_norm")
                        x2 = cbind(pt_df[,-t],prt_m_norm[,-c(1,num_time_point)])
                        
                        xte_vec = prt_m_norm[,t]
                        xtr_vec = prt_m_norm[,-t]
                    } else {
                        print("artificial_prt_tensor")
                        x2 = cbind(pt_df[,-t],artificial_prt_tensor[[v]][,-c(1,num_time_point)])
                        
                        xte_vec = artificial_prt_tensor[[v]][,t]
                        xtr_vec = artificial_prt_tensor[[v]][,-t]
                    }

                    if (obs_only) {
                        r_v = r_vlist[[v]][,-t]
                    } else {
                        r_v = matrix(TRUE, num_pt, num_time_point-1)
                    }

                    w_dir = sprintf("%s/em_param/imp_%s_iter_%s",out_cdn,i,k)

                    dir.create(w_dir,recursive=TRUE)
                    w_fn = sprintf("%s/val%s_tp%s",w_dir,v,t)

                    if (model_type == "rrg" || model_type == "both") {
                        if (gpmodel_dir == "") {
                            gpmodel_dir = sprintf("%s/GP_models",out_cdn)
                            dir.create(gpmodel_dir,recursive=TRUE)
                        }

                        GPmodel_fn = sprintf("%s/val%s_tp%s",gpmodel_dir,v,t)
                        if (!file.exists(GPmodel_fn)) {
                            print("fitting GPmodel")
                            GPmodel_vec = mclapply(1:num_pt, function(pt) fit_gp(xtr_vec[pt,][r_v[pt,]],pt_df[pt,-t][r_v[pt,]]),mc.cores=num_cores)
                            dump("GPmodel_vec",GPmodel_fn)
                        } else {
                            GPmodel_vec = source(GPmodel_fn)$value
                            print("loaded GPmodel_vec")
                        }

                        lvec = sapply(GPmodel_vec, function(x) { if (is.list(x)) {x$beta} else {NA}})
                        lvec = lvec[!is.na(lvec)]
                        l = median(lvec)
                        print(l)
                    } else {
                        l = 0
                    }
                    

                    if (sum(!ry) > 0) {
                        imp_res <- impute_em_rrg_obs_only(model_type,i,num_time_point,v,y,ry,x1,x2,pt_df,ori_y,xtr_vec,xte_vec,t,r_v,mix_model_num,mix_model_1_param,mix_model_2_param,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn,ridge,observation)

                        if (model_type == "rrg" || model_type == "both") {
                            mix_model_1_param$pi_1_m[v,t,i] = (imp_res$rrg_em_param)$pi1
                            mix_model_1_param$pi_2_m[v,t,i] = (imp_res$rrg_em_param)$pi2
                            mix_model_1_param$pi_3_m[v,t,i] = (imp_res$rrg_em_param)$pi3
                            mix_model_1_param$w1_lst[[i]][[(v-1)*num_time_point+t]] = (imp_res$rrg_em_param)$w1
                            mix_model_1_param$w2_lst[[i]][[(v-1)*num_time_point+t]] = (imp_res$rrg_em_param)$w2
                            mix_model_1_param$w3_lst[[i]][[(v-1)*num_time_point+t]] = (imp_res$rrg_em_param)$w3
                        }

                        if (model_type == "rr" || model_type == "both") {
                            mix_model_2_param$pi_1_m[v,t,i] = (imp_res$rr_em_param)$pi1
                            mix_model_2_param$pi_2_m[v,t,i] = (imp_res$rr_em_param)$pi2
                            mix_model_2_param$w1_lst[[i]][[(v-1)*num_time_point+t]] = (imp_res$rr_em_param)$w1
                            mix_model_2_param$w2_lst[[i]][[(v-1)*num_time_point+t]] = (imp_res$rr_em_param)$w2
                        }

                        imp_tensor[[t]][[v]][,i] = imp_res$pred
                        
                    }
                    
                    pv_tensor[[t]][!r[,v],v] <- imp_tensor[[t]][[v]][,i]
                    pt_df[!r[,v],t] <- imp_tensor[[t]][[v]][,i]

                }
            }

            for (t in 1:length(pv_tensor)) {
                imputed <- pv_tensor[[t]]
                for (v in visit_col_sequence_list[[t]]) {    
                    imputed[!r_list[[t]][,v],v] <- imp_tensor[[t]][[v]][,i]
                }
                write.csv(imputed, file=sprintf("%s/imp_%d_attime_%s_miceiter_%d.csv",out_cdn,i,t,k), row.names=FALSE)
            }
        }
    }
    return(imp_tensor)
}

impute_em_rrg_obs_only <- function(model_type,impi,num_time_point,v,y,ry,x1,x2,pt_df,ori_y,xtr_vec,xte_vec,t,r_v,mix_model_num,mix_model_1_param,mix_model_2_param,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn,ridge=1e-5,observation=FALSE)
{
    print("impute_em_rrg_obs_only")

    rr_param = 0
    rrg_param = 0
    res <- list(0,rr_param,rrg_param)
    names(res) <- c("pred","rr_em_param","rrg_em_param")

    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)

    if (!is.null(ori_y)) {
        ori_y_te = ori_y[!ry]
        test_r = !is.na(ori_y_te)
    }
    Nstar = sum(!ry)

    if (!file.exists(sprintf("%s_pred_error.txt",w_fn))) {
        # Train with selected records
        # that have at least one measurement
        # 
        sy = ry
        for (i in 1:length(sy)) {
            if (sy[i] == TRUE) {
                ts = pt_df[i,-t][r_v[i,]]
                if (sum(r_v[i,]) == 0) {
                    sy[i] = FALSE
                }
            }
        }

        # Train with all records
        # 
        # N = sum(ry)
        # S = y[ry]
        # Z = x1[ry,]
        # Yreg = x2[ry,]
        # Ygp = pt_df[ry,]

        # xtr_vec_tr = xtr_vec[ry,]
        # xte_vec_tr = xte_vec[ry]

        # r_v_tr = r_v[ry,]

        # lr_param1 <- norm_fix(y, ry, x1)
        # lr_param2 <- norm_fix(y, ry, x2)

        if (model_type == "rrg" || model_type == "both") {
            # Train rrg model
            print("training EM rrg")

            lr_param1 <- norm_fix(y, sy, x1)
            lr_param2 <- norm_fix(y, sy, x2)
            T = dim(pt_df)[2]
            N = sum(sy)
            S = y[sy]
            Z = x1[sy,]
            Yreg = x2[sy,]
            Ygp = pt_df[sy,]

            xtr_vec_tr = xtr_vec[sy,]
            xte_vec_tr = xte_vec[sy]

            r_v_tr = r_v[sy,]

            pi1 = mix_model_1_param$pi_1_m[v,t,impi]
            pi2 = mix_model_1_param$pi_2_m[v,t,impi]
            pi3 = mix_model_1_param$pi_3_m[v,t,impi]
            w1 = mix_model_1_param$w1_lst[[impi]][[(v-1)*num_time_point+t]]
            w2 = mix_model_1_param$w2_lst[[impi]][[(v-1)*num_time_point+t]]
            w3 = mix_model_1_param$w3_lst[[impi]][[(v-1)*num_time_point+t]]
            if (length(w1) != N) {
                w1 = rep(pi1,N)
                w2 = rep(pi2,N)
                w3 = rep(pi3,N)
            }

            X = cbind(Z,Yreg)
            U1 = apply(X,2,mean)
            U2 = U1; U3 = U1

            S1 = Reduce('+',lapply(split(X,1:nrow(X)),function(row) {(row-U1)%*%t(row-U1)})) / N
            S2 = S1; S3 = S1

            sink(sprintf("%s_rrg_em_params.txt",w_fn))

            rrg_param = em_rrg_obs_only(S,Z,Yreg,Ygp,xte_vec_tr,xtr_vec_tr,t,r_v_tr,mix_model_num,w1,w2,w3,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3,lr_param1$beta,lr_param1$sigma,lr_param2$beta,lr_param2$sigma,l,em_max_iter,tolerance,step,gd_miter,gd_precision,ridge)
            
            sink()

            # save rrg params
            pi1 = rrg_param$pi1
            pi2 = rrg_param$pi2
            pi3 = rrg_param$pi3
            w1 = rrg_param$w1
            w2 = rrg_param$w2
            w3 = rrg_param$w3
            U1 = rrg_param$U1
            U2 = rrg_param$U2
            U3 = rrg_param$U3
            S1 = rrg_param$S1
            S2 = rrg_param$S2
            S3 = rrg_param$S3
            lr_beta1 = rrg_param$lr_beta1
            lr_sigma1 = rrg_param$lr_sigma1
            lr_beta2 = rrg_param$lr_beta2
            lr_sigma2 = rrg_param$lr_sigma2
            ll = rrg_param$ll
            dump("ll",sprintf("%s.ll",w_fn))
            dump("pi1", sprintf("%s_rrg.pi1",w_fn))
            dump("pi2", sprintf("%s_rrg.pi2",w_fn))
            dump("pi3", sprintf("%s_rrg.pi3",w_fn))
            dump("w1", sprintf("%s_rrg.w1",w_fn))
            dump("w2", sprintf("%s_rrg.w2",w_fn))
            dump("w3", sprintf("%s_rrg.w3",w_fn))
            dump("U1", sprintf("%s_rrg.U1",w_fn))
            dump("U2", sprintf("%s_rrg.U2",w_fn))
            dump("U3", sprintf("%s_rrg.U3",w_fn))
            dump("S1", sprintf("%s_rrg.S1",w_fn))
            dump("S2", sprintf("%s_rrg.S2",w_fn))
            dump("S3", sprintf("%s_rrg.S3",w_fn))
            dump("lr_beta1", sprintf("%s_rrg.lr_beta1",w_fn))
            dump("lr_sigma1", sprintf("%s_rrg.lr_sigma1",w_fn))
            dump("lr_beta2", sprintf("%s_rrg.lr_beta2",w_fn))
            dump("lr_sigma2", sprintf("%s_rrg.lr_sigma2",w_fn)) 
            
            lr_prediction1 = x1[!ry,  ] %*% lr_beta1
            lr_prediction2 = x2[!ry,  ] %*% lr_beta2
            Ystar = pt_df[!ry,]
        
            xtr_vec_star = xtr_vec[!ry,]
            xte_vec_star = xte_vec[!ry]
            r_v_star = r_v[!ry,]

            if (Nstar == 1) {
                Ystar = t(as.matrix(Ystar))
            }

            GPprediction_res = mclapply(1:Nstar, function(i) gp_predict_one_rt(ll,xtr_vec_star[i,][r_v_star[i,]],Ystar[i,-t][r_v_star[i,]],xte_vec_star[i]), mc.cores=num_cores)
            gp_prediction = sapply(GPprediction_res, function(x) x$pred)
            
            wws = get_ww(Nstar,t,Ystar,x1[!ry,],x2[!ry,],pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
            ww1 = wws$w1; ww2 = wws$w2; ww3 = wws$w3

            dump("ww1", sprintf("%s_rrg.ww1",w_fn))
            dump("ww2", sprintf("%s_rrg.ww2",w_fn))
            dump("ww3", sprintf("%s_rrg.ww3",w_fn))

            if (rrg_param$mix_model_num == 1) {
                # RRG
                rrg_prediction = ww1 * lr_prediction1 + ww2 * lr_prediction2 + ww3 * gp_prediction
                rrg_pi_prediction = pi1 * lr_prediction1 + pi2 * lr_prediction2 + pi3 * gp_prediction
            } else {
                # GP
                rrg_prediction = gp_prediction
                rrg_pi_prediction = gp_prediction
            }
            rrg_rescale_rr_prediction = (ww1 * lr_prediction1 + ww2 * lr_prediction2) / (1 - ww3)
        }

        if (model_type == "rr" || model_type == "both") {
            # Train rr model
            print("training EM rr")
 
            lr_param1 <- norm_fix(y, sy, x1)
            lr_param2 <- norm_fix(y, sy, x2)
            S = y[sy]
            Z = x1[sy,]
            Y = x2[sy,]
            N = length(S)

            pi1 = mix_model_2_param$pi_1_m[v,t,impi]
            pi2 = mix_model_2_param$pi_2_m[v,t,impi]
            w1 = mix_model_2_param$w1_lst[[impi]][[(v-1)*num_time_point+t]]
            w2 = mix_model_2_param$w2_lst[[impi]][[(v-1)*num_time_point+t]]
            if (length(w1) != N) {
                w1 = rep(pi1,N)
                w2 = rep(pi2,N)
            }

            X = cbind(Z,Y)
            U1 = apply(X,2,mean)
            U2 = U1

            S1 = Reduce('+',lapply(split(X,1:nrow(X)),function(row) {(row-U1)%*%t(row-U1)})) / N
            S2 = S1

            sink(sprintf("%s_rr_em_params.txt",w_fn))

            rr_param = em_double_reg(S,Z,Y,T,t,w1,w2,pi1,pi2,U1,U2,S1,S2,lr_param1$beta,lr_param1$sigma,lr_param2$beta,lr_param2$sigma,em_max_iter,tolerance,ridge)
            
            sink()

            # save rr params
            pi1 = rr_param$pi1
            pi2 = rr_param$pi2
            w1 = rr_param$w1
            w2 = rr_param$w2
            U1 = rr_param$U1
            U2 = rr_param$U2
            S1 = rr_param$S1
            S2 = rr_param$S2
            lr_beta1 = rr_param$lr_beta1
            lr_sigma1 = rr_param$lr_sigma1
            lr_beta2 = rr_param$lr_beta2
            lr_sigma2 = rr_param$lr_sigma2
            dump("pi1", sprintf("%s_rr.pi1",w_fn))
            dump("pi2", sprintf("%s_rr.pi2",w_fn))
            dump("w1", sprintf("%s_rr.w1",w_fn))
            dump("w2", sprintf("%s_rr.w2",w_fn))
            dump("U1", sprintf("%s_rr.U1",w_fn))
            dump("U2", sprintf("%s_rr.U2",w_fn))
            dump("S1", sprintf("%s_rr.S1",w_fn))
            dump("S2", sprintf("%s_rr.S2",w_fn))
            dump("lr_beta1", sprintf("%s_rr.lr_beta1",w_fn))
            dump("lr_sigma1", sprintf("%s_rr.lr_sigma1",w_fn))
            dump("lr_beta2", sprintf("%s_rr.lr_beta2",w_fn))
            dump("lr_sigma2", sprintf("%s_rr.lr_sigma2",w_fn)) 

            rr_lr_prediction1 = x1[!ry,  ] %*% lr_beta1
            rr_lr_prediction2 = x2[!ry,  ] %*% lr_beta2
            
            wws = get_ww_rr(Nstar,x1[!ry,],x2[!ry,],pi1,pi2,U1,U2,S1,S2)
            ww1 = wws$w1; ww2 = wws$w2

            dump("ww1", sprintf("%s_rr.ww1",w_fn))
            dump("ww2", sprintf("%s_rr.ww2",w_fn))

            rr_prediction = ww1 * rr_lr_prediction1 + ww2 * rr_lr_prediction2
            rr_pi_prediction = pi1 * rr_lr_prediction1 + pi2 * rr_lr_prediction2
        }

        if (model_type == "both") {
            if (rr_param$abs_error < rrg_param$abs_error) {
                mix_model_num = rr_param$mix_model_num
                prediction = rr_prediction
            } else {
                mix_model_num = rrg_param$mix_model_num
                prediction = rrg_prediction
            }
        }

        if (model_type == "rrg") {
            mix_model_num = rrg_param$mix_model_num
            prediction = rrg_prediction
        }

        if (model_type == "rr") {
            mix_model_num = rr_param$mix_model_num
            prediction = rr_prediction
        }
        
        dump("mix_model_num", sprintf("%s.mix_model_num",w_fn))

        sink(sprintf("%s_pred_error.txt",w_fn))
        print(sprintf("mix_model_num: %s", mix_model_num))

        if (!is.null(ori_y)) {
            print(sprintf("num masked: %s", sum(test_r)))
            print(sprintf("pred_error: %s", sum(abs(ori_y_te[test_r] - prediction[test_r]))))

            if (model_type == "both" || model_type == "rr") {
                print(sprintf("rr_pred_error: %s", sum(abs(ori_y_te[test_r] - rr_prediction[test_r]))))
                print(sprintf("rr_pi_pred_error: %s", sum(abs(ori_y_te[test_r] - rr_pi_prediction[test_r]))))
                print(sprintf("rr_reg1_pred_error: %s", sum(abs(ori_y_te[test_r] - rr_lr_prediction1[test_r]))))
                print(sprintf("rr_reg2_pred_error: %s", sum(abs(ori_y_te[test_r] - rr_lr_prediction2[test_r]))))
            }

            if (model_type == "both" || model_type == "rrg") {
                print(sprintf("rrg_pred_error: %s", sum(abs(ori_y_te[test_r] - rrg_prediction[test_r]))))
                print(sprintf("rrg_pi_pred_error: %s", sum(abs(ori_y_te[test_r] - rrg_pi_prediction[test_r]))))
                print(sprintf("rrg_reg1_pred_error: %s", sum(abs(ori_y_te[test_r] - lr_prediction1[test_r]))))
                print(sprintf("rrg_reg2_pred_error: %s", sum(abs(ori_y_te[test_r] - lr_prediction2[test_r]))))
                print(sprintf("rrg_rescale_rr_pred_error: %s", sum(abs(ori_y_te[test_r] - rrg_rescale_rr_prediction[test_r]))))
                if (mix_model_num != 2) {
                    print(sprintf("GP_pred_error: %s", sum(abs(ori_y_te[test_r] - gp_prediction[test_r]))))
                }
            }
        }

        sink()

    } else {

        if (observation) {
            sy = ry
            for (i in 1:length(sy)) {
                if (sy[i] == TRUE) {
                    ts = pt_df[i,-t][r_v[i,]]
                    if (sum(r_v[i,]) == 0) {
                        sy[i] = FALSE
                    }
                }
            }

            T = dim(pt_df)[2]
            N = sum(sy)
            S = y[sy]
            Z = x1[sy,]
            Yreg = x2[sy,]
            Ygp = pt_df[sy,]
        }

        mix_model_num = source(sprintf("%s.mix_model_num",w_fn))$value
        
        if (model_type == "rr" || model_type == "both") {
            print("loading rr EM params")

            pi1 = source(sprintf("%s_rr.pi1",w_fn))$value
            pi2 = source(sprintf("%s_rr.pi2",w_fn))$value
            w1 = source(sprintf("%s_rr.w1",w_fn))$value
            w2 = source(sprintf("%s_rr.w2",w_fn))$value
            U1 = source(sprintf("%s_rr.U1",w_fn))$value
            U2 = source(sprintf("%s_rr.U2",w_fn))$value
            S1 = source(sprintf("%s_rr.S1",w_fn))$value
            S2 = source(sprintf("%s_rr.S2",w_fn))$value
            ww1 = source(sprintf("%s_rr.ww1",w_fn))$value
            ww2 = source(sprintf("%s_rr.ww2",w_fn))$value
            lr_beta1 = source(sprintf("%s_rr.lr_beta1",w_fn))$value
            lr_sigma1 = source(sprintf("%s_rr.lr_sigma1",w_fn))$value
            lr_beta2 = source(sprintf("%s_rr.lr_beta2",w_fn))$value
            lr_sigma2 = source(sprintf("%s_rr.lr_sigma2",w_fn))$value

            if (observation) {
                print("saving ww_tr")
                wws = get_ww_rr(N,Z,Yreg,pi1,pi2,U1,U2,S1,S2)
                ww1_tr = wws$w1; ww2_tr = wws$w2
                dump("ww1_tr", sprintf("%s_rr.ww1_tr",w_fn))
                dump("ww2_tr", sprintf("%s_rr.ww2_tr",w_fn))
            }

            rr_param = list(pi1,pi2,w1,w2)
            names(rr_param) = c('pi1','pi2','w1','w2')

            if (mix_model_num == 2) {
                rr_lr_prediction1 = x1[!ry,  ] %*% lr_beta1
                rr_lr_prediction2 = x2[!ry,  ] %*% lr_beta2
                prediction = ww1 * rr_lr_prediction1 + ww2 * rr_lr_prediction2
            }
        }
        if (model_type == "rrg" || model_type == "both") {
            print("loading rrg EM params")

            ll = source(sprintf("%s.ll",w_fn))$value
            pi1 = source(sprintf("%s_rrg.pi1",w_fn))$value
            pi2 = source(sprintf("%s_rrg.pi2",w_fn))$value
            pi3 = source(sprintf("%s_rrg.pi3",w_fn))$value
            w1 = source(sprintf("%s_rrg.w1",w_fn))$value
            w2 = source(sprintf("%s_rrg.w2",w_fn))$value
            w3 = source(sprintf("%s_rrg.w3",w_fn))$value
            U1 = source(sprintf("%s_rrg.U1",w_fn))$value
            U2 = source(sprintf("%s_rrg.U2",w_fn))$value
            U3 = source(sprintf("%s_rrg.U3",w_fn))$value
            S1 = source(sprintf("%s_rrg.S1",w_fn))$value
            S2 = source(sprintf("%s_rrg.S2",w_fn))$value
            S3 = source(sprintf("%s_rrg.S3",w_fn))$value
            ww1 = source(sprintf("%s_rrg.ww1",w_fn))$value
            ww2 = source(sprintf("%s_rrg.ww2",w_fn))$value
            ww3 = source(sprintf("%s_rrg.ww3",w_fn))$value
            lr_beta1 = source(sprintf("%s_rrg.lr_beta1",w_fn))$value
            lr_sigma1 = source(sprintf("%s_rrg.lr_sigma1",w_fn))$value
            lr_beta2 = source(sprintf("%s_rrg.lr_beta2",w_fn))$value
            lr_sigma2 = source(sprintf("%s_rrg.lr_sigma2",w_fn))$value

            if (observation) {
                print("saving ww_tr")
                wws = get_ww(N,t,Ygp,Z,Yreg,pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
                ww1_tr = wws$w1; ww2_tr = wws$w2; ww3_tr = wws$w3
                dump("ww1_tr", sprintf("%s_rrg.ww1_tr",w_fn))
                dump("ww2_tr", sprintf("%s_rrg.ww2_tr",w_fn))
                dump("ww3_tr", sprintf("%s_rrg.ww3_tr",w_fn))
            }

            rrg_param <- list(pi1,pi2,pi3,w1,w2,w3)
            names(rrg_param) <- c('pi1','pi2','pi3','w1','w2','w3')

            lr_prediction1 = x1[!ry,  ] %*% lr_beta1
            lr_prediction2 = x2[!ry,  ] %*% lr_beta2

            Ystar = pt_df[!ry,]
            
            xtr_vec_star = xtr_vec[!ry,]
            xte_vec_star = xte_vec[!ry]
            r_v_star = r_v[!ry,]

            if (Nstar == 1) {
                Ystar = t(as.matrix(Ystar))
            }

            if (mix_model_num != 2) {
                lr_prediction1 = x1[!ry,  ] %*% lr_beta1
                lr_prediction2 = x2[!ry,  ] %*% lr_beta2

                Ystar = pt_df[!ry,]
                
                xtr_vec_star = xtr_vec[!ry,]
                xte_vec_star = xte_vec[!ry]
                r_v_star = r_v[!ry,]

                if (Nstar == 1) {
                    Ystar = t(as.matrix(Ystar))
                }

                GPprediction_res = mclapply(1:Nstar, function(i) gp_predict_one_rt(ll,xtr_vec_star[i,][r_v_star[i,]],Ystar[i,-t][r_v_star[i,]],xte_vec_star[i]), mc.cores=num_cores)
                gp_prediction = sapply(GPprediction_res, function(x) x$pred)

                if (mix_model_num == 1) {
                    prediction = ww1 * lr_prediction1 + ww2 * lr_prediction2 + ww3 * gp_prediction
                } else {
                    prediction = gp_prediction
                }
            }
        }
        
        sink(sprintf("%s_loadEM_pred_error.txt",w_fn))
        
        print(sprintf("mix_model_num: %s", mix_model_num))

        if (!is.null(ori_y)) {
            print(sprintf("num masked: %s", sum(test_r)))
            print(sprintf("pred_error: %s", sum(abs(ori_y_te[test_r] - prediction[test_r]))))
        }
        sink()
    }

    res$pred = as.matrix(prediction)
    res$rr_em_param = rr_param
    res$rrg_em_param = rrg_param
    return(res)
}

mixtureMITemporal <- function(pv_tensor, prt_m=NULL,
    artificial_prt_tensor = NULL,
    ori_tensor = NULL,
    model_type = "rrg",
    m = 5,
    exclude = "",
    maxit = 5,
    obs_only = FALSE,
    rg = FALSE,
    em_max_iter = 10,
    tolerance = 1,
    step = 0.01,
    gd_miter = 15,
    gd_precision = 0.1,
    ridge = 1e-5,
    out_cdn = "",
    gpmodel_dir = "",
    imp_tensor = NA,
    printFlag = TRUE,
    seed = NA,
    observation = FALSE,
    ...)
{
    set.seed(seed)

    dir.create(out_cdn,recursive=TRUE)

    # real (observed) values indicator, per time point
    r_list <- lapply(pv_tensor, function(x) !is.na(x))

    # real (observed) values indicator, per variable
    r_vlist = NULL
    if (obs_only) {
        r_vlist = list()
        num_time_point = length(pv_tensor)
        num_var <- dim(pv_tensor[[1]])[2]
        for (v in 2:num_var) {
            # assume the first variable (time) is excluded
            pt_list <- vector('list',num_time_point)
            for (t in 1:num_time_point) {
                pt_list[[t]] <- pv_tensor[[t]][,v]
            }
            pt_df <- do.call(cbind,pt_list)
            r_vlist[[v]] = !is.na(pt_df)
        }
    }

    visit_col_sequence_list <- lapply(pv_tensor, function(x) (1:ncol(x))[apply(is.na(x),2,any)])
    
    # create predictor_matrix
    predictor_matrix_list <- lapply(pv_tensor, function(x) (1 - diag(1, ncol(x))))
    
    # exclude the "exclude" columns from predictor_matrix
    predictor_matrix_list <- lapply(predictor_matrix_list, function(x) {x[,pmatch(exclude, names(pv_tensor[[1]]))] <- 0; x})

    imp_tensor <- initialize_imp_tensor(pv_tensor,m)

    sink(sprintf("%s/em_gd_params.txt",out_cdn))
    print(sprintf("em_max_iter = %s",em_max_iter))
    print(sprintf("tolerance = %s",tolerance))
    print(sprintf("step = %s",step))
    print(sprintf("gd_miter = %s",gd_miter))
    print(sprintf("gd_precision = %s",gd_precision))
    sink()

    imp_tensor = sampler(pv_tensor, prt_m, artificial_prt_tensor, ori_tensor, model_type, out_cdn, gpmodel_dir, m, maxit, obs_only, imp_tensor, r_list, r_vlist, predictor_matrix_list, visit_col_sequence_list, em_max_iter, tolerance, step, gd_miter, gd_precision, ridge, observation, printFlag, ...)

    for (t in 1:length(pv_tensor)) {
        for (i in 1:m) {
            for (v in visit_col_sequence_list[[t]]) {   
                pv_tensor[[t]][!r_list[[t]][,v],v] <- imp_tensor[[t]][[v]][,i]
                write.csv(pv_tensor[[t]], file=sprintf("%s/imp_%d_attime_%s.csv",out_cdn,i,t), row.names=FALSE)
            }
        }
    }
}


