# Impute missing values for examples in test set
# by the imputation model trained on training set

source('em.R')
source('regression.R')
source('GP.R')
library(parallel)

initial_impute_sample <- function(y, tr_ry, ry, x=NULL, ...)
{   
    return(sample(y[tr_ry], size=sum(!ry), replace=TRUE))
}

initialize_imp <- function(data, tr_data, m, visit_col_sequence, r, tr_r, prev_data=NA)
{
  imp <- vector("list", ncol(data))
  nmis <- apply(is.na(data),2,sum)

  if(m > 0) {
    for(j in visit_col_sequence) {
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(!r[,j]), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[r[,j] == FALSE], 1:m)
      ry <- r[,j]
      
      for(i in 1:m) {
        if (nmis[j]<nrow(data)) {
          if (is.data.frame(prev_data)) {

            imp[[j]][,i] <- prev_data[[j]][!ry]
          }
          else {
            imp[[j]][,i] <- initial_impute_sample(tr_data[,j], tr_r[,j], ry)
          }
        }
        else imp[[j]][,i] <- rnorm(nrow(data))
        
      }

    }
  }

  return(imp)
}

initialize_imp_tensor <- function(pv_tensor,m,tr_pv_tensor) {
  imp_tensor <- vector('list', length(pv_tensor))
  for (i in 1:length(pv_tensor)) {
    visit_col_sequence <- (1:ncol(pv_tensor[[i]]))[apply(is.na(pv_tensor[[i]]),2,any)]
    r <- (!is.na(pv_tensor[[i]]))
    tr_r <- (!is.na(tr_pv_tensor[[i]]))
    imp_tensor[[i]] <- initialize_imp(pv_tensor[[i]],tr_pv_tensor[[i]],m,visit_col_sequence,r,tr_r)

  }
  return(imp_tensor)
}

sampler <- function(pv_tensor, tr_imputed_dir, prt_m, artificial_prt_tensor, model_type, out_cdn, gpmodel_dir, m, maxit, obs_only, imp_tensor, r_list, r_vlist, predictor_matrix_list, visit_col_sequence_list, em_max_iter, tolerance, step, gd_miter, gd_precision, ridge, observation, printFlag, ...)
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

                    w_dir = sprintf("%s/em_param/imp_%s_iter_%s",tr_imputed_dir,i,k)
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

                        imp_res <- impute_em_rrg_obs_only(model_type,i,num_time_point,v,y,ry,x1,x2,pt_df,xtr_vec,xte_vec,t,r_v,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn,ridge,observation)
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

impute_em_rrg_obs_only <- function(model_type,impi,num_time_point,v,y,ry,x1,x2,pt_df,xtr_vec,xte_vec,t,r_v,l,em_max_iter,tolerance,step,gd_miter,gd_precision,w_fn,ridge=1e-5,observation=FALSE)
{
    print("impute_em_rrg_obs_only")

    rr_param = 0
    rrg_param = 0
    res <- list(0,rr_param,rrg_param)
    names(res) <- c("pred","rr_em_param","rrg_em_param")

    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)

    Nstar = sum(!ry)

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
        lr_beta1 = source(sprintf("%s_rr.lr_beta1",w_fn))$value
        lr_sigma1 = source(sprintf("%s_rr.lr_sigma1",w_fn))$value
        lr_beta2 = source(sprintf("%s_rr.lr_beta2",w_fn))$value
        lr_sigma2 = source(sprintf("%s_rr.lr_sigma2",w_fn))$value

        rr_param = list(pi1,pi2,w1,w2)
        names(rr_param) = c('pi1','pi2','w1','w2')

        if (mix_model_num == 2) {
            rr_lr_prediction1 = x1[!ry,  ] %*% lr_beta1
            rr_lr_prediction2 = x2[!ry,  ] %*% lr_beta2

            wws = get_ww_rr(Nstar,x1[!ry,],x2[!ry,],pi1,pi2,U1,U2,S1,S2)
            ww1 = wws$w1; ww2 = wws$w2
   
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
        lr_beta1 = source(sprintf("%s_rrg.lr_beta1",w_fn))$value
        lr_sigma1 = source(sprintf("%s_rrg.lr_sigma1",w_fn))$value
        lr_beta2 = source(sprintf("%s_rrg.lr_beta2",w_fn))$value
        lr_sigma2 = source(sprintf("%s_rrg.lr_sigma2",w_fn))$value

        rrg_param <- list(pi1,pi2,pi3,w1,w2,w3)
        names(rrg_param) <- c('pi1','pi2','pi3','w1','w2','w3')

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
                wws = get_ww(Nstar,t,Ystar,x1[!ry,],x2[!ry,],pi1,pi2,pi3,U1,U2,U3,S1,S2,S3)
                ww1 = wws$w1; ww2 = wws$w2; ww3 = wws$w3
                prediction = ww1 * lr_prediction1 + ww2 * lr_prediction2 + ww3 * gp_prediction
            } else {
                prediction = gp_prediction
            }
        }
    }

    res$pred = as.matrix(prediction)
    res$rr_em_param = rr_param
    res$rrg_em_param = rrg_param
    return(res)
}

mixtureMITemporalImputeTest <- function(pv_tensor, 
    tr_pv_tensor,
    tr_imputed_dir,
    prt_m=NULL,
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

    r_list <- lapply(pv_tensor, function(x) !is.na(x))

    r_vlist = NULL
    if (obs_only) {
        r_vlist = list()
        num_time_point = length(pv_tensor)
        num_var <- dim(pv_tensor[[1]])[2]
        for (v in 2:num_var) {
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

    imp_tensor <- initialize_imp_tensor(pv_tensor,m,tr_pv_tensor)

    sink(sprintf("%s/em_gd_params.txt",out_cdn))
    print(sprintf("em_max_iter = %s",em_max_iter))
    print(sprintf("tolerance = %s",tolerance))
    print(sprintf("step = %s",step))
    print(sprintf("gd_miter = %s",gd_miter))
    print(sprintf("gd_precision = %s",gd_precision))
    sink()

    imp_tensor = sampler(pv_tensor, tr_imputed_dir, prt_m, artificial_prt_tensor, model_type, out_cdn, gpmodel_dir, m, maxit, obs_only, imp_tensor, r_list, r_vlist, predictor_matrix_list, visit_col_sequence_list, em_max_iter, tolerance, step, gd_miter, gd_precision, ridge, observation, printFlag, ...)

    for (t in 1:length(pv_tensor)) {
        for (i in 1:m) {
            for (v in visit_col_sequence_list[[t]]) {   
                pv_tensor[[t]][!r_list[[t]][,v],v] <- imp_tensor[[t]][[v]][,i]
                write.csv(pv_tensor[[t]], file=sprintf("%s/imp_%d_attime_%s.csv",out_cdn,i,t), row.names=FALSE)
            }
        }
    }
}