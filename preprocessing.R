tv_tensor_generator <- function(cdn) {
	cat('tv_tensor_generator\n')

	file_list = sort(list.files(cdn))
	num_pt = length(file_list)
	df_list = vector('list',num_pt)
	
	for (i in 1:length(file_list)) {
		df_list[[i]] = read.csv(sprintf("%s/%s",cdn,file_list[i]))
	}
	return(df_list)
}

write_tv_tensor <- function(tv_tensor, out_cdn, file_list) {
	dir.create(out_cdn,recursive=TRUE)

	print(length(tv_tensor))
	print(length(file_list))

	for (i in 1:length(tv_tensor)) {
		data = tv_tensor[[i]]
		write.csv(data, file=sprintf("%s/%s",out_cdn,file_list[[i]]), row.names=FALSE)
	}
}

unfold_tensor_by_patients <- function(nimp,dir,out_dir,filename_list,tp,impi=1,iter=0,normpred=FALSE) {
    for (imp_i in impi:nimp) {
		print(imp_i)

		w_dir = sprintf("%s/%s",out_dir,imp_i)

		dir.create(w_dir, recursive=TRUE)
	    tensor_list = vector('list',tp)
	    for (i in 1:tp) {
	    	if (normpred) {
	    		tensor_list[[i]] = read.csv(sprintf("%s/imp_%s_attime_%s_normpred.csv",dir,imp_i,i))
	    	}
	    	else {
	    		if (iter == 0) {
	    			tensor_list[[i]] = read.csv(sprintf("%s/imp_%s_attime_%s.csv",dir,imp_i,i))
	    		}
	    		else {
	    			tensor_list[[i]] = read.csv(sprintf("%s/imp_%s_attime_%s_miceiter_%s.csv",dir,imp_i,i,iter))
	    		}
	    	}
	    }

	    for (i in 1:nrow(tensor_list[[1]])) {
	        l = lapply(tensor_list, function(x) unlist(x[i,]))
	        df = do.call(rbind,l)

	        write.csv(df,file=sprintf("%s/%s",w_dir,filename_list[i]),row.names=FALSE)
	    }

    }
}

combine_m_imp_ptview <- function(in_dir,out_dir,nimp) {
	dir.create(out_dir,recursive=TRUE)
	filename_list = list.files(sprintf("%s/1",in_dir))
	filename_list = sort(filename_list)
	for (file in filename_list) {
		imp = read.csv(sprintf('%s/%s/%s',in_dir,1,file))
		cols = colnames(imp)
		xt = array(NA, dim=c(dim(imp)[1], dim(imp)[2], nimp))
		xt[,,1] = as.matrix(imp)
		for (i in 2:nimp) {
			imp = read.csv(sprintf('%s/%s/%s',in_dir,i,file))
			xt[,,i] = as.matrix(imp)
		}
		xm = apply(xt, c(1,2), mean)
		df = data.frame(xm)
		colnames(df) = cols
		write.csv(df,file=sprintf('%s/%s',out_dir,file),row.names=FALSE)
	}
}

self_norm <- function(Y, center=F, logt=F, fncf='mixtureMITemporalConfig.R') {
    source(fncf)
    V = Y; n = length(V)
	hmax=list()
	hmin=list()
    
    cat(sprintf('self normalizing\n'))
    # print(colnames(V[[1]]))
    for (ipt in 1:length(V)) {
        for (col in tests) {
            # print(col)
            k = sprintf('%s|%s', ipt, col)
            pv = V[[ipt]][,col]
            kmax = max(pv, na.rm=T);
            hmax[[k]] = kmax
            
            kmin = min(pv, na.rm=T);
            hmin[[k]] = kmin
            if (center) {
                offset = (kmin+kmax)/2
            }else {
                offset = kmin
            }
            if (kmax == kmin) {
            	V[[ipt]][!is.na(V[[ipt]][,col]),col] = 1
            }
            else {
            	V[[ipt]][,col] = (pv - offset) / (kmax - kmin)
            }
        }
    }
    if (logt) {
        for (ipt in 1:length(V)) {
            V[[ipt]][,cols] = log(V[[ipt]][,cols]+1)
        }
    }
    cat(sprintf('done\n'))

    res = list(V,hmax,hmin)
    names(res) = c('V','hmax','hmin')
    return(res)
}

self_denorm <- function(V, hmax, hmin, center=F, logt=F, fncf='mixtureMITemporalConfig.R') {
    source(fncf)
    V.raw = V; n = length(V)
    cat(sprintf('self denormalizing\n'))
    for (ipt in 1:length(V)) {
        for (col in tests) {
            if (logt) {
                V.raw[[ipt]][,col] = exp(V[[ipt]][,col]) - 1
            }
            k = sprintf('%s|%s', ipt, col)
            kmax = hmax[[k]]
            kmin = hmin[[k]]
            if (center) {
                offset = (kmin+kmax)/2
            }else {
                offset = kmin
            }
            V.raw[[ipt]][,col] = V.raw[[ipt]][,col] * (kmax - kmin) + offset
            
        }
    }
    cat(sprintf('done\n'))
    return (V.raw)
}

denorm_dir <- function(cdn, out_cdn, hmax, hmin, center=F, logt=F) {
	tv_tensor = tv_tensor_generator(cdn)

	imp_tv_tensor = self_denorm(tv_tensor, hmax, hmin)

	file_list <- list.files(cdn)
	file_list <- sort(file_list)
	write_tv_tensor(imp_tv_tensor, out_cdn, file_list)
}