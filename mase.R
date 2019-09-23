MASE <- function(ori_df, marking_df, imp_df, tests, ts=TRUE) {
	mis_col = rep(0,length(tests))
	mase_vec = rep(0,length(tests))
	names(mis_col) = tests # colnames(ori_df)
	names(mase_vec) = tests # colnames(ori_df)
	res = list(mase_vec, mis_col)
	names(res) = c('mase_vec', 'mis_col')
	# tp = dim(ori_df)[1]

	for (val in tests) {
	# for (val in 2:14) {
		# print(val)
		denomi = 0
		marked_mis = marking_df[,val]
		tr = ori_df[!is.na(ori_df[,val]),val]
		tp = length(tr)

		if (length(tr) < 2) {
			res$mase_vec[val] = 0
			# res$failed_cnt[val] = res$failed_cnt[val] + 1
		}
		else {
			for (i in 2:length(tr)) {
				denomi = denomi + abs(tr[i] - tr[i-1])
			}
			te_ori = ori_df[marked_mis,val]
			te_imp = imp_df[marked_mis,val]

			e = sum(abs(te_ori - te_imp))

			if (denomi == 0) {
				res$mase_vec[val] = 0
			}
			else {
				res$mase_vec[val] = e/denomi * (tp-1)/tp
			}
		}
		if (ts) {
			res$mis_col[val] = sum(marked_mis) > 0
		}
		else {
			res$mis_col[val] = sum(marked_mis)
		}
	}

	return(res)
}

cal_MASE_file <- function(ori_dir, marking_dir, imp_dir, out_dir, tests, fn="", ts=TRUE, wr=TRUE) {
	dir.create(out_dir,recursive=TRUE)
	
	file_list = sort(list.files(ori_dir))

	mase_vec = rep(0,length(tests))
	mis_col_cnt = rep(0,length(tests))
	names(mis_col_cnt) = tests
	names(mase_vec) = tests

	# for (i in 1:length(imp_tv_tensor)) {
	for (file in file_list) {
		ori_df = read.csv(sprintf("%s/%s",ori_dir,file))
		marking_df = read.csv(sprintf("%s/%s",marking_dir,file))
		imp_df = read.csv(sprintf("%s/%s",imp_dir,file))
		res = MASE(ori_df, marking_df, imp_df, tests, ts=ts)
		
		if (sum(is.na(res$mase_vec)) > 0) {
			print(file)
		}

		mase_vec = mase_vec + res$mase_vec
		mis_col_cnt = mis_col_cnt + res$mis_col
	}

	mase_one = sum(mase_vec) / sum(mis_col_cnt)
	
	print(mis_col_cnt)

	print("MASE_one")
	print(mase_one)

	mase_vec = mase_vec / mis_col_cnt
	print("MASE_ave")
	print(mase_vec)

	if (wr) {
		if (ts) {
			# sink(sprintf("%s/MASE_divide_by_missing_timeseries_scores.txt",out_dir))
			# mase_vec = c(unlist(mase_vec),mase_one)
			# names(mase_vec) = c(names(mase_vec),"ONE")
			
			mase_vec["ONE"] = mase_one
			write.csv(data.frame(mase_vec),file=sprintf("%s/mase.csv",out_dir))
		}
		else {
			mase_vec["ONE"] = mase_one
			if (fn == "") {
				write.csv(data.frame(mase_vec),file=sprintf("%s/mase_nmis.csv",out_dir))
			}
			else {
				write.csv(data.frame(mase_vec),file=sprintf("%s/%s",out_dir,fn))
			}
			# sink(sprintf("%s/MASE_divide_by_missing_value_scores.txt",out_dir))
		}
		# sink()
	}

	return(mase_one)
}