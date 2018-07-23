run <- function(out_cdn, data_prefix, tp, missing_pcnt, gpmodel_dir="", obs_only = TRUE)
{
	source('mixtureMITemporalConfig.R')
	source('mixtureMITemporal.R')
	library(parallel)

	em_max_iter = 30
	tolerance = 1
	step = 0.02
	gd_miter = 30
	gd_precision = 0.1
	exclude = c(tcolname)

	dir.create(out_cdn,recursive=TRUE)

	# norm_marked_prt_m
	load(sprintf("../../data/%s_norm_%smispcnt_test_marked.prt_m",data_prefix,missing_pcnt*100))
	
	# norm_marked_pv_tensor
	load(sprintf("../../data/%s_norm_%smispcnt_test_marked.pv_tensor",data_prefix,missing_pcnt*100))
	
	# file_list
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked.file_list",data_prefix,missing_pcnt*100))
	
	# ori_norm_pv_tensor
	load(sprintf("../../data/%s_norm.pv_tensor",data_prefix))

	mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, ori_tensor=ori_norm_pv_tensor, m = 3, exclude = exclude, maxit = 2, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)

}

tp = 11
missing_pcnt = 0.2
data_prefix = "LabViewCase_11tp_1measure"
out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_em30",tp,missing_pcnt*100)
gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir)
