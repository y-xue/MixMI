run <- function(out_cdn, data_prefix, tp, missing_pcnt, gpmodel_dir="", model_type="rrg", ridge=1e-5, obs_only = TRUE)
{
	source('mixtureMITemporalConfig.R')
	source('mixtureMITemporal.R')
	library(parallel)
	library(mvtnorm)

	em_max_iter = 30
	tolerance = 1
	step = 0.02
	gd_miter = 10
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

	# artificial_prt_tensor
	load(sprintf("../../data/%s_norm_%smispcnt_test_marked_artificial_time_per_v.prt_tensor",data_prefix,missing_pcnt*100))

	# # # tuned_artificial_prt_tensor
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked_tuned_artificial_time_per_v_a0.5.prt_tensor",data_prefix,missing_pcnt*100))


	# artificial experiment
	mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, artificial_prt_tensor=artificial_prt_tensor, ori_tensor=ori_norm_pv_tensor, model_type=model_type, m = 3, exclude = exclude, maxit = 5, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)
	
	# # real data experiment
	# mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, ori_tensor=ori_norm_pv_tensor, m = 1, exclude = exclude, maxit = 5, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, ridge = ridge, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)

}

tp = 11
missing_pcnt = 0.2
data_prefix = "LabViewCase_11tp_1measure"
model_type = "rrg"
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/tuned_artificial_0.5_joint_rrg_sameXweight31_em30/GP_models",tp,missing_pcnt*100)
# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/joint_rr_sameXweight31_em30_3imp",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)
# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/artificial_joint_rrg_sameXweight_em30",tp,missing_pcnt*100)
gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/artificial_rrg_equalpi_TregAT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/artificial_joint_rrg_sameXweight31_em30_3imp",tp,missing_pcnt*100)
run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir,model_type=model_type)

# for (ridge in c(1e-4, 1e-3, 0.01, 0.1, 1, 10)) {
# 	out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/joint_rr_sameXweight31_em30_one_ridge_%s",tp,missing_pcnt*100,ridge)
# 	run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir,ridge=ridge)
# }
