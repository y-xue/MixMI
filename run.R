run <- function(out_cdn, data_prefix, tp, missing_pcnt, gpmodel_dir="", model_type="rrg", ridge=1e-5, obs_only = TRUE, observation = FALSE)
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

	# MIMIC-III

	# norm_marked_prt_m
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked.prt_m",data_prefix,missing_pcnt*100))
	
	# prt_m
	load(sprintf("../../data/%s_tr_seed5346.prt_m",data_prefix))
	
	# norm_marked_pv_tensor
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked.pv_tensor",data_prefix,missing_pcnt*100))
	
	# file_list
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked.file_list",data_prefix,missing_pcnt*100))
	
	# ori_norm_pv_tensor
	load(sprintf("../../data/%s_tr_seed5346_norm.pv_tensor",data_prefix))

	# # artificial_prt_tensor
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked_artificial_time_per_v.prt_tensor",data_prefix,missing_pcnt*100))

	# # # tuned_artificial_prt_tensor
	# load(sprintf("../../data/%s_norm_%smispcnt_test_marked_tuned_artificial_time_per_v_a0.5.prt_tensor",data_prefix,missing_pcnt*100))


	# # artificial experiment
	# mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, artificial_prt_tensor=artificial_prt_tensor, ori_tensor=ori_norm_pv_tensor, model_type=model_type, m = 3, exclude = exclude, maxit = 3, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)
	
	# real data experiment
	# mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, ori_tensor=ori_norm_pv_tensor, model_type=model_type, m = 3, exclude = exclude, maxit = 2, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, ridge = ridge, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed, observation=observation)

	# impute original data
	mixtureMITemporal(ori_norm_pv_tensor, prt_m=prt_m, ori_tensor=ori_norm_pv_tensor, model_type=model_type, m = 3, exclude = exclude, maxit = 2, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, ridge = ridge, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed, observation=observation)


	# # EDW

	# # norm_marked_prt_m
	# load(sprintf("../../data/edw_time_series/%s_norm_%smispcnt_test_marked.prt_m",data_prefix,missing_pcnt*100))
	
	# # norm_marked_pv_tensor
	# load(sprintf("../../data/edw_time_series/%s_norm_%smispcnt_test_marked.pv_tensor",data_prefix,missing_pcnt*100))

	# # norm_res
	# load(sprintf("../../data/edw_time_series/%s.norm_res",data_prefix))

	# if (!file.exists(sprintf("../../data/edw_time_series/%s_norm.pv_tensor",data_prefix))) {
	# 	norm_tv_tensor = norm_res$V
	# 	num_pt = length(norm_tv_tensor)
	# 	num_time_point = nrow(norm_tv_tensor[[1]])
	# 	ori_norm_pv_tensor <- vector('list',num_time_point)
	# 	for (t in 1:num_time_point) {
	# 		pv_list <- vector('list',num_pt)
	# 		for (i in 1:num_pt) {
	# 			pv_list[[i]] <- norm_tv_tensor[[i]][t,]
	# 		}
	# 		pv_df <- do.call(rbind,pv_list)
	# 		ori_norm_pv_tensor[[t]] <- pv_df
	# 	}
	# 	save("ori_norm_pv_tensor",file=sprintf("../../data/edw_time_series/%s_norm.pv_tensor",data_prefix))
	# }

	# # ori_norm_pv_tensor
	# load(sprintf("../../data/edw_time_series/%s_norm.pv_tensor",data_prefix))

	# # # norm_marked_tv_tensor
	# # load("../../data/edw_time_series/all_data_7tp_1measure_norm_20mispcnt_test_marked.tv_tensor")
	
	# # # file_list
	# # load("../../data/edw_time_series/all_data_7tp_1measure_norm_20mispcnt_test_marked.file_list")

	# mixtureMITemporal(norm_marked_pv_tensor, prt_m=norm_marked_prt_m, ori_tensor=ori_norm_pv_tensor, model_type=model_type, m = 3, exclude = exclude, maxit = 5, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, ridge = ridge, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)

}

tp = 11
# tp=7
missing_pcnt = 0.2
data_prefix = "LabViewCase_11tp_1measure"
# data_prefix = "all_data_7tp_1measure"
model_type = "both"

# model_type = "rr"
# out_cdn = sprintf("../../non-equidistant_experiments/edw_%stp_1measure_norm_%smispcnt/real/joint_rr_sameXweight31_em30_3imp",tp,missing_pcnt*100)
# gpmodel_dir = ""
# run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir,model_type=model_type)

# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/joint_both_sameXweight31_em10_pi3_0.1_3imp_overfittingEM",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/joint_both_sameXweight31_em30_pi3_0.1_3imp",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/artificial_joint_rrg_sameXweight_em30",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/artificial_rrg_equalpi_TregAT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/tuned_artificial_0.5_joint_rrg_sameXweight31_em30/GP_models",tp,missing_pcnt*100)
# out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/tuned_artificial_0.5_joint_both_sameXweight31_em30_3imp",tp,missing_pcnt*100)


# model_type = "both"
# out_cdn = sprintf("../../non-equidistant_experiments/edw_%stp_1measure_norm_%smispcnt/real/joint_both_sameXweight31_em30_pi3_0.1_3imp",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/edw_%stp_1measure_norm_%smispcnt/real/joint_both_sameXweight31_em30_pi3_0.1_3imp/GP_models",tp,missing_pcnt*100)

# out_cdn = sprintf("../../impute_original_data_joint_both_sameXweight31_em30_pi3_0.33_3imp",tp,missing_pcnt*100)
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)

out_cdn = "../../impute_training_data_seed5346_joint_both_sameXweight31_em30_pi3_0.1_3imp"
gpmodel_dir = ""

run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir,model_type=model_type,obs_only=TRUE,observation=FALSE)

# for (ridge in c(1e-4, 1e-3, 0.01, 0.1, 1, 10)) {
# 	out_cdn = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/joint_rr_sameXweight31_em30_one_ridge_%s",tp,missing_pcnt*100,ridge)
# 	run(out_cdn,data_prefix,tp,missing_pcnt,gpmodel_dir=gpmodel_dir,ridge=ridge)
# }


impute_test <- function(out_cdn,tr_imputed_dir,data_prefix,gpmodel_dir="",model_type="both",obs_only=TRUE,ridge=1e-5) {
	source('mixtureMITemporalConfig.R')
	source('mixtureMITemporalImputeTest.R')
	library(parallel)
	library(mvtnorm)

	em_max_iter = 30
	tolerance = 1
	step = 0.02
	gd_miter = 10
	gd_precision = 0.1
	exclude = c(tcolname)

	dir.create(out_cdn,recursive=TRUE)

	load(sprintf('../../data/%s_te_seed5346_norm.pv_tensor',data_prefix))
	pv_tensor = ori_norm_pv_tensor

	load(sprintf('../../data/%s_tr_seed5346_norm.pv_tensor',data_prefix))
	tr_pv_tensor = ori_norm_pv_tensor

	load(sprintf("../../data/%s_te_seed5346.prt_m",data_prefix))

	mixtureMITemporalImputeTest(pv_tensor, tr_pv_tensor, tr_imputed_dir, prt_m=prt_m,model_type=model_type, m = 3, exclude = exclude, maxit = 2, obs_only = obs_only, em_max_iter = em_max_iter, tolerance = tolerance, step = step, gd_miter = gd_miter, gd_precision = gd_precision, ridge = ridge, out_cdn=out_cdn,  gpmodel_dir = gpmodel_dir, imp_tensor=NA, seed=seed)


}

# out_cdn = '../../impute_test'
# gpmodel_dir = sprintf("../../non-equidistant_experiments/mimic_%stp_1measure_norm_%smispcnt_test/real/rrg_equalpi_TregRT_GPObsOnly_gd30_em10/GP_models",tp,missing_pcnt*100)
# tr_imputed_dir = '../../non-equidistant_experiments/mimic_11tp_1measure_norm_20mispcnt_test/real/joint_both_sameXweight31_em30_pi3_0.1_3imp'

# impute_test(out_cdn,data_prefix,gpmodel_dir=gpmodel_dir,model_type=model_type)


out_cdn = '../../impute_test_seed5346'
gpmodel_dir = '../../impute_training_data_seed5346_joint_both_sameXweight31_em30_pi3_0.1_3imp/GP_models'
tr_imputed_dir = '../../impute_training_data_seed5346_joint_both_sameXweight31_em30_pi3_0.1_3imp'

impute_test(out_cdn,tr_imputed_dir,data_prefix,gpmodel_dir=gpmodel_dir,model_type=model_type)