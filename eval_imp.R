eval_imp <- function(cdn, nimp, iteri, norm_res_fn, file_list_fn, tp, ori_dir, marking_dir,denorm=FALSE) {
	source('preprocessing.R')
	source('mase.R')
	source('mixtureMITemporalConfig.R')

	load(norm_res_fn)
	load(file_list_fn)

	hmax = norm_res$hmax
	hmin = norm_res$hmin

	unfold_dir = sprintf("%s/ptview_nimp_iter%s",cdn,iteri)
	denorm_imp_dir = sprintf("%s/ptview_%simp_iter%s_denorm",cdn,nimp,iteri)
	combined_dir = sprintf("%s/ptview_%simp_iter%s",cdn,nimp,iteri)

	if(denorm){
		unfold_tensor_by_patients(nimp,cdn,unfold_dir,file_list,tp,iter=iteri)
		combine_m_imp_ptview(unfold_dir,combined_dir,nimp)
		denorm_dir(combined_dir, denorm_imp_dir, hmax, hmin)
	}

	cal_MASE_file(ori_dir, marking_dir, denorm_imp_dir, cdn, tests, fn=sprintf("mase_nmis_%simp_iter%s.csv",nimp,iteri), ts=FALSE)
}

cdn = "/home/yxe836/micegp/non-equidistant_experiments/mimic_11tp_1measure_norm_20mispcnt_test/real/joint_both_sameXweight31_em30_pi3_0.1_3imp_exclude_HGB"
norm_res_fn = "/home/yxe836/micegp/data/LabViewCase_11tp_1measure.norm_res"
file_list_fn = "/home/yxe836/micegp/data/LabViewCase_11tp_1measure_norm_20mispcnt_test_marked.file_list"
ori_dir = "/home/yxe836/micegp/data/LabViewCase_11tp_1measure"
marking_dir = "/home/yxe836/micegp/data/LabViewCase_11tp_1measure_norm_20mispcnt_test_markings"

nimp=3
iteri=2
tp=11
eval_imp(cdn, nimp, iteri, norm_res_fn, file_list_fn, tp, ori_dir, marking_dir, denorm=TRUE)