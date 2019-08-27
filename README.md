Source code of MixMI proposed in "[Mixture-based Multiple Imputation Models for Clinical Data with a Temporal Dimension](https://arxiv.org/abs/1908.04209)".

Example of running imputation model:

```
mixtureMITemporal(pv_tensor=patient_variable_view_tensor, prt_m=patient_time_matrix, model_type="both", 
m = 3, exclude = "time", maxit = 2, obs_only = True, em_max_iter = 10, tolerance = 1, step = 0.02, 
gd_miter = 10, gd_precision = 0.1, ridge = 1e-5, out_cdn="", seed=8888)
```

pv_tensor is a list (indexed by time indices) of patient\*variable matrices of measurements, containing missing values that to be imputed. 
prt_m is a patient\*time matrix that stores times when patients are measured.

Configurations can be set in mixtureMITemporalConfig.R



