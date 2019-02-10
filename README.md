Example of running imputation model:

```
mixtureMITemporal(pv_tensor=patient_variable_view_tensor, prt_m=patient_time_matrix, model_type="both", m = 3, exclude = "time", maxit = 2, obs_only = True, em_max_iter = 10, tolerance = 1, step = 0.02, gd_miter = 10, gd_precision = 0.1, ridge = 1e-5, out_cdn="", seed=8888)
```

Configurations can be set in mixtureMITemporalConfig.R



