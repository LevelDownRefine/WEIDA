# WEiDA Python migration

Four functions in `weida.py` correspond to the four numbered sections of [StateSpace_modelling.m](https://github.com/chenfei-ye/WEIDA-DBS/blob/main/StateSpace_modelling.m). Install `requirements.txt`, then run:

```
python run_weida.py /path/to/PD_DBS_zcc/derivatives/fmripost output --tr 2.1
```

The loader selects HC and STN on/off AAL1ctx CSVs, removes the index column, checks region names, and excludes STN18 in both conditions because its ON scan has 53 zero regions. Results and subject ordering are stored in NPZ and JSON. State labels are zero-based. Input is never modified. Confirm TR and inclusion criteria before interpreting results.

This is a preliminary translation, not a validated reproduction. MATLAB eigs and NumPy eigenvectors have arbitrary individual signs; eigenvector degeneracy and k-means initialization can materially change results. MATLAB/Python numerical equivalence and scientific reproducibility remain to be checked. Unlike the original, absent transition rows are NaN, and variable scan lengths are supported. The MATLAB silhouette statistic is not currently calculated.
