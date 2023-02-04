#!/usr/bin/bash
gfortran Frac_VAL_CCA.f90 Ctes.f90 CCA_module.f90 PCA_cca.f90 PCA_Subclusters_module.f90 RAND_SAMPLE.f90 random.f90 a_Random_PP.f90 Save_results_CC.f90 -o go.out
start=`date +%s`
./go.out
end=`date +%s`
runtime=$((end-start))
echo $runtime
