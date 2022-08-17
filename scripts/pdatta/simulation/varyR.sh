#!/bin/bash

imin=$1
imax=$2

for (( i=$imin; i<$imax; i++ ))
do
    # We need to do some floating point arithmetic, 
    # let's use bc (basic/bench calculator)
    R=$(echo "scale=2; (80+$i)/100" | bc)
    analyzer -b -q 'simu_data_fit.cpp("sbs4_sbs50p_hcaldxdy_simu.cfg",'$R',"siout/simu_data_fit_SBS50_.82T_R'$R'.root")'
done

# R=$1
# analyzer -b -q 'simu_data_fit.cpp("sbs4_sbs50p_hcaldxdy_simu.cfg",'$R',"simu_data_fit_SBS50_.82T_R'$R'.root")'

