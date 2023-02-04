#!/usr/bin/bash
start=`date +%s`
Num_Loops=5
#Wcb=(215 265 315 365 415)
Wcb=(315)
Wcb_length=${#Wcb[@]}
Q_vals=(100000)
#Q_vals=(20 40 60 80 100)
Q_length=${#Q_vals[@]}
t_arr=(200)
#t_arr=(50 63 75 88 100)
Num_t=${#t_arr[@]}
fill_arr=(25 50 75 100)
Num_fill=${#fill_arr[@]}
XY=4000
for ((k=0; k<$Wcb_length; k++))
do
    for ((z=0; z<$Q_length; z++))
    do
        for ((i=0; i<$Num_t; i++))
        do
            for ((j=0; j<$Num_Loops; j++))
            do
            mpirun -hostfile machines Generate_Training_Data.sh ${Q_vals[$z]} ${Q_vals[$z]} ${t_arr[$i]} $XY ${Wcb[$k]} ${fill_arr[$m]}
            mpirun -hostfile compression compression_MPI.out
            python Move_Data.py
            mid=`date +%s`
            runtime_mid=$((mid-start))
            echo $runtime_mid ${Wcb[$k]}
            done
        done
    done
done
end=`date +%s`
runtime=$((end-start))
echo $runtime

sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Loop Complete

Total Time = $runtime
EOF
