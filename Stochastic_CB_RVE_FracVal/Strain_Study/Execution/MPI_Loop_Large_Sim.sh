#!/usr/bin/bash
start=`date +%s`
Num_Loops=5
#Wcb=(25 30 35)
Wcb=(25)
Wcb_length=${#Wcb[@]}
#Q_vals=(10000)
Q_vals=(250 500 1000)
#Q_vals=(4000 5000 6000)
#Q_denominater is 10 million, (Q_vals/Q_denom) * # particles = # rand
Q_length=${#Q_vals[@]}
#t_arr=(500 600 700)
#t_arr=(200 300 400)
t_arr=(2000)
Num_t=${#t_arr[@]}
#fill_arr=(80 85 90 95)
fill_arr=(20 40 60 80)
Num_fill=${#fill_arr[@]}
XY=13333
for ((k=0; k<$Wcb_length; k++))
do
    for ((z=0; z<$Q_length; z++))
    do
        for ((i=0; i<$Num_t; i++))
        do
            for ((m=0; m<$Num_fill; m++))
            do
                for ((j=0; j<$Num_Loops; j++))
                do
                mpirun -hostfile largemachines Generate_Training_Data.sh ${Q_vals[$z]} ${Q_vals[$z]} ${t_arr[$i]} $XY ${Wcb[$k]} ${fill_arr[$m]}
                echo "Compressing"
		mpirun -hostfile compressionlarge compression_MPI.out
                python Move_Data.py
                mid=`date +%s`
                runtime_mid=$((mid-start))
                echo $runtime_mid ${Wcb[$k]}
                done
            done
        done
    done
done
end=`date +%s`
runtime=$((end-start))
echo $runtime

sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Loop_Large Complete

Total Time = $runtime
EOF
