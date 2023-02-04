#!/usr/bin/bash
start=`date +%s`
Num_Loops=1
#Wcb=(25 30 35)
Wcb=(315)
Wcb_length=${#Wcb[@]}
#Q_vals=(10000)
Q_vals=(28000)
#Q_vals=(4000 12000 20000 28000)
#Q_denominater is 10 million, (Q_vals/Q_denom) * # particles = # rand
Q_length=${#Q_vals[@]}
#t_arr=(550 600)
t_arr=(3320)
Num_t=${#t_arr[@]}
#fill_arr=(80 85 90 95)
fill_arr=(80)
Num_fill=${#fill_arr[@]}
XY=3320
Rc=20000
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
                mpirun -hostfile machines --bind-to none Generate_Training_Data.sh ${Q_vals[$z]} ${Q_vals[$z]} ${t_arr[$i]} $XY ${Wcb[$k]} ${fill_arr[$m]} $Rc
                echo "Compressing"
		mpirun -hostfile compression compression_MPI.out
                python Move_Data.py
		python Move_Output.py
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
Subject: MPI_Loop Complete

Total Time = $runtime
EOF
