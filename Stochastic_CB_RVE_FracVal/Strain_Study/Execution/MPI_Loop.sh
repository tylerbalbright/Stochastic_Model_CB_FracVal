#!/usr/bin/bash
#####################################################################
#This bash script uses MPI to distribute simulations across Crunch.
#The MPIrun called does not actually run an MPI enabled script,
#rather it commands each node (delineated in the "machines" files)
#to run a bash script that calls an instance of the stochastic model.
#This results in each node running multiple instances of the model
#executable. The number of executables run per node is equal to the
#number of slots set in the machines file.
#####################################################################

#####################################################################
#Variable Declarations
#####################################################################
#Code timing "start stopwatch" variable
start=`date +%s`

#Number of times to run the MPI command per set of sim parameters
Num_Loops=1

#Weight percentage desired in simulations. Set as an array to be
#looped over (for running multiple desired simulation sets)
#Wcb=(25 30 35)
Wcb=(315)

#Length of weight percentage array
Wcb_length=${#Wcb[@]}

#Dispersion quality desired in simulations. Set as an array to be
#looped over (for running multiple desired simulation sets).
#Q_denominater is 10 million, (Q_vals/Q_denom) * # particles = # rand
Q_vals=(28000)
#Q_vals=(4000 12000 20000 28000)

#Length of q value array
Q_length=${#Q_vals[@]}

#RVE thickness desired in simulations. Set as an array to be
#looped over (for running multiple desired simulation sets)
#t_arr=(550 600)
t_arr=(5000)

#Length of thickness array
Num_t=${#t_arr[@]}

#Fill percentage desired in simulations. Set as an array to be
#looped over (for running multiple desired simulation sets).
#fill_arr=(80 85 90 95)
fill_arr=(80)

#Length of fill percentage array
Num_fill=${#fill_arr[@]}

#XY length of the simulated RVE. X direction is the direction of
#potential difference for the conductivity analysis.
XY=5000
#####################################################################

#####################################################################
# Main Loop
#####################################################################

# Loop through each index in the weight percentage array
for ((k=0; k<$Wcb_length; k++))
do
# Loop through each index in the q val array
    for ((z=0; z<$Q_length; z++))
    do
# Loop through each index in the thickness array
        for ((i=0; i<$Num_t; i++))
        do
# Loop through each index in the fill pct array
            for ((m=0; m<$Num_fill; m++))
            do
# Loop through the Num_Loops
                for ((j=0; j<$Num_Loops; j++))
                do
# Run the simulations by calling the shell script on each node
                mpirun -hostfile machines --bind-to none Generate_Training_Data.sh ${Q_vals[$z]} ${Q_vals[$z]} ${t_arr[$i]} $XY ${Wcb[$k]} ${fill_arr[$m]}
                echo "Compressing"
# Run the compression script utilizing the same machines and
# a new number of slots per node
		mpirun -hostfile compression compression_MPI.out
                python Move_Data.py
#Move the data files from the temporary folder to a permenant
#storage folder in the partition
		python Move_Output.py
#Output simulation runtime
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

#Send a notification via "sendmail" to notify the user the script
#has successfully finished.
sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Loop Complete

Total Time = $runtime
EOF
