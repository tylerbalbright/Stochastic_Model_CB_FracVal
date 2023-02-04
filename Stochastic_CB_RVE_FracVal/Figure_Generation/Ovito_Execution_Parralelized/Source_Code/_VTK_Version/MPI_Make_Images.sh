#!/usr/bin/bash
#$1 is the folder name w/out a forward slash, $2 is pixel count (i.e. 200x200)
start=`date +%s`
mpirun -hostfile machines a.out $1 $2
end=`date +%s`
runtime=$((end-start))
echo $runtime
sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Images Complete

Total Time = $runtime
EOF
