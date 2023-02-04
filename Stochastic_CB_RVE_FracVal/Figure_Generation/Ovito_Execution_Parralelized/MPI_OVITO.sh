#!/usr/bin/bash
#$1 is the folder name w/out a forward slash, all the way from origin (i.e. /hom/talbright/.../Working_Folder), $2 is pixel count (i.e. 200x200), $3 is the XY size in nm (i.e. if working folder is ".../5000_3.15", you enter 5000 for $3 and $4 is 3.15
#$5 is the number of threads to give ovito to use
start=`date +%s`
mpirun -hostfile machines OVITO_Image_Gen.out $1 $2 $3 $4 $5
end=`date +%s`
runtime=$((end-start))
echo $runtime
sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Images Complete

Total Time = $runtime

EOF
