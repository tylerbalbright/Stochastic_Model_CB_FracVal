#!/usr/bin/bash
uid=$(lsblk -nro SERIAL)
bid=$$
Lower_Q=100000
Upper_Q=100000
thickness=200
XY=13333
WCB=315
fill=50
start=`date +%s`
./a.out $Lower_Q $Upper_Q $XY $thickness $uid $bid $WCB $fill
end=`date +%s`
runtime=$((end-start))
echo $runtime
