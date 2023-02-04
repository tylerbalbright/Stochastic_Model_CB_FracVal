#!/usr/bin/bash
start=`date +%s`
uid=12345
bid=12345
Lower_Q=2000
Upper_Q=2000
thickness=300
XY=13333
WCB=30
fill=20
./a.out $Lower_Q $Upper_Q $XY $thickness $uid $bid $WCB $fill
end=`date +%s`
runtime=$((end-start))
echo $runtime
