#!/usr/bin/bash
uid=$(lsblk -nro SERIAL)
bid=$$
Lower_Q=$1
Upper_Q=$2
thickness=$3
XY=$4
WCB=$5
fill=$6
./a.out $Lower_Q $Upper_Q $XY $thickness $uid $bid $WCB $fill
