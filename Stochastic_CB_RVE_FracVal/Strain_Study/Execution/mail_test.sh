#!/usr/bin/bash
sendmail -F STATUS_UPDATE -t talbright@ksu.edu,tylerbalbright@gmail.com << EOF
Subject: MPI_Loop Complete

Total Time = taco
EOF
