#!/bin/bash

if [ "x$1" == "x" ];
then
  echo "Must specify run number!"
  exit
fi
runnum=$1

lastEntry=-1
if [ "x$2" != "x" ];
then
  lastEntry=$2
fi


analyzer -b -q -l replay_hcal.C\(${runnum},${lastEntry}\)
analyzer -b -q -l make_cosmics_tdc_rotation.C\(${runnum}\)
