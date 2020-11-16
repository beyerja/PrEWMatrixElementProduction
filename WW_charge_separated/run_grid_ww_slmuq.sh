#!/bin/bash
#

. /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v01-19-02/init_ilcsoft.sh
cd /nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/WW_charge_separated
pwd
currentdate=$(date +"%Y-%m-%d %T %z")
echo $currentdate
echo "ww_slmuq-"$1"-"$2"-"$3"-"
/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/second_copy_GF/omega/grid_ww_sl0muq.exe $1 $2 $3
currentdate=$(date +"%Y-%m-%d %T %z")
echo $currentdate
