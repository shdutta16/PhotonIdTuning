#!/bin/sh

DIR=/afs/cern.ch/user/s/shdutta/public/ServiceTask/PhotonIdTuning/merged/AllPV/PF/Trainner/bar

echo "Working directory: "$DIR
cd $DIR
echo "Moved to working directory."

export SCRAM_ARCH=slc7_amd64_gcc820
export CPATH=$CPATH:$DIR
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$DIR

export X509_USER_PROXY=/afs/cern.ch/user/s/shdutta/proxies/x509up_u128040

source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "sourced"
pwd

eval cmsenv

cd ${_CONDOR_SCRATCH_DIR}
pwd

echo "Proxy info:"
voms-proxy-info -all
echo ""

hostname

root -l -b Reg.C
cp -r weights/*.xml .
cp -r weights/*.C .
