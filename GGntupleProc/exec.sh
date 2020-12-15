#!/bin/sh

DIR=/afs/cern.ch/user/s/shdutta/public/ServiceTask/PhotonIdTuning/GGntupleProc/

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

#cd ${_CONDOR_SCRATCH_DIR}
#pwd

#hostname
./Run3.exe /eos/cms/store/group/phys_egamma/shdutta/PhotonIdTuning/Input/in_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2021/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2021.root /eos/cms/store/group/phys_egamma/shdutta/PhotonIdTuning/Output/out_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2021/ntuple_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2021.root 10000 100 1.0

./Run3.exe /eos/cms/store/group/phys_egamma/shdutta/PhotonIdTuning/Input/in_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2021/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2021.root /eos/cms/store/group/phys_egamma/shdutta/PhotonIdTuning/Output/out_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2021/ntuple_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2021.root 10000 100 1.0
