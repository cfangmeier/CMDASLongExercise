#!/usr/bin/env bash
function start_job {
    ./ZTT_XSection.exe ${1}.root  root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/${1}.root  > ${1}.log &
}

start_job DYJetsToLL
start_job TTJets
start_job WJetsToLNu
#start_job SingleMu




#  ./ZTT_XSection.exe DYJetsToLL.root  root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/DYJetsToLL.root
#  ./ZTT_XSection.exe TTJets.root      root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/TTJets.root 
#  ./ZTT_XSection.exe WJetsToLNu.root  root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/WJetsToLNu.root 
#  ./ZTT_XSection.exe SingleMu.root    root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/SingleMu.root  
