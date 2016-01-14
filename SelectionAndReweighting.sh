#!/usr/bin/env bash
function start_job {
    ./ZTT_XSection.exe ${1}.root  root://cmseos.fnal.gov://store/user/abdollah/ROOTHadd/${1}.root  > ${1}.log &
}

start_job DYJetsToLL
start_job TTJets
start_job WJetsToLNu
start_job SingleEle

