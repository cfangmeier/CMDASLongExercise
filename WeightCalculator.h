#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>
using namespace std;

float XSection(std::string OutName) {
    
    
//    https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591
    
    if (OutName.compare("WJets") == 0)           return 61526;
    else if (OutName.compare("DYJetsToLL") == 0) return 6025.2;
    else if (OutName.compare("TTJets") == 0)     return 831.76 ;
    else if (OutName.compare("WW") == 0)         return 309.7;
    else if (OutName.compare("WZ") == 0)         return 30400;
    else if (OutName.compare("ZZ") == 0)         return 5400;
    else                                         return 0;
}


float weightCalc(TH1F *Histo,std::string outputName) {
    
//    cout<< "outputName is "<<outputName << "  and histoname is " <<Histo->GetName()<<  " Histo->GetBinContent(1)="<<Histo->GetBinContent(1)<< " XSection(wjet)=" <<XSection("WJets")<<"\n";
    
//    float luminosity=1264-55;
    float luminosity=1560;
    
    bool isSingleMu  = outputName.find("SingleMu")   != string::npos;
    bool isSingleEle = outputName.find("SingleEle")  != string::npos;
    bool isWjet      = outputName.find("WJets")      != string::npos;
    bool isZJet      = outputName.find("DYJetsToLL") != string::npos;
    bool isTTbar     = outputName.find("TTJets")     != string::npos;
    bool isWW        = outputName.find("WW")         != string::npos;
    bool isWZ        = outputName.find("WZ")         != string::npos;
    bool isZZ        = outputName.find("ZZ")         != string::npos;
    
    if (isSingleMu || isSingleEle)   return 1;
    else if (isWjet)  return (luminosity * XSection("WJets")*1.0) / Histo->GetBinContent(2);
    else if (isZJet)  return luminosity * XSection("DYJetsToLL")*1.0 / Histo->GetBinContent(2);
    else if (isTTbar) return luminosity * XSection("TTJets")*1.0 / Histo->GetBinContent(2);
    else if (isWW)    return luminosity * XSection("WW")*1.0 / Histo->GetBinContent(2);
    else if (isWZ)    return luminosity * XSection("WZ")*1.0 / Histo->GetBinContent(2);
    else if (isZZ)    return luminosity * XSection("ZZ")*1.0 / Histo->GetBinContent(2);
    else                                return 0;
    
}


