////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TreeReader.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "WeightCalculator.h"
#include <string>
#include <ostream>

/* electron_selection
 * Returns the index of the selected electron.
 * If not only a single electron is found, return -1.
 * Cuts out all di-photon events
 */
int electronSelection(TH1F* histo){
    int selection = -1;
//    cout << "====================" << endl;
//    cout << "\n\nElectrons in event: " << nEle << endl;
    for(int iele=0; iele<nEle; iele++){
        //Trigger cuts
        bool passTrigger =((HLTEleMuX >> 6 & 1) == 1  &&  isData ) ||
            ((HLTEleMuX >> 11 & 1) == 1 && !isData ); 
        //Isolation measurement
        float isoEle=elePFChIso->at(iele)/elePt->at(iele);
        if ((elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5*elePFPUIso->at(iele)) > 0.0)
            isoEle = (elePFChIso->at(iele)/elePt->at(iele) +
                    elePFNeuIso->at(iele)                +
                    elePFPhoIso->at(iele)                -
                    0.5*elePFPUIso->at(iele))/elePt->at(iele);

        bool eleMVAId; // Cuts out less-sensitive ECAL eta regions
        if (fabs(eleSCEta->at(iele)) < 0.8 && eleIDMVANonTrg->at(iele) > 0.967083)
            eleMVAId= true;
        else if (fabs(eleSCEta->at(iele)) > 0.8 && fabs(eleSCEta->at(iele)) < 1.5 && 
                eleIDMVANonTrg->at(iele) > 0.929117)
            eleMVAId= true;
        else if (fabs(eleSCEta->at(iele)) > 1.5 && eleIDMVANonTrg->at(iele) > 0.726311 )
            eleMVAId= true;
        else
            eleMVAId= false; 
        bool passCuts = elePt->at(iele) > 15    &&
            fabs(eleEta->at(iele))      < 2.4   && 
//            isoEle                      >= 0.30  &&
            isoEle                      >= 0.01  &&
            fabs(eleD0->at(iele))       < 0.045 &&
            fabs(eleDz->at(iele))       < 0.200;
//        cout << "Electron info:" << endl;
//        cout << "isoEle:" << isoEle << endl;
//        cout << "eleMVAId:" << eleMVAId << endl;
//        cout << "passTrigger:" << passTrigger << endl;
//        cout << "passCuts:" << passCuts << endl;

        if(passTrigger && eleMVAId && passCuts){
            histo->Fill(isoEle);
            if(selection == -1){
 //               cout <<"found first electron" << endl;
                selection = iele; //Finds first electron
            } else{
  //              cout << "found second electron" << endl;
                selection == -1; //Finds second electron, so fail!
                break;
            }
        }
    }
    return selection;
}

int tauIsoID()
{
    int idx = -1;
    float max_iso = -999.0;
    for(int itau = 0; itau < nTau; itau++)
    {
        bool TauPtCut = tauPt->at(itau) > 20  && fabs(tauEta->at(itau)) < 2.3 ;
        bool TauPreSelection = tauByLooseMuonRejection3->at(itau) > 0;
        //TauPreSelection = TauPreSelection && fabs(tauZImpact->at(itau)) < 0.2;
        TauPreSelection =  TauPreSelection && tauByMVA5TightElectronRejection->at(itau) > 0;
        TauPreSelection =  TauPreSelection && fabs(tauDxy->at(itau)) < 0.05 ;
        if(!(TauPtCut && TauPreSelection)) continue;
        //bool iso_pass = !tauByLooseCombinedIsolationDeltaBetaCorr3Hits->at(itau);
        float iso = tauCombinedIsolationDeltaBetaCorrRaw3Hits->at(itau);
        if(max_iso < iso)
        {
            max_iso = iso;
            idx = itau;
        }

    }
    return idx;
}

int main(int argc, char** argv) {
    using namespace std;

    std::string out = *(argv + 1);

    cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;
    TFile *fout = TFile::Open(out.c_str(), "RECREATE");

    std::string input = *(argv + 2);
    cout << "\n\n\n INPUT NAME IS:    " << input << endl;
    TFile * myFile = TFile::Open(input.c_str());
    TH1F * HistoTot = (TH1F*) myFile->Get("hcount");

    //add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
    TH1F * visibleMassOS = new TH1F ("visibleMassOS","visibleMassOS", 300, 0, 300);
    TH1F * visibleMassSS = new TH1F ("visibleMassSS","visibleMassSS", 300, 0, 300);
    TH1F * visibleMassRatio = new TH1F ("visibleMassRatio","visibleMassRatio", 300, 0, 300);
    TH1F * electronIsolation = new TH1F ("electronIsolation","electronIsolation", 300, 0, .5);

    TTree *Run_Tree = (TTree*) myFile->Get("EventTree"); //Associate branches w/ predeclared variables
    associateTree(Run_Tree);
    cout.setf(ios::fixed, ios::floatfield);
    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    cout << "nentries_wtn====" << nentries_wtn << "\n";
    int eleCount = 0;
    int tauCount = 0;
    for (Int_t i = 0; i < nentries_wtn; i++) {
//        if(i > 5000) break;
        Run_Tree->GetEntry(i);

        if (i % 1000 == 0){
            fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
            fflush(stdout);
        }
        int eleI = electronSelection(electronIsolation);
        int tauI = tauIsoID();
        if(tauI == -1) continue; else tauCount++;
        if(eleI == -1) continue; else eleCount++;
        TLorentzVector tau4Vector, ele4Vector;
        ele4Vector.SetPtEtaPhiE(elePt->at(eleI), eleEta->at(eleI),
                                elePhi->at(eleI), eleEn->at(eleI));
        tau4Vector.SetPtEtaPhiE(tauPt->at(tauI), tauEta->at(tauI),
                                tauPhi->at(tauI), tauEnergy->at(tauI));


        if(tauCharge->at(tauI)*eleCharge->at(eleI) < 0){
            visibleMassOS->Fill((ele4Vector+tau4Vector).M()); 
        }else{
            visibleMassSS->Fill((ele4Vector+tau4Vector).M()); 
        }

    }
    visibleMassRatio->Divide(visibleMassOS, visibleMassSS);
    double fullRatio = visibleMassOS->Integral()/visibleMassSS->Integral();
    cout << "\nTau passing selection" << tauCount << endl;
    cout << "Ele passing selection" << eleCount << endl;
    cout << "OS/SS=" << fullRatio << endl;


    //end of analysis code, close and write histograms/file
    fout->cd();
    visibleMassOS->Write(); 
    visibleMassSS->Write(); 
    visibleMassRatio->Write(); 
    electronIsolation->Write(); 
    fout->Close();

}
