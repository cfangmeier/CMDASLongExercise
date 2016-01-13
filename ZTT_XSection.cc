////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TreeReader.h"
#include "WeightCalculator.h"
#include <string>
#include <ostream>

/* electron_selection
 * Returns the index of the selected electron.
 * If not only a single electron is found, return -1.
 * Cuts out all di-photon events
 */
int electronSelection(){
    int selection = -1;
    cout << "====================" << endl;
    cout << "\n\nElectrons in event: " << nEle << endl;
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
        bool passCuts = elePt->at(iele)       > 15    &&
                        eleEta->at(iele)      > 2.4   && 
                        isoEle                < 0.30  &&
                        fabs(eleD0->at(iele)) < 0.045 &&
                        fabs(eleDz->at(iele)) < 0.200;
        cout << "Electron info:" << endl;
        cout << "isoEle:" << isoEle << endl;
        cout << "eleMVAId:" << eleMVAId << endl;
        cout << "passTrigger:" << passTrigger << endl;
        cout << "passCuts:" << passCuts << endl;

        if(passTrigger && eleMVAId && passCuts){
            if(selection == -1){
                cout <<"found first electron" << endl;
                selection = iele; //Finds first electron
            } else{
                cout << "found second electron" << endl;
                selection == -1; //Finds second electron, so fail!
                break;
            }
        }
    }
    return selection;
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
    
    TTree *Run_Tree = (TTree*) myFile->Get("EventTree"); //Associate branches w/ predeclared variables
    associateTree(Run_Tree);
    cout.setf(ios::fixed, ios::floatfield);
    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    cout << "nentries_wtn====" << nentries_wtn << "\n";
    for (Int_t i = 0; i < nentries_wtn; i++) {
        if(i > 3000) break;
        Run_Tree->GetEntry(i);

        if (i % 1000 == 0){
            fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
            fflush(stdout);
        }
        int selected_electron = electronSelection();
        if(selected_electron != -1){
            fprintf(stdout, "electron index: %d\n", selected_electron);
        }
        //Tau Loop

    }


    //end of analysis code, close and write histograms/file
    fout->cd();
    visibleMassOS->Write();
    visibleMassSS->Write();
    //visibleMassOSRelaxedTauIso->Write();
    //visibleMassSSRelaxedTauIso->Write();
    fout->Close();
    
}
