////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TreeReader.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "WeightCalculator.h"
#include "selection.h"
#include <string>
#include <ostream>

//#define OS_SS_VERBOSE
bool transverseMassOk(int eleI){
        float MtCut = 40;
        TLorentzVector ele4Vec;
        ele4Vec.SetPtEtaPhiE(elePt->at(eleI), eleEta->at(eleI),
                             elePhi->at(eleI), eleEn->at(eleI));
        float Mt = TMass_F(ele4Vec.Pt(), ele4Vec.Px(), ele4Vec.Py(),
                           pfMET, pfMETPhi);
#ifdef OS_SS_VERBOSE
        cout << "Mt: " << ((Mt < MtCut)?"PASS":"FAIL")
             << " "    << Mt << endl;
#endif
        return Mt < MtCut;
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
    associateTree(Run_Tree, false);
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
        int eleI = electronSelectionInverse();
        int tauI = tauSelectionInverse();
        if(tauI == -1) continue; else tauCount++;
        if(eleI == -1) continue; else eleCount++;
        if(!transverseMassOk(eleI)) continue;
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
    double a = visibleMassOS->Integral();
    double b = visibleMassSS->Integral();
    double a_b = a/b;
    double err = sqrt(1./a+1./b)*a_b;
    cout << "\nTau passing selection " << tauCount << endl;
    cout << "Ele passing selection " << eleCount << endl;
    cout << "OS/SS=" << a_b << "+-" << err << endl;


    //end of analysis code, close and write histograms/file
    fout->cd();
    visibleMassOS->Write(); 
    visibleMassSS->Write(); 
    visibleMassRatio->Write(); 
    electronIsolation->Write(); 
    fout->Close();

}
