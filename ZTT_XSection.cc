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

int main(int argc, char** argv) {

    string out = *(argv + 1);

    cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;
    TFile *fout = TFile::Open(out.c_str(), "RECREATE");

    string input = *(argv + 2);
    cout << "\n\n\n INPUT NAME IS:    " << input << endl;
    TFile * myFile = TFile::Open(input.c_str());
    TH1F * HistoTot = (TH1F*) myFile->Get("hcount");

    //add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
    TH1F * visibleMassOS = new TH1F ("visibleMassOS","visibleMassOS", 300, 0, 300);
    TH1F * visibleMassSS = new TH1F ("visibleMassSS","visibleMassSS", 300, 0, 300);
    TH1F * Mt_h = new TH1F ("Mt_h","Mt_h", 300, 0, 300);
    TH1F * Mt_hSS = new TH1F ("Mt_hSS","Mt_hSS", 300, 0, 300);
    TH1F * tauPt_h = new TH1F ("tauPt_h","tauPt_h", 300, 0, 300);
    TH1F * tauPt_hSS = new TH1F ("tauPt_hSS","tauPt_hSS", 300, 0, 300);
    TH1F * elePt_h = new TH1F ("elePt_h","elePt_h", 300, 0, 300);
    TH1F * elePt_hSS = new TH1F ("elePt_hSS","elePt_hSS", 300, 0, 300);
    TH1F * tauEta_h = new TH1F ("tauEta_h","tauEta_h", 100, -2.5, 2.5);
    TH1F * tauZI_h = new TH1F ("tauZI_h","tauZI_h", 200, -100.0, 100.0);
    TH1F * NPV_h = new TH1F ("NPV_h","NPV_h", 50, 0, 50);

    TFile * PUData= new TFile("Data_Pileup_2015D_1p56fb.root");
    TH1F * HistoPUData= (TH1F *) PUData->Get("pileup");
    HistoPUData->Scale(1.0/HistoPUData->Integral());

    TFile * PUMC= new TFile("MC_Spring15_PU25_Startup.root");
    TH1F * HistoPUMC= (TH1F *) PUMC->Get("pileup");
    HistoPUMC->Scale(1.0/HistoPUMC->Integral());

    float lumiW = weightCalc(HistoTot,out);
    if(lumiW == 1) genWeight = 1.0;

    TTree *Run_Tree = (TTree*) myFile->Get("EventTree"); //Associate branches w/ predeclared variables
    associateTree(Run_Tree);
    cout.setf(ios::fixed, ios::floatfield);
    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    cout << "nentries_wtn====" << nentries_wtn << "\n";

    for (Int_t i = 0; i < nentries_wtn; i++) {
//        if(i > 5000) break;
        Run_Tree->GetEntry(i);

        if (i % 1000 == 0){
            fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
            fflush(stdout);
        }

        float PUWeight = 1;
        if(lumiW != 1) 
        {
            int puNUmmc=int(puTrue->at(0)*10);
            int puNUmdata=int(puTrue->at(0)*10);
            float PUMC_=HistoPUMC->GetBinContent(puNUmmc+1);
            float PUData_=HistoPUData->GetBinContent(puNUmdata+1);
            PUWeight= PUData_/PUMC_;
        }
        int eleI = electronSelection();
        int tauI = tauSelection();
        if(tauI == -1) continue;
        if(eleI == -1) continue;
        if(!transverseMassOk(eleI)) continue;
        TLorentzVector tau4Vector, ele4Vector;
        ele4Vector.SetPtEtaPhiE(elePt->at(eleI), eleEta->at(eleI),
                elePhi->at(eleI), eleEn->at(eleI));
        tau4Vector.SetPtEtaPhiE(tauPt->at(tauI), tauEta->at(tauI),
                tauPhi->at(tauI), tauEnergy->at(tauI));
        float MT = TMass_F(ele4Vector.Pt(),ele4Vector.Px(),ele4Vector.Py(),pfMET,pfMETPhi); 

        if(MT > 60.0) continue;
        

        if(tauCharge->at(tauI)*eleCharge->at(eleI) < 0){
            visibleMassOS->Fill((ele4Vector+tau4Vector).M(),lumiW*PUWeight*genWeight); 
            Mt_h->Fill(MT,lumiW*PUWeight*genWeight);       
            tauPt_h->Fill(tauPt->at(tauI),lumiW*PUWeight*genWeight);       
            elePt_h->Fill(elePt->at(eleI),lumiW*PUWeight*genWeight);       
            tauEta_h->Fill(tauEta->at(tauI),lumiW*PUWeight*genWeight);      
            tauZI_h->Fill(tauZImpact->at(tauI),lumiW*PUWeight*genWeight);       
            NPV_h->Fill(nVtx,lumiW*PUWeight*genWeight); 	      
        }else{
            visibleMassSS->Fill((ele4Vector+tau4Vector).M(),lumiW*PUWeight*genWeight); 
            Mt_hSS->Fill(MT,lumiW*PUWeight*genWeight);       
            tauPt_hSS->Fill(tauPt->at(tauI),lumiW*PUWeight*genWeight);       
            elePt_hSS->Fill(elePt->at(eleI),lumiW*PUWeight*genWeight);       
        }

    }


    //end of analysis code, close and write histograms/file
    fout->cd();
    visibleMassOS->Write(); 
    visibleMassSS->Write(); 
    tauPt_h->Write();       
    Mt_h->Write();       
    Mt_hSS->Write();       
    tauEta_h->Write();      
    tauZI_h->Write();       
    NPV_h->Write(); 	      
    //visibleMassOSRelaxedTauIso->Write();
    //visibleMassSSRelaxedTauIso->Write();
    fout->Close();

}
