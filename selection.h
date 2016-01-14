#include "TreeReader.h"

//#define SEL_VERBOSE
// Selects proper electron trigger
bool eleTriggerOk(){
    return ((HLTEleMuX >> 6 & 1) == 1  &&  isData ) ||
           ((HLTEleMuX >> 11 & 1) == 1 && !isData ); 
}

// Calculates electron isolation
float eleIsoCalculation(int iele){
    float isoEle=elePFChIso->at(iele)/elePt->at(iele);
    if ((elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5*elePFPUIso->at(iele)) > 0.0)
        isoEle = (elePFChIso->at(iele)/elePt->at(iele) + elePFNeuIso->at(iele)                +
                  elePFPhoIso->at(iele)                - 0.5*elePFPUIso->at(iele))/elePt->at(iele);
    return isoEle;
}

// Cuts out less-sensitive ECAL eta regions
bool eleMVAIdSelection(int iele){
        bool eleMVAId;
        if (fabs(eleSCEta->at(iele)) < 0.8 && eleIDMVANonTrg->at(iele) > 0.967083)
            eleMVAId= true;
        else if (fabs(eleSCEta->at(iele)) > 0.8 && fabs(eleSCEta->at(iele)) < 1.5 && 
                eleIDMVANonTrg->at(iele) > 0.929117)
            eleMVAId= true;
        else if (fabs(eleSCEta->at(iele)) > 1.5 && eleIDMVANonTrg->at(iele) > 0.726311 )
            eleMVAId= true;
        else
            eleMVAId= false; 
        return eleMVAId;
}


int electronSelectionInverse(){
    /* electron_selection
     * Returns the index of the selected electron.
     * If not only a single electron is found, return -1.
     * Cuts out all di-photon events
     */
    int selection = -1;
#ifdef SEL_VERBOSE
    cout << "====================" << endl;
    cout << "\n\nElectrons in event: " << nEle << endl;
#endif
    for(int iele=0; iele<nEle; iele++){
        //Trigger cuts
        bool trigOk = eleTriggerOk();
        //Isolation measurement
        float eleIso = eleIsoCalculation(iele);
        //Region selection
        bool eleMVAId = eleMVAIdSelection(iele);

        bool pTCut  = elePt->at(iele)        > 15.000;
        bool etaCut = fabs(eleEta->at(iele)) <  2.400;
        bool isoCut = eleIso                 >  0.100;
        bool d0Cut  = fabs(eleD0->at(iele))  <  0.045;
        bool dzCut  = fabs(eleDz->at(iele))  < 0.2000;

#ifdef SEL_VERBOSE
        cout << "trg : " << (trigOk    ?"PASS":"FAIL") << endl;
        cout << "mva : " << (eleMVAId  ?"PASS":"FAIL") << endl;
        cout << "Pt  : " << (pTCut     ?"PASS":"FAIL") << "  " << elePt->at(iele)  << endl;
        cout << "Eta : " << (etaCut    ?"PASS":"FAIL") << "  " << eleEta->at(iele) << endl;
        cout << "iso : " << (isoCut    ?"PASS":"FAIL") << "  " << eleIso << endl;
        cout << "d0  : " << (d0Cut     ?"PASS":"FAIL") << "  " << eleD0->at(iele)  << endl;
        cout << "dz  : " << (dzCut     ?"PASS":"FAIL") << "  " << eleDz->at(iele)  << endl;
#endif

        if(trigOk && eleMVAId && pTCut && etaCut && isoCut && d0Cut && dzCut){
            if(selection == -1){
                selection = iele; //Finds first electron
            } else{
                selection == -1; //Finds second electron, so fail!
                break;
            }
        }
    }
#ifdef SEL_VERBOSE
    if(selection != -1) cout <<"SELECTION: " << selection << endl<<endl;
#endif
    return selection;
}

int tauSelectionInverse()
{
    int idx = -1;
//    float max_iso = -999.0;
    for(int itau = 0; itau < nTau; itau++)
    {
        bool TauPtCut = tauPt->at(itau) > 20  && fabs(tauEta->at(itau)) < 2.3 ;
        bool TauPreSelection = tauByLooseMuonRejection3->at(itau)        > 0  && 
                               tauByMVA5TightElectronRejection->at(itau) > 0  &&
                               fabs(tauDxy->at(itau))                    < 0.05;
        //TauPreSelection = TauPreSelection && fabs(tauZImpact->at(itau)) < 0.2;
        if(!(TauPtCut && TauPreSelection)) continue;
        if(tauByLooseCombinedIsolationDeltaBetaCorr3Hits->at(itau)) continue;
//        float iso = tauCombinedIsolationDeltaBetaCorrRaw3Hits->at(itau);
//        if(max_iso < iso)
//        {
//            max_iso = iso;
            idx = itau;
//        }

    }
    return idx;
}

int electronSelection(){
    /* electron_selection
     * Returns the index of the selected electron.
     * If not only a single electron is found, return -1.
     * Cuts out all di-photon events
     */
    int selection = -1;
#ifdef SEL_VERBOSE
    cout << "====================" << endl;
    cout << "\n\nElectrons in event: " << nEle << endl;
#endif
    for(int iele=0; iele<nEle; iele++){
        //Trigger cuts
        bool trigOk = eleTriggerOk();
        //Isolation measurement
        float eleIso = eleIsoCalculation(iele);
        //Region selection
        bool eleMVAId = eleMVAIdSelection(iele);

        bool pTCut  = elePt->at(iele)        > 15.000;
        bool etaCut = fabs(eleEta->at(iele)) <  2.400;
        bool isoCut = eleIso                 <  0.100;
        bool d0Cut  = fabs(eleD0->at(iele))  <  0.045;
        bool dzCut  = fabs(eleDz->at(iele))  < 0.2000;

#ifdef SEL_VERBOSE
        cout << "trg : " << (trigOk    ?"PASS":"FAIL") << endl;
        cout << "mva : " << (eleMVAId  ?"PASS":"FAIL") << endl;
        cout << "Pt  : " << (pTCut     ?"PASS":"FAIL") << "  " << elePt->at(iele)  << endl;
        cout << "Eta : " << (etaCut    ?"PASS":"FAIL") << "  " << eleEta->at(iele) << endl;
        cout << "iso : " << (isoCut    ?"PASS":"FAIL") << "  " << eleIso << endl;
        cout << "d0  : " << (d0Cut     ?"PASS":"FAIL") << "  " << eleD0->at(iele)  << endl;
        cout << "dz  : " << (dzCut     ?"PASS":"FAIL") << "  " << eleDz->at(iele)  << endl;
#endif

        if(trigOk && eleMVAId && pTCut && etaCut && isoCut && d0Cut && dzCut){
            if(selection == -1){
                selection = iele; //Finds first electron
            } else{
                selection == -1; //Finds second electron, so fail!
                break;
            }
        }
    }
#ifdef SEL_VERBOSE
    if(selection != -1) cout <<"SELECTION: " << selection << endl<<endl;
#endif
    return selection;
}

int tauSelection() {
    int idx = -1;
    float min_iso = 999.0;
    for(int itau = 0; itau < nTau; itau++)
    {
        bool TauPtCut = tauPt->at(itau) > 20  && fabs(tauEta->at(itau)) < 2.3 ;
        bool TauPreSelection = tauByLooseMuonRejection3->at(itau) > 0;
        //TauPreSelection = TauPreSelection && fabs(tauZImpact->at(itau)) < 0.2;
        TauPreSelection =  TauPreSelection && tauByMVA5TightElectronRejection->at(itau) > 0;
        TauPreSelection =  TauPreSelection && fabs(tauDxy->at(itau)) < 0.05 ;
        if(!(TauPtCut && TauPreSelection)) continue;
        bool iso_pass = tauByLooseCombinedIsolationDeltaBetaCorr3Hits->at(itau);
        float iso = tauCombinedIsolationDeltaBetaCorrRaw3Hits->at(itau);
        if(iso_pass && min_iso > iso)
        {
            min_iso = iso;
            idx = itau;
        }
    }
    return idx;
}
