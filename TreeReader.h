#ifndef TREE_READER_H
#define	TREE_READER_H

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "math.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
#include "TLorentzVector.h"

using namespace std;

float TMass_F(float pt3lep, float px3lep, float py3lep, float met, float metPhi) {
    return sqrt(pow(pt3lep + met, 2) - pow(px3lep + met * cos(metPhi), 2) - pow(py3lep + met * sin(metPhi), 2));
}

// Declaration of leaf types
Int_t           run;
Long64_t        event;
Int_t           lumis;
Bool_t          isData;
Int_t           nVtx;
Int_t           nTrksPV;
Float_t         vtx;
Float_t         vty;
Float_t         vtz;
Float_t         rho;
Float_t         rhoCentral;
ULong64_t       HLTEleMuX;
ULong64_t       HLTPho;
ULong64_t       HLTJet;
ULong64_t       HLTEleMuXIsPrescaled;
ULong64_t       HLTPhoIsPrescaled;
ULong64_t       HLTJetIsPrescaled;
Float_t         genWeight;  //MC Only!
vector<float>   *puTrue;    //MC Only!
Int_t           metFilters;
Float_t         pfMET;
Float_t         pfMETPhi;
Float_t         pfMETsumEt;
Float_t         pfMETmEtSig;
Float_t         pfMETSig;
Float_t         pfMET_T1JERUp;
Float_t         pfMET_T1JERDo;
Float_t         pfMET_T1JESUp;
Float_t         pfMET_T1JESDo;
Float_t         pfMET_T1MESUp;
Float_t         pfMET_T1MESDo;
Float_t         pfMET_T1EESUp;
Float_t         pfMET_T1EESDo;
Float_t         pfMET_T1PESUp;
Float_t         pfMET_T1PESDo;
Float_t         pfMET_T1TESUp;
Float_t         pfMET_T1TESDo;
Float_t         pfMET_T1UESUp;
Float_t         pfMET_T1UESDo;
Int_t           nPho;
vector<float>   *phoE;
vector<float>   *phoEt;
vector<float>   *phoEta;
vector<float>   *phoPhi;
vector<float>   *phoSCE;
vector<float>   *phoSCRawE;
vector<float>   *phoESEn;
vector<float>   *phoESEnP1;
vector<float>   *phoESEnP2;
vector<float>   *phoSCEta;
vector<float>   *phoSCPhi;
vector<float>   *phoSCEtaWidth;
vector<float>   *phoSCPhiWidth;
vector<float>   *phoSCBrem;
vector<int>     *phohasPixelSeed;
vector<int>     *phoEleVeto;
vector<float>   *phoR9;
vector<float>   *phoHoverE;
vector<float>   *phoSigmaIEtaIEta;
vector<float>   *phoSigmaIEtaIPhi;
vector<float>   *phoSigmaIPhiIPhi;
vector<float>   *phoE1x3;
vector<float>   *phoE2x2;
vector<float>   *phoE2x5Max;
vector<float>   *phoE5x5;
vector<float>   *phoESEffSigmaRR;
vector<float>   *phoSigmaIEtaIEtaFull5x5;
vector<float>   *phoSigmaIEtaIPhiFull5x5;
vector<float>   *phoSigmaIPhiIPhiFull5x5;
vector<float>   *phoE1x3Full5x5;
vector<float>   *phoE2x2Full5x5;
vector<float>   *phoE2x5MaxFull5x5;
vector<float>   *phoE5x5Full5x5;
vector<float>   *phoR9Full5x5;
vector<float>   *phoSeedBCE;
vector<float>   *phoSeedBCEta;
vector<float>   *phoPFChIso;
vector<float>   *phoPFPhoIso;
vector<float>   *phoPFNeuIso;
vector<float>   *phoPFChWorstIso;
vector<float>   *phoPFChIsoFrix1;
vector<float>   *phoPFChIsoFrix2;
vector<float>   *phoPFChIsoFrix3;
vector<float>   *phoPFChIsoFrix4;
vector<float>   *phoPFChIsoFrix5;
vector<float>   *phoPFChIsoFrix6;
vector<float>   *phoPFChIsoFrix7;
vector<float>   *phoPFChIsoFrix8;
vector<float>   *phoPFPhoIsoFrix1;
vector<float>   *phoPFPhoIsoFrix2;
vector<float>   *phoPFPhoIsoFrix3;
vector<float>   *phoPFPhoIsoFrix4;
vector<float>   *phoPFPhoIsoFrix5;
vector<float>   *phoPFPhoIsoFrix6;
vector<float>   *phoPFPhoIsoFrix7;
vector<float>   *phoPFPhoIsoFrix8;
vector<float>   *phoPFNeuIsoFrix1;
vector<float>   *phoPFNeuIsoFrix2;
vector<float>   *phoPFNeuIsoFrix3;
vector<float>   *phoPFNeuIsoFrix4;
vector<float>   *phoPFNeuIsoFrix5;
vector<float>   *phoPFNeuIsoFrix6;
vector<float>   *phoPFNeuIsoFrix7;
vector<float>   *phoPFNeuIsoFrix8;
vector<float>   *phoEcalRecHitSumEtConeDR03;
vector<float>   *phohcalDepth1TowerSumEtConeDR03;
vector<float>   *phohcalDepth2TowerSumEtConeDR03;
vector<float>   *phohcalTowerSumEtConeDR03;
vector<float>   *photrkSumPtHollowConeDR03;
vector<float>   *phoIDMVA;
vector<int>     *phoFiredSingleTrgs;
vector<int>     *phoFiredDoubleTrgs;
vector<unsigned short> *phoIDbit;
Int_t           nEle;
vector<int>     *eleCharge;
vector<int>     *eleChargeConsistent;
vector<float>   *eleEn;
vector<float>   *eleSCEn;
vector<float>   *eleESEn;
vector<float>   *eleESEnP1;
vector<float>   *eleESEnP2;
vector<float>   *eleD0;
vector<float>   *eleDz;
vector<float>   *elePt;
vector<float>   *eleEta;
vector<float>   *elePhi;
vector<float>   *eleR9;
vector<float>   *eleSCEta;
vector<float>   *eleSCPhi;
vector<float>   *eleSCRawEn;
vector<float>   *eleSCEtaWidth;
vector<float>   *eleSCPhiWidth;
vector<float>   *eleHoverE;
vector<float>   *eleEoverP;
vector<float>   *eleEoverPout;
vector<float>   *eleEoverPInv;
vector<float>   *eleBrem;
vector<float>   *eledEtaAtVtx;
vector<float>   *eledPhiAtVtx;
vector<float>   *eledEtaAtCalo;
vector<float>   *eleSigmaIEtaIEta;
vector<float>   *eleSigmaIEtaIPhi;
vector<float>   *eleSigmaIPhiIPhi;
vector<float>   *eleSigmaIEtaIEtaFull5x5;
vector<float>   *eleSigmaIPhiIPhiFull5x5;
vector<int>     *eleConvVeto;
vector<int>     *eleMissHits;
vector<float>   *eleESEffSigmaRR;
vector<float>   *elePFChIso;
vector<float>   *elePFPhoIso;
vector<float>   *elePFNeuIso;
vector<float>   *elePFPUIso;
vector<float>   *eleIDMVANonTrg;
vector<float>   *eleIDMVATrg;
vector<float>   *eledEtaseedAtVtx;
vector<float>   *eleE1x5;
vector<float>   *eleE2x5;
vector<float>   *eleE5x5;
vector<float>   *eleE1x5Full5x5;
vector<float>   *eleE2x5Full5x5;
vector<float>   *eleE5x5Full5x5;
vector<float>   *eleR9Full5x5;
vector<int>     *eleEcalDrivenSeed;
vector<float>   *eleDr03EcalRecHitSumEt;
vector<float>   *eleDr03HcalDepth1TowerSumEt;
vector<float>   *eleDr03HcalDepth2TowerSumEt;
vector<float>   *eleDr03HcalTowerSumEt;
vector<float>   *eleDr03TkSumPt;
vector<float>   *elecaloEnergy;
vector<float>   *eleTrkdxy;
vector<float>   *eleKFHits;
vector<float>   *eleKFChi2;
vector<vector<float> > *eleGSFPt;
vector<vector<float> > *eleGSFEta;
vector<vector<float> > *eleGSFPhi;
vector<vector<float> > *eleGSFCharge;
vector<vector<int> > *eleGSFHits;
vector<vector<int> > *eleGSFMissHits;
vector<vector<int> > *eleGSFNHitsMax;
vector<vector<float> > *eleGSFVtxProb;
vector<vector<float> > *eleGSFlxyPV;
vector<vector<float> > *eleGSFlxyBS;
vector<vector<float> > *eleBCEn;
vector<vector<float> > *eleBCEta;
vector<vector<float> > *eleBCPhi;
vector<vector<float> > *eleBCS25;
vector<vector<float> > *eleBCS15;
vector<vector<float> > *eleBCSieie;
vector<vector<float> > *eleBCSieip;
vector<vector<float> > *eleBCSipip;
vector<int>     *eleFiredTrgs;
vector<unsigned short> *eleIDbit;
Int_t           nMu;
vector<float>   *muPt;
vector<float>   *muEn;
vector<float>   *muEta;
vector<float>   *muPhi;
vector<int>     *muCharge;
vector<int>     *muType;
vector<bool>    *muIsLooseID;
vector<bool>    *muIsMediumID;
vector<bool>    *muIsTightID;
vector<bool>    *muIsSoftID;
vector<bool>    *muIsHighPtID;
vector<float>   *muD0;
vector<float>   *muDz;
vector<float>   *muChi2NDF;
vector<float>   *muInnerD0;
vector<float>   *muInnerDz;
vector<int>     *muTrkLayers;
vector<int>     *muPixelLayers;
vector<int>     *muPixelHits;
vector<int>     *muMuonHits;
vector<int>     *muStations;
vector<int>     *muTrkQuality;
vector<float>   *muIsoTrk;
vector<float>   *muPFChIso;
vector<float>   *muPFPhoIso;
vector<float>   *muPFNeuIso;
vector<float>   *muPFPUIso;
vector<int>     *muFiredTrgs;
vector<float>   *muInnervalidFraction;
vector<float>   *musegmentCompatibility;
vector<float>   *muchi2LocalPosition;
vector<float>   *mutrkKink;
vector<float>   *muBestTrkPtError;
vector<float>   *muBestTrkPt;
Int_t           nTau;
vector<bool>    *pfTausDiscriminationByDecayModeFinding;
vector<bool>    *pfTausDiscriminationByDecayModeFindingNewDMs;
vector<bool>    *tauByMVA5LooseElectronRejection;
vector<bool>    *tauByMVA5MediumElectronRejection;
vector<bool>    *tauByMVA5TightElectronRejection;
vector<bool>    *tauByMVA5VTightElectronRejection;
vector<bool>    *tauByLooseMuonRejection3;
vector<bool>    *tauByTightMuonRejection3;
vector<bool>    *tauByLooseCombinedIsolationDeltaBetaCorr3Hits;
vector<bool>    *tauByMediumCombinedIsolationDeltaBetaCorr3Hits;
vector<bool>    *tauByTightCombinedIsolationDeltaBetaCorr3Hits;
vector<float>   *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
vector<bool>    *tauByVLooseIsolationMVA3oldDMwLT;
vector<bool>    *tauByLooseIsolationMVA3oldDMwLT;
vector<bool>    *tauByMediumIsolationMVA3oldDMwLT;
vector<bool>    *tauByTightIsolationMVA3oldDMwLT;
vector<bool>    *tauByVTightIsolationMVA3oldDMwLT;
vector<bool>    *tauByVVTightIsolationMVA3oldDMwLT;
vector<float>   *tauByIsolationMVA3oldDMwLTraw;
vector<bool>    *tauByLooseIsolationMVA3newDMwLT;
vector<bool>    *tauByVLooseIsolationMVA3newDMwLT;
vector<bool>    *tauByMediumIsolationMVA3newDMwLT;
vector<bool>    *tauByTightIsolationMVA3newDMwLT;
vector<bool>    *tauByVTightIsolationMVA3newDMwLT;
vector<bool>    *tauByVVTightIsolationMVA3newDMwLT;
vector<float>   *tauByIsolationMVA3newDMwLTraw;
vector<float>   *tauEta;
vector<float>   *tauPhi;
vector<float>   *tauPt;
vector<float>   *tauEt;
vector<float>   *tauCharge;
vector<float>   *tauP;
vector<float>   *tauPx;
vector<float>   *tauPy;
vector<float>   *tauPz;
vector<float>   *tauVz;
vector<float>   *tauEnergy;
vector<float>   *tauMass;
vector<float>   *tauDxy;
vector<float>   *tauZImpact;
vector<int>     *tauDecayMode;
vector<bool>    *tauLeadChargedHadronExists;
vector<float>   *tauLeadChargedHadronEta;
vector<float>   *tauLeadChargedHadronPhi;
vector<float>   *tauLeadChargedHadronPt;
vector<float>   *tauChargedIsoPtSum;
vector<float>   *tauNeutralIsoPtSum;
vector<float>   *tauPuCorrPtSum;
vector<float>   *tauNumSignalPFChargedHadrCands;
vector<float>   *tauNumSignalPFNeutrHadrCands;
vector<float>   *tauNumSignalPFGammaCands;
vector<float>   *tauNumSignalPFCands;
vector<float>   *tauNumIsolationPFChargedHadrCands;
vector<float>   *tauNumIsolationPFNeutrHadrCands;
vector<float>   *tauNumIsolationPFGammaCands;
vector<float>   *tauNumIsolationPFCands;
Int_t           nJet;
vector<float>   *jetPt;
vector<float>   *jetEn;
vector<float>   *jetEta;
vector<float>   *jetPhi;
vector<float>   *jetRawPt;
vector<float>   *jetRawEn;
vector<float>   *jetArea;
vector<float>   *jetpfCombinedInclusiveSecondaryVertexV2BJetTags;
vector<float>   *jetJetProbabilityBJetTags;
vector<float>   *jetpfCombinedMVABJetTags;
vector<bool>    *jetPFLooseId;
vector<float>   *jetPUidFullDiscriminant;
vector<float>   *jetJECUnc;
vector<int>     *jetFiredTrgs;
Int_t           nAK8Jet;
vector<float>   *AK8JetPt;
vector<float>   *AK8JetEn;
vector<float>   *AK8JetRawPt;
vector<float>   *AK8JetRawEn;
vector<float>   *AK8JetEta;
vector<float>   *AK8JetPhi;
vector<float>   *AK8JetMass;
vector<float>   *AK8Jet_tau1;
vector<float>   *AK8Jet_tau2;
vector<float>   *AK8Jet_tau3;
vector<float>   *AK8JetCHF;
vector<float>   *AK8JetNHF;
vector<float>   *AK8JetCEF;
vector<float>   *AK8JetNEF;
vector<int>     *AK8JetNCH;
vector<int>     *AK8Jetnconstituents;
vector<bool>    *AK8JetPFLooseId;
vector<float>   *AK8CHSSoftDropJetMass;
vector<float>   *AK8JetpfBoostedDSVBTag;
vector<float>   *AK8JetJECUnc;
vector<int>     *nAK8softdropSubjet;
vector<vector<float> > *AK8softdropSubjetPt;
vector<vector<float> > *AK8softdropSubjetEta;
vector<vector<float> > *AK8softdropSubjetPhi;
vector<vector<float> > *AK8softdropSubjetMass;
vector<vector<float> > *AK8softdropSubjetE;
vector<vector<int> > *AK8softdropSubjetCharge;
vector<vector<int> > *AK8softdropSubjetFlavour;
vector<vector<float> > *AK8softdropSubjetCSV;

   // List of branches
TBranch        *b_run;   //!
TBranch        *b_event;   //!
TBranch        *b_lumis;   //!
TBranch        *b_isData;   //!
TBranch        *b_nVtx;   //!
TBranch        *b_nTrksPV;   //!
TBranch        *b_vtx;   //!
TBranch        *b_vty;   //!
TBranch        *b_vtz;   //!
TBranch        *b_rho;   //!
TBranch        *b_rhoCentral;   //!
TBranch        *b_HLTEleMuX;   //!
TBranch        *b_HLTPho;   //!
TBranch        *b_HLTJet;   //!
TBranch        *b_HLTEleMuXIsPrescaled;   //!
TBranch        *b_HLTPhoIsPrescaled;   //!
TBranch        *b_HLTJetIsPrescaled;   //!
TBranch        *b_genWeight;  //MC Only!
TBranch        *b_puTrue;    //MC Only!
TBranch        *b_metFilters;   //!
TBranch        *b_pfMET;   //!
TBranch        *b_pfMETPhi;   //!
TBranch        *b_pfMETsumEt;   //!
TBranch        *b_pfMETmEtSig;   //!
TBranch        *b_pfMETSig;   //!
TBranch        *b_pfMET_T1JERUp;   //!
TBranch        *b_pfMET_T1JERDo;   //!
TBranch        *b_pfMET_T1JESUp;   //!
TBranch        *b_pfMET_T1JESDo;   //!
TBranch        *b_pfMET_T1MESUp;   //!
TBranch        *b_pfMET_T1MESDo;   //!
TBranch        *b_pfMET_T1EESUp;   //!
TBranch        *b_pfMET_T1EESDo;   //!
TBranch        *b_pfMET_T1PESUp;   //!
TBranch        *b_pfMET_T1PESDo;   //!
TBranch        *b_pfMET_T1TESUp;   //!
TBranch        *b_pfMET_T1TESDo;   //!
TBranch        *b_pfMET_T1UESUp;   //!
TBranch        *b_pfMET_T1UESDo;   //!
TBranch        *b_nPho;   //!
TBranch        *b_phoE;   //!
TBranch        *b_phoEt;   //!
TBranch        *b_phoEta;   //!
TBranch        *b_phoPhi;   //!
TBranch        *b_phoSCE;   //!
TBranch        *b_phoSCRawE;   //!
TBranch        *b_phoESEn;   //!
TBranch        *b_phoESEnP1;   //!
TBranch        *b_phoESEnP2;   //!
TBranch        *b_phoSCEta;   //!
TBranch        *b_phoSCPhi;   //!
TBranch        *b_phoSCEtaWidth;   //!
TBranch        *b_phoSCPhiWidth;   //!
TBranch        *b_phoSCBrem;   //!
TBranch        *b_phohasPixelSeed;   //!
TBranch        *b_phoEleVeto;   //!
TBranch        *b_phoR9;   //!
TBranch        *b_phoHoverE;   //!
TBranch        *b_phoSigmaIEtaIEta;   //!
TBranch        *b_phoSigmaIEtaIPhi;   //!
TBranch        *b_phoSigmaIPhiIPhi;   //!
TBranch        *b_phoE1x3;   //!
TBranch        *b_phoE2x2;   //!
TBranch        *b_phoE2x5Max;   //!
TBranch        *b_phoE5x5;   //!
TBranch        *b_phoESEffSigmaRR;   //!
TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
TBranch        *b_phoE1x3Full5x5;   //!
TBranch        *b_phoE2x2Full5x5;   //!
TBranch        *b_phoE2x5MaxFull5x5;   //!
TBranch        *b_phoE5x5Full5x5;   //!
TBranch        *b_phoR9Full5x5;   //!
TBranch        *b_phoSeedBCE;   //!
TBranch        *b_phoSeedBCEta;   //!
TBranch        *b_phoPFChIso;   //!
TBranch        *b_phoPFPhoIso;   //!
TBranch        *b_phoPFNeuIso;   //!
TBranch        *b_phoPFChWorstIso;   //!
TBranch        *b_phoPFChIsoFrix1;   //!
TBranch        *b_phoPFChIsoFrix2;   //!
TBranch        *b_phoPFChIsoFrix3;   //!
TBranch        *b_phoPFChIsoFrix4;   //!
TBranch        *b_phoPFChIsoFrix5;   //!
TBranch        *b_phoPFChIsoFrix6;   //!
TBranch        *b_phoPFChIsoFrix7;   //!
TBranch        *b_phoPFChIsoFrix8;   //!
TBranch        *b_phoPFPhoIsoFrix1;   //!
TBranch        *b_phoPFPhoIsoFrix2;   //!
TBranch        *b_phoPFPhoIsoFrix3;   //!
TBranch        *b_phoPFPhoIsoFrix4;   //!
TBranch        *b_phoPFPhoIsoFrix5;   //!
TBranch        *b_phoPFPhoIsoFrix6;   //!
TBranch        *b_phoPFPhoIsoFrix7;   //!
TBranch        *b_phoPFPhoIsoFrix8;   //!
TBranch        *b_phoPFNeuIsoFrix1;   //!
TBranch        *b_phoPFNeuIsoFrix2;   //!
TBranch        *b_phoPFNeuIsoFrix3;   //!
TBranch        *b_phoPFNeuIsoFrix4;   //!
TBranch        *b_phoPFNeuIsoFrix5;   //!
TBranch        *b_phoPFNeuIsoFrix6;   //!
TBranch        *b_phoPFNeuIsoFrix7;   //!
TBranch        *b_phoPFNeuIsoFrix8;   //!
TBranch        *b_phoEcalRecHitSumEtConeDR03;   //!
TBranch        *b_phohcalDepth1TowerSumEtConeDR03;   //!
TBranch        *b_phohcalDepth2TowerSumEtConeDR03;   //!
TBranch        *b_phohcalTowerSumEtConeDR03;   //!
TBranch        *b_photrkSumPtHollowConeDR03;   //!
TBranch        *b_phoIDMVA;   //!
TBranch        *b_phoFiredSingleTrgs;   //!
TBranch        *b_phoFiredDoubleTrgs;   //!
TBranch        *b_phoIDbit;   //!
TBranch        *b_nEle;   //!
TBranch        *b_eleCharge;   //!
TBranch        *b_eleChargeConsistent;   //!
TBranch        *b_eleEn;   //!
TBranch        *b_eleSCEn;   //!
TBranch        *b_eleESEn;   //!
TBranch        *b_eleESEnP1;   //!
TBranch        *b_eleESEnP2;   //!
TBranch        *b_eleD0;   //!
TBranch        *b_eleDz;   //!
TBranch        *b_elePt;   //!
TBranch        *b_eleEta;   //!
TBranch        *b_elePhi;   //!
TBranch        *b_eleR9;   //!
TBranch        *b_eleSCEta;   //!
TBranch        *b_eleSCPhi;   //!
TBranch        *b_eleSCRawEn;   //!
TBranch        *b_eleSCEtaWidth;   //!
TBranch        *b_eleSCPhiWidth;   //!
TBranch        *b_eleHoverE;   //!
TBranch        *b_eleEoverP;   //!
TBranch        *b_eleEoverPout;   //!
TBranch        *b_eleEoverPInv;   //!
TBranch        *b_eleBrem;   //!
TBranch        *b_eledEtaAtVtx;   //!
TBranch        *b_eledPhiAtVtx;   //!
TBranch        *b_eledEtaAtCalo;   //!
TBranch        *b_eleSigmaIEtaIEta;   //!
TBranch        *b_eleSigmaIEtaIPhi;   //!
TBranch        *b_eleSigmaIPhiIPhi;   //!
TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
TBranch        *b_eleConvVeto;   //!
TBranch        *b_eleMissHits;   //!
TBranch        *b_eleESEffSigmaRR;   //!
TBranch        *b_elePFChIso;   //!
TBranch        *b_elePFPhoIso;   //!
TBranch        *b_elePFNeuIso;   //!
TBranch        *b_elePFPUIso;   //!
TBranch        *b_eleIDMVANonTrg;   //!
TBranch        *b_eleIDMVATrg;   //!
TBranch        *b_eledEtaseedAtVtx;   //!
TBranch        *b_eleE1x5;   //!
TBranch        *b_eleE2x5;   //!
TBranch        *b_eleE5x5;   //!
TBranch        *b_eleE1x5Full5x5;   //!
TBranch        *b_eleE2x5Full5x5;   //!
TBranch        *b_eleE5x5Full5x5;   //!
TBranch        *b_eleR9Full5x5;   //!
TBranch        *b_eleEcalDrivenSeed;   //!
TBranch        *b_eleDr03EcalRecHitSumEt;   //!
TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
TBranch        *b_eleDr03HcalTowerSumEt;   //!
TBranch        *b_eleDr03TkSumPt;   //!
TBranch        *b_elecaloEnergy;   //!
TBranch        *b_eleTrkdxy;   //!
TBranch        *b_eleKFHits;   //!
TBranch        *b_eleKFChi2;   //!
TBranch        *b_eleGSFPt;   //!
TBranch        *b_eleGSFEta;   //!
TBranch        *b_eleGSFPhi;   //!
TBranch        *b_eleGSFCharge;   //!
TBranch        *b_eleGSFHits;   //!
TBranch        *b_eleGSFMissHits;   //!
TBranch        *b_eleGSFNHitsMax;   //!
TBranch        *b_eleGSFVtxProb;   //!
TBranch        *b_eleGSFlxyPV;   //!
TBranch        *b_eleGSFlxyBS;   //!
TBranch        *b_eleBCEn;   //!
TBranch        *b_eleBCEta;   //!
TBranch        *b_eleBCPhi;   //!
TBranch        *b_eleBCS25;   //!
TBranch        *b_eleBCS15;   //!
TBranch        *b_eleBCSieie;   //!
TBranch        *b_eleBCSieip;   //!
TBranch        *b_eleBCSipip;   //!
TBranch        *b_eleFiredTrgs;   //!
TBranch        *b_eleIDbit;   //!
TBranch        *b_nMu;   //!
TBranch        *b_muPt;   //!
TBranch        *b_muEn;   //!
TBranch        *b_muEta;   //!
TBranch        *b_muPhi;   //!
TBranch        *b_muCharge;   //!
TBranch        *b_muType;   //!
TBranch        *b_muIsLooseID;   //!
TBranch        *b_muIsMediumID;   //!
TBranch        *b_muIsTightID;   //!
TBranch        *b_muIsSoftID;   //!
TBranch        *b_muIsHighPtID;   //!
TBranch        *b_muD0;   //!
TBranch        *b_muDz;   //!
TBranch        *b_muChi2NDF;   //!
TBranch        *b_muInnerD0;   //!
TBranch        *b_muInnerDz;   //!
TBranch        *b_muTrkLayers;   //!
TBranch        *b_muPixelLayers;   //!
TBranch        *b_muPixelHits;   //!
TBranch        *b_muMuonHits;   //!
TBranch        *b_muStations;   //!
TBranch        *b_muTrkQuality;   //!
TBranch        *b_muIsoTrk;   //!
TBranch        *b_muPFChIso;   //!
TBranch        *b_muPFPhoIso;   //!
TBranch        *b_muPFNeuIso;   //!
TBranch        *b_muPFPUIso;   //!
TBranch        *b_muFiredTrgs;   //!
TBranch        *b_muInnervalidFraction;   //!
TBranch        *b_musegmentCompatibility;   //!
TBranch        *b_muchi2LocalPosition;   //!
TBranch        *b_mutrkKink;   //!
TBranch        *b_muBestTrkPtError;   //!
TBranch        *b_muBestTrkPt;   //!
TBranch        *b_nTau;   //!
TBranch        *b_pfTausDiscriminationByDecayModeFinding;   //!
TBranch        *b_pfTausDiscriminationByDecayModeFindingNewDMs;   //!
TBranch        *b_tauByMVA5LooseElectronRejection;   //!
TBranch        *b_tauByMVA5MediumElectronRejection;   //!
TBranch        *b_tauByMVA5TightElectronRejection;   //!
TBranch        *b_tauByMVA5VTightElectronRejection;   //!
TBranch        *b_tauByLooseMuonRejection3;   //!
TBranch        *b_tauByTightMuonRejection3;   //!
TBranch        *b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
TBranch        *b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
TBranch        *b_tauByTightCombinedIsolationDeltaBetaCorr3Hits;   //!
TBranch        *b_tauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
TBranch        *b_tauByVLooseIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByLooseIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByMediumIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByTightIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByVTightIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByVVTightIsolationMVA3oldDMwLT;   //!
TBranch        *b_tauByIsolationMVA3oldDMwLTraw;   //!
TBranch        *b_tauByLooseIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByVLooseIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByMediumIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByTightIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByVTightIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByVVTightIsolationMVA3newDMwLT;   //!
TBranch        *b_tauByIsolationMVA3newDMwLTraw;   //!
TBranch        *b_tauEta;   //!
TBranch        *b_tauPhi;   //!
TBranch        *b_tauPt;   //!
TBranch        *b_tauEt;   //!
TBranch        *b_tauCharge;   //!
TBranch        *b_tauP;   //!
TBranch        *b_tauPx;   //!
TBranch        *b_tauPy;   //!
TBranch        *b_tauPz;   //!
TBranch        *b_tauVz;   //!
TBranch        *b_tauEnergy;   //!
TBranch        *b_tauMass;   //!
TBranch        *b_tauDxy;   //!
TBranch        *b_tauZImpact;   //!
TBranch        *b_tauDecayMode;   //!
TBranch        *b_tauLeadChargedHadronExists;   //!
TBranch        *b_tauLeadChargedHadronEta;   //!
TBranch        *b_tauLeadChargedHadronPhi;   //!
TBranch        *b_tauLeadChargedHadronPt;   //!
TBranch        *b_tauChargedIsoPtSum;   //!
TBranch        *b_tauNeutralIsoPtSum;   //!
TBranch        *b_tauPuCorrPtSum;   //!
TBranch        *b_tauNumSignalPFChargedHadrCands;   //!
TBranch        *b_tauNumSignalPFNeutrHadrCands;   //!
TBranch        *b_tauNumSignalPFGammaCands;   //!
TBranch        *b_tauNumSignalPFCands;   //!
TBranch        *b_tauNumIsolationPFChargedHadrCands;   //!
TBranch        *b_tauNumIsolationPFNeutrHadrCands;   //!
TBranch        *b_tauNumIsolationPFGammaCands;   //!
TBranch        *b_tauNumIsolationPFCands;   //!
TBranch        *b_nJet;   //!
TBranch        *b_jetPt;   //!
TBranch        *b_jetEn;   //!
TBranch        *b_jetEta;   //!
TBranch        *b_jetPhi;   //!
TBranch        *b_jetRawPt;   //!
TBranch        *b_jetRawEn;   //!
TBranch        *b_jetArea;   //!
TBranch        *b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
TBranch        *b_jetJetProbabilityBJetTags;   //!
TBranch        *b_jetpfCombinedMVABJetTags;   //!
TBranch        *b_jetPFLooseId;   //!
TBranch        *b_jetPUidFullDiscriminant;   //!
TBranch        *b_jetJECUnc;   //!
TBranch        *b_jetFiredTrgs;   //!
TBranch        *b_nAK8Jet;   //!
TBranch        *b_AK8JetPt;   //!
TBranch        *b_AK8JetEn;   //!
TBranch        *b_AK8JetRawPt;   //!
TBranch        *b_AK8JetRawEn;   //!
TBranch        *b_AK8JetEta;   //!
TBranch        *b_AK8JetPhi;   //!
TBranch        *b_AK8JetMass;   //!
TBranch        *b_AK8Jet_tau1;   //!
TBranch        *b_AK8Jet_tau2;   //!
TBranch        *b_AK8Jet_tau3;   //!
TBranch        *b_AK8JetCHF;   //!
TBranch        *b_AK8JetNHF;   //!
TBranch        *b_AK8JetCEF;   //!
TBranch        *b_AK8JetNEF;   //!
TBranch        *b_AK8JetNCH;   //!
TBranch        *b_AK8Jetnconstituents;   //!
TBranch        *b_AK8JetPFLooseId;   //!
TBranch        *b_AK8CHSSoftDropJetMass;   //!
TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
TBranch        *b_AK8JetJECUnc;   //!
TBranch        *b_nAK8softdropSubjet;   //!
TBranch        *b_AK8softdropSubjetPt;   //!
TBranch        *b_AK8softdropSubjetEta;   //!
TBranch        *b_AK8softdropSubjetPhi;   //!
TBranch        *b_AK8softdropSubjetMass;   //!
TBranch        *b_AK8softdropSubjetE;   //!
TBranch        *b_AK8softdropSubjetCharge;   //!
TBranch        *b_AK8softdropSubjetFlavour;   //!
TBranch        *b_AK8softdropSubjetCSV;   //!

void associateTree(TTree *tree, bool isMonteCarlo = true){
   tree->SetBranchAddress("run", &run, &b_run);
   tree->SetBranchAddress("event", &event, &b_event);
   tree->SetBranchAddress("lumis", &lumis, &b_lumis);
   tree->SetBranchAddress("isData", &isData, &b_isData);
   tree->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   tree->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   tree->SetBranchAddress("vtx", &vtx, &b_vtx);
   tree->SetBranchAddress("vty", &vty, &b_vty);
   tree->SetBranchAddress("vtz", &vtz, &b_vtz);
   tree->SetBranchAddress("rho", &rho, &b_rho);
   tree->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   tree->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   tree->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   tree->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   tree->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   tree->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   if(isMonteCarlo){
       tree->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
       tree->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   }
   tree->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   tree->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   tree->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   tree->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   tree->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   tree->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   tree->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   tree->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   tree->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   tree->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   tree->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   tree->SetBranchAddress("pfMET_T1MESUp", &pfMET_T1MESUp, &b_pfMET_T1MESUp);
   tree->SetBranchAddress("pfMET_T1MESDo", &pfMET_T1MESDo, &b_pfMET_T1MESDo);
   tree->SetBranchAddress("pfMET_T1EESUp", &pfMET_T1EESUp, &b_pfMET_T1EESUp);
   tree->SetBranchAddress("pfMET_T1EESDo", &pfMET_T1EESDo, &b_pfMET_T1EESDo);
   tree->SetBranchAddress("pfMET_T1PESUp", &pfMET_T1PESUp, &b_pfMET_T1PESUp);
   tree->SetBranchAddress("pfMET_T1PESDo", &pfMET_T1PESDo, &b_pfMET_T1PESDo);
   tree->SetBranchAddress("pfMET_T1TESUp", &pfMET_T1TESUp, &b_pfMET_T1TESUp);
   tree->SetBranchAddress("pfMET_T1TESDo", &pfMET_T1TESDo, &b_pfMET_T1TESDo);
   tree->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   tree->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   tree->SetBranchAddress("nPho", &nPho, &b_nPho);
   tree->SetBranchAddress("phoE", &phoE, &b_phoE);
   tree->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   tree->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   tree->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   tree->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   tree->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   tree->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   tree->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   tree->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   tree->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   tree->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   tree->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   tree->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   tree->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   tree->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   tree->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   tree->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   tree->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   tree->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   tree->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   tree->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   tree->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   tree->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   tree->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   tree->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   tree->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   tree->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   tree->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   tree->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   tree->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   tree->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   tree->SetBranchAddress("phoE2x5MaxFull5x5", &phoE2x5MaxFull5x5, &b_phoE2x5MaxFull5x5);
   tree->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   tree->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   tree->SetBranchAddress("phoSeedBCE", &phoSeedBCE, &b_phoSeedBCE);
   tree->SetBranchAddress("phoSeedBCEta", &phoSeedBCEta, &b_phoSeedBCEta);
   tree->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   tree->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   tree->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   tree->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   tree->SetBranchAddress("phoPFChIsoFrix1", &phoPFChIsoFrix1, &b_phoPFChIsoFrix1);
   tree->SetBranchAddress("phoPFChIsoFrix2", &phoPFChIsoFrix2, &b_phoPFChIsoFrix2);
   tree->SetBranchAddress("phoPFChIsoFrix3", &phoPFChIsoFrix3, &b_phoPFChIsoFrix3);
   tree->SetBranchAddress("phoPFChIsoFrix4", &phoPFChIsoFrix4, &b_phoPFChIsoFrix4);
   tree->SetBranchAddress("phoPFChIsoFrix5", &phoPFChIsoFrix5, &b_phoPFChIsoFrix5);
   tree->SetBranchAddress("phoPFChIsoFrix6", &phoPFChIsoFrix6, &b_phoPFChIsoFrix6);
   tree->SetBranchAddress("phoPFChIsoFrix7", &phoPFChIsoFrix7, &b_phoPFChIsoFrix7);
   tree->SetBranchAddress("phoPFChIsoFrix8", &phoPFChIsoFrix8, &b_phoPFChIsoFrix8);
   tree->SetBranchAddress("phoPFPhoIsoFrix1", &phoPFPhoIsoFrix1, &b_phoPFPhoIsoFrix1);
   tree->SetBranchAddress("phoPFPhoIsoFrix2", &phoPFPhoIsoFrix2, &b_phoPFPhoIsoFrix2);
   tree->SetBranchAddress("phoPFPhoIsoFrix3", &phoPFPhoIsoFrix3, &b_phoPFPhoIsoFrix3);
   tree->SetBranchAddress("phoPFPhoIsoFrix4", &phoPFPhoIsoFrix4, &b_phoPFPhoIsoFrix4);
   tree->SetBranchAddress("phoPFPhoIsoFrix5", &phoPFPhoIsoFrix5, &b_phoPFPhoIsoFrix5);
   tree->SetBranchAddress("phoPFPhoIsoFrix6", &phoPFPhoIsoFrix6, &b_phoPFPhoIsoFrix6);
   tree->SetBranchAddress("phoPFPhoIsoFrix7", &phoPFPhoIsoFrix7, &b_phoPFPhoIsoFrix7);
   tree->SetBranchAddress("phoPFPhoIsoFrix8", &phoPFPhoIsoFrix8, &b_phoPFPhoIsoFrix8);
   tree->SetBranchAddress("phoPFNeuIsoFrix1", &phoPFNeuIsoFrix1, &b_phoPFNeuIsoFrix1);
   tree->SetBranchAddress("phoPFNeuIsoFrix2", &phoPFNeuIsoFrix2, &b_phoPFNeuIsoFrix2);
   tree->SetBranchAddress("phoPFNeuIsoFrix3", &phoPFNeuIsoFrix3, &b_phoPFNeuIsoFrix3);
   tree->SetBranchAddress("phoPFNeuIsoFrix4", &phoPFNeuIsoFrix4, &b_phoPFNeuIsoFrix4);
   tree->SetBranchAddress("phoPFNeuIsoFrix5", &phoPFNeuIsoFrix5, &b_phoPFNeuIsoFrix5);
   tree->SetBranchAddress("phoPFNeuIsoFrix6", &phoPFNeuIsoFrix6, &b_phoPFNeuIsoFrix6);
   tree->SetBranchAddress("phoPFNeuIsoFrix7", &phoPFNeuIsoFrix7, &b_phoPFNeuIsoFrix7);
   tree->SetBranchAddress("phoPFNeuIsoFrix8", &phoPFNeuIsoFrix8, &b_phoPFNeuIsoFrix8);
   tree->SetBranchAddress("phoEcalRecHitSumEtConeDR03", &phoEcalRecHitSumEtConeDR03, &b_phoEcalRecHitSumEtConeDR03);
   tree->SetBranchAddress("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03, &b_phohcalDepth1TowerSumEtConeDR03);
   tree->SetBranchAddress("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03, &b_phohcalDepth2TowerSumEtConeDR03);
   tree->SetBranchAddress("phohcalTowerSumEtConeDR03", &phohcalTowerSumEtConeDR03, &b_phohcalTowerSumEtConeDR03);
   tree->SetBranchAddress("photrkSumPtHollowConeDR03", &photrkSumPtHollowConeDR03, &b_photrkSumPtHollowConeDR03);
   tree->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   tree->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   tree->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   tree->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   tree->SetBranchAddress("nEle", &nEle, &b_nEle);
   tree->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   tree->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   tree->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   tree->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   tree->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   tree->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   tree->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   tree->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   tree->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   tree->SetBranchAddress("elePt", &elePt, &b_elePt);
   tree->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   tree->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   tree->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   tree->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   tree->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   tree->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   tree->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   tree->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   tree->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   tree->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   tree->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   tree->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   tree->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   tree->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   tree->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   tree->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   tree->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   tree->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   tree->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   tree->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   tree->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   tree->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   tree->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   tree->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   tree->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   tree->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   tree->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   tree->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   tree->SetBranchAddress("eleIDMVANonTrg", &eleIDMVANonTrg, &b_eleIDMVANonTrg);
   tree->SetBranchAddress("eleIDMVATrg", &eleIDMVATrg, &b_eleIDMVATrg);
   tree->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   tree->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   tree->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   tree->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   tree->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   tree->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   tree->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   tree->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   tree->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   tree->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   tree->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   tree->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   tree->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   tree->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   tree->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   tree->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   tree->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   tree->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   tree->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   tree->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   tree->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   tree->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   tree->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   tree->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   tree->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   tree->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   tree->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   tree->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   tree->SetBranchAddress("eleBCEn", &eleBCEn, &b_eleBCEn);
   tree->SetBranchAddress("eleBCEta", &eleBCEta, &b_eleBCEta);
   tree->SetBranchAddress("eleBCPhi", &eleBCPhi, &b_eleBCPhi);
   tree->SetBranchAddress("eleBCS25", &eleBCS25, &b_eleBCS25);
   tree->SetBranchAddress("eleBCS15", &eleBCS15, &b_eleBCS15);
   tree->SetBranchAddress("eleBCSieie", &eleBCSieie, &b_eleBCSieie);
   tree->SetBranchAddress("eleBCSieip", &eleBCSieip, &b_eleBCSieip);
   tree->SetBranchAddress("eleBCSipip", &eleBCSipip, &b_eleBCSipip);
   tree->SetBranchAddress("eleFiredTrgs", &eleFiredTrgs, &b_eleFiredTrgs);
   tree->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   tree->SetBranchAddress("nMu", &nMu, &b_nMu);
   tree->SetBranchAddress("muPt", &muPt, &b_muPt);
   tree->SetBranchAddress("muEn", &muEn, &b_muEn);
   tree->SetBranchAddress("muEta", &muEta, &b_muEta);
   tree->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   tree->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   tree->SetBranchAddress("muType", &muType, &b_muType);
   tree->SetBranchAddress("muIsLooseID", &muIsLooseID, &b_muIsLooseID);
   tree->SetBranchAddress("muIsMediumID", &muIsMediumID, &b_muIsMediumID);
   tree->SetBranchAddress("muIsTightID", &muIsTightID, &b_muIsTightID);
   tree->SetBranchAddress("muIsSoftID", &muIsSoftID, &b_muIsSoftID);
   tree->SetBranchAddress("muIsHighPtID", &muIsHighPtID, &b_muIsHighPtID);
   tree->SetBranchAddress("muD0", &muD0, &b_muD0);
   tree->SetBranchAddress("muDz", &muDz, &b_muDz);
   tree->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   tree->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   tree->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   tree->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   tree->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   tree->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   tree->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   tree->SetBranchAddress("muStations", &muStations, &b_muStations);
   tree->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   tree->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   tree->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   tree->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   tree->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   tree->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   tree->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   tree->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   tree->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   tree->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   tree->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   tree->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   tree->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   tree->SetBranchAddress("nTau", &nTau, &b_nTau);
   tree->SetBranchAddress("pfTausDiscriminationByDecayModeFinding", &pfTausDiscriminationByDecayModeFinding, &b_pfTausDiscriminationByDecayModeFinding);
   tree->SetBranchAddress("pfTausDiscriminationByDecayModeFindingNewDMs", &pfTausDiscriminationByDecayModeFindingNewDMs, &b_pfTausDiscriminationByDecayModeFindingNewDMs);
   tree->SetBranchAddress("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection, &b_tauByMVA5LooseElectronRejection);
   tree->SetBranchAddress("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection, &b_tauByMVA5MediumElectronRejection);
   tree->SetBranchAddress("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection, &b_tauByMVA5TightElectronRejection);
   tree->SetBranchAddress("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection, &b_tauByMVA5VTightElectronRejection);
   tree->SetBranchAddress("tauByLooseMuonRejection3", &tauByLooseMuonRejection3, &b_tauByLooseMuonRejection3);
   tree->SetBranchAddress("tauByTightMuonRejection3", &tauByTightMuonRejection3, &b_tauByTightMuonRejection3);
   tree->SetBranchAddress("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits);
   tree->SetBranchAddress("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits);
   tree->SetBranchAddress("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits, &b_tauByTightCombinedIsolationDeltaBetaCorr3Hits);
   tree->SetBranchAddress("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tauCombinedIsolationDeltaBetaCorrRaw3Hits);
   tree->SetBranchAddress("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT, &b_tauByVLooseIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT, &b_tauByLooseIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT, &b_tauByMediumIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT, &b_tauByTightIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT, &b_tauByVTightIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT, &b_tauByVVTightIsolationMVA3oldDMwLT);
   tree->SetBranchAddress("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw, &b_tauByIsolationMVA3oldDMwLTraw);
   tree->SetBranchAddress("tauByLooseIsolationMVA3newDMwLT", &tauByLooseIsolationMVA3newDMwLT, &b_tauByLooseIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByVLooseIsolationMVA3newDMwLT", &tauByVLooseIsolationMVA3newDMwLT, &b_tauByVLooseIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByMediumIsolationMVA3newDMwLT", &tauByMediumIsolationMVA3newDMwLT, &b_tauByMediumIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByTightIsolationMVA3newDMwLT", &tauByTightIsolationMVA3newDMwLT, &b_tauByTightIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByVTightIsolationMVA3newDMwLT", &tauByVTightIsolationMVA3newDMwLT, &b_tauByVTightIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByVVTightIsolationMVA3newDMwLT", &tauByVVTightIsolationMVA3newDMwLT, &b_tauByVVTightIsolationMVA3newDMwLT);
   tree->SetBranchAddress("tauByIsolationMVA3newDMwLTraw", &tauByIsolationMVA3newDMwLTraw, &b_tauByIsolationMVA3newDMwLTraw);
   tree->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   tree->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   tree->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   tree->SetBranchAddress("tauEt", &tauEt, &b_tauEt);
   tree->SetBranchAddress("tauCharge", &tauCharge, &b_tauCharge);
   tree->SetBranchAddress("tauP", &tauP, &b_tauP);
   tree->SetBranchAddress("tauPx", &tauPx, &b_tauPx);
   tree->SetBranchAddress("tauPy", &tauPy, &b_tauPy);
   tree->SetBranchAddress("tauPz", &tauPz, &b_tauPz);
   tree->SetBranchAddress("tauVz", &tauVz, &b_tauVz);
   tree->SetBranchAddress("tauEnergy", &tauEnergy, &b_tauEnergy);
   tree->SetBranchAddress("tauMass", &tauMass, &b_tauMass);
   tree->SetBranchAddress("tauDxy", &tauDxy, &b_tauDxy);
   tree->SetBranchAddress("tauZImpact", &tauZImpact, &b_tauZImpact);
   tree->SetBranchAddress("tauDecayMode", &tauDecayMode, &b_tauDecayMode);
   tree->SetBranchAddress("tauLeadChargedHadronExists", &tauLeadChargedHadronExists, &b_tauLeadChargedHadronExists);
   tree->SetBranchAddress("tauLeadChargedHadronEta", &tauLeadChargedHadronEta, &b_tauLeadChargedHadronEta);
   tree->SetBranchAddress("tauLeadChargedHadronPhi", &tauLeadChargedHadronPhi, &b_tauLeadChargedHadronPhi);
   tree->SetBranchAddress("tauLeadChargedHadronPt", &tauLeadChargedHadronPt, &b_tauLeadChargedHadronPt);
   tree->SetBranchAddress("tauChargedIsoPtSum", &tauChargedIsoPtSum, &b_tauChargedIsoPtSum);
   tree->SetBranchAddress("tauNeutralIsoPtSum", &tauNeutralIsoPtSum, &b_tauNeutralIsoPtSum);
   tree->SetBranchAddress("tauPuCorrPtSum", &tauPuCorrPtSum, &b_tauPuCorrPtSum);
   tree->SetBranchAddress("tauNumSignalPFChargedHadrCands", &tauNumSignalPFChargedHadrCands, &b_tauNumSignalPFChargedHadrCands);
   tree->SetBranchAddress("tauNumSignalPFNeutrHadrCands", &tauNumSignalPFNeutrHadrCands, &b_tauNumSignalPFNeutrHadrCands);
   tree->SetBranchAddress("tauNumSignalPFGammaCands", &tauNumSignalPFGammaCands, &b_tauNumSignalPFGammaCands);
   tree->SetBranchAddress("tauNumSignalPFCands", &tauNumSignalPFCands, &b_tauNumSignalPFCands);
   tree->SetBranchAddress("tauNumIsolationPFChargedHadrCands", &tauNumIsolationPFChargedHadrCands, &b_tauNumIsolationPFChargedHadrCands);
   tree->SetBranchAddress("tauNumIsolationPFNeutrHadrCands", &tauNumIsolationPFNeutrHadrCands, &b_tauNumIsolationPFNeutrHadrCands);
   tree->SetBranchAddress("tauNumIsolationPFGammaCands", &tauNumIsolationPFGammaCands, &b_tauNumIsolationPFGammaCands);
   tree->SetBranchAddress("tauNumIsolationPFCands", &tauNumIsolationPFCands, &b_tauNumIsolationPFCands);
   tree->SetBranchAddress("nJet", &nJet, &b_nJet);
   tree->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   tree->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   tree->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   tree->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   tree->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   tree->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   tree->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   tree->SetBranchAddress("jetpfCombinedInclusiveSecondaryVertexV2BJetTags", &jetpfCombinedInclusiveSecondaryVertexV2BJetTags, &b_jetpfCombinedInclusiveSecondaryVertexV2BJetTags);
   tree->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   tree->SetBranchAddress("jetpfCombinedMVABJetTags", &jetpfCombinedMVABJetTags, &b_jetpfCombinedMVABJetTags);
   tree->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   tree->SetBranchAddress("jetPUidFullDiscriminant", &jetPUidFullDiscriminant, &b_jetPUidFullDiscriminant);
   tree->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   tree->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   tree->SetBranchAddress("nAK8Jet", &nAK8Jet, &b_nAK8Jet);
   tree->SetBranchAddress("AK8JetPt", &AK8JetPt, &b_AK8JetPt);
   tree->SetBranchAddress("AK8JetEn", &AK8JetEn, &b_AK8JetEn);
   tree->SetBranchAddress("AK8JetRawPt", &AK8JetRawPt, &b_AK8JetRawPt);
   tree->SetBranchAddress("AK8JetRawEn", &AK8JetRawEn, &b_AK8JetRawEn);
   tree->SetBranchAddress("AK8JetEta", &AK8JetEta, &b_AK8JetEta);
   tree->SetBranchAddress("AK8JetPhi", &AK8JetPhi, &b_AK8JetPhi);
   tree->SetBranchAddress("AK8JetMass", &AK8JetMass, &b_AK8JetMass);
   tree->SetBranchAddress("AK8Jet_tau1", &AK8Jet_tau1, &b_AK8Jet_tau1);
   tree->SetBranchAddress("AK8Jet_tau2", &AK8Jet_tau2, &b_AK8Jet_tau2);
   tree->SetBranchAddress("AK8Jet_tau3", &AK8Jet_tau3, &b_AK8Jet_tau3);
   tree->SetBranchAddress("AK8JetCHF", &AK8JetCHF, &b_AK8JetCHF);
   tree->SetBranchAddress("AK8JetNHF", &AK8JetNHF, &b_AK8JetNHF);
   tree->SetBranchAddress("AK8JetCEF", &AK8JetCEF, &b_AK8JetCEF);
   tree->SetBranchAddress("AK8JetNEF", &AK8JetNEF, &b_AK8JetNEF);
   tree->SetBranchAddress("AK8JetNCH", &AK8JetNCH, &b_AK8JetNCH);
   tree->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
   tree->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
   tree->SetBranchAddress("AK8CHSSoftDropJetMass", &AK8CHSSoftDropJetMass, &b_AK8CHSSoftDropJetMass);
   tree->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
   tree->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
   tree->SetBranchAddress("nAK8softdropSubjet", &nAK8softdropSubjet, &b_nAK8softdropSubjet);
   tree->SetBranchAddress("AK8softdropSubjetPt", &AK8softdropSubjetPt, &b_AK8softdropSubjetPt);
   tree->SetBranchAddress("AK8softdropSubjetEta", &AK8softdropSubjetEta, &b_AK8softdropSubjetEta);
   tree->SetBranchAddress("AK8softdropSubjetPhi", &AK8softdropSubjetPhi, &b_AK8softdropSubjetPhi);
   tree->SetBranchAddress("AK8softdropSubjetMass", &AK8softdropSubjetMass, &b_AK8softdropSubjetMass);
   tree->SetBranchAddress("AK8softdropSubjetE", &AK8softdropSubjetE, &b_AK8softdropSubjetE);
   tree->SetBranchAddress("AK8softdropSubjetCharge", &AK8softdropSubjetCharge, &b_AK8softdropSubjetCharge);
   tree->SetBranchAddress("AK8softdropSubjetFlavour", &AK8softdropSubjetFlavour, &b_AK8softdropSubjetFlavour);
   tree->SetBranchAddress("AK8softdropSubjetCSV", &AK8softdropSubjetCSV, &b_AK8softdropSubjetCSV);
}

#endif	/* TREE_READER_H */

