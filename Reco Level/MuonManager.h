#ifndef SRC_MuonManager_H_
#define SRC_MuonManager_H_
#include <iostream>
#include "RoccoR.h"
#include <vector>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include "Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "AuxiliaryFunctions.h"
#include "cmdparser.hpp"

using namespace std;

const int LESS_THAN_2_MUONS= 1;
const int TOO_MANY_MUONS = 2;
const int NO_MATCHING_CHARGE =3;

//class RoccoR;

class MuonManager
    {
public:
    MuonManager(TTree* BOOM, bool IsRealData, cli::Parser& program_options);
    friend std::ostream &operator<<( std::ostream  &output, MuonManager &manager );
    double maxRelIso;
    double minMuonPt;
    bool noGenMatching;
    bool noRochesterCorrections;
    RoccoR* rc;
    int nentries;
    bool IsRealData;
    void setBranchAddress();
    virtual ~MuonManager();
    int size();
    std::vector<Muon> selectMuons(bool useTrackIso, bool useFsrIso);
    std::vector<Muon> selectMuonsPFIso();
    std::vector<Muon> selectMuonsTrackIso();
    std::vector<Muon>  selectMuonsPFIsoMinusGamma();
    double getGenPt(int muonIndex);
    Muon makeMuon(int index);
    TTree* BOOM;
    std::vector<bool> *Muon_isTriggerMatched;
    std::vector<double> *Muon_pt;
    std::vector<double> *Muon_eta;
    std::vector<double> *Muon_phi;
    std::vector<double> *Muon_energy;
    std::vector<double> *Muon_charge;
    std::vector<bool> *Muon_tight;
    std::vector<bool> *Muon_loose;
    std::vector<bool> *Muon_medium;
    std::vector<bool> *Muon_soft;
    std::vector<bool> *Muon_pf;
    std::vector<bool> *Muon_isGlobal;
    std::vector<double> *Muon_isoCharged;
    std::vector<double> *Muon_isoSum;
    std::vector<double> *Muon_isoCharParPt;
    std::vector<double> *Muon_isTrackerMuon;
    std::vector<double> *Muon_chi2;
    std::vector<double> *Muon_validHits;
    std::vector<double> *Muon_validHitsInner;
    std::vector<double> *Muon_matchedStat;
    std::vector<double> *Muon_dxy_pv;
    std::vector<double> *Muon_TLayers;
    std::vector<double> *Muon_dz_bs;
    std::vector<double> *Muon_isoNeutralHadron;
    std::vector<double> *Muon_isoPhoton;
    std::vector<double> *Muon_isoPU;
    std::vector<double> *Muon_combinedIso;
    std::vector<double> *Muon_trackRe_iso;



    Double_t Gen_Met;
    std::vector<double> *Gen_pt;
    std::vector<double> *Gen_eta;
    std::vector<double> *Gen_phi;
    std::vector<double> *Gen_energy;
    std::vector<double> *Gen_pdg_id;
    std::vector<double> *Gen_motherpdg_id;
    std::vector<double> *Gen_status;
    std::vector<int> *Gen_BmotherIndex;
    std::vector<double> *Gen_charge;

    // List of branches
    TBranch *b_Muon_isTriggerMatched;
    TBranch *b_Muon_pt;   //!
    TBranch *b_Muon_eta;   //!
    TBranch *b_Muon_phi;   //!
    TBranch *b_Muon_energy;   //!
    TBranch *b_Muon_charge;   //!
    TBranch *b_Muon_tight;   //!
    TBranch *b_Muon_loose;
    TBranch *b_Muon_medium;
    TBranch *b_Muon_soft;   //!
    TBranch *b_Muon_pf;   //!
    TBranch *b_Muon_isoCharged;   //!
    TBranch *b_Muon_isoSum;   //!
    TBranch *b_Muon_isoCharParPt;   //!
    TBranch *b_Muon_isTrackerMuon;
    TBranch *b_Muon_chi2;   //!
    TBranch *b_Muon_validHits;   //!
    TBranch *b_Muon_validHitsInner;   //!
    TBranch *b_Muon_matchedStat;   //!
    TBranch *b_Muon_dxy_pv;   //!
    TBranch *b_Muon_TLayers;   //!
    TBranch *b_Muon_dz_bs;   //!
    TBranch *b_Muon_isoNeutralHadron;   //!
    TBranch *b_Muon_isoPhoton;   //!
    TBranch *b_Muon_isoPU;   //!
    TBranch *b_Muon_combinedIso;
    TBranch *b_Muon_trackRe_iso;
    TBranch *b_Muon_isGlobal;


    TBranch *b_Gen_pt;   //!
    TBranch *b_Gen_eta;   //!
    TBranch *b_Gen_phi;   //!
    TBranch *b_Gen_energy;   //!
    TBranch *b_Gen_pdg_id;   //!
    TBranch *b_Gen_motherpdg_id;   //!
    TBranch *b_Gen_status;   //!
    TBranch *b_Gen_BmotherIndex;   //!
    TBranch *b_Gen_Met;   //!
    TBranch *b_Gen_charge;
    TBranch *b_Met_type1PF_shiftedPtUp;   //!
    TBranch *b_Met_type1PF_shiftedPtDown;   //!
    TBranch *b_pvertex_z;
    //TBranch        *b_Met_shiftedPtUp;   //!
    //TBranch        *b_Met_shiftedPtDown;   //!
    };

#endif /* SRC_MuonManager_H_ */
