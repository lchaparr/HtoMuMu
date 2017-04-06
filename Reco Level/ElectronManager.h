#ifndef SRC_ElectronManager_H_
#define SRC_ElectronManager_H_
#include <iostream>

#include <vector>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include "Muon.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

const int LESS_THAN_2_ElectronS = 1;
const int TOO_MANY_ElectronS = 2;

class ElectronManager
    {
public:
    ElectronManager(TTree* BOOM);
    int nentries;
    friend std::ostream &operator<<( std::ostream  &output, ElectronManager &manager );

    void setBranchAddress();
    virtual ~ElectronManager();
    int size();
    bool electronVeto(std::vector<Muon>& muons);
    TTree* BOOM;
    vector<double> *patElectron_pt;
    vector<double> *patElectron_eta;
    vector<double> *patElectron_phi;
    vector<double> *patElectron_energy;
    vector<int> *patElectron_isPassVeto;
    vector<int> *patElectron_isPassLoose;
    vector<int> *patElectron_isPassMedium;
    vector<int> *patElectron_isPassTight;
    vector<int> *patElectron_isPassHEEPId;
    vector<double> *patElectron_isoChargedHadrons;
    vector<double> *patElectron_isoNeutralHadrons;
    vector<double> *patElectron_isoPhotons;
    vector<double> *patElectron_isoPU;
    vector<double> *patElectron_charge;

    TBranch *b_patElectron_pt;   //!
    TBranch *b_patElectron_eta;   //!
    TBranch *b_patElectron_phi;   //!
    TBranch *b_patElectron_energy;   //!
    TBranch *b_patElectron_isPassVeto;   //!
    TBranch *b_patElectron_isPassLoose;   //!
    TBranch *b_patElectron_isPassMedium;   //!
    TBranch *b_patElectron_isPassTight;   //!
    TBranch *b_patElectron_isPassHEEPId;   //!
    TBranch *b_patElectron_isoChargedHadrons;   //!
    TBranch *b_patElectron_isoNeutralHadrons;   //!
    TBranch *b_patElectron_isoPhotons;   //!
    TBranch *b_patElectron_isoPU;   //!
    TBranch *b_patElectron_charge;   //!

    };

#endif /* SRC_ElectronManager_H_ */
