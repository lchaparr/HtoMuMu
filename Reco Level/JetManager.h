#ifndef SRC_JetManager_H_
#define SRC_JetManager_H_
#include <iostream>

#include <vector>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include "Jet.h"
#include "Muon.h"
#include <algorithm>
#include "AuxiliaryFunctions.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "cmdparser.hpp"



using namespace std;

const int TOO_FEW_JETS= 1;
const int NO_JET_ABOVE_40_PT = 2;
class JetManager
    {
public:
    JetManager(TTree* BOOM, cli::Parser& program_options);
    friend std::ostream &operator<<( std::ostream  &output, JetManager &manager );
    double minJetPt;
    int nentries;
    void setBranchAddress();
    virtual ~JetManager();
    int size();
    std::vector<Jet> selectJets(std::vector<Muon>& muons);
    bool bJetVeto();
    Jet makeJet(int index);
    TTree* BOOM;
    std::vector<double> *Jet_pt;
    std::vector<double> *Jet_eta;
    std::vector<double> *Jet_phi;
    std::vector<double> *Jet_energy;
    std::vector<double> *Jet_mass;
    std::vector<double> *Jet_neutralHadEnergyFraction;
    std::vector<double> *Jet_neutralEmEmEnergyFraction;
    std::vector<double> *Jet_chargedHadronEnergyFraction;
    std::vector<double> *Jet_chargedEmEnergyFraction;
    std::vector<double> *Jet_muonEnergyFraction;
    std::vector<double> *Jet_electronEnergy;
    std::vector<double> *Jet_photonEnergy;
    std::vector<double> *UncorrJet_pt;
    std::vector<double> * Jet_bDiscriminator_pfCISVV2;

    // List of branches
    TBranch *b_Jet_pt;   //!
    TBranch *b_Jet_eta;   //!
    TBranch *b_Jet_phi;   //!
    TBranch *b_Jet_energy;   //!
    TBranch *b_Jet_mass;   //!
    TBranch *b_Jet_neutralHadEnergyFraction;   //!
    TBranch *b_Jet_neutralEmEmEnergyFraction;   //!
    TBranch *b_Jet_chargedHadronEnergyFraction;   //!
    TBranch *b_Jet_chargedEmEnergyFraction;   //!
    TBranch *b_Jet_muonEnergyFraction;   //!
    TBranch *b_Jet_electronEnergy;   //!
    TBranch *b_Jet_photonEnergy;   //!
    TBranch *b_UncorrJet_pt;   //!
    TBranch *b_Jet_puppi_pt;   //!
    TBranch *b_Jet_bDiscriminator_pfCISVV2;
    };

#endif /* SRC_JetManager_H_ */
