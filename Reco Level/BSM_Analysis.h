////////////////////////////////////////////////////////////////////
// Developer: Andres Florez, Universidad de los Andes, Colombia. //
// Developer: Denis Rathjens, University of Hamburg, Germany     //
//////////////////////////////////////////////////////////////////


#ifndef BSM_Analysis_h
#define BSM_Analysis_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
using namespace std;

class BSM_Analysis {
public :
   BSM_Analysis(TFile*, TDirectory* dir[], int nDir, char*);
   ~BSM_Analysis();

   // create Histo maps
   void crateHistoMasps (int);
   // Define maps for histograms
   // For muons
   std::map<unsigned int, TH1*> _hmap_lead_muon_pT;
   std::map<unsigned int, TH1*> _hmap_slead_muon_pT;
   std::map<unsigned int, TH1*> _hmap_lead_muon_eta;
   std::map<unsigned int, TH1*> _hmap_slead_muon_eta;
   std::map<unsigned int, TH1*> _hmap_lead_muon_phi;
   std::map<unsigned int, TH1*> _hmap_slead_muon_phi;
   std::map<unsigned int, TH1*> _hmap_diMuon_mass;
   std::map<unsigned int, TH1*> _hmap_diMuon_mass_dmass;
   std::map<unsigned int, TH1*> _hmap_lead_muon_pT_dmass;
   std::map<unsigned int, TH1*> _hmap_slead_muon_pT_dmass;
   std::map<unsigned int, TH1*> _hmap_slead_muon_eta_dmass;
   std::map<unsigned int, TH1*> _hmap_lead_muon_eta_dmass;
   std::map<unsigned int, TH1*> _hmap_slead_muon_phi_dmass;
   std::map<unsigned int, TH1*> _hmap_lead_muon_phi_dmass;
   
 	std::map<unsigned int, TH1*> _hmap_DeltaR_Genp_Recp;
	std::map<unsigned int, TH1*> _hmap_DeltaR_Genm_Recm;
	std::map<unsigned int, TH1*> _hmap_DeltaR_Gen2_Rec1;
	std::map<unsigned int, TH1*> _hmap_DeltaR_Gen1_Rec2; 
   
   // For Photons
   std::map<unsigned int, TH1*> _hmap_threebody_light;
   std::map<unsigned int, TH1*> _hmap_threebody_heavy;
   std::map<unsigned int, TH1*> _hmap_threebody_light_dmass;
   std::map<unsigned int, TH1*> _hmap_threebody_heavy_dmass;
   std::map<unsigned int, TH1*> _hmap_photon_pT;
   std::map<unsigned int, TH1*> _hmap_photon_eta;
   std::map<unsigned int, TH1*> _hmap_photon_phi;
   std::map<unsigned int, TH1*> _hmap_DeltaR_lep_light;
   std::map<unsigned int, TH1*> _hmap_DeltaR_sublep_light;
   std::map<unsigned int, TH1*> _hmap_DeltaR_lep_high;
   std::map<unsigned int, TH1*> _hmap_DeltaR_sublep_high;

   // Define Branches
   void setBranchAddress(TTree* BOOM);
   vector<string>  *Trigger_names;
   vector<int>     *Trigger_decision;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_energy;
   vector<double>  *Muon_charge;
   vector<bool>    *Muon_tight;
   vector<bool>    *Muon_loose;
   vector<bool>    *Muon_medium;
   vector<bool>    *Muon_soft;
   vector<double>  *Muon_isoCharged;
   vector<double>  *Muon_isoSum;
   vector<double>  *Muon_isoCharParPt;
   vector<double>  *Muon_isTrackerMuon;
   vector<double>  *Muon_chi2;
   vector<double>  *Muon_validHits;
   vector<double>  *Muon_validHitsInner;
   vector<double>  *Muon_matchedStat;
   vector<double>  *Muon_dxy_pv;
   vector<double>  *Muon_TLayers;
   vector<double>  *Muon_dz_bs;
   vector<double>  *Muon_isoNeutralHadron;
   vector<double>  *Muon_isoPhoton;
   vector<double>  *Muon_isoPU;
   
   vector<double>  *patElectron_pt;
   vector<double>  *patElectron_eta;
   vector<double>  *patElectron_phi;
   vector<double>  *patElectron_energy;
   vector<int>     *patElectron_isPassVeto;
   vector<int>     *patElectron_isPassLoose;
   vector<int>     *patElectron_isPassMedium;
   vector<int>     *patElectron_isPassTight;
   vector<int>     *patElectron_isPassHEEPId;
   vector<double>  *patElectron_isoChargedHadrons;
   vector<double>  *patElectron_isoNeutralHadrons;
   vector<double>  *patElectron_isoPhotons;
   vector<double>  *patElectron_isoPU;
   vector<double>  *patElectron_charge;
   
   vector<float>   *Photon_eta;
   vector<float>   *Photon_phi;
   vector<float>   *Photon_pt;
   vector<float>   *Photon_energy;
   vector<float>   *Photon_et;
   vector<float>   *Photon_HoverE;
   vector<float>   *Photon_phoR9;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_PFChIso;
   vector<float>   *Photon_PFPhoIso;
   vector<float>   *Photon_PFNeuIso;
   vector<int>     *Photon_EleVeto;
   
   vector<double>  *Jet_pt;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_bDiscriminator;
   vector<double>  *Jet_mass;
   vector<double>  *Jet_neutralHadEnergyFraction;
   vector<double>  *Jet_neutralEmEmEnergyFraction;
   vector<double>  *Jet_chargedHadronEnergyFraction;
   vector<double>  *Jet_chargedEmEnergyFraction;
   vector<double>  *Jet_muonEnergyFraction;
   vector<double>  *Jet_electronEnergy;
   vector<double>  *Jet_photonEnergy;
   vector<double>  *UncorrJet_pt;
   //Int_t           npuVertices;
   Float_t         ntruePUInteractions;
   //Int_t           ootnpuVertices;
   //Int_t           npuVerticesp1;
   Int_t           bestVertices;
   Double_t        Met_puppi_pt;
   Double_t        Met_puppi_sumEt;
   Double_t        Met_puppi_phi;
   Double_t        Met_puppi_px;
   Double_t        Met_puppi_py;
   Double_t        Met_puppi_pz;
   
   Double_t        Gen_Met;
   vector<double>  *Gen_pt;
   vector<double>  *Gen_eta;
   vector<double>  *Gen_phi;
   vector<double>  *Gen_energy;
   vector<double>  *Gen_pdg_id;
   vector<double>  *Gen_motherpdg_id;
   vector<double>  *Gen_status;
   vector<int>  *Gen_BmotherIndex;
   vector <double> *Gen_charge;
   
   Double_t        Met_type1PF_shiftedPtUp;
   Double_t        Met_type1PF_shiftedPtDown;
     
      // List of branches
   TBranch        *b_Trigger_names;
   TBranch        *b_Trigger_decision;
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tight;   //!
   TBranch        *b_Muon_loose;
   TBranch        *b_Muon_medium;
   TBranch        *b_Muon_soft;   //!
   TBranch        *b_Muon_pf;   //!
   TBranch        *b_Muon_isoCharged;   //!
   TBranch        *b_Muon_isoSum;   //!
   TBranch        *b_Muon_isoCharParPt;   //!
   TBranch		  *b_Muon_isTrackerMuon;
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_validHitsInner;   //!
   TBranch        *b_Muon_matchedStat;   //!
   TBranch        *b_Muon_dxy_pv;   //!
   TBranch        *b_Muon_TLayers;   //!
   TBranch        *b_Muon_dz_bs;   //!
   TBranch        *b_Muon_isoNeutralHadron;   //!
   TBranch        *b_Muon_isoPhoton;   //!
   TBranch        *b_Muon_isoPU;   //!
   
   TBranch        *b_patElectron_pt;   //!
   TBranch        *b_patElectron_eta;   //!
   TBranch        *b_patElectron_phi;   //!
   TBranch        *b_patElectron_energy;   //!
   TBranch        *b_patElectron_isPassVeto;   //!
   TBranch        *b_patElectron_isPassLoose;   //!
   TBranch        *b_patElectron_isPassMedium;   //!
   TBranch        *b_patElectron_isPassTight;   //!
   TBranch        *b_patElectron_isPassHEEPId;   //!
   TBranch        *b_patElectron_isoChargedHadrons;   //!
   TBranch        *b_patElectron_isoNeutralHadrons;   //!
   TBranch        *b_patElectron_isoPhotons;   //!
   TBranch        *b_patElectron_isoPU;   //!
   TBranch        *b_patElectron_charge;   //!
	 
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_et;
   TBranch        *b_Photon_HoverE;
   TBranch        *b_Photon_phoR9;
   TBranch        *b_Photon_SigmaIEtaIEta;
   TBranch        *b_Photon_SigmaIPhiIPhi;
   TBranch        *b_Photon_PFChIso;
   TBranch        *b_Photon_PFPhoIso;
   TBranch        *b_Photon_PFNeuIso;
   TBranch        *b_Photon_EleVeto;
   
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_bDiscriminator;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_UncorrJet_pt;   //!
   //TBranch        *b_npuVertices;   //!
   TBranch        *b_ntruePUInteractions;   //!
   //TBranch        *b_ootnpuVertices;   //!
   //TBranch        *b_npuVerticesp1;   //!
   TBranch        *b_bestVertices;   //!
   TBranch        *b_Met_puppi_pt;   //!
   TBranch        *b_Met_puppi_sumEt;   //!
   TBranch        *b_Met_puppi_phi;   //!
   TBranch        *b_Met_puppi_px;   //!
   TBranch        *b_Met_puppi_py;   //!
   TBranch        *b_Met_puppi_pz;   //!
   
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   TBranch        *b_Gen_energy;   //!
   TBranch        *b_Gen_pdg_id;   //!
   TBranch        *b_Gen_motherpdg_id;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_BmotherIndex;   //!  
   TBranch        *b_Gen_Met;   //!
   TBranch        *b_Gen_charge;
   TBranch        *b_Met_type1PF_shiftedPtUp;   //!
   TBranch        *b_Met_type1PF_shiftedPtDown;   //!
   //TBranch        *b_Met_shiftedPtUp;   //!
   //TBranch        *b_Met_shiftedPtDown;   //!
};
#endif
