#include "ElectronManager.h"

ElectronManager::ElectronManager(TTree* BOOM)
    {
    this->BOOM = BOOM;
    this->setBranchAddress();
    }

ElectronManager::~ElectronManager()
    {
    }

int ElectronManager::size()
    {
    return patElectron_pt->size();
    }

void ElectronManager::setBranchAddress()
    {
    // Set object pointer

	patElectron_eta = 0;
	patElectron_phi = 0;
	patElectron_pt = 0;
	patElectron_energy = 0;
	patElectron_isPassVeto = 0;
	patElectron_isPassLoose = 0;
	patElectron_isPassMedium = 0;
	patElectron_isPassTight = 0;
	patElectron_isPassHEEPId = 0;
	patElectron_isoChargedHadrons = 0;
	patElectron_isoNeutralHadrons = 0;
	patElectron_isoPhotons = 0;
	patElectron_isoPU = 0;
	patElectron_charge = 0;


    // Set branch addresses and branch pointers
	BOOM->SetBranchAddress("patElectron_eta", &patElectron_eta, &b_patElectron_eta);
	BOOM->SetBranchAddress("patElectron_phi", &patElectron_phi, &b_patElectron_phi);
	BOOM->SetBranchAddress("patElectron_pt", &patElectron_pt, &b_patElectron_pt);
	BOOM->SetBranchAddress("patElectron_energy", &patElectron_energy, &b_patElectron_energy);
	BOOM->SetBranchAddress("patElectron_isPassVeto", &patElectron_isPassVeto, &b_patElectron_isPassVeto);
	BOOM->SetBranchAddress("patElectron_isPassLoose", &patElectron_isPassLoose, &b_patElectron_isPassLoose);
	BOOM->SetBranchAddress("patElectron_isPassMedium", &patElectron_isPassMedium, &b_patElectron_isPassMedium);
	BOOM->SetBranchAddress("patElectron_isPassTight", &patElectron_isPassTight, &b_patElectron_isPassTight);
	BOOM->SetBranchAddress("patElectron_isPassHEEPId", &patElectron_isPassHEEPId, &b_patElectron_isPassHEEPId);
	BOOM->SetBranchAddress("patElectron_isoChargedHadrons", &patElectron_isoChargedHadrons, &b_patElectron_isoChargedHadrons);
	BOOM->SetBranchAddress("patElectron_isoNeutralHadrons", &patElectron_isoNeutralHadrons, &b_patElectron_isoNeutralHadrons);
	BOOM->SetBranchAddress("patElectron_isoPhotons", &patElectron_isoPhotons, &b_patElectron_isoPhotons);
	BOOM->SetBranchAddress("patElectron_isoPU", &patElectron_isoPU, &b_patElectron_isoPU);
	BOOM->SetBranchAddress("patElectron_charge", &patElectron_charge, &b_patElectron_charge);

    }

// True rechazar evento
bool ElectronManager::electronVeto(std::vector<Muon>& muons)
    {
    //cout << "In electron Veto" << endl;
    //cout << "No. electrons = " << this->size() << endl;
    for (int electronIndex = 0; electronIndex < this->size(); electronIndex++)
	{
	//electron veto: no electrons with cutBasedElectronID-Summer16-80X-V1-medium
	//and pT > 10 GeV
	// |eta| < 1.4442 ||
	// 1.566 <|eta| <2.5
	//dR(e,mu) > 0.4
	TLorentzVector tempMuon1 = muons.at(0).getTLorentz();
	TLorentzVector tempMuon2 = muons.at(1).getTLorentz();
	TLorentzVector electron;
	electron.SetPtEtaPhiE(patElectron_pt->at(electronIndex),patElectron_eta->at(electronIndex),patElectron_phi->at(electronIndex), patElectron_energy->at(electronIndex));
	double deltaR1 = tempMuon1.DeltaR(electron);
	double deltaR2 = tempMuon2.DeltaR(electron);

	//double deltaR1 = deltaR(patElectron_eta->at(electronIndex), patElectron_phi->at(electronIndex), tempMuon1.eta(), tempMuon1.phi());
	//double deltaR2 = deltaR(patElectron_eta->at(electronIndex), patElectron_phi->at(electronIndex), tempMuon2.eta(), tempMuon2.phi());
	bool cutBasedCondition = patElectron_isPassMedium->at(electronIndex) != 0;
	bool ptCond = this->patElectron_pt->at(electronIndex) > 10;
	bool etaCond = abs(this->patElectron_eta->at(electronIndex)) < 1.4442;
	bool etaCond2 = (abs(this->patElectron_eta->at(electronIndex)) > 1.566) && (abs(this->patElectron_eta->at(electronIndex)) < 2.5);
	//cout << "Electron: " << "cutBasedCondition: " << patElectron_isPassMedium->at(electronIndex) << ", pt: "<< patElectron_pt->at(electronIndex) << "|eta|: "<< abs(this->patElectron_eta->at(electronIndex)) << " deltaR M1: "<< deltaR1 << "deltaR M2: "<< deltaR2 <<endl;
	bool deltaRCond = (deltaR1 > 0.4 && deltaR2 > 0.4);
	if (cutBasedCondition && ptCond && (etaCond || etaCond2) && deltaRCond)
	    {
	    return true;
	    }
	}
    return false;
    }




;
