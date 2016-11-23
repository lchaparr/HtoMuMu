#include "BSM_Analysis.h"

string myTrigger1 = "HLT_IsoMu24_v";
string myTrigger2 = "HLT_IsoMu24_eta2p1_v";

int main (int argc, char *argv[])
{
	//-----------------------histogram folders--------------------------
	TFile * MassHisto = new TFile(argv[2], "RECREATE");
	int nDir = 1;
	TDirectory *theDirectory[nDir];
	theDirectory[0] = MassHisto->mkdir("AfterDimuonSelection");
	cout <<argv[1]<<endl;
	BSM_Analysis BSM_Analysis_(MassHisto, theDirectory, nDir, argv[1]);	
}

BSM_Analysis::BSM_Analysis(TFile* theFile, TDirectory *cdDir[], int nDir, char* fname)
{
	crateHistoMasps(nDir);
	//-----------------------load PU weights----------------------------
	
	TFile file_PUdata("PU2016data_15p9ifb.root","read");
	TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
	PUweights->Scale(1/PUweights->Integral());
	TFile file_PUsim("PU2016MC.root","read");
	TH1F *PUsim = (TH1F*)file_PUsim.Get("analyzeHiMassTau/NVertices_0");
	PUsim->Scale(1/PUsim->Integral());
  
	PUweights->Divide(PUsim);
	

	//----------------------configure input file--------------------------
	TFile *f = TFile::Open(fname);
	f->cd("TNT");
	TTree* BOOM = (TTree*)f->Get("TNT/BOOM");
  
	int nentries = (int) BOOM->GetEntries();
	setBranchAddress(BOOM);
  
	cout <<"Number of entries in the sample:  "<<nentries<<endl;
	//int prueba = 100;
	//cout <<"Number of entries in the sample:  "<<prueba<<endl;
	
	//---------------------counters---------------------------------------
		
	int trigger_evt =0;
	int dimuon_evt =0;
	/*int match_evt = 0;
	int dimuon_match_evt = 0;
	
	int muonid_evt = 0;
	int muon_iso_evt = 0;
	
	int genmuon_counter = 0;
	int genantimuon_counter = 0;
	int gen_dimuons = 0;
	int tight_id_evt = 0;
	int Isolation_evt = 0;*/

	
	//---------------------Loop over the events----------------------------
	
	//for (int i = 0; i < prueba; ++i)
	for (int i = 0; i < nentries; ++i)
	{
		BOOM->GetEntry(i);
      
		//------------------------define global event weight------------------------
		double pu_weight = 1.;
		pu_weight=PUweights->GetBinContent(PUweights->FindBin(ntruePUInteractions));

		//-------------------------- TLorentz vector -------------------------------

		TLorentzVector dalitz_first_muon_vec(0., 0., 0., 0.);
		TLorentzVector dalitz_Subfirst_muon_vec(0., 0., 0., 0.);
		TLorentzVector Reco_lepton1(0.,0.,0.,0.);
		TLorentzVector Reco_lepton2(0.,0.,0.,0.);
		TLorentzVector Photon_TL(0., 0., 0., 0.);

		int lmuon_counter = 0;
		int smuon_counter = 0;
		double RelIso1 = 0, RelIso2 = 0;
		
	
		int pass_dalitz_id[nDir] = {0};
		
		//----------------------control booleans------------------------
		bool lepton1_tight = false;
		bool lepton2_tight = false;
		bool pass_dimuon = false;
		bool gen_events = false;
		bool pass_match = false;
		bool pass_uno_uno= false;
		bool pass_uno_dos = false;
        
		//-------------------For Trigger-------------------------------
		if(passRecoTrigger(myTrigger1, myTrigger2))
		{
			trigger_evt++;
		}
		
		//--------------------Dimuons (1st selection)------------------
		MuonsVectors(Reco_lepton1, Reco_lepton2);
		if((Reco_lepton1 + Reco_lepton2).M()!= 0.)
		{
			pass_dalitz_id[0] = 1;
			dimuon_evt++;
			

		//	RelIso(RelIso1, RelIso2); Falta hacer seleccion de eta, pt y reliso.

		}

		
		//----------------------Fill the histograms-----------------------
		for (int i = 0; i < nDir; i++)
		{
			if (pass_dalitz_id[i] == 1)
			{			
				_hmap_diMuon_mass[i]->Fill((Reco_lepton1 + Reco_lepton2).M(),pu_weight);
				_hmap_lead_muon_pT[i]->Fill(Reco_lepton1.Pt(),pu_weight);
				_hmap_slead_muon_pT[i]->Fill(Reco_lepton2.Pt(),pu_weight);
				_hmap_lead_muon_eta[i]->Fill(Reco_lepton1.Eta(),pu_weight);
				_hmap_slead_muon_eta[i]->Fill(Reco_lepton2.Eta(),pu_weight);
				_hmap_lead_muon_phi[i]->Fill(Reco_lepton1.Phi(),pu_weight);
				_hmap_slead_muon_phi[i]->Fill(Reco_lepton2.Phi(),pu_weight);
				_hmap_photon_E[i]->Fill(Photon_TL.E(),pu_weight);
			}     
		}
	}

	cout<<"==============================="<<endl;
	cout<<"Events that pass the trigger:     "<<trigger_evt<<endl;
	cout<<"Events that pass the dimuon selection:     "<<dimuon_evt<<endl;
	cout<<"==============================="<<endl;
	/*cout<<"events that pass match and dimuon selection:  "<<dimuon_match_evt<<endl;
	cout<<"Events that pass the Muon ID(match+dimuon):    "<<tight_id_evt<<endl;
	cout<<"Events that pass the Isolation(MuonId+match+dimuon):     "<<Isolation_evt<<endl;
	cout<<"==============================="<<endl;
	cout<<"Events that pass the Muon ID(+dimuon):    "<<muonid_evt<<endl;
	cout<<"Events that pass the Isolation(MuonId+dimuon):     "<<muon_iso_evt<<endl;*/
	
	
	//-------------------------Write the histograms----------------------------	
	theFile->cd();
	
	for (int d = 0; d < nDir; d++)
	{
		cdDir[d]->cd();	
		_hmap_diMuon_mass[d]->Write();
		_hmap_lead_muon_pT[d]->Write();
		_hmap_slead_muon_pT[d]->Write();
		_hmap_lead_muon_eta[d]->Write();
		_hmap_slead_muon_eta[d]->Write();
		_hmap_lead_muon_phi[d]->Write();
		_hmap_slead_muon_phi[d]->Write();
		_hmap_photon_E[d]->Write();
	}
	theFile->Close();
}
BSM_Analysis::~BSM_Analysis()
{
	// do anything here that needs to be done at desctruction time
}


bool BSM_Analysis::passRecoTrigger(string myTrigger1, string myTrigger2) 
{
	for(int nTrig = 0 ; nTrig < Trigger_names->size(); ++nTrig)
	{
		string trigName = Trigger_names->at(nTrig);
		if( ((trigName.find(myTrigger1) != string::npos) && (Trigger_decision->at(nTrig) == 1)) ||
			((trigName.find(myTrigger2) != string::npos) && (Trigger_decision->at(nTrig) == 1)))
		{
			return true;
		}
	}
	return false;
}

TLorentzVector BSM_Analysis::MuonsVectors(TLorentzVector& Reco_lepton1,TLorentzVector& Reco_lepton2)
{
	float dimuon_mass_int = 999999.;
	TLorentzVector first_muon_vec(0., 0., 0., 0.);
	TLorentzVector Subfirst_muon_vec(0., 0., 0., 0.);
			
	for (int m1 = 0; m1 < Muon_pt->size(); m1++)
	{
		if ((Muon_pt->size() > 1) && (passRecoTrigger(myTrigger1, myTrigger2)))
		{			
			for (int m2 = 0; m2 < Muon_pt->size(); m2++)
			{
					if (Muon_charge->at(m1)*Muon_charge->at(m2) < 0)
					{
						first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
						Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
						float dimuon_mass = (first_muon_vec + Subfirst_muon_vec).M();	
								    
						if (dimuon_mass < dimuon_mass_int )
						{
							first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
							Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));							
							Reco_lepton1 = first_muon_vec;
							Reco_lepton2 = Subfirst_muon_vec;
							dimuon_mass_int = dimuon_mass;							
						}													
					}								
			}				
		}
	}	
}

double BSM_Analysis::RelIso(double& RelIso1, double& RelIso2)
{
	float dimuon_mass_int = 999999.;
	TLorentzVector first_muon_vec(0., 0., 0., 0.);
	TLorentzVector Subfirst_muon_vec(0., 0., 0., 0.);

			
	for (int m1 = 0; m1 < Muon_pt->size(); m1++)
	{
		if ((Muon_pt->size() > 1) && (passRecoTrigger(myTrigger1, myTrigger2)))
		{			
			for (int m2 = 0; m2 < Muon_pt->size(); m2++)
			{
					if (Muon_charge->at(m1)*Muon_charge->at(m2) < 0)
					{
						first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
						Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
						float dimuon_mass = (first_muon_vec + Subfirst_muon_vec).M();	
								    
						if (dimuon_mass < dimuon_mass_int )
						{							
							dimuon_mass_int = dimuon_mass;	
							RelIso1 = Muon_isoSum->at(m1)/Muon_pt->at(m1);
							RelIso2 = Muon_isoSum->at(m2)/Muon_pt->at(m2);					
						}													
					}								
			}				
		}
	}	
}


TLorentzVector BSM_Analysis::PhotonsVectors(TLorentzVector& Photon_TL)
{
	vector<TLorentzVector> Photon_vec;
	Photon_vec.erase (Photon_vec.begin(),Photon_vec.end());
	
	for(int ph = 0; ph < Photon_pt->size(); ph++)
	{
		if(Photon_pt->size()>0)
		{
			Photon_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));	
			Photon_vec.push_back(Photon_TL);
		}		
	}
}

void BSM_Analysis::crateHistoMasps (int directories)
{
	for (int i = 0; i < directories; i++)
	{
		// Muon distributions
		_hmap_diMuon_mass[i]      = new TH1F("diMuonMass",      "m_{#mu, #mu}", 600., 0., 300.);
		_hmap_lead_muon_pT[i]     = new TH1F("lead_muon_pT",    "#mu p_{T}",    600, 0., 300.);
		_hmap_slead_muon_pT[i]    = new TH1F("slead_muon_pT",   "#mu p_{T}",    600, 0., 300.);
		_hmap_lead_muon_eta[i]    = new TH1F("lead_muon_eta",   "#mu #eta",     100, -2.5, 2.5);
		_hmap_slead_muon_eta[i]   = new TH1F("slead_muon_eta",  "#mu #eta",     100, -2.5, 2.5);
		_hmap_lead_muon_phi[i]    = new TH1F("lead_muon_phi",   "#mu #phi",     140, -3.5, 3.5);
		_hmap_slead_muon_phi[i]    = new TH1F("slead_muon_phi", "#mu #phi",     140, -3.5, 3.5);    
		_hmap_photon_E[i] = new TH1F("PhotonEnergy","E#gamma", 200., 0., 50.);
	}
}

void BSM_Analysis::setBranchAddress(TTree* BOOM)
{
  
	// Set object pointer
  
	Muon_pt = 0;
	Muon_eta = 0;
	Muon_phi = 0;
	Muon_energy = 0;
	Muon_charge = 0;
	Muon_tight = 0;
	Muon_loose = 0;
	Muon_medium = 0;
	Muon_pf = 0;
	Muon_soft = 0;
	Muon_isoCharged = 0;
	Muon_isoSum = 0;
	Muon_isoCharParPt = 0;
	Muon_isTrackerMuon = 0;
	Muon_chi2 = 0;
	Muon_validHits = 0;
	Muon_validHitsInner = 0;
	Muon_matchedStat = 0;
	Muon_dxy_pv = 0;
	Muon_TLayers = 0;
	Muon_dz_bs = 0;
	Muon_isoNeutralHadron = 0;
	Muon_isoPhoton = 0;
	Muon_isoPU = 0;
	
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
	
	Photon_eta = 0;
	Photon_phi = 0;
	Photon_pt = 0;
	Photon_energy = 0;
	Photon_et = 0;
	Photon_HoverE = 0;
	Photon_phoR9 = 0;
	Photon_SigmaIEtaIEta = 0;
	Photon_SigmaIPhiIPhi = 0;
	Photon_PFChIso = 0;
	Photon_PFPhoIso = 0;
	Photon_PFNeuIso = 0;
	Photon_EleVeto = 0;
	
	Jet_pt = 0;
	Jet_eta = 0;
	Jet_phi = 0;
	Jet_energy = 0;
	Jet_bDiscriminator = 0;
	Jet_mass = 0;
	Jet_neutralHadEnergyFraction = 0;
	Jet_neutralEmEmEnergyFraction = 0;
	Jet_chargedHadronEnergyFraction = 0;
	Jet_chargedEmEnergyFraction = 0;
	Jet_muonEnergyFraction = 0;
	Jet_electronEnergy = 0;
	Jet_photonEnergy = 0;
	
	Gen_eta = 0;
	Gen_phi = 0;
	Gen_pt = 0;
	Gen_energy = 0;
	Gen_pdg_id = 0;
	Gen_motherpdg_id = 0;
	Gen_status = 0;
	Gen_BmotherIndex = 0;
	Gen_charge =0;
	  
	UncorrJet_pt = 0;
	Trigger_names = 0;
	Trigger_decision = 0;
	ntruePUInteractions = 0;
	
	// Set branch addresses and branch pointers
	if(!BOOM) return;
	BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
	BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
	BOOM->SetBranchAddress("nTruePUInteractions", &ntruePUInteractions, &b_ntruePUInteractions);
	
	BOOM->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
	BOOM->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
	BOOM->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
	BOOM->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
	BOOM->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
	BOOM->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
	BOOM->SetBranchAddress("Muon_loose", &Muon_loose, &b_Muon_loose);
	BOOM->SetBranchAddress("Muon_medium", &Muon_medium, &b_Muon_medium);
	BOOM->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
	BOOM->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
	BOOM->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
	BOOM->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
	BOOM->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
	BOOM->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
	BOOM->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
	BOOM->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
	BOOM->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
	BOOM->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
	BOOM->SetBranchAddress("Muon_dxy_pv", &Muon_dxy_pv, &b_Muon_dxy_pv);
	BOOM->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
	BOOM->SetBranchAddress("Muon_dz_bs", &Muon_dz_bs, &b_Muon_dz_bs);
	BOOM->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
	BOOM->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
	BOOM->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
	
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
	
	BOOM->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
	BOOM->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
	BOOM->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
	BOOM->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
	BOOM->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
	BOOM->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
	BOOM->SetBranchAddress("Photon_phoR9", &Photon_phoR9, &b_Photon_phoR9);
	BOOM->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
	BOOM->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
	BOOM->SetBranchAddress("Photon_PFChIso", &Photon_PFChIso, &b_Photon_PFChIso);
	BOOM->SetBranchAddress("Photon_PFPhoIso", &Photon_PFPhoIso, &b_Photon_PFPhoIso);
	BOOM->SetBranchAddress("Photon_PFNeuIso", &Photon_PFNeuIso, &b_Photon_PFNeuIso);
	BOOM->SetBranchAddress("Photon_EleVeto", &Photon_EleVeto, &b_Photon_EleVeto);
	
	BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
	BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
	BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
	BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
	BOOM->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
	BOOM->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
	BOOM->SetBranchAddress("Jet_neutralHadEnergyFraction", &Jet_neutralHadEnergyFraction, &b_Jet_neutralHadEnergyFraction);
	BOOM->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &Jet_neutralEmEmEnergyFraction, &b_Jet_neutralEmEmEnergyFraction);
	BOOM->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
	BOOM->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
	BOOM->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
	BOOM->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
	BOOM->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
	BOOM->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
	BOOM->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
	BOOM->SetBranchAddress("Met_puppi_pt", &Met_puppi_pt, &b_Met_puppi_pt);
	BOOM->SetBranchAddress("Met_puppi_sumEt", &Met_puppi_sumEt, &b_Met_puppi_sumEt);
	BOOM->SetBranchAddress("Met_puppi_phi", &Met_puppi_phi, &b_Met_puppi_phi);
	BOOM->SetBranchAddress("Met_puppi_px", &Met_puppi_px, &b_Met_puppi_px);
	BOOM->SetBranchAddress("Met_puppi_py", &Met_puppi_py, &b_Met_puppi_py);
	BOOM->SetBranchAddress("Met_puppi_pz", &Met_puppi_pz, &b_Met_puppi_pz);
	
	BOOM->SetBranchAddress("Gen_pt", &Gen_pt, &b_Gen_pt);
	BOOM->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
	BOOM->SetBranchAddress("Gen_phi", &Gen_phi, &b_Gen_phi);
	BOOM->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
	BOOM->SetBranchAddress("Gen_pdg_id", &Gen_pdg_id, &b_Gen_pdg_id);
	BOOM->SetBranchAddress("Gen_motherpdg_id", &Gen_motherpdg_id, &b_Gen_motherpdg_id);
	BOOM->SetBranchAddress("Gen_status", &Gen_status, &b_Gen_status);
	BOOM->SetBranchAddress("Gen_BmotherIndex", &Gen_BmotherIndex, &b_Gen_BmotherIndex);
	BOOM->SetBranchAddress("Gen_Met", &Gen_Met, &b_Gen_Met);
	BOOM->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
	
	BOOM->SetBranchAddress("Met_type1PF_shiftedPtUp", &Met_type1PF_shiftedPtUp, &b_Met_type1PF_shiftedPtUp);
	BOOM->SetBranchAddress("Met_type1PF_shiftedPtDown", &Met_type1PF_shiftedPtDown, &b_Met_type1PF_shiftedPtDown);
};
