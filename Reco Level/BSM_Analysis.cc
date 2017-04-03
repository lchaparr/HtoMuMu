//*********************************************************************************************
// BSM Analysis: Main analysis code for H2Mu. The code includes the main selection for decay 
// 				 channel, also includes FSR recovery. 
//*********************************************************************************************

#include "BSM_Analysis.h"
#include "HistogramManager.h"

string myTrigger1 = "HLT_IsoMu24_";
string myTrigger2 = "HLT_IsoTkMu24_";

bool Gluon_Fusion = false;
bool Vector_Boson = false;
bool W_PlusH = false;
bool W_MinusH = false;
bool ZH = false;
bool ttH = false;

std::map<int, std::string> eventCategories;
std::map<std::string, int> counters;
std::map<std::string, std::string> countersDescription;

int main(int argc, char *argv[])
{
	cout << argv[1] << endl;
	BSM_Analysis BSM_Analysis_(argv[2], NULL, 0, argv[1], argv[3], argv[4]);
}

BSM_Analysis::BSM_Analysis(std::string theFile, TDirectory *cdDir[], int nDir, char* fname, char* dataType, char* RelIsotag)
{
	data = false;
	RoccoR rc("rcdata.2016.v3");
	std::string RelIsotagstr(RelIsotag);
	double RelIsomax = 0.25;
	bool usetrackiso = false;
	if (RelIsotagstr == "PFIsoLoose")
	{
		RelIsomax = 0.25;
		usetrackiso = false;
	}
	if (RelIsotagstr == "PFIsoTight")
	{
		RelIsomax = 0.12;
		usetrackiso = false;
	}
	if (RelIsotagstr == "TkrIsoLoose")
	{
		RelIsomax = 0.10;
		usetrackiso = true;
	}
	if (RelIsotagstr == "TkrIsoTight")
	{
		RelIsomax = 0.05;
		usetrackiso = true;
	}
	if (dataType == std::string("DATA"))
	{
		cout << "is data" << endl;
		data = true;
	}

	eventCategories[0] = "Unclassified";
	eventCategories[1] = "VBF/Tight";
	eventCategories[2] = "VBF/Loose";
	eventCategories[3] = "ggF/Tight";
	eventCategories[4] = "0,1_Jet/Tight";
	eventCategories[5] = "0,1_Jet/Loose";
	// SubCarpetas
	string strRelIso = "/RelIso";
	string strNonRelIso = "/NonRelIso";
	string strAllEvents = "/AllEvents"; //all the selected events, including rel iso and no rel iso.
	string strFSREvents = "/FSREvents";
	string strNOFSREvents = "/NOFSREvents";
	// Sub Sub Carpetas
	string tagOnlyDimuon = "/OnlyDimuon";
	string tagDimuonPlusPH = "/DimuonPlusPH";
	string tagInvariantMass = "OnlyDimuon DimuonPlusPH";
	string strAfterDimuonSelection = "AfterDimuonSelection";
	string strAfterDimuonRelIso = "/AfterDimuonRelIso";
	string oneGamma = "/oneGamma";
	string twoGamma = "/twoGamma";
	string passTrigger = "PassTrigger";
	string passDimuon = "passDimuon";

	HistogramManager histos(theFile + RelIsotagstr + ".root");
	// Crear la carpeta AfterDimuonSelection
	histos.createMuonHistograms(strAfterDimuonSelection);
	// Recorrer todas las categorias, crear las subcarpetas
	counters[strAllEvents] = 0;
	countersDescription[strAllEvents] = "0. All Events: ";
	counters[passTrigger] = 0;
	countersDescription[passTrigger] = "0. Events that pass the trigger: ";
	counters[passDimuon] = 0;
	countersDescription[passDimuon] = "0. Events that pass the dimuon selection: ";
	counters[strRelIso] = 0;
	countersDescription[strRelIso] = "0. Events that pass the RelIso: ";
	counters[strNonRelIso] = 0;
	countersDescription[strNonRelIso] = "0. Events that do NOT pass the RelIso: ";

	for (auto cat : eventCategories)
	{
		std::string categoryName = cat.second;
		std::string categoryNumber = std::to_string(cat.first + 1) + ". ";
		//Create counters
		counters[categoryName + strAllEvents] = 0;
		counters[categoryName + strRelIso] = 0;
		counters[categoryName + strRelIso + strNOFSREvents] = 0;
		counters[categoryName + strRelIso + strFSREvents] = 0;
		counters[categoryName + strRelIso + strFSREvents + oneGamma] = 0;
		counters[categoryName + strRelIso + strFSREvents + twoGamma] = 0;
		counters[categoryName + strNonRelIso] = 0;
		counters[categoryName + strNonRelIso + strNOFSREvents] = 0;
		counters[categoryName + strNonRelIso + strFSREvents] = 0;
		counters[categoryName + strNonRelIso + strFSREvents + oneGamma] = 0;
		counters[categoryName + strNonRelIso + strFSREvents + twoGamma] = 0;
		// Create Counter Description
		countersDescription[categoryName + strAllEvents] = categoryNumber + categoryName + ", All Events :";
		countersDescription[categoryName + strRelIso] = categoryNumber + categoryName + ", Events that pass the RelIso: ";
		;
		countersDescription[categoryName + strRelIso + strNOFSREvents] = categoryNumber + categoryName + ", Events pass RelIso, closer to H or without FSR Ph: ";
		countersDescription[categoryName + strRelIso + strFSREvents] = categoryNumber + categoryName + ", Events pass RelIso, with * FSR ph: ";
		countersDescription[categoryName + strRelIso + strFSREvents + oneGamma] = categoryNumber + categoryName + ", Events pass RelIso, with 1 FSR ph: ";
		countersDescription[categoryName + strRelIso + strFSREvents + twoGamma] = categoryNumber + categoryName + ", Events pass RelIso, with 2 FSR ph: ";
		countersDescription[categoryName + strNonRelIso] = categoryNumber + categoryName + ", Events that do not pass the RelIso: ";
		countersDescription[categoryName + strNonRelIso + strNOFSREvents] = categoryNumber + categoryName
				+ ", Events that do NOT pass RelIso, closer to H or without FSR Ph: ";
		countersDescription[categoryName + strNonRelIso + strFSREvents] = categoryNumber + categoryName + ", Events that do NOT pass RelIso, with * FSR ph: ";
		countersDescription[categoryName + strNonRelIso + strFSREvents + oneGamma] = categoryNumber + categoryName
				+ ", Events that do not pass iso, with 1 FSR ph: ";
		countersDescription[categoryName + strNonRelIso + strFSREvents + twoGamma] = categoryNumber + categoryName
				+ ", Events that do not pass iso, with 2 FSR ph: ";

		// All Events histograms
		histos.createDimuonHistograms(categoryName + strAllEvents);
		//RelIso Histograms
		histos.createDimuonHistograms(categoryName + strRelIso + strNOFSREvents, "OnlyDimuon");
		histos.createDimuonHistograms(categoryName + strRelIso + strFSREvents, tagInvariantMass);
		//Non RelIso Histograms
		histos.createDimuonHistograms(categoryName + strNonRelIso + strNOFSREvents, "OnlyDimuon");
		histos.createDimuonHistograms(categoryName + strNonRelIso + strFSREvents, tagInvariantMass);
	}

	//-----------------------load PU weights----------------------------

	TFile file_PUdata("PU2016data_15p9ifb.root", "read");
	TH1F *PUweights = (TH1F*) file_PUdata.Get("analyzeHiMassTau/NVertices_0");
	PUweights->Scale(1 / PUweights->Integral());
	TFile file_PUsim("PU2016MC.root", "read");
	TH1F *PUsim = (TH1F*) file_PUsim.Get("analyzeHiMassTau/NVertices_0");
	PUsim->Scale(1 / PUsim->Integral());
	PUweights->Divide(PUsim);

	//----------------------configure input file--------------------------

	TFile *f = TFile::Open(fname);
	f->cd("TNT");
	TTree* BOOM = (TTree*) f->Get("TNT/BOOM");

	int nentries = (int) BOOM->GetEntries();
	setBranchAddress(BOOM);

	cout << "No. of entries in the sample:  " << nentries << endl;
	int prueba = 10000;
	int trigger_evt = 0;
	int electronveto = 0;
	int btagveto = 0;
	//---------------------Loop over the events----------------------------
	//for (int i = 0; i < prueba; ++i)
	for (int i = 0; i < nentries; ++i)
	{
		//if(i!=193554)continue;
		BOOM->GetEntry(i);


		//------------------------define global event weight------------------------
		double pu_weight = 1.;
		pu_weight = PUweights->GetBinContent(PUweights->FindBin(ntruePUInteractions));
		//-------------------------- TLorentz vector -------------------------------
		TLorentzVector Reco_lepton1(0., 0., 0., 0.);
		TLorentzVector Reco_lepton2(0., 0., 0., 0.);
		TLorentzVector Dimuon;
		TLorentzVector Gen_lepton1_vec(0., 0., 0., 0.);
		TLorentzVector Gen_lepton2_vec(0., 0., 0., 0.);
		TLorentzVector Photon_TL(0., 0., 0., 0.);
		TLorentzVector *pEmpty = NULL;
		TLorentzVector Empty;
		pEmpty = &Empty;

		vector<TLorentzVector> Photon_vec;
		Photon_vec.erase(Photon_vec.begin(), Photon_vec.end());

		TLorentzVector FSR_Photon_TL(0., 0., 0., 0.);

		vector<TLorentzVector> FSRPhoton_vec;
		FSRPhoton_vec.erase(FSRPhoton_vec.begin(), FSRPhoton_vec.end());

		TLorentzVector Photon_Gen_TL(0., 0., 0., 0.);

		vector<TLorentzVector> Photon_Gen_vec;
		Photon_Gen_vec.erase(Photon_Gen_vec.begin(), Photon_Gen_vec.end());

		vector<float> FSRPhotonISO_vec;
		FSRPhotonISO_vec.erase(FSRPhotonISO_vec.begin(), FSRPhotonISO_vec.end());
		TLorentzVector FSR_PhotonToLep1(0., 0., 0., 0.);
		TLorentzVector FSR_PhotonToLep2(0., 0., 0., 0.);

		double RelIso1 = 0, RelIso2 = 0, RelIso_NoFSR1 = 0, RelIso_NoFSR2 = 0;
		float FSRIso = 0;
		int eventCategory = 0, GF_cat = 0;

		//----------------------control booleans------------------------
		bool pass_dimuon = false;
		bool pass_trigger = false;
		bool pass_dimuon_RelIso = false;
		//bool pass_noFsr_RelIso = false;

		//-------------------For Trigger-------------------------------
		if (passRecoTrigger(myTrigger1, myTrigger2))
			pass_trigger = true;

		//--------------------Dimuons (1st selection)------------------
		pass_dimuon_RelIso = MuonsVectors(Reco_lepton1, Reco_lepton2, RelIso_NoFSR1, RelIso_NoFSR2, rc);
		//Reco_lepton1.Print();
		//cout<<"lepton2"<<endl;
		//Reco_lepton2.Print();
		Dimuon = Reco_lepton1 + Reco_lepton2;
		bool passDimuonMass = false;

		if (pass_trigger)
		{
			trigger_evt++;
			if (Electron_veto(Reco_lepton1, Reco_lepton2) && pass_dimuon_RelIso)
			{
				electronveto++;
			}

			if (BtagVeto())
			{
				btagveto++;
			}
		}

		if (Dimuon.M() < 12 || Electron_veto(Reco_lepton1, Reco_lepton2) || !pass_trigger || BtagVeto())
		{
			continue;
		}

		if (Dimuon.M() > 115.0 && Dimuon.M() < 135.0)
		{
			passDimuonMass = true;
		}

		//--------------------FSR Recovery algorithm---------------------
		double DR_FSR_lep1, DR_FSR_lep2, DROverET_1, DROverET_2, minDrOEt1 = 999.9, minDrOEt2 = 999.9;
		bool fsr_photonlep1 = false, fsr_photonlep2 = false;

		for (int ph = 0; ph < FSRPhoton_pt->size(); ph++)
		{
			FSR_Photon_TL.SetPtEtaPhiE(FSRPhoton_pt->at(ph), FSRPhoton_eta->at(ph), FSRPhoton_phi->at(ph), FSRPhoton_energy->at(ph));
			FSRIso = (FSRPhoton_isoCHPU->at(ph) + FSRPhoton_isoNHPhot->at(ph)) / FSRPhoton_pt->at(ph);
			if (FSR_Photon_TL.Pt() > 2.0 && abs(FSR_Photon_TL.Eta()) < 2.4 && FSRIso < 1.8)
			{
				DR_FSR_lep1 = FSRdeltaR(FSR_Photon_TL, Reco_lepton1);
				DR_FSR_lep2 = FSRdeltaR(FSR_Photon_TL, Reco_lepton2);
				DROverET_1 = FSRDROverET2(FSR_Photon_TL, DR_FSR_lep1);
				DROverET_2 = FSRDROverET2(FSR_Photon_TL, DR_FSR_lep2);

				//For lepton 1 Only;
				//if (DR_FSR_lep1 >= 0.5 && DROverET_1 < 0.012)
				if (DROverET_1 < 0.012)
				{
					if (FSR_PhotonToLep2 != FSR_Photon_TL)
					{
						if (DROverET_1 < minDrOEt1)
						{
							minDrOEt1 = DROverET_1;
							FSR_PhotonToLep1 = FSR_Photon_TL;
							fsr_photonlep1 = true;
						}
					}
				}
				//For lepton 2 Only;
				//if (DR_FSR_lep2 >= 0.5 && DROverET_2 < 0.012)
				if (DROverET_2 < 0.012)
				{
					if (FSR_PhotonToLep1 != FSR_Photon_TL)
					{
						if (DROverET_2 < minDrOEt2)
						{
							minDrOEt2 = DROverET_1;
							FSR_PhotonToLep2 = FSR_Photon_TL;
							fsr_photonlep2 = true;
						}
					}
				}
			}
		}
		// Is the Mass better with FSR?
		bool noFSR = false;
		if ((abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M())) || (!fsr_photonlep1 && !fsr_photonlep2))
		{
			noFSR = true;
		}

		// Check if the event passes RelIso or RelIsoNoFSR -- pt>26 and trigger match included in MuonsVector function.
		pass_dimuon = true;
		//------------Rel iso Selection-------------------------------------------
		//------- If the events didn`t pass the rel iso, but pass reliso-gamma----
	/*	if ((RelIso1 >= RelIsomax && RelIso2 >= RelIsomax) || (RelIso1 >= RelIsomax && RelIso2 < RelIsomax) || (RelIso1 < RelIsomax && RelIso2 >= RelIsomax))
		{
			if (!usetrackiso && RelIso_NoFSR1 < RelIsomax && RelIso_NoFSR2 < RelIsomax)
			{
				pass_dimuon_RelIso = false;
			}
			if (usetrackiso && TrackIsoGamma(RelIso1, RelIso2, Reco_lepton1, Reco_lepton2, FSR_PhotonToLep1, FSR_PhotonToLep2))
			{
				pass_dimuon_RelIso = false;
			}
		}

		//-----------If the events pass the rel iso-----------------------------------
		if (RelIso1 < RelIsomax && RelIso2 < RelIsomax)
		{
			pass_dimuon_RelIso = true;
		}*/

		counters[strAllEvents] = counters[strAllEvents] + 1;
		if (pass_trigger)
		{
			counters[passTrigger] = counters[passTrigger] + 1;
		}
		if (pass_dimuon)
		{
			// Add Muon Histograms to AfterDimuonSelection
			histos.fillMuonHist(strAfterDimuonSelection, Reco_lepton1, Reco_lepton2, pu_weight);
			counters[passDimuon] = counters[passDimuon] + 1;
			// Classify the event!
			eventCategory = EventCategory(Reco_lepton1, Reco_lepton2);
			// Write Out Histograms according to category.
			string categoryEventTag = eventCategories[eventCategory];
			// Fill counter for All Events in category
			counters[categoryEventTag + strAllEvents] = counters[categoryEventTag + strAllEvents] + 1;
			//Counters if pass RelIso or Not
			string relIsoEventTag = pass_dimuon_RelIso ? strRelIso : strNonRelIso;

			if (pass_dimuon_RelIso || (!pass_dimuon_RelIso && passDimuonMass))
			{
				counters[relIsoEventTag] = counters[relIsoEventTag] + 1;
				counters[categoryEventTag + relIsoEventTag] = counters[categoryEventTag + relIsoEventTag] + 1;

				string fsrEventTag = noFSR ? strNOFSREvents : strFSREvents;
				// Counter if FSR used or not
				counters[categoryEventTag + relIsoEventTag + fsrEventTag] = counters[categoryEventTag + relIsoEventTag + fsrEventTag] + 1;

				//Fill All Events in Category
				histos.fillHiggsHist(categoryEventTag + strAllEvents, Dimuon, pu_weight);
				// Fill all RelIso and non RelIso Events
				if (noFSR)
				{
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
				}
				else if (fsr_photonlep1 && !fsr_photonlep2)
				{
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
					counters[categoryEventTag + relIsoEventTag + strFSREvents + oneGamma] = counters[categoryEventTag + relIsoEventTag + strFSREvents + oneGamma] + 1;

				}
				else if (!fsr_photonlep1 && fsr_photonlep2)
				{
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
					counters[categoryEventTag + relIsoEventTag + strFSREvents + oneGamma] = counters[categoryEventTag + relIsoEventTag + strFSREvents + oneGamma] + 1;
				}
				else if (fsr_photonlep1 && fsr_photonlep2)
				{
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
					histos.fillHiggsHist(categoryEventTag + relIsoEventTag + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
					counters[categoryEventTag + relIsoEventTag + strFSREvents + twoGamma] = counters[categoryEventTag + relIsoEventTag + strFSREvents + twoGamma] + 1;
				}
				else
				{
					cout << "Error check FSR" << endl;
				}
			}
			// Write AfterDimuonRelIso, when events pass RelIso, per category
		}
		//End of Event loop
	}
	//print counters
	for (auto counter : counters)
	{
		cout << RelIsotagstr << ";" << countersDescription[counter.first] << "|" << counter.second << endl;
	}
	histos.writeTFile();

cout<<"trigger pass:   "<< trigger_evt<<endl;
cout<<"Electron veto:   "<<electronveto<<endl;
cout<<"btag veto:   "<<btagveto<<endl;
}
BSM_Analysis::~BSM_Analysis()
{
	// do anything here that needs to be done at destruction time
}

//*********************************************************************************************
//  passRecoTrigger: This function give trigger decision, "true" or "false" for each event. 
//                   Receive the name of the trigger.
//*********************************************************************************************

bool BSM_Analysis::passRecoTrigger(string myTrigger1, string myTrigger2)
{
	for (int nTrig = 0; nTrig < Trigger_names->size(); ++nTrig)
	{
		string trigName = Trigger_names->at(nTrig);
		if (((trigName.find(myTrigger1) != string::npos) && (Trigger_decision->at(nTrig) == 1))
				|| ((trigName.find(myTrigger2) != string::npos) && (Trigger_decision->at(nTrig) == 1)))
		{
			return true;
		}
	}
	return false;
}

//*********************************************************************************************
//  MuonsVectors: Give a vector for each lepton (1 and 2). The selection taking into account 
//                the charge of the leptons in order to make pairs with opposite signs.
//*********************************************************************************************

double BSM_Analysis::Rhocorfactor(RoccoR& rc, int lepton_index)
{
	//RoccoR rc("rcdata.2016.v3");
	float u1 = gRandom->Rndm();
	float u2 = gRandom->Rndm();
	double SF = 1.0;
	int Gen_index = 0;
	bool Gen_muon = false;

	if (data == false)
	{
		for (int g = 0; g < Gen_pt->size(); g++)
		{

			if (Gen_status->at(g) == 1 && abs(Gen_pdg_id->at(g)) == 13 && Gen_charge->at(g) == Muon_charge->at(lepton_index))
			{
				if (deltaR(Gen_eta->at(g), Gen_phi->at(g), Muon_eta->at(lepton_index), Muon_phi->at(lepton_index)) < 0.005)
				{
					Gen_index = g;
					Gen_muon = true;
					break;
				}
			}
		}
		if (Gen_muon && Gen_pt->at(Gen_index) > 0)
		{
			//for MC, if matched gen-level muon (genPt) is available, use this function
			SF = rc.kScaleFromGenMC(Muon_charge->at(lepton_index), Muon_pt->at(lepton_index), Muon_eta->at(lepton_index), Muon_phi->at(lepton_index),
					Muon_TLayers->at(lepton_index), Gen_pt->at(Gen_index), u1, 0, 0);
		}
		else
		{
			//if not, then:
			//cout<<Muon_TLayers->at(lepton_index);
			SF = rc.kScaleAndSmearMC(Muon_charge->at(lepton_index), Muon_pt->at(lepton_index), Muon_eta->at(lepton_index), Muon_phi->at(lepton_index),
					Muon_TLayers->at(lepton_index), u1, u2, 0, 0);
		}
	}
	else
	{
//for each data muon in the loop, use this function to get a scale factor for its momentum:
		SF = rc.kScaleDT(Muon_charge->at(lepton_index), Muon_pt->at(lepton_index), Muon_eta->at(lepton_index), Muon_phi->at(lepton_index), 0, 0);
	}
	//cout<<"RoccoR:"<<SF<<endl;
	return SF;
}

bool BSM_Analysis::MuonsVectors(TLorentzVector& Reco_lepton1, TLorentzVector& Reco_lepton2, double& RelIso_NoFSR1, double& RelIso_NoFSR2, RoccoR& rc)
//, double& ChargedOverall1, double& ChargedOverall2,double& NeutralOverall1, double& NeutralOverall2)
{
	//float dimuon_pt_int = 0.;
	TLorentzVector Muon(0., 0., 0., 0.);
	TLorentzVector first_muon_vec(0., 0., 0., 0.);
	TLorentzVector Subfirst_muon_vec(0., 0., 0., 0.);
	vector<int> Muons;
	bool RelIso = false;

	for (int m = 0; m < Muon_pt->size(); m++)
	{
		double corr = Rhocorfactor(rc, m);
		if (Muon_medium->at(m) > 0 && Muon_isGlobal->at(m) > 0 && Muon_isTrackerMuon->at(m) > 0 && (Muon_pt->at(m)) * corr > 10 && abs(Muon_eta->at(m)) < 2.4
				&& Muon_combinedIso->at(m) < 0.25)
		{
			Muons.push_back(m);
		}
	}


	if (Muons.size() == 2)
	{
		double corr1 = Rhocorfactor(rc, Muons.at(0));
		double corr2 = Rhocorfactor(rc, Muons.at(1));
		if (Muon_charge->at(Muons.at(0)) * Muon_charge->at(Muons.at(1)) < 0)
		{
			if (Muon_pt->at(Muons.at(0)) > Muon_pt->at(Muons.at(1)))
			{
				if ((Muon_pt->at(Muons.at(0))  > 26 && Muon_isTriggerMatched->at(Muons.at(0)))||(Muon_pt->at(Muons.at(1))  > 26 && Muon_isTriggerMatched->at(Muons.at(1))))
				{
					first_muon_vec.SetPtEtaPhiE((Muon_pt->at(Muons.at(0)) ), Muon_eta->at(Muons.at(0)), Muon_phi->at(Muons.at(0)), Muon_energy->at(Muons.at(0)));
					Subfirst_muon_vec.SetPtEtaPhiE((Muon_pt->at(Muons.at(1))), Muon_eta->at(Muons.at(1)), Muon_phi->at(Muons.at(1)),
							Muon_energy->at(Muons.at(1)));
					Reco_lepton1 = first_muon_vec;
					Reco_lepton2 = Subfirst_muon_vec;
					RelIso = true;
				}
			}
			else
			{
				if ((Muon_pt->at(Muons.at(1))  > 26 && Muon_isTriggerMatched->at(Muons.at(1)))||(Muon_pt->at(Muons.at(0))  > 26 && Muon_isTriggerMatched->at(Muons.at(0))))
				{
					first_muon_vec.SetPtEtaPhiE((Muon_pt->at(Muons.at(1)) ), Muon_eta->at(Muons.at(1)), Muon_phi->at(Muons.at(1)), Muon_energy->at(Muons.at(1)));
					Subfirst_muon_vec.SetPtEtaPhiE((Muon_pt->at(Muons.at(0)) ), Muon_eta->at(Muons.at(0)), Muon_phi->at(Muons.at(0)),
							Muon_energy->at(Muons.at(0)));
					Reco_lepton1 = first_muon_vec;
					Reco_lepton2 = Subfirst_muon_vec;
					RelIso = true;
				}
			}

			RelIso_NoFSR1 = (Muon_isoCharged->at(Muons.at(0)) + max(0., Muon_isoNeutralHadron->at(Muons.at(0)) - 0.5 * (Muon_isoPU->at(Muons.at(0)))))
					/ Muon_pt->at(Muons.at(0));
			RelIso_NoFSR2 = (Muon_isoCharged->at(Muons.at(1)) + max(0., Muon_isoNeutralHadron->at(Muons.at(1)) - 0.5 * (Muon_isoPU->at(Muons.at(1)))))
					/ Muon_pt->at(Muons.at(1));
		}
	}
	return RelIso;
/*
	for (int m1 = 0; m1 < Muon_pt->size(); m1++)
	{
		if (Muon_pt->size() > 1)
		{
			for (int m2 = 0; m2 < Muon_pt->size(); m2++)
			{
				double corr1 = Rhocorfactor(rc, m1);
				double corr2 = Rhocorfactor(rc, m2);
				if (Muon_medium->at(m1) > 0 && Muon_isGlobal->at(m1) > 0 && Muon_isTrackerMuon->at(m1) > 0 && Muon_pt->at(m1) > 10 && abs(Muon_eta->at(m1)) < 2.4
						&& Muon_medium->at(m2) > 0 && Muon_isGlobal->at(m2) > 0 && Muon_isTrackerMuon->at(m2) > 0 && Muon_pt->at(m2) > 10 && abs(Muon_eta->at(m2)) < 2.4)
				{
					if (Muon_charge->at(m1) * Muon_charge->at(m2) < 0)
					{
						if (Muon_pt->at(m1) > Muon_pt->at(m2))
						{
							first_muon_vec.SetPtEtaPhiE((Muon_pt->at(m1)) * corr1, Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
							Subfirst_muon_vec.SetPtEtaPhiE((Muon_pt->at(m2)) * corr2, Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
						}
						else
						{
							Subfirst_muon_vec.SetPtEtaPhiE((Muon_pt->at(m1)) * corr1, Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
							first_muon_vec.SetPtEtaPhiE((Muon_pt->at(m2)) * corr2, Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
						}
						float dimuon_pt = first_muon_vec.Pt() + Subfirst_muon_vec.Pt();

						if (dimuon_pt > dimuon_pt_int
								&& (((Muon_pt->at(m1)) * corr1 > 26 && Muon_isTriggerMatched->at(m1)) || ((Muon_pt->at(m2)) * corr2 > 26 && Muon_isTriggerMatched->at(m2))))
						{
							first_muon_vec.SetPtEtaPhiE((Muon_pt->at(m1)) * corr1, Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
							Subfirst_muon_vec.SetPtEtaPhiE((Muon_pt->at(m2)) * corr2, Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
							Reco_lepton1 = first_muon_vec;
							Reco_lepton2 = Subfirst_muon_vec;
							dimuon_pt_int = dimuon_pt;
							if (usetrackiso == false)
							{
								RelIso1 = Muon_combinedIso->at(m1);				//Muon_isoSum->at(m1) / Muon_pt->at(m1);
								RelIso2 = Muon_combinedIso->at(m2);				//Muon_isoSum->at(m2) / Muon_pt->at(m2);
							}
							if (usetrackiso == true)
							{
								RelIso1 = Muon_trackRe_iso->at(m1);			//track iso
								RelIso2 = Muon_trackRe_iso->at(m2);	    //track iso
							}

							RelIso_NoFSR1 = (Muon_isoCharged->at(m1) + max(0., Muon_isoNeutralHadron->at(m1) - 0.5 * (Muon_isoPU->at(m1)))) / Muon_pt->at(m1);
							RelIso_NoFSR2 = (Muon_isoCharged->at(m2) + max(0., Muon_isoNeutralHadron->at(m2) - 0.5 * (Muon_isoPU->at(m2)))) / Muon_pt->at(m2);

							ChargedOverall1 = Muon_isoCharged->at(m1) / Muon_pt->at(m1);
							 ChargedOverall2 = Muon_isoCharged->at(m2) / Muon_pt->at(m2);
							 NeutralOverall1 = Muon_isoNeutralHadron->at(m1) / Muon_pt->at(m1);
							 NeutralOverall2 = Muon_isoNeutralHadron->at(m2) / Muon_pt->at(m2);

						}
					}
				}
			}
		}
	}*/
}

bool BSM_Analysis::TrackIsoGamma(double& RelIso1, double& RelIso2, TLorentzVector& Muon1, TLorentzVector& Muon2, TLorentzVector& PFPhoton1,
		TLorentzVector& PFPhoton2)
{
	bool noelectron = false;
	bool nojet = false;
	bool onlyfsr = false;

	TLorentzVector Electron(0., 0., 0., 0);
	for (int e1 = 0; e1 < patElectron_pt->size(); e1++)
	{
		Electron.SetPtEtaPhiE(patElectron_pt->at(e1), patElectron_eta->at(e1), patElectron_phi->at(e1), patElectron_energy->at(e1));
		if (Electron.DeltaR(Muon1) >= 0.3 && Electron.DeltaR(Muon2) >= 0.3)
			noelectron = true;
	}
	if (patElectron_pt->size() == 0)
		noelectron = true;

	TLorentzVector Jet(0., 0., 0., 0);
	for (int j1 = 0; j1 < Jet_pt->size(); j1++)
	{
		Jet.SetPtEtaPhiE(Jet_pt->at(j1), Jet_eta->at(j1), Jet_phi->at(j1), Jet_energy->at(j1));
		if (Jet.DeltaR(Muon1) >= 0.3 && Jet.DeltaR(Muon2) >= 0.3)
			nojet = true;
	}
	if (Jet_pt->size() == 0)
		nojet = true;

	if ((PFPhoton1.DeltaR(Muon1) < 0.3 && PFPhoton2.DeltaR(Muon2) < 0.3) || (PFPhoton2.DeltaR(Muon1) < 0.3 && PFPhoton1.DeltaR(Muon2) < 0.3))
		onlyfsr = true;

	if (noelectron && nojet && onlyfsr)
		return true;
	return false;
}

//*********************************************************************************************
// PhotonsVectors: this function gives the vector that contains all the vectors of the photons
//                 produced in each event
//*********************************************************************************************

void BSM_Analysis::PhotonsVectors(vector<TLorentzVector>& Photon_vec)
{
	TLorentzVector Photon_TL(0., 0., 0., 0.);
	for (int ph = 0; ph < Photon_pt->size(); ph++)
	{
		if (Photon_pt->size() > 0)
		{
			Photon_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
			Photon_vec.push_back(Photon_TL);
		}
	}
}

//*********************************************************************************************
// GenleptonVectors: this function gives the vector for each generated lepton
//*********************************************************************************************

void BSM_Analysis::GenleptonVector(TLorentzVector& Gen_lepton1_vec, TLorentzVector& Gen_lepton2_vec)
{
	float dimuon_mass_int = 999999.;
	double Gen_lepton1_charge;
	double Gen_lepton2_charge;
	TLorentzVector Gen_lepton_vec(0., 0., 0., 0.);
	TLorentzVector Gen_antilepton_vec(0., 0., 0., 0.);

	for (int g1 = 0; g1 < Gen_pt->size(); g1++)
	{
		if (Gen_pt->size() > 1 && (Gen_status->at(g1) == 1 || Gen_status->at(g1) == 2))
		{
			if (Gen_pdg_id->at(g1) == 13)
			{
				Gen_lepton_vec.SetPtEtaPhiE(Gen_pt->at(g1), Gen_eta->at(g1), Gen_phi->at(g1), Gen_energy->at(g1));
				Gen_lepton1_charge = Gen_charge->at(g1);
			}

			if (Gen_pdg_id->at(g1) == -13)
			{
				Gen_antilepton_vec.SetPtEtaPhiE(Gen_pt->at(g1), Gen_eta->at(g1), Gen_phi->at(g1), Gen_energy->at(g1));
				Gen_lepton2_charge = Gen_charge->at(g1);
			}
			if (Gen_lepton1_charge * Gen_lepton2_charge < 0)
			{
				float dimuon_gen_mass = (Gen_lepton_vec + Gen_antilepton_vec).M();
				if (dimuon_gen_mass < dimuon_mass_int)
				{
					Gen_lepton1_vec = Gen_lepton_vec;
					Gen_lepton2_vec = Gen_antilepton_vec;
					dimuon_mass_int = dimuon_gen_mass;
				}
			}
		}
	}
}

//*********************************************************************************************
// GenPhotonVectors: this function gives the vector for each generated photon
//*********************************************************************************************

void BSM_Analysis::PhotonsGenVectors(vector<TLorentzVector>& Photon_Gen_vec, TLorentzVector& Photon_Gen_TL)
{
	for (int ph = 0; ph < Gen_pt->size(); ph++)
	{
		if (Gen_pt->size() > 0 && (Gen_status->at(ph) == 1 || Gen_status->at(ph) == 2))
		{
			if (Gen_pdg_id->at(ph) == 22)
			{
				Photon_Gen_TL.SetPtEtaPhiE(Gen_pt->at(ph), Gen_eta->at(ph), Gen_phi->at(ph), Gen_energy->at(ph));
				Photon_Gen_vec.push_back(Photon_Gen_TL);
			}
		}
	}
}

//*********************************************************************************************
//FSRdeltaR: Gives the Delta R between a muon and the PF photon
//*********************************************************************************************

double BSM_Analysis::FSRdeltaR(TLorentzVector& FSRPhotonVec, TLorentzVector& LeptonVec)
{
	double FSRDR;
	FSRDR = deltaR(LeptonVec.Eta(), LeptonVec.Phi(), FSRPhotonVec.Eta(), FSRPhotonVec.Phi());
	return FSRDR;
}

//*********************************************************************************************
//FSRDROverET2: Gives the radio between DeltaR and Transversal Energy
//*********************************************************************************************

double BSM_Analysis::FSRDROverET2(TLorentzVector& FSRPhotonVec, double fsrdr)
{
	double fsrDrOEt2;
	fsrDrOEt2 = fsrdr / (FSRPhotonVec.Pt() * FSRPhotonVec.Pt());
	return fsrDrOEt2;
}

bool BSM_Analysis::BtagVeto()
{
	bool bTaggedJetVeto = false;
	for (int m1 = 0; m1 < Jet_pt->size(); m1++)
	{
		if (Jet_pt->at(m1) > 30.0 && abs(Jet_eta->at(m1)) < 2.4 && Jet_bDiscriminator_pfCMVAV2->at(m1) > 0.8484)
		{
			bTaggedJetVeto = true;
		}
	}
	return bTaggedJetVeto;
}
int BSM_Analysis::EventCategory(TLorentzVector& Muon1, TLorentzVector& Muon2)
{
	TLorentzVector Jet_temp(0., 0., 0., 0.);
	TLorentzVector Jet1(0., 0., 0., 0.);
	TLorentzVector Jet2(0., 0., 0., 0.);
	vector<TLorentzVector> Jets_passSel_vec;
	Jets_passSel_vec.erase(Jets_passSel_vec.begin(), Jets_passSel_vec.end());
	bool VBF_Tight = false;
	bool VBF_Loose = false;
	bool ggF_Tight = false;
	bool b0_1_jet_Tight = false;
	bool b0_1_jet_Loose = false;
	float dijet_mass = 0.0;
	//Jet selection:

	if (Jet_pt->size() > 1)
	{
		for (int m1 = 0; m1 < Jet_pt->size(); m1++)
		{
			Jet_temp.SetPtEtaPhiE(Jet_pt->at(m1), Jet_eta->at(m1), Jet_phi->at(m1), Jet_energy->at(m1));
			if (Jet_pt->at(m1) > 30.0 && abs(Jet_eta->at(m1)) < 4.7 && (Jet_temp.DeltaR(Muon1) > 0.4 && Jet_temp.DeltaR(Muon2) > 0.4))
			{
				Jets_passSel_vec.push_back(Jet_temp);
			}
		}
	}

	if (Jets_passSel_vec.size() > 1)
	{
		for (int j1 = 0; j1 < Jets_passSel_vec.size(); j1++)
		{
			for (int j2 = 0; j2 < Jets_passSel_vec.size(); j2++)
			{
				if (j2 != j1)	    //&& (Jets_passSel_vec[j1].Eta() * Jets_passSel_vec[j2].Eta() < 0))
				{
					if (Jets_passSel_vec[j1].Pt() > 40.0 || Jets_passSel_vec[j2].Pt() > 40.0)
					{
						// Paso la preseleccion
						Jet1 = Jets_passSel_vec[j1];
						Jet2 = Jets_passSel_vec[j2];
						dijet_mass = (Jet1 + Jet2).M();
						// VBF Tight selection
						if (dijet_mass > 650.0 && abs(Jet1.Eta() - Jet2.Eta()) > 3.5)
						{
							VBF_Tight = true;
						}
						// ggF Tight selection
						else if (dijet_mass > 250.0 && (Muon1 + Muon2).Pt() > 50.0)
						{
							ggF_Tight = true;
						}
						else
						{
							VBF_Loose = true;
						}
					}
				}
			}
		}
	}
	// No pasa preselection!

	if ((Muon1 + Muon2).Pt() >= 25.0)
	{
		b0_1_jet_Tight = true;
	}
	else
	{
		b0_1_jet_Loose = true;
	}

	if (VBF_Tight == true)
		return 1;
	if (ggF_Tight == true)
		return 3;
	if (VBF_Loose == true)
		return 2;
	if (b0_1_jet_Tight == true)
		return 4;
	if (b0_1_jet_Loose == true)
		return 5;
	return 0;
}

//*********************************************************************************************
// setBranchAdress: initialization of each object pointer. Also Set branch addresses and branch 
//                  pointers.
//*********************************************************************************************

bool BSM_Analysis::Electron_veto(TLorentzVector& Muon1, TLorentzVector& Muon2)
{
	TLorentzVector electron(0., 0., 0., 0.);
	bool passElectronVeto = false;
	for (int e1 = 0; e1 < patElectron_pt->size(); e1++)
	{
		electron.SetPtEtaPhiE(patElectron_pt->at(e1), patElectron_eta->at(e1), patElectron_phi->at(e1), patElectron_energy->at(e1));
		//cout<< patElectron_isPassMedium->at(e1)<<","<<electron.Pt()<<","<<electron.Eta()<<","<<Muon1.DeltaR(electron)<<","<<Muon2.DeltaR(electron)<<endl;
		if (patElectron_isPassMedium->at(e1) > 0 && electron.Pt() > 10)
		{
			if (abs(electron.Eta()) < 1.4442 || (abs(electron.Eta()) > 1.566 && abs(electron.Eta()) < 2.5))
			{
				if (Muon1.DeltaR(electron) > 0.4 && Muon2.DeltaR(electron) > 0.4)
				{
					passElectronVeto = true;
					break;
				}
			}
		}
	}
	return passElectronVeto;
}

void BSM_Analysis::setBranchAddress(TTree* BOOM)
{
	// Set object pointer
	Muon_isTriggerMatched = 0;
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
	Muon_isGlobal = 0;
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
	Muon_combinedIso = 0;
	Muon_trackRe_iso = 0;

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

	FSRPhoton_pt = 0;
	FSRPhoton_eta = 0;
	FSRPhoton_phi = 0;
	FSRPhoton_energy = 0;
	FSRPhoton_pX = 0;
	FSRPhoton_pY = 0;
	FSRPhoton_pZ = 0;
	FSRPhoton_isoNH = 0;
	FSRPhoton_isoCH = 0;
	FSRPhoton_isoCHPU = 0;
	FSRPhoton_isoPhot = 0;
	FSRPhoton_et = 0;
	FSRPhoton_isoNHPhot = 0;

	Jet_pt = 0;
	Jet_eta = 0;
	Jet_phi = 0;
	Jet_energy = 0;

	Jet_mass = 0;
	Jet_neutralHadEnergyFraction = 0;
	Jet_neutralEmEmEnergyFraction = 0;
	Jet_chargedHadronEnergyFraction = 0;
	Jet_chargedEmEnergyFraction = 0;
	Jet_muonEnergyFraction = 0;
	Jet_electronEnergy = 0;
	Jet_photonEnergy = 0;
	Jet_bDiscriminator_pfCMVAV2 = 0;
	Jet_puppi_bDiscriminator_pfCMVAV2 = 0;
	Jet_puppi_pt = 0;

	Gen_eta = 0;
	Gen_phi = 0;
	Gen_pt = 0;
	Gen_energy = 0;
	Gen_pdg_id = 0;
	Gen_motherpdg_id = 0;
	Gen_status = 0;
	Gen_BmotherIndex = 0;
	Gen_charge = 0;

	UncorrJet_pt = 0;
	Trigger_names = 0;
	Trigger_decision = 0;
	ntruePUInteractions = 0;

	Met_type1PF_sumEt = 0;
	Met_type1PF_pt = 0;

	pvertex_z = 0;
	eventNumber = 0;

	// Set branch addresses and branch pointers
	if (!BOOM)
		return;
	BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
	BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
	BOOM->SetBranchAddress("nTruePUInteractions", &ntruePUInteractions, &b_ntruePUInteractions);
	BOOM->SetBranchAddress("Muon_isTriggerMatched", &Muon_isTriggerMatched, &b_Muon_isTriggerMatched);
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
	BOOM->SetBranchAddress("Muon_combinedIso", &Muon_combinedIso, &b_Muon_combinedIso);
	BOOM->SetBranchAddress("Muon_trackRe_iso", &Muon_trackRe_iso, &b_Muon_trackRe_iso);
	BOOM->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal, &b_Muon_isGlobal);

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

	BOOM->SetBranchAddress("FSRPhoton_pt", &FSRPhoton_pt, &b_FSRPhoton_pt);
	BOOM->SetBranchAddress("FSRPhoton_eta", &FSRPhoton_eta, &b_FSRPhoton_eta);
	BOOM->SetBranchAddress("FSRPhoton_phi", &FSRPhoton_phi, &b_FSRPhoton_phi);
	BOOM->SetBranchAddress("FSRPhoton_energy", &FSRPhoton_energy, &b_FSRPhoton_energy);
	BOOM->SetBranchAddress("FSRPhoton_pX", &FSRPhoton_pX, &b_FSRPhoton_pX);
	BOOM->SetBranchAddress("FSRPhoton_pY", &FSRPhoton_pY, &b_FSRPhoton_pY);
	BOOM->SetBranchAddress("FSRPhoton_pZ", &FSRPhoton_pZ, &b_FSRPhoton_pZ);
	BOOM->SetBranchAddress("FSRPhoton_isoNH", &FSRPhoton_isoNH, &b_FSRPhoton_isoNH);
	BOOM->SetBranchAddress("FSRPhoton_isoCH", &FSRPhoton_isoCH, &b_FSRPhoton_isoCH);
	BOOM->SetBranchAddress("FSRPhoton_isoCHPU", &FSRPhoton_isoCHPU, &b_FSRPhoton_isoCHPU);
	BOOM->SetBranchAddress("FSRPhoton_isoPhot", &FSRPhoton_isoPhot, &b_FSRPhoton_isoPhot);
	BOOM->SetBranchAddress("FSRPhoton_et", &FSRPhoton_et, &b_FSRPhoton_et);
	BOOM->SetBranchAddress("FSRPhoton_isoNHPhot", &FSRPhoton_isoNHPhot, &b_FSRPhoton_isoNHPhot);

	BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
	BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
	BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
	BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
	BOOM->SetBranchAddress("Jet_bDiscriminator_pfCMVAV2", &Jet_bDiscriminator_pfCMVAV2, &b_Jet_bDiscriminator_pfCMVAV2);
	BOOM->SetBranchAddress("Jet_puppi_bDiscriminator_pfCMVAV2", &Jet_puppi_bDiscriminator_pfCMVAV2, &b_Jet_puppi_bDiscriminator_pfCMVAV2);
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
	BOOM->SetBranchAddress("Met_type1PF_pt", &Met_type1PF_pt, &b_Met_type1PF_pt);
	BOOM->SetBranchAddress("Met_type1PF_sumEt", &Met_type1PF_sumEt, &b_Met_type1PF_sumEt);
	BOOM->SetBranchAddress("Met_puppi_pt", &Met_puppi_pt, &b_Met_puppi_pt);
	BOOM->SetBranchAddress("Met_puppi_sumEt", &Met_puppi_sumEt, &b_Met_puppi_sumEt);
	BOOM->SetBranchAddress("Met_puppi_phi", &Met_puppi_phi, &b_Met_puppi_phi);
	BOOM->SetBranchAddress("Met_puppi_px", &Met_puppi_px, &b_Met_puppi_px);
	BOOM->SetBranchAddress("Met_puppi_py", &Met_puppi_py, &b_Met_puppi_py);
	BOOM->SetBranchAddress("Met_puppi_pz", &Met_puppi_pz, &b_Met_puppi_pz);
	BOOM->SetBranchAddress("Jet_puppi_pt", &Jet_puppi_pt, &b_Jet_puppi_pt);

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
	BOOM->SetBranchAddress("pvertex_z", &pvertex_z, &b_pvertex_z);
	BOOM->SetBranchAddress("eventNumber", &eventNumber , &b_eventNumber);
}
;
