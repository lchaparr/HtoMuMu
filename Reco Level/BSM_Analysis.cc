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
bool data = false;

int main(int argc, char *argv[])
{
	cout << argv[1] << endl;
	BSM_Analysis BSM_Analysis_(NULL, NULL, 0, argv[1], argv[3]);
}

BSM_Analysis::BSM_Analysis(TFile* theFile, TDirectory *cdDir[], int nDir, char* fname, char* dataType)
{
	if (dataType == std::string("GF"))
	{
		Gluon_Fusion = true;
		cout << "is Gluon Fusion" << endl;
	}
	if (dataType == std::string("VBF"))
	{
		Vector_Boson = true;
		cout << "is VBF" << endl;
	}
	if (dataType == std::string("WP"))
		W_PlusH = true;
	if (dataType == std::string("WM"))
		W_MinusH = true;
	if (dataType == std::string("ZH"))
		ZH = true;
	if (dataType == std::string("ttH"))
		ttH = true;
	if (dataType == std::string("DATA"))
	{
		cout << "is data or BKG" << endl;
		data = true;
	}

	HistogramManager histos("H2Mu.root");
	string strAfterDimuonSelection = "AfterDimuonSelection";
	string strGluonFusion = "GluonFusion";
	string strGFBBTight = "/GFBBTight";
	string strGFBEBOTight = "/GFBEBOTight";
	string strGFEEOOTight = "/GFEEOOTight";
	string strGFBBLoose = "/GFBBLoose";
	string strGFBEBOLoose = "/GFBEBOLoose";
	string strGFEEOOLoose = "/GFEEOOLoose";
	string strGFTwoJets = "/GFTwoJets";
	string strVectorBoson = "VectorBoson";
	string strVBFTight = "/VBFTight";
	string strVBFLoose = "/VBFLoose";
	string strAfterDimuonRelIso = "/AfterDimuonRelIso";

	string strRelIso = "/RelIso";
	string strNonRelIso = "/NonRelIso";
	string strAllEvents = "/AllEvents"; //all the selected events, including rel iso and no rel iso.
	string strFSREvents = "/FSREvents";
	string strNOFSREvents = "/NOFSREvents";
	string tagOnlyDimuon = "/OnlyDimuon";
	string tagDimuonPlusPH = "/DimuonPlusPH";

	string tagInvariantMass = "OnlyDimuon DimuonPlusPH";
	string tagAllEvents = "AllEvents";
	string tagPhotons = "/Photons";

	string strtest = "testDirectory";
	string tagothers = "/others";

	histos.createMuonHistograms(strAfterDimuonSelection);
	if (Gluon_Fusion || data)
	{
		histos.createMuonHistograms(strGluonFusion + strAfterDimuonRelIso);
		histos.createDimuonHistograms(strGluonFusion + strGFBBTight, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFBBTight + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBTight + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBTight + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOTight, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOTight + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOTight + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOTight, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOTight + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOTight + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBLoose, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFBBLoose + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBLoose + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOLoose, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOLoose + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOLoose + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOLoose, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOLoose + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOLoose + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFTwoJets, tagAllEvents);
		histos.createDimuonHistograms(strGluonFusion + strGFTwoJets + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFTwoJets + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents, tagInvariantMass);
	}
	if (Vector_Boson || data)
	{
		histos.createMuonHistograms(strVectorBoson + strAfterDimuonRelIso);
		histos.createDimuonHistograms(strVectorBoson + strVBFTight, tagAllEvents);
		histos.createDimuonHistograms(strVectorBoson + strVBFTight + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFTight + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFTight + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFLoose, tagAllEvents);
		histos.createDimuonHistograms(strVectorBoson + strVBFLoose + strRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFLoose + strRelIso + strFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFLoose + strNonRelIso + strNOFSREvents, tagInvariantMass);
		histos.createDimuonHistograms(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents, tagInvariantMass);
	}

	string strIsolation1 = "Isolation1";
	string strIsolation2 = "Isolation2";
	histos.addHistogram(strtest + tagothers, strIsolation1, "RelIso1", 100, 0, 2);
	histos.addHistogram(strtest + tagothers, strIsolation2, "RelIso2", 100, 0, 2);

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
	//int prueba = 100;

	//-------------------General counters---------------------------------------

	int trigger_evt = 0, dimuon_evt = 0, dimuon_initialsel_evt = 0, dimuon_passRelIso = 0;
	int AllFSRPhotons = 0;
	int fsr_photons = 0;
	int BothNoreliso = 0, OneNoreliso = 0, RelIsoReal = 0;
	int PassRelisoNoFSR = 0;

	//·················Gluon Fusion or data counters...........................
	int GF_BB_Tight_cat = 0, GF_BBTight_Iso = 0, GF_BBTightIso_NOFSR = 0, GF_BBTightIso_OneFSR = 0, GF_BBTightIso_TwoFSR = 0, GF_BBTightNoIso_NOFSR = 0,
			GF_BBTightNoIso_OneFSR = 0, GF_BBTightNoIso_TwoFSR = 0;
	int GF_BEBO_Tight_cat = 0, GF_BEBOTight_Iso = 0, GF_BEBOTightIso_NOFSR = 0, GF_BEBOTightIso_OneFSR = 0, GF_BEBOTightIso_TwoFSR = 0,
			GF_BEBOTightNoIso_NOFSR = 0, GF_BEBOTightNoIso_OneFSR = 0, GF_BEBOTightNoIso_TwoFSR = 0;
	int GF_EEOO_Tight_cat = 0, GF_EEOOTight_Iso = 0, GF_EEOOTightIso_NOFSR = 0, GF_EEOOTightIso_OneFSR = 0, GF_EEOOTightIso_TwoFSR = 0,
			GF_EEOOTightNoIso_NOFSR = 0, GF_EEOOTightNoIso_OneFSR = 0, GF_EEOOTightNoIso_TwoFSR = 0;
	int GF_BB_Loose_cat = 0, GF_BBLoose_Iso = 0, GF_BBLooseIso_NOFSR = 0, GF_BBLooseIso_OneFSR = 0, GF_BBLooseIso_TwoFSR = 0, GF_BBLooseNoIso_NOFSR = 0,
			GF_BBLooseNoIso_OneFSR = 0, GF_BBLooseNoIso_TwoFSR = 0;
	int GF_BEBO_Loose_cat = 0, GF_BEBOLoose_Iso = 0, GF_BEBOLooseIso_NOFSR = 0, GF_BEBOLooseIso_OneFSR = 0, GF_BEBOLooseIso_TwoFSR = 0,
			GF_BEBOLooseNoIso_NOFSR = 0, GF_BEBOLooseNoIso_OneFSR = 0, GF_BEBOLooseNoIso_TwoFSR = 0;
	int GF_EEOO_Loose_cat = 0, GF_EEOOLoose_Iso = 0, GF_EEOOLooseIso_NOFSR = 0, GF_EEOOLooseIso_OneFSR = 0, GF_EEOOLooseIso_TwoFSR = 0,
			GF_EEOOLooseNoIso_NOFSR = 0, GF_EEOOLooseNoIso_OneFSR = 0, GF_EEOOLooseNoIso_TwoFSR = 0;
	int GF_TwoJets_cat = 0, GF_TwoJets_Iso = 0, GF_TwoJetsIso_NOFSR = 0, GF_TwoJetsIso_OneFSR = 0, GF_TwoJetsIso_TwoFSR = 0, GF_TwoJetsNoIso_NOFSR = 0,
			GF_TwoJetsNoIso_OneFSR = 0, GF_TwoJetsNoIso_TwoFSR = 0;

	//·················Vector_Boson or data counters···························
	int VBF_Tight_cat = 0, VBF_Tight_Iso = 0, VBF_TightIso_NOFSR = 0, VBF_TightIso_OneFSR = 0, VBF_TightIso_TwoFSR = 0, VBF_TightNoIso_NOFSR = 0,
			VBF_TightNoIso_OneFSR = 0, VBF_TightNoIso_TwoFSR = 0;
	int VBF_Loose_cat = 0, VBF_Loose_Iso = 0, VBF_LooseIso_NOFSR = 0, VBF_LooseIso_OneFSR = 0, VBF_LooseIso_TwoFSR = 0, VBF_LooseNoIso_NOFSR = 0,
			VBF_LooseNoIso_OneFSR = 0, VBF_LooseNoIso_TwoFSR = 0;

	//---------------------Loop over the events----------------------------
	//for (int i = 0; i < prueba; ++i)
	for (int i = 0; i < nentries; ++i)
	{
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
		int VBF_cat = 0, GF_cat = 0;

		//----------------------control booleans------------------------
		bool pass_dimuon = false, pass_dimuon_RelIso = false;

		//-------------------For Trigger-------------------------------
		if (passRecoTrigger(myTrigger1, myTrigger2))
			trigger_evt++;

		//--------------------Dimuons (1st selection)------------------
		MuonsVectors(Reco_lepton1, Reco_lepton2, RelIso1, RelIso2, RelIso_NoFSR1, RelIso_NoFSR2);

		//--------------------FSR Recovery algorithm---------------------
		double DR_FSR_lep1, DR_FSR_lep2, DROverET_1, DROverET_2, minDrOEt1 = 999.9, minDrOEt2 = 999.9;
		bool fsr_photonlep1 = false, fsr_photonlep2 = false;

		for (int ph = 0; ph < FSRPhoton_pt->size(); ph++)
		{
			FSR_Photon_TL.SetPtEtaPhiE(FSRPhoton_pt->at(ph), FSRPhoton_eta->at(ph), FSRPhoton_phi->at(ph), FSRPhoton_energy->at(ph));
			FSRIso = (FSRPhoton_isoCHPU->at(ph) + FSRPhoton_isoNHPhot->at(ph)) / FSRPhoton_pt->at(ph);
			AllFSRPhotons++;
			if (FSR_Photon_TL.Pt() > 2.0 && abs(FSR_Photon_TL.Eta()) < 2.4 && FSRIso < 1.8)
			{
				fsr_photons++;
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

		//---------------------Muons selections--------------------------
		if ((Reco_lepton1 + Reco_lepton2).M() != 0. && Reco_lepton1.M() != 0 && Reco_lepton2.M() != 0)
		{
			Dimuon = Reco_lepton1 + Reco_lepton2;
			dimuon_evt++;
			//Only for GF categories or data
			if (Gluon_Fusion || data)
			{
				GF_cat = GFCategories(Reco_lepton1, Reco_lepton2);
				if (GF_cat == 1)
					GF_BB_Tight_cat++;
				if (GF_cat == 2)
					GF_BEBO_Tight_cat++;
				if (GF_cat == 3)
					GF_EEOO_Tight_cat++;
				if (GF_cat == 4)
					GF_BB_Loose_cat++;
				if (GF_cat == 5)
					GF_BEBO_Loose_cat++;
				if (GF_cat == 6)
					GF_EEOO_Loose_cat++;
				if (GF_cat == 7)
					GF_TwoJets_cat++;
			}
			//Only for VBF or data
			if (Vector_Boson || data)
			{
				VBF_cat = JetsVector(Reco_lepton1, Reco_lepton2);
				if (VBF_cat == 1)
					VBF_Tight_cat++;
				if (VBF_cat == 2)
					VBF_Loose_cat++;
			}

			if (Reco_lepton1.Pt() > 26 && abs(Reco_lepton1.Eta()) < 2.4 && passRecoTrigger(myTrigger1, myTrigger2))
			{
				if (Reco_lepton2.Pt() > 10 && abs(Reco_lepton2.Eta()) < 2.4)
				{
					histos.fillMuonHist(strAfterDimuonSelection, Reco_lepton1, Reco_lepton2, pu_weight);
					dimuon_initialsel_evt++;
					pass_dimuon = true;

					//------------Rel iso Selection-------------------------------------------
					//------- If the events didn`t pass the rel iso, but pass reliso-gamma----
					if (RelIso1 >= 0.05 && RelIso2 >= 0.05)
					{
						histos.fillHistogram(strtest + tagothers, strIsolation1, RelIso1, pu_weight);
						histos.fillHistogram(strtest + tagothers, strIsolation2, RelIso2, pu_weight);
						BothNoreliso++;
						if (TrackIsoGamma(RelIso1, RelIso2, Reco_lepton1, Reco_lepton2, FSR_PhotonToLep1, FSR_PhotonToLep2))
						{
							pass_dimuon_RelIso = true;
							PassRelisoNoFSR++;
						}
						/*if (RelIso_NoFSR1 < 0.25 && RelIso_NoFSR2 < 0.25)
						 {
						 pass_dimuon_RelIso = true;
						 PassRelisoNoFSR++;
						 }*/
					}

					if ((RelIso1 >= 0.05 && RelIso2 < 0.05) || (RelIso1 < 0.05 && RelIso2 >= 0.05))
					{
						OneNoreliso++;
						if (TrackIsoGamma(RelIso1, RelIso2, Reco_lepton1, Reco_lepton2, FSR_PhotonToLep1, FSR_PhotonToLep2))
						{
							pass_dimuon_RelIso = true;
							PassRelisoNoFSR++;
						}
						/*if (RelIso_NoFSR1 < 0.25 && RelIso_NoFSR2 < 0.25)
						 {
						 pass_dimuon_RelIso = true;
						 PassRelisoNoFSR++;
						 }*/
					}
					//-----------If the events pass the rel iso-----------------------------------
					if (RelIso1 < 0.05 && RelIso2 < 0.05)
					{
						RelIsoReal++;
						pass_dimuon_RelIso = true;
					}

					//---------------For all the events that pass relIso-----------------------
					if (pass_dimuon_RelIso)
					{
						dimuon_passRelIso++;

						if (Gluon_Fusion || data)
						{
							histos.fillMuonHist(strGluonFusion + strAfterDimuonRelIso, Reco_lepton1, Reco_lepton2, pu_weight);
							{
								if (GF_cat == 1)
								{
									GF_BBTight_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_BBTightIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_BBTightIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_BBTightIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_BBTightIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 2)
								{
									GF_BEBOTight_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_BEBOTightIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_BEBOTightIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_BEBOTightIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_BEBOTightIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 3)
								{
									GF_EEOOTight_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_EEOOTightIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_EEOOTightIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_EEOOTightIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_EEOOTightIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 4)
								{
									GF_BBLoose_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_BBLooseIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_BBLooseIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_BBLooseIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_BBLooseIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 5)
								{
									GF_BEBOLoose_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_BEBOLooseIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_BEBOLooseIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_BEBOLooseIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_BEBOLooseIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 6)
								{
									GF_EEOOLoose_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_EEOOLooseIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_EEOOLooseIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_EEOOLooseIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_EEOOLooseIso_TwoFSR++;
										}
									}
								}
								if (GF_cat == 7)
								{
									GF_TwoJets_Iso++;
									if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
									{
										histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										GF_TwoJetsIso_NOFSR++;
									}
									else
									{
										if (fsr_photonlep1 && !fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
											GF_TwoJetsIso_OneFSR++;
										}
										if (fsr_photonlep2 && !fsr_photonlep1)
										{
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
											GF_TwoJetsIso_OneFSR++;
										}
										if (fsr_photonlep1 && fsr_photonlep2)
										{
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
											histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
													FSR_PhotonToLep2, pu_weight);
											GF_TwoJetsIso_TwoFSR++;
										}
									}
								}
							}
						}

						if (Vector_Boson || data)
						{
							histos.fillMuonHist(strVectorBoson + strAfterDimuonRelIso, Reco_lepton1, Reco_lepton2, pu_weight);
							if (VBF_cat == 1)
							{
								VBF_Tight_Iso++;
								if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
								{
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									VBF_TightIso_NOFSR++;
								}
								else
								{
									if (fsr_photonlep1 && !fsr_photonlep2)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
										VBF_TightIso_OneFSR++;
									}
									if (fsr_photonlep2 && !fsr_photonlep1)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
										VBF_TightIso_OneFSR++;
									}
									if (fsr_photonlep1 && fsr_photonlep2)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFTight + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										VBF_TightIso_TwoFSR++;
									}
								}
							}

							if (VBF_cat == 2)
							{
								VBF_Loose_Iso++;
								if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
								{
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									VBF_LooseIso_NOFSR++;
								}
								else
								{
									if (fsr_photonlep1 && !fsr_photonlep2)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
										VBF_LooseIso_OneFSR++;
									}
									if (fsr_photonlep2 && !fsr_photonlep1)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
										VBF_LooseIso_OneFSR++;
									}
									if (fsr_photonlep1 && fsr_photonlep2)
									{
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
										histos.fillHiggsHist(strVectorBoson + strVBFLoose + strRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
												FSR_PhotonToLep2, pu_weight);
										VBF_LooseIso_TwoFSR++;
									}
								}
							}
						}
					}
				}
			}
			//------------If the muons didn't pass any reliso--------------------
			if (pass_dimuon && !pass_dimuon_RelIso)
			{
				if (Dimuon.M() > 115.0 && Dimuon.M() < 135.0)
				{
					if (Gluon_Fusion || data)
					{
						if (GF_cat == 1)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_BBTightNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_BBTightNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_BBTightNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_BBTightNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 2)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_BEBOTightNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_BEBOTightNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_BEBOTightNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_BEBOTightNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 3)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_EEOOTightNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_EEOOTightNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_EEOOTightNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_EEOOTightNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 4)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_BBLooseNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_BBLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_BBLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBBLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_BBLooseNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 5)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_BEBOLooseNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_BEBOLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_BEBOLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFBEBOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_BEBOLooseNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 6)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_EEOOLooseNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_EEOOLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_EEOOLooseNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFEEOOLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_EEOOLooseNoIso_TwoFSR++;
								}
							}
						}
						if (GF_cat == 7)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								GF_TwoJetsNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									GF_TwoJetsNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									GF_TwoJetsNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strGluonFusion + strGFTwoJets + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									GF_TwoJetsNoIso_TwoFSR++;
								}
							}
						}
					}

					if (Vector_Boson || data)
					{
						if (VBF_cat == 1)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								VBF_TightNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									VBF_TightNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									VBF_TightNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFTight + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									VBF_TightNoIso_TwoFSR++;
								}
							}
						}
						if (VBF_cat == 2)
						{
							if (abs(125.0 - Dimuon.M()) <= abs(125.0 - (Dimuon + FSR_PhotonToLep1 + FSR_PhotonToLep2).M()))
							{
								histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, pu_weight);
								histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
								histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strNOFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
										FSR_PhotonToLep2, pu_weight);
								VBF_LooseNoIso_NOFSR++;
							}
							else
							{
								if (fsr_photonlep1 && !fsr_photonlep2)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, FSR_PhotonToLep1, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1, pu_weight);
									VBF_LooseNoIso_OneFSR++;
								}
								if (fsr_photonlep2 && !fsr_photonlep1)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strAllEvents, Dimuon, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep2, pu_weight);
									VBF_LooseNoIso_OneFSR++;
								}
								if (fsr_photonlep1 && fsr_photonlep2)
								{
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + tagAllEvents, Dimuon, FSR_PhotonToLep1, FSR_PhotonToLep2, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
									histos.fillHiggsHist(strVectorBoson + strVBFLoose + strNonRelIso + strFSREvents + tagDimuonPlusPH, Dimuon, FSR_PhotonToLep1,
											FSR_PhotonToLep2, pu_weight);
									VBF_LooseNoIso_TwoFSR++;
								}
							}
						}
					}
				}
			}
		}
	}

	cout << "==========================================================" << endl;
	cout << "Events that pass the trigger:     " << trigger_evt << endl;
	cout << "Events that pass the dimuon preselection:     " << dimuon_evt << endl;
	cout << "Events that pass the dimuon selection:     " << dimuon_initialsel_evt << endl;
	cout << "Events that pass the isolation (reliso or relisoFSR):  " << dimuon_passRelIso << endl;
	cout << "==========================================================" << endl;
	cout << "················ Relative to Isolation ···················" << endl;
	cout << "Events that both muons pass the RelIso :  " << RelIsoReal << endl;
	cout << "Events in which at least one muon do not pass the RelIso:  " << OneNoreliso << endl;
	cout << "Events in which both muons do not pass the RelIso:  " << BothNoreliso << endl;
	cout << "Events that no pass RelIso, but pass PassRelIso-gamma:  " << PassRelisoNoFSR << endl;
	cout << "==========================================================" << endl;

	if (Gluon_Fusion || data)
	{
		cout << "····················Gluon Fusion Categories···············" << endl;
		cout << "------------------------GF BB Tight-----------------------" << endl;
		cout << "Events in category GF BB Tight:         " << GF_BB_Tight_cat << endl;
		cout << "Events that pass isolation:         " << GF_BBTight_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_BBTightIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_BBTightIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_BBTightIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_BB_Tight_cat - GF_BBTight_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_BBTightNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_BBTightNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_BBTightNoIso_TwoFSR << endl;
		cout << "------------------------GF BEBO Tight-----------------------" << endl;
		cout << "Events in category GF BB Tight:         " << GF_BEBO_Tight_cat << endl;
		cout << "Events that pass isolation:         " << GF_BEBOTight_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_BEBOTightIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_BEBOTightIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_BEBOTightIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_BEBO_Tight_cat - GF_BEBOTight_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_BEBOTightNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_BEBOTightNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_BEBOTightNoIso_TwoFSR << endl;
		cout << "------------------------GF EEOO Tight-----------------------" << endl;
		cout << "Events in category GF BB Tight:         " << GF_EEOO_Tight_cat << endl;
		cout << "Events that pass isolation:         " << GF_EEOOTight_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_EEOOTightIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_EEOOTightIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_EEOOTightIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_EEOO_Tight_cat - GF_EEOOTight_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_EEOOTightNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_EEOOTightNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_EEOOTightNoIso_TwoFSR << endl;
		cout << "------------------------GF BB Loose-----------------------" << endl;
		cout << "Events in category GF BB Loose:         " << GF_BB_Loose_cat << endl;
		cout << "Events that pass isolation:         " << GF_BBLoose_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_BBLooseIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_BBLooseIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_BBLooseIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_BB_Loose_cat - GF_BBLoose_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_BBLooseNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_BBLooseNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_BBLooseNoIso_TwoFSR << endl;
		cout << "------------------------GF BEBO Loose-----------------------" << endl;
		cout << "Events in category GF BB Loose:         " << GF_BEBO_Loose_cat << endl;
		cout << "Events that pass isolation:         " << GF_BEBOLoose_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_BEBOLooseIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_BEBOLooseIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_BEBOLooseIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_BEBO_Loose_cat - GF_BEBOLoose_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_BEBOLooseNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_BEBOLooseNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_BEBOLooseNoIso_TwoFSR << endl;
		cout << "------------------------GF EEOO Loose-----------------------" << endl;
		cout << "Events in category GF BB Loose:         " << GF_EEOO_Loose_cat << endl;
		cout << "Events that pass isolation:         " << GF_EEOOLoose_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_EEOOLooseIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_EEOOLooseIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_EEOOLooseIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_EEOO_Loose_cat - GF_EEOOLoose_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_EEOOLooseNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_EEOOLooseNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_EEOOLooseNoIso_TwoFSR << endl;
		cout << "------------------------GF Two Jets-----------------------" << endl;
		cout << "Events in category GF BB Loose:         " << GF_TwoJets_cat << endl;
		cout << "Events that pass isolation:         " << GF_TwoJets_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << GF_TwoJetsIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << GF_TwoJetsIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << GF_TwoJetsIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << GF_TwoJets_cat - GF_TwoJets_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << GF_TwoJetsNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << GF_TwoJetsNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << GF_TwoJetsNoIso_TwoFSR << endl;
		cout << "==========================================================" << endl;
	}

	if (Vector_Boson || data)
	{
		cout << "················ VBF Categories ························" << endl;
		cout << "----------------VBF Tight category----------------------" << endl;
		cout << "Events in category VBF Tight:       " << VBF_Tight_cat << endl;
		cout << "Events that pass isolation:         " << VBF_Tight_Iso << endl;
		cout << "Events pass iso, closer to H or without FSR Ph:    " << VBF_TightIso_NOFSR << endl;
		cout << "Events pass iso with one FSR Ph:        " << VBF_TightIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << VBF_TightIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << VBF_Tight_cat - VBF_Tight_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << VBF_TightNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << VBF_TightNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << VBF_TightNoIso_TwoFSR << endl;

		cout << "----------------VBF Loose category----------------------" << endl;
		cout << "Events in category VBF Loose:       " << VBF_Loose_cat << endl;
		cout << "Events that pass isolation:         " << VBF_Loose_Iso << endl;
		cout << "Events with iso, closer to H or without FSR Ph:    " << VBF_LooseIso_NOFSR << endl;
		cout << "Events with iso with one FSR Ph:        " << VBF_LooseIso_OneFSR << endl;
		cout << "Events pass iso with two FSR Ph:        " << VBF_LooseIso_TwoFSR << endl;
		cout << "Events that not pass isolation:         " << VBF_Loose_cat - VBF_Loose_Iso << endl;
		cout << "Events not pass iso, closer to H or without FSR Ph:   " << VBF_LooseNoIso_NOFSR << endl;
		cout << "Events not pass iso with one FSR Ph:        " << VBF_LooseNoIso_OneFSR << endl;
		cout << "Events not pass iso with two FSR Ph:        " << VBF_LooseNoIso_TwoFSR << endl;
		cout << "==========================================================" << endl;
	}

	//-------------------------Write the histograms----------------------------
	histos.writeTFile();
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

void BSM_Analysis::MuonsVectors(TLorentzVector& Reco_lepton1, TLorentzVector& Reco_lepton2, double& RelIso1, double& RelIso2, double& RelIso_NoFSR1,
		double& RelIso_NoFSR2)
//, double& ChargedOverall1, double& ChargedOverall2,double& NeutralOverall1, double& NeutralOverall2)
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
				if (Muon_charge->at(m1) * Muon_charge->at(m2) < 0)
				{
					first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
					Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
					float dimuon_mass = (first_muon_vec + Subfirst_muon_vec).M();

					if (dimuon_mass < dimuon_mass_int)
					{
						first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
						Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
						Reco_lepton1 = first_muon_vec;
						Reco_lepton2 = Subfirst_muon_vec;
						dimuon_mass_int = dimuon_mass;
						//RelIso1 = Muon_combinedIso->at(m1);				//Muon_isoSum->at(m1) / Muon_pt->at(m1);
						//RelIso2 = Muon_combinedIso->at(m2);				//Muon_isoSum->at(m2) / Muon_pt->at(m2);
						RelIso1 = Muon_trackRe_iso->at(m1);			//track iso
						RelIso2 = Muon_trackRe_iso->at(m2);	    //track iso

						RelIso_NoFSR1 = (Muon_isoCharged->at(m1) + max(0., Muon_isoNeutralHadron->at(m1) - 0.5 * (Muon_isoPU->at(m1)))) / Muon_pt->at(m1);
						RelIso_NoFSR2 = (Muon_isoCharged->at(m2) + max(0., Muon_isoNeutralHadron->at(m2) - 0.5 * (Muon_isoPU->at(m2)))) / Muon_pt->at(m2);

						/*ChargedOverall1 = Muon_isoCharged->at(m1) / Muon_pt->at(m1);
						 ChargedOverall2 = Muon_isoCharged->at(m2) / Muon_pt->at(m2);
						 NeutralOverall1 = Muon_isoNeutralHadron->at(m1) / Muon_pt->at(m1);
						 NeutralOverall2 = Muon_isoNeutralHadron->at(m2) / Muon_pt->at(m2);*/

					}
				}
			}
		}
	}
}

bool BSM_Analysis::TrackIsoGamma(double& RelIso1, double& RelIso2, TLorentzVector& Muon1, TLorentzVector& Muon2, TLorentzVector& PFPhoton1, TLorentzVector& PFPhoton2)
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
	if (patElectron_pt->size() == 0) noelectron = true;

	TLorentzVector Jet(0., 0., 0., 0);
	for (int j1 = 0; j1 < Jet_pt->size(); j1++)
	{
		Jet.SetPtEtaPhiE(Jet_pt->at(j1), Jet_eta->at(j1), Jet_phi->at(j1), Jet_energy->at(j1));
		if (Jet.DeltaR(Muon1) >= 0.3 && Jet.DeltaR(Muon2) >= 0.3)
			nojet = true;
	}
	if (Jet_pt->size()==0) nojet = true;

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

int BSM_Analysis::JetsVector(TLorentzVector& Muon1, TLorentzVector& Muon2)
{
	TLorentzVector Jet_temp(0., 0., 0., 0.);
	TLorentzVector Jet1(0., 0., 0., 0.);
	TLorentzVector Jet2(0., 0., 0., 0.);
	vector<TLorentzVector> Jets_passpresel_vec;
	Jets_passpresel_vec.erase(Jets_passpresel_vec.begin(), Jets_passpresel_vec.end());
	bool VBF_Tight = false;
	bool VBF_Loose = false;
	float dijet_mass = 0.0;

	if (Jet_pt->size() > 1 && Met_type1PF_pt < 40.0)
	{
		for (int m1 = 0; m1 < Jet_pt->size(); m1++)
		{
			if (Jet_pt->at(m1) > 30.0 && Jet_eta->at(m1) < 2.4 && Jet_bDiscriminator_cMVAv2->at(m1) > 0.8484)
				continue;
			if (Jet_pt->at(m1) > 30.0 && Jet_eta->at(m1) < 4.7)
			{
				Jet_temp.SetPtEtaPhiE(Jet_pt->at(m1), Jet_eta->at(m1), Jet_phi->at(m1), Jet_energy->at(m1));
				Jets_passpresel_vec.push_back(Jet_temp);
			}

		}
	}

	if (Jets_passpresel_vec.size() > 1)
	{
		for (int j1 = 0; j1 < Jets_passpresel_vec.size(); j1++)
		{
			for (int j2 = j1 + 1; j2 < Jets_passpresel_vec.size(); j2++)
			{
				if (j2 != j1 && (Jets_passpresel_vec[j1].Eta() * Jets_passpresel_vec[j2].Eta() < 0))
				{
					if (Jets_passpresel_vec[j1].Pt() > 40.0 || Jets_passpresel_vec[j2].Pt() > 40.0)
					{
						Jet1 = Jets_passpresel_vec[j1];
						Jet2 = Jets_passpresel_vec[j2];

						dijet_mass = (Jet1 + Jet2).M();

						if (Jet1.DeltaR(Muon1) > 0.4 && Jet1.DeltaR(Muon2) > 0.4 && Jet2.DeltaR(Muon1) > 0.4 && Jet2.DeltaR(Muon2) > 0.4)
						{
							if (dijet_mass > 650.0 && abs(Jet1.Eta() - Jet2.Eta()) > 2.5)
							{
								VBF_Tight = true;
							}
							else if (dijet_mass > 250.0 && abs(Jet1.Eta() - Jet2.Eta()) > 2.5)
							{
								VBF_Loose = true;
							}
						}
					}
				}
			}
		}
	}
	if (VBF_Tight == true)
	{
		return 1;
	}
	if (VBF_Loose == true)
	{
		return 2;
	}
	return 0;
}

int BSM_Analysis::GFCategories(TLorentzVector& Muon1, TLorentzVector& Muon2)
{
	TLorentzVector Dimuon = Muon1 + Muon2;
	double MuEta1 = abs(Muon1.Eta());
	double MuEta2 = abs(Muon2.Eta());
	bool GF_BB_Tight = false, GF_BB_Loose = false, GF_BEBO_Tight = false, GF_BEBO_Loose = false, GF_EEOO_Tight = false, GF_EEOO_Loose = false;
	bool GF_TwoJets = false;

	if (Jet_pt->size() < 2 || Jet_puppi_pt->size() < 2)
	{
		if (Dimuon.Pt() < 25.0)
		{
			if (MuEta1 < 0.8 && MuEta2 < 0.8)
			{
				GF_BB_Tight = true;
			}
			if ((MuEta1 < 0.8 && MuEta2 > 0.8 && MuEta2 < 1.6) || (MuEta2 < 0.8 && MuEta1 > 0.8 && MuEta1 < 1.6) || (MuEta1 < 0.8 && MuEta2 > 1.6 && MuEta2 < 2.4)
					|| (MuEta2 < 0.8 && MuEta1 > 1.6 && MuEta1 < 2.4))
			{
				GF_BEBO_Tight = true;
			}
			if (MuEta1 > 1.6 && MuEta1 < 2.4 && MuEta2 > 1.6 && MuEta2 < 2.4)
			{
				GF_EEOO_Tight = true;
			}
		}
		else
		{
			if (MuEta1 < 0.8 && MuEta2 < 0.8)
			{
				GF_BB_Loose = true;
			}
			if ((MuEta1 < 0.8 && MuEta2 > 0.8 && MuEta2 < 1.6) || (MuEta2 < 0.8 && MuEta1 > 0.8 && MuEta1 < 1.6) || (MuEta1 < 0.8 && MuEta2 > 1.6 && MuEta2 < 2.4)
					|| (MuEta2 < 0.8 && MuEta1 > 1.6 && MuEta1 < 2.4))
			{
				GF_BEBO_Loose = true;
			}
			if (MuEta1 > 1.6 && MuEta1 < 2.4 && MuEta2 > 1.6 && MuEta2 < 2.4)
			{
				GF_EEOO_Loose = true;
			}
		}
	}
	if (Jet_pt->size() >= 2 || Jet_puppi_pt->size() >= 2)
	{
		GF_TwoJets = true;
	}

	if (GF_BB_Tight == true)
		return 1;
	if (GF_BEBO_Tight == true)
		return 2;
	if (GF_EEOO_Tight == true)
		return 3;
	if (GF_BB_Loose == true)
		return 4;
	if (GF_BEBO_Loose == true)
		return 5;
	if (GF_EEOO_Loose == true)
		return 6;
	if (GF_TwoJets == true)
		return 7;
	return 0;
}

//*********************************************************************************************
// setBranchAdress: initialization of each object pointer. Also Set branch addresses and branch 
//                  pointers.
//*********************************************************************************************

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
	Jet_bDiscriminator = 0;
	Jet_mass = 0;
	Jet_neutralHadEnergyFraction = 0;
	Jet_neutralEmEmEnergyFraction = 0;
	Jet_chargedHadronEnergyFraction = 0;
	Jet_chargedEmEnergyFraction = 0;
	Jet_muonEnergyFraction = 0;
	Jet_electronEnergy = 0;
	Jet_photonEnergy = 0;
	Jet_bDiscriminator_cMVAv2 = 0;
	Jet_puppi_bDiscriminator_cMVAv2 = 0;
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

	// Set branch addresses and branch pointers
	if (!BOOM)
		return;
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
	BOOM->SetBranchAddress("Muon_combinedIso", &Muon_combinedIso, &b_Muon_combinedIso);
	BOOM->SetBranchAddress("Muon_trackRe_iso", &Muon_trackRe_iso, &b_Muon_trackRe_iso);

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
	BOOM->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
	BOOM->SetBranchAddress("Jet_bDiscriminator_cMVAv2", &Jet_bDiscriminator_cMVAv2, &b_Jet_bDiscriminator_cMVAv2);
	BOOM->SetBranchAddress("Jet_puppi_bDiscriminator_cMVAv2", &Jet_puppi_bDiscriminator_cMVAv2, &b_Jet_puppi_bDiscriminator_cMVAv2);
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
}
;
