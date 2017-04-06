#include "MuonManager.h"


MuonManager::MuonManager(TTree* BOOM, bool IsRealData, cli::Parser& program_options)
    {
    this->BOOM = BOOM;
    this->IsRealData = IsRealData;
    this->setBranchAddress();
    this->rc = new RoccoR();
    this->rc->init("rcdata.2016.v3");
    this->maxRelIso = program_options.getAlternative<double>("maxRelIso");;
    this->minMuonPt = program_options.getAlternative<double>("minMuonPt");
    this->noGenMatching = program_options.getAlternative<bool>("noGenMatching");
    this->noRochesterCorrections = program_options.getAlternative<bool>("noRochesterCorrections");
    }

MuonManager::~MuonManager()
    {
    delete rc;
    }

int MuonManager::size()
    {
    return Muon_pt->size();
    }

void MuonManager::setBranchAddress()
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
    Gen_pt = 0;
    Gen_eta = 0;
    Gen_phi = 0;
    Gen_energy = 0;
    Gen_pdg_id = 0;
    Gen_motherpdg_id = 0;
    Gen_status = 0;
    Gen_BmotherIndex = 0;
    Gen_Met = 0;
    Gen_charge = 0;

    // Set branch addresses and branch pointers
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
    }

std::vector<Muon> MuonManager::selectMuons(bool useTrackIso, bool useFSRIso)
    {
    std::vector<Muon> returnVector;
    std::vector<Muon> returnVector2;

    if (this->size() < 2)
	{
	throw LESS_THAN_2_MUONS;
	}

    for (int i = 0; i < this->size(); i++)
	{
	Muon muon = this->makeMuon(i);
	double muonrelIsoNoFSR =  (Muon_isoCharged->at(i) + std::max(0., Muon_isoNeutralHadron->at(i) - 0.5 * (Muon_isoPU->at(i)))) / Muon_pt->at(i);
	double muonPFIso = Muon_combinedIso->at(i);
	double muonTrackIso = Muon_trackRe_iso->at(i);

	double muonRelIso;
	if (useTrackIso == true){ muonRelIso = muonTrackIso;}
	if (useTrackIso == false){ muonRelIso = muonPFIso;}
	if (useFSRIso == true) { muonRelIso = muonrelIsoNoFSR;}

	if (muon.pt() > this->minMuonPt
		&& muon.isMuonIsGlobal() == true
		&& muon.getMuonIsTrackerMuon() == true
		&& abs(muon.eta()) < 2.4
		&& muon.isMuonMedium() == true
		&& muonRelIso < this->maxRelIso)
		{
		returnVector.push_back(muon);
		}

	}
    if (returnVector.size() == 2)
	{
	std::sort(returnVector.begin(), returnVector.end(), aux::greaterPt<Muon>);
	double dimuon_pt_inicial = 0.0;
	double dimuon_pt;
	for (auto muon1 : returnVector){
	    for (auto muon2: returnVector){
		if (!aux::isSameParticle<Muon>(muon1, muon2)){
		    dimuon_pt = muon1.pt() + muon2.pt();
		    if (muon1.getMuonCharge() * muon2.getMuonCharge() < 0
		    		&& ((muon1.pt()> 26.0 &&  muon1.isMuonIsTriggerMatched()   )
		    		|| (muon2.pt()> 26.0 &&  muon2.isMuonIsTriggerMatched() ))
				&& dimuon_pt > dimuon_pt_inicial){
			returnVector2.clear();
			returnVector2.push_back(muon1);
			returnVector2.push_back(muon2);
			dimuon_pt_inicial = dimuon_pt;
		    }
		}
	    }

	}
	std::sort(returnVector2.begin(), returnVector2.end(), aux::greaterPt<Muon>);
	if (returnVector2.size() == 2){
	return returnVector2;}
	else {
	    throw LESS_THAN_2_MUONS;
	}
	}
	else {
	    throw NO_MATCHING_CHARGE;
	}

    }
std::vector<Muon> MuonManager::selectMuonsPFIso(){
    return this->selectMuons(false, false);
}

std::vector<Muon> MuonManager::selectMuonsPFIsoMinusGamma(){
    return this->selectMuons(false, true);
}

std::vector<Muon> MuonManager::selectMuonsTrackIso(){
    return this->selectMuons(true, false);
}

double MuonManager::getGenPt(int muonIndex){
    int selectedGen = 66666666;
    double minDeltaR = 0.005;
    if (this->noGenMatching) return -1.0;
    for (int genIndex = 0; genIndex < Gen_pt->size(); genIndex++)
	{
	//cout <<"Decay||"<< Gen_pdg_id->at(genIndex) << "||"<< Gen_status->at(genIndex)<< endl;
	if (abs(Gen_pdg_id->at(genIndex)) == 13 && Gen_charge->at(genIndex) == Muon_charge->at(muonIndex) && (Gen_status->at(genIndex) == 1  /*|| Gen_status->at(muonIndex) == 2*/))
	    {
	    double currentDeltaR = deltaR(Muon_eta->at(muonIndex), Muon_phi->at(muonIndex), Gen_eta->at(genIndex), Gen_phi->at(genIndex));
	    //cout << "DeltaR " << currentDeltaR << "\n";
	    if (currentDeltaR < minDeltaR)
		{
		selectedGen = genIndex;
		minDeltaR = currentDeltaR;
		}
	    }
	}
    if (selectedGen < 66666666)
	{
	return Gen_pt->at(selectedGen);
	}
    else
	{
	return -1.0;
	}
    }

Muon MuonManager::makeMuon(int index){
    Muon newMuon(this, index);
    return newMuon;
}

std::ostream &operator<<( std::ostream  &output, MuonManager &manager ){
    for (int i = 0; i < manager.size(); i++){
	Muon tempMuon = manager.makeMuon(i);
	output << tempMuon << std::endl;
    }
	return output;

}

;
