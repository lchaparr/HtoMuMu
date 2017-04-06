#include "JetManager.h"


JetManager::JetManager(TTree* BOOM, cli::Parser& program_options)
    {
    this->BOOM = BOOM;
    this->setBranchAddress();
    //this->program_options = NULL;
    this->minJetPt = program_options.getAlternative<double>("minJetPt");
    }

JetManager::~JetManager()
    {
    }

int JetManager::size()
    {
    return Jet_pt->size();
    }

void JetManager::setBranchAddress()
    {
    // Set object pointer
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
    UncorrJet_pt = 0;
    Jet_bDiscriminator_pfCISVV2 = 0;

    // Set branch addresses and branch pointers
    BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
    BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
    BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
    BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
    BOOM->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
    BOOM->SetBranchAddress("Jet_neutralHadEnergyFraction", &Jet_neutralHadEnergyFraction, &b_Jet_neutralHadEnergyFraction);
    BOOM->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &Jet_neutralEmEmEnergyFraction, &b_Jet_neutralEmEmEnergyFraction);
    BOOM->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
    BOOM->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
    BOOM->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
    BOOM->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
    BOOM->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
    BOOM->SetBranchAddress("Jet_bDiscriminator_pfCISVV2", &Jet_bDiscriminator_pfCISVV2, &b_Jet_bDiscriminator_pfCISVV2);
    //BOOM->SetBranchAddress("Jet_bDiscriminator_pfCMVAV2", &Jet_bDiscriminator_pfCMVAV2, &b_Jet_bDiscriminator_pfCMVAV2);
    //BOOM->SetBranchAddress("Jet_puppi_bDiscriminator_pfCMVAV2", &Jet_bDiscriminator_pfCMVAV2, &b_Jet_bDiscriminator_pfCMVAV2);
    }

/*Preselection:
 at least 2 jets passing jet selection
 leading jet pT > 40 subleading jet pT > 30 GeV
 Jet selection:
 dR(jet,mu) > 0.4
 pT > 30 GeV
 |eta| < 4.7*/
std::vector<Jet> JetManager::selectJets(std::vector<Muon>& muons)
    {
    std::vector<Jet> returnVector;
    for (int jetIndex = 0; jetIndex < this->size(); jetIndex++)
	{
	Jet tempJet = this->makeJet(jetIndex);
	//cout << tempJet.pt()  << ", abs eta" << abs(tempJet.eta());
	if (tempJet.pt() > 30.0 && abs(tempJet.eta()) < 4.7)
	    {
	    Muon tempMuon1 = muons.at(0);
	    Muon tempMuon2 = muons.at(1);
	    double deltaR1 = deltaR(tempJet.eta(), tempJet.phi(), tempMuon1.eta(), tempMuon1.phi());
	    double deltaR2 = deltaR(tempJet.eta(), tempJet.phi(), tempMuon2.eta(), tempMuon2.phi());
	    //cout  << tempJet << "deltaR1 " <<deltaR1 << "deltaR2 " <<deltaR2 << endl;

	    if (deltaR1 > 0.4 && deltaR2 > 0.4)
		{
		returnVector.push_back(tempJet);
		}
	    }
	}

    std::sort(returnVector.begin(), returnVector.end(), aux::greaterPt<Jet>);
    if (returnVector.size() > 1)
	{
	for (auto jet : returnVector)
	    {
	    if (jet.pt() > 40)
		{
		return returnVector;
		}
	    }
	throw NO_JET_ABOVE_40_PT;
	}
    else
	{
	throw TOO_FEW_JETS;
	}
    }

bool JetManager::bJetVeto()
    {
    //cout << "Entre en el BJetVeto, size =" << this->size() << endl;
    for (int jetIndex = 0; jetIndex < this->size(); jetIndex++)
	{
	//cout << "entre en el for" << endl;
	//cout  << this->size() <<"\t"<<  this->Jet_pt->at(jetIndex) <<"\t"<<abs(this->Jet_eta->at(jetIndex)) << "\t" << this->Jet_bDiscriminator_pfCMVAV2->at(jetIndex) << endl;
	//no jet with pT > 30 GeV, |eta| < 2.4 and CSVv2 > 0.8484 (Medium wp)
	if (this->Jet_pt->at(jetIndex) > 30 && abs(this->Jet_eta->at(jetIndex)) < 2.4 && this->Jet_bDiscriminator_pfCISVV2->at(jetIndex) > 0.8484)
	    {
	    return true;
	    }
	}
    return false;
    }
//
Jet JetManager::makeJet(int index)
    {
    Jet newJet(this, index);
    return newJet;
    }


std::ostream &operator<<( std::ostream  &output, JetManager &manager ){
    for (int i = 0; i < manager.size(); i++){
	Jet tempJet = manager.makeJet(i);
	output << tempJet << std::endl;
    }
	return output;

}


;
