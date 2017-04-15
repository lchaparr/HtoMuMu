#include "FSRManager.h"

FSRManager::FSRManager(TTree* BOOM)
    {
    this->BOOM = BOOM;
    this->setBranchAddress();
    }

FSRManager::~FSRManager()
    {
    }

int FSRManager::size()
    {
    return FSRPhoton_pt->size();
    }

void FSRManager::setBranchAddress()
    {
    // Set object pointer
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

    // Set branch addresses and branch pointers
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

    }

std::ostream &operator<<(std::ostream &output, FSRManager &manager)
    {
    for (int i = 0; i < manager.size(); i++)
	{
	output << "Not implemented" << std::endl;
	}
    return output;

    }

std::vector<FSRPhoton> FSRManager::selectFSRPhotons(std::vector<Muon>& muons)
    {
    std::vector<FSRPhoton> returnVector;
    std::vector<std::pair<double, FSRPhoton>> selectedPhotonsM1;
    std::vector<std::pair<double, FSRPhoton>> selectedPhotonsM2;
    if (muons.size() != 2)
	{
	return returnVector;
	}
    for (int i = 0; i < this->size(); i++)
	{
	FSRPhoton tempPhoton = this->makeFSRPhoton(i);
	double fsrIso = (tempPhoton.getFsrPhotonIsoChpu() + tempPhoton.getFsrPhotonIsoNhPhot()) / tempPhoton.pt();
	if (tempPhoton.pt() > 2.0 && abs(tempPhoton.eta()) < 2.4 && fsrIso < 1.8)
	    {
	    double deltaR1 = aux::calcDeltaR<Muon, FSRPhoton>(muons.at(0), tempPhoton);
	    double deltaR2 = aux::calcDeltaR<Muon, FSRPhoton>(muons.at(1), tempPhoton);
	    double dROverEt1 = deltaR1 / (tempPhoton.pt() * tempPhoton.pt());
	    double dROverEt2 = deltaR2 / (tempPhoton.pt() * tempPhoton.pt());
	    if (dROverEt1 < 0.012)
		selectedPhotonsM1.push_back(make_pair(dROverEt1, tempPhoton));
	    if (dROverEt2 < 0.012)
		selectedPhotonsM2.push_back(make_pair(dROverEt2, tempPhoton));
	    }
	}
	std::sort(selectedPhotonsM1.begin(), selectedPhotonsM1.end(), aux::lesserFirstElement<std::pair<double, FSRPhoton>>);
	std::sort(selectedPhotonsM2.begin(), selectedPhotonsM2.end(), aux::lesserFirstElement<std::pair<double, FSRPhoton>>);
	if (selectedPhotonsM1.size() > 0)
	    returnVector.push_back(selectedPhotonsM1.at(0).second);
	if (selectedPhotonsM2.size() > 0)
	    {
	    if (returnVector.size() == 0)
		{
		returnVector.push_back(selectedPhotonsM2.at(0).second);
		}
	    else
		{
		for (auto selectedPhoton : selectedPhotonsM2)
		    {
		    if (!aux::isSameParticle<FSRPhoton>(returnVector.at(0), selectedPhoton.second))
			{
			returnVector.push_back(selectedPhoton.second);
			break;
			}
		    }
		}
	    }
    return returnVector;
    }

FSRPhoton FSRManager::makeFSRPhoton(int index){
    FSRPhoton newFSRPhoton(this, index);
    return newFSRPhoton;
}

;
