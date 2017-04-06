#include "Muon.h"
#include "MuonManager.h"

Muon::Muon(MuonManager* manager, int index){
    this->manager = manager;
    this->index = index;
    this->rochCorr = 1.0;
    if (!(this->manager->noRochesterCorrections)){
    this->setRochCorr();}
}
double Muon::phi() const
    {
    return manager->Muon_phi->at(index);
    }
double Muon::eta() const
    {
    return manager->Muon_eta->at(index);
    }

double Muon::pt() const
    {
    return manager->Muon_pt->at(index)*this->rochCorr;
    }

double Muon::getMuonCharge() const
    {
    return manager->Muon_charge->at(index);
    }

double Muon::getMuonChi2() const
    {
    return manager->Muon_chi2->at(index);
    }

double Muon::getMuonCombinedIso() const
    {
    return manager->Muon_combinedIso->at(index);
    }

double Muon::getMuonDxyPv() const
    {
    return manager->Muon_dxy_pv->at(index);
    }

double Muon::getMuonDzBs() const
    {
    return manager->Muon_dz_bs->at(index);
    }

double Muon::getMuonEnergy() const
    {
    return manager->Muon_energy->at(index);
    }

double Muon::getMuonEta() const
    {
    return manager->Muon_eta->at(index);
    }

bool Muon::isMuonIsGlobal() const
    {
    return manager->Muon_isGlobal->at(index);
    }

double Muon::getMuonIsoCharged() const
    {
    return manager->Muon_isoCharged->at(index);
    }

double Muon::getMuonIsoCharParPt() const
    {
    return manager->Muon_isoCharParPt->at(index);
    }

double Muon::getMuonIsoNeutralHadron() const
    {
    return manager->Muon_isoNeutralHadron->at(index);
    }

double Muon::getMuonIsoPhoton() const
    {
    return manager->Muon_isoPhoton->at(index);
    }

double Muon::getMuonIsoPu() const
    {
    return manager->Muon_isoPU->at(index);
    }

double Muon::getMuonIsoSum() const
    {
    return manager->Muon_isoSum->at(index);
    }

double Muon::getMuonIsTrackerMuon() const
    {
    return manager->Muon_isTrackerMuon->at(index);
    }

bool Muon::isMuonIsTriggerMatched() const
    {
    return manager->Muon_isTriggerMatched->at(index);
    }

bool Muon::isMuonLoose() const
    {
    return manager->Muon_loose->at(index);
    }

double Muon::getMuonMatchedStat() const
    {
    return manager->Muon_matchedStat->at(index);
    }

bool Muon::isMuonMedium() const
    {
    return manager->Muon_medium->at(index);
    }

bool Muon::isMuonPf() const
    {
    return manager->Muon_pf->at(index);
    }

double Muon::getMuonPhi() const
    {
    return manager->Muon_phi->at(index);
    }

double Muon::getMuonPt() const
    {
    return manager->Muon_pt->at(index);
    }

bool Muon::isMuonSoft() const
    {
    return manager->Muon_soft->at(index);
    }

bool Muon::isMuonTight() const
    {
    return manager->Muon_tight->at(index);
    }

double Muon::getMuonTLayers() const
    {
    return manager->Muon_TLayers->at(index);
    }

double Muon::getMuonTrackReIso() const
    {
    return manager->Muon_trackRe_iso->at(index);
    }

double Muon::getMuonValidHits() const
    {
    return manager->Muon_validHits->at(index);
    }

double Muon::getMuonValidHitsInner() const
    {
    return manager->Muon_validHitsInner->at(index);
    }

bool Muon::getIsRealData() const
    {
    return manager->IsRealData;
    }
TLorentzVector Muon::getTLorentz(){
    TLorentzVector returnTLorentz;
    returnTLorentz.SetPtEtaPhiE(this->pt(), this->eta(), this->phi(), this->getMuonEnergy());
    return returnTLorentz;
}

void Muon::setRochCorr(){
    float u1 = gRandom->Rndm();
    float u2 = gRandom->Rndm();

    if (getIsRealData()){
	this->rochCorr = manager->rc->kScaleDT((int) this->getMuonCharge(), this->pt(),this->eta(), this->phi(), 0, 0);
    }
    else {
	double genPt = this->manager->getGenPt(this->index);
	if (genPt > 0.0){
	    this->rochCorr =  manager->rc->kScaleFromGenMC((int) this->getMuonCharge(), this->pt(), this->eta(), this->phi(), this->manager->Muon_TLayers->at(this->index), genPt, u1, 0, 0);
	}
	else{
	    this->rochCorr = manager->rc->kScaleAndSmearMC((int) this->getMuonCharge(), this->pt(), this->eta(), this->phi(), this->manager->Muon_TLayers->at(this->index), u1, u2, 0, 0);
	}
    }

}
std::ostream &operator<<( std::ostream  &output, Muon &muon ){
    output << "Muon particle, pt=" << muon.pt() << ", eta=" << muon.eta() << ", phi=" << muon.phi() <<", charge= " << muon.getMuonCharge() << ", isTriggerMatched=" << muon.isMuonIsTriggerMatched();
    return output;
}
Muon::~Muon()
    {
    // TODO Auto-generated destructor stub
    }

