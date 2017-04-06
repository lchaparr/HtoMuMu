#include "Jet.h"
#include "JetManager.h"

Jet::Jet(JetManager* manager, int index)
    {
    this->manager = manager;
    this->index = index;
    }

int Jet::getIndex() const
    {
    return index;
    }

double Jet::getJetChargedEmEnergyFraction() const
    {
    return manager->Jet_chargedEmEnergyFraction->at(index);
    }

double Jet::getJetChargedHadronEnergyFraction() const
    {
    return manager->Jet_chargedHadronEnergyFraction->at(index);
    }

double Jet::getJetElectronEnergy() const
    {
    return manager->Jet_electronEnergy->at(index);
    }

double Jet::getJetEnergy() const
    {
    return manager->Jet_energy->at(index);
    }

double Jet::eta() const
    {
    return manager->Jet_eta->at(index);
    }

double Jet::getJetMass() const
    {
    return manager->Jet_mass->at(index);
    }

double Jet::getJetMuonEnergyFraction() const
    {
    return manager->Jet_muonEnergyFraction->at(index);
    }

double Jet::getJetNeutralEmEmEnergyFraction() const
    {
    return manager->Jet_neutralEmEmEnergyFraction->at(index);
    }

double Jet::getJetNeutralHadEnergyFraction() const
    {
    return manager->Jet_neutralHadEnergyFraction->at(index);
    }

double Jet::phi() const
    {
    return manager->Jet_phi->at(index);
    }

double Jet::getJetPhotonEnergy() const
    {
    return manager->Jet_photonEnergy->at(index);
    }

double Jet::pt() const
    {
    return manager->Jet_pt->at(index);
    }

double Jet::getUncorrJetPt() const
    {
    return manager->UncorrJet_pt->at(index);
    }

void Jet::setUncorrJetPt(double uncorrJetPt)
    {
    UncorrJet_pt = uncorrJetPt;
    }

TLorentzVector Jet::getTLorentz(){
    TLorentzVector returnTLorentz;
    returnTLorentz.SetPtEtaPhiE(this->pt(), this->eta(), this->phi(), this->getJetEnergy());
    return returnTLorentz;
}
std::ostream &operator<<(std::ostream &output, Jet &Jet)
    {
    output << "Jet particle, pt=" << Jet.pt() << ", eta=" << Jet.eta() << ", phi=" << Jet.phi();
    return output;
    }
Jet::~Jet()
    {
    // TODO Auto-generated destructor stub
    }

