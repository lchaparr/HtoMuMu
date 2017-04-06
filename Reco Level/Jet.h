#ifndef Jet_H_
#define Jet_H_

#include <iostream>
#include "Muon.h"
#include <TLorentzVector.h>

class JetManager;

class Jet
    {
public:
    JetManager* manager;
    int index;
    double Jet_pt;
    Jet(JetManager* manager, int index);
    virtual ~Jet();
    double Jet_eta;
    double Jet_phi;
    double Jet_energy;

    double Jet_mass;
    double Jet_neutralHadEnergyFraction;
    double Jet_neutralEmEmEnergyFraction;
    double Jet_chargedHadronEnergyFraction;
    double Jet_chargedEmEnergyFraction;
    double Jet_muonEnergyFraction;
    double Jet_electronEnergy;
    double Jet_photonEnergy;
    double UncorrJet_pt;

    int getIndex() const;
    double getJetChargedEmEnergyFraction() const;
    double getJetChargedHadronEnergyFraction() const;
    double getJetElectronEnergy() const;
    double getJetEnergy() const;
    double eta() const;
    double getJetMass() const;
    double getJetMuonEnergyFraction() const;
    double getJetNeutralEmEmEnergyFraction() const;
    double getJetNeutralHadEnergyFraction() const;
    double phi() const;
    double getJetPhotonEnergy() const;
    double pt() const;
    double getUncorrJetPt() const;
    void setUncorrJetPt(double uncorrJetPt);
    TLorentzVector getTLorentz();
    friend std::ostream &operator<<(std::ostream &output, Jet &Jet);

    };

#endif /* Jet_H_ */
