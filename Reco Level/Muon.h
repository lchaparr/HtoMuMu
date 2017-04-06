#ifndef MUON_H_
#define MUON_H_

#include <iostream>
#include "TRandom3.h"
#include <TLorentzVector.h>
#include "cmdparser.hpp"


class MuonManager;

class Muon
    {
public:
    MuonManager* manager;
    int index;
    double rochCorr;
    double phi() const;
    double eta() const;
    double pt() const;
    Muon(MuonManager* manager, int index);
    virtual ~Muon();
    double getMuonCharge() const;
    double getMuonChi2() const;
    double getMuonCombinedIso() const;
    double getMuonDxyPv() const;
    double getMuonDzBs() const;
    double getMuonEnergy() const;
    double getMuonEta() const;
    bool isMuonIsGlobal() const;
    double getMuonIsoCharged() const;
    double getMuonIsoCharParPt() const;
    double getMuonIsoNeutralHadron() const;
    double getMuonIsoPhoton() const;
    double getMuonIsoPu() const;
    double getMuonIsoSum() const;
    double getMuonIsTrackerMuon() const;
    bool isMuonIsTriggerMatched() const;
    bool isMuonLoose() const;
    double getMuonMatchedStat() const;
    bool isMuonMedium() const;
    bool isMuonPf() const;
    double getMuonPhi() const;
    double getMuonPt() const;
    bool isMuonSoft() const;
    bool isMuonTight() const;
    double getMuonTLayers() const;
    double getMuonTrackReIso() const;
    double getMuonValidHits() const;
    double getMuonValidHitsInner() const;
    bool getIsRealData() const;
    void setRochCorr();
    TLorentzVector getTLorentz();
    friend std::ostream &operator<<( std::ostream  &output, Muon &muon );
    };

#endif /* MUON_H_ */
