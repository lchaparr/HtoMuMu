/*
 * FSRPhoton.cpp
 *
 *  Created on: Apr 6, 2017
 *      Author: gastonlp
 */

#include "FSRPhoton.h"
#include "FSRManager.h"

FSRPhoton::FSRPhoton(FSRManager* manager, int index)
    {
    this->manager = manager;
    this->index = index;
    }

FSRPhoton::~FSRPhoton()
    {
    // TODO Auto-generated destructor stub
    }

float FSRPhoton::getFsrPhotonEnergy() const
    {
    return this->manager->FSRPhoton_energy->at(index);
    }

float FSRPhoton::getFsrPhotonEt() const
    {
    return this->manager->FSRPhoton_et->at(index);
    }

float FSRPhoton::eta() const
    {
    return this->manager->FSRPhoton_eta->at(index);
    }

float FSRPhoton::getFsrPhotonIsoCh() const
    {
    return this->manager->FSRPhoton_isoCH->at(index);
    }

float FSRPhoton::getFsrPhotonIsoChpu() const
    {
    return this->manager->FSRPhoton_isoCHPU->at(index);
    }

float FSRPhoton::getFsrPhotonIsoNh() const
    {
    return this->manager->FSRPhoton_isoNH->at(index);
    }

float FSRPhoton::getFsrPhotonIsoNhPhot() const
    {
    return this->manager->FSRPhoton_isoNHPhot->at(index);
    }

float FSRPhoton::getFsrPhotonIsoPhot() const
    {
    return this->manager->FSRPhoton_isoPhot->at(index);
    }

float FSRPhoton::phi() const
    {
    return this->manager->FSRPhoton_phi->at(index);
    }

float FSRPhoton::pt() const
    {
    return this->manager->FSRPhoton_pt->at(index);
    }

float FSRPhoton::getFsrPhotonPX() const
    {
    return this->manager->FSRPhoton_pX->at(index);
    }

float FSRPhoton::getFsrPhotonPY() const
    {
    return this->manager->FSRPhoton_pY->at(index);
    }

float FSRPhoton::getFsrPhotonPZ() const
    {
    return this->manager->FSRPhoton_pZ->at(index);
    }
TLorentzVector FSRPhoton::getTLorentz(){
    TLorentzVector returnTLorentz;
    returnTLorentz.SetPtEtaPhiE(this->pt(), this->eta(), this->phi(), this->getFsrPhotonEnergy());
    return returnTLorentz;
}
