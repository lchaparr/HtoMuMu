/*
 * FSRPhoton.h
 *
 *  Created on: Apr 6, 2017
 *      Author: gastonlp
 */

#ifndef FSRPHOTON_H_
#define FSRPHOTON_H_
#include <TLorentzVector.h>


class FSRManager;

class FSRPhoton
    {
public:
    FSRPhoton(FSRManager* manager, int index);
    int index;
    FSRManager* manager;
    virtual ~FSRPhoton();
    TLorentzVector getTLorentz();
    float getFsrPhotonEnergy() const;
    float getFsrPhotonEt() const;
    float eta() const;
    float getFsrPhotonIsoCh() const;
    float getFsrPhotonIsoChpu() const;
    float getFsrPhotonIsoNh() const;
    float getFsrPhotonIsoNhPhot() const;
    float getFsrPhotonIsoPhot() const;
    float phi() const;
    float pt() const;
    float getFsrPhotonPX() const;
    float getFsrPhotonPY() const;
    float getFsrPhotonPZ() const;
    };

#endif /* FSRPHOTON_H_ */
