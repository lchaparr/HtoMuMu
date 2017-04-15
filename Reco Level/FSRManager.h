#ifndef SRC_FSRManager_H_
#define SRC_FSRManager_H_
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TBranch.h>
#include <TTree.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "AuxiliaryFunctions.h"
#include "FSRPhoton.h"
#include "Muon.h"
#include <algorithm>

using namespace std;

class FSRManager
    {
public:
    FSRManager(TTree* BOOM);
    friend std::ostream &operator<<(std::ostream &output, FSRManager &manager);
    void setBranchAddress();
    virtual ~FSRManager();
    int size();
    std::vector<FSRPhoton> selectFSRPhotons(std::vector<Muon>& muons);
    FSRPhoton makeFSRPhoton(int index);
    TTree* BOOM;
    vector<float> *FSRPhoton_pt;
    vector<float> *FSRPhoton_eta;
    vector<float> *FSRPhoton_phi;
    vector<float> *FSRPhoton_energy;
    vector<float> *FSRPhoton_pX;
    vector<float> *FSRPhoton_pY;
    vector<float> *FSRPhoton_pZ;
    vector<float> *FSRPhoton_isoNH;
    vector<float> *FSRPhoton_isoCH;
    vector<float> *FSRPhoton_isoCHPU;
    vector<float> *FSRPhoton_isoPhot;
    vector<float> *FSRPhoton_et;
    vector<float> *FSRPhoton_isoNHPhot;
    // List of branches
    TBranch *b_FSRPhoton_pt;
    TBranch *b_FSRPhoton_eta;
    TBranch *b_FSRPhoton_phi;
    TBranch *b_FSRPhoton_energy;
    TBranch *b_FSRPhoton_pX;
    TBranch *b_FSRPhoton_pY;
    TBranch *b_FSRPhoton_pZ;
    TBranch *b_FSRPhoton_isoNH;
    TBranch *b_FSRPhoton_isoCH;
    TBranch *b_FSRPhoton_isoCHPU;
    TBranch *b_FSRPhoton_isoPhot;
    TBranch *b_FSRPhoton_et;
    TBranch *b_FSRPhoton_isoNHPhot;
    };

#endif /* SRC_FSRManager_H_ */
