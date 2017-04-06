#ifndef EVENTMANAGER_H_
#define EVENTMANAGER_H_
#include <TBranch.h>
#include <vector>
#include <TTree.h>
#include "Muon.h"
#include "Jet.h"
#include <TLorentzVector.h>

using namespace std;

class EventManager
    {
public:
    EventManager(TTree* BOOM);
    TTree* BOOM;
    std::map<bool, std::string> relIsoTag;
    std::map<int, std::string> fsrTag;
    std::vector<Muon> muons;
    std::vector<Jet> jets;
    TLorentzVector Dimuon;
    void setBranchAddress();
    virtual ~EventManager();
    bool passEventSelection;
    bool passDimuonSelection;
    bool isTriggerMatched();
    bool getIsRealData();
    int getRunNumber();
    int getEventNumber();
    int getLumiBlock();

    std::vector<std::string> *Trigger_names;
    std::vector<int> *Trigger_decision;
    Float_t ntruePUInteractions;

    bool isRealData;
    int runNumber;
    int eventNumber;
    int lumiBlock;

    TBranch *b_Trigger_names;
    TBranch *b_Trigger_decision;
    TBranch *b_isRealData;
    TBranch *b_ntruePUInteractions;
    TBranch *b_runNumber;
    TBranch *b_eventNumber;
    TBranch *b_lumiBlock;

    };

#endif /* EVENTMANAGER_H_ */
