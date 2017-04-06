#include "EventManager.h"

EventManager::EventManager(TTree* BOOM)
    {
    this->BOOM = BOOM;
    this->setBranchAddress();
    relIsoTag[true] = "RelIso";
    relIsoTag[false] ="NonRelIso";
    fsrTag[0] = "noFSR";
    fsrTag[1] = "oneGamma";
    fsrTag[2] = "TwoGamma";
    }

EventManager::~EventManager()
    {
    // TODO Auto-generated destructor stub
    }

void EventManager::setBranchAddress(){
    Trigger_names = 0;
    Trigger_decision = 0;
    ntruePUInteractions = 0;
    isRealData = 0;
    runNumber=0;
    eventNumber=0;
    lumiBlock=0;
    BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
    BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
    BOOM->SetBranchAddress("nTruePUInteractions", &ntruePUInteractions, &b_ntruePUInteractions);
    BOOM->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
    BOOM->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    BOOM->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    BOOM->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

}

bool EventManager::isTriggerMatched()
    {
    static std::string myTrigger1 = "HLT_IsoMu24_";
    static std::string myTrigger2 = "HLT_IsoTkMu24_";
    for (int nTrig = 0; nTrig < Trigger_names->size(); ++nTrig)
	{
	if (Trigger_decision->at(nTrig) == 1)
	    {
	    string trigName = Trigger_names->at(nTrig);
		if ((trigName.find(myTrigger1) != string::npos) || (trigName.find(myTrigger2) != string::npos))
		    {
		    return true;
		    }
	    }
	}
    return false;
    }
bool EventManager::getIsRealData(){
    return this->isRealData;
}

int EventManager::getRunNumber(){
    return this->runNumber;
}

int EventManager::getEventNumber(){
    return this->eventNumber;
}
int EventManager::getLumiBlock(){
    return this->lumiBlock;
}

