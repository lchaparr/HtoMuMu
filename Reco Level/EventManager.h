#ifndef EVENTMANAGER_H_
#define EVENTMANAGER_H_
#include <TBranch.h>
#include <vector>
#include <TTree.h>
#include "Muon.h"
#include "Jet.h"
#include "FSRPhoton.h"
#include <TLorentzVector.h>
#include "AuxiliaryFunctions.h"
#include "HistogramManager.h"
using namespace std;

const std::string vbfTightStr = "VBF_Tight";
const std::string vbfLooseStr = "VBF_Loose";
const std::string ggfTightStr = "ggF_Tight";
const std::string zeroOr1JetTightStr = "0_or_1_Jet_Tight";
const std::string zeroOr1JetLooseStr = "0_or_1_Jet_Loose";
const std::string unclassifiedStr = "Unclassified";
// SubCarpetas
const string strRelIso = "/RelIso";
const std::string strNonRelIso = "/NonRelIso";
const std::string strAllEvents = "/AllEvents"; //all the selected events, including rel iso and no rel iso.
const std::string strFSREvents = "/FSREvents";
const std::string strNOFSREvents = "/NOFSREvents";
// Sub Sub Carpetas
const std::string tagOnlyDimuon = "/OnlyDimuon";
const std::string tagDimuonPlusPH = "/DimuonPlusPH";
const std::string tagInvariantMass = "OnlyDimuon DimuonPlusPH";
const std::string strAfterDimuonSelection = "AfterDimuonSelection";
const std::string strAfterDimuonRelIso = "/AfterDimuonRelIso";
const std::string oneGamma = "/oneGamma";
const std::string twoGamma = "/twoGamma";
const std::string passTrigger = "PassTrigger";
const std::string passDimuon = "passDimuon";

const std::vector<std::string> categoryNames = {vbfTightStr,vbfLooseStr,ggfTightStr, zeroOr1JetTightStr,zeroOr1JetLooseStr, unclassifiedStr};

class EventManager
    {
public:
    EventManager(TTree* BOOM, cli::Parser& program_options, std::string& running_configuration);
    TH1F *PUweights;
    TTree* BOOM;
    TFile* file_PUdata;
    TFile* file_PUsim;
    TH1F *PUsim;
    HistogramManager histos;
    std::map<bool, std::string> relIsoTag;
    std::map<int, std::string> fsrTag;
    std::vector<Muon> muons;
    std::vector<Jet> jets;
    std::vector<FSRPhoton> fsrPhotons;
    TLorentzVector Dimuon;
    std::string classifyEvent();
    std::string eventCategory;
    void setBranchAddress();
    virtual ~EventManager();
    bool passEventSelection;
    bool passDimuonSelection;
    bool isTriggerMatched();
    bool getIsRealData();
    int getRunNumber();
    int getEventNumber();
    int getLumiBlock();
    bool eventPassesTrigger;
    bool eventPassesMuonSelection;
    bool eventPassesRelIso;
    bool eventisNotElectronVeto;
    bool eventisNotBJetVetoed;
    bool eventPassesEventSelection();
    bool eventPassesJetPreselection;
    void saveEvent();
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
