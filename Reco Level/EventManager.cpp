#include "EventManager.h"

EventManager::EventManager(TTree* BOOM, cli::Parser& program_options, std::string& running_configuration)
    {
    this->BOOM = BOOM;
    this->setBranchAddress();
    relIsoTag[true] = "RelIso";
    relIsoTag[false] ="NonRelIso";
    fsrTag[0] = "noFSR";
    fsrTag[1] = "oneGamma";
    fsrTag[2] = "TwoGamma";
    std::string input_file = program_options.getAlternative<std::string>("input-file");
    std::string output_file= "H2Mu_"+input_file.substr(input_file.find_last_of("/")+1, input_file.find_last_of(".") - input_file.find_last_of("/") -1) + "__" +  running_configuration + ".root";
    this->histos.init(output_file);
    this->histos.createMuonHistograms(strAfterDimuonSelection);

    for (auto categoryName : categoryNames){
	// All Events histograms
	this->histos.createDimuonHistograms(categoryName + strAllEvents);
	//RelIso Histograms
	histos.createDimuonHistograms(categoryName + strRelIso + strNOFSREvents, "OnlyDimuon");
	histos.createDimuonHistograms(categoryName + strRelIso + strFSREvents, tagInvariantMass);
	//Non RelIso Histograms
	histos.createDimuonHistograms(categoryName + strNonRelIso + strNOFSREvents, "OnlyDimuon");
	histos.createDimuonHistograms(categoryName + strNonRelIso + strFSREvents, tagInvariantMass);
    }

	//-----------------------load PU weights----------------------------

	this->file_PUdata = new TFile("PU2016data_15p9ifb.root", "read");
	this->PUweights = (TH1F*) file_PUdata->Get("analyzeHiMassTau/NVertices_0");
	PUweights->Scale(1 / PUweights->Integral());
	this->file_PUsim= new TFile("PU2016MC.root", "read");
	TH1F *PUsim = (TH1F*) file_PUsim->Get("analyzeHiMassTau/NVertices_0");
	PUsim->Scale(1 / PUsim->Integral());
	this->PUweights->Divide(PUsim);
    }

EventManager::~EventManager()
    {
    this->histos.writeTFile();
    delete file_PUdata;
    delete file_PUsim;
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

std::string EventManager::classifyEvent(){
    bool VbfTight = false;
    	bool ggFTight = false;
    	bool VbfLoose = false;
    	bool zero1jetTight = false;
    	bool zero1jetLoose = false;
    	if (eventPassesEventSelection() && eventisNotBJetVetoed)
    	    {
            TLorentzVector dimuon = muons.at(0).getTLorentz() + muons.at(1).getTLorentz();

    	    if (eventPassesJetPreselection)
    		{
    		VbfLoose = true;
    		for (auto jet1 : this->jets)
    		    {
    		    if (VbfTight == true)
    			break;
    		    for (auto jet2 : this->jets)
    			{
    			if (!aux::isSameParticle<Jet>(jet1, jet2) && (jet1.pt() > 40 || jet2.pt() > 40))
    			    {
    			    TLorentzVector dijet = jet1.getTLorentz() + jet2.getTLorentz();
    			    if (dijet.M() > 650.0 && abs(jet1.eta() - jet2.eta()) > 3.5)
    				{
    				VbfTight = true;
    				break;
    				}
    			    if (dijet.M() > 250.0 && dimuon.Pt() > 50.0)
    				{
    				ggFTight = true;
    				}
    			    }
    			}
    		    }
    		}
    	    else if (!eventPassesJetPreselection)
    		{
    		zero1jetLoose = true;
    		if (dimuon.Pt() >= 25.0)
    		    {
    		    zero1jetTight = true;
    		    }
    		}

    	    if (VbfTight)
    		{
    		this->eventCategory = vbfTightStr;
    		return this->eventCategory;

    		}
    	    else if (ggFTight)
    		{
    		this->eventCategory =  ggfTightStr;
    		return this->eventCategory;


    		}
    	    else if (VbfLoose)
    		{
    		this->eventCategory =  vbfLooseStr;
    		return this->eventCategory;


    		}
    	    else if (zero1jetTight)
    		{
    		this->eventCategory =  zeroOr1JetTightStr;
    		return this->eventCategory;


    		}
    	    else if (zero1jetLoose)
    		{
    		this->eventCategory =  zeroOr1JetLooseStr;
    		return this->eventCategory;

    		}

    	    }

this->eventCategory = unclassifiedStr;
return this->eventCategory;
}
bool EventManager::eventPassesEventSelection(){
    bool ret = this->eventPassesTrigger && this->eventPassesMuonSelection && this->eventisNotElectronVeto && this->eventisNotBJetVetoed;
    return ret;
}
void EventManager::saveEvent()
    {
    double pu_weight = 1.0;

    pu_weight = PUweights->GetBinContent(PUweights->FindBin(ntruePUInteractions));

    std::string relIsoEventTag = this->eventPassesRelIso ? strRelIso : strNonRelIso;
    std::string fsrEventTag = this->fsrPhotons.size() == 0 ? strNOFSREvents : strFSREvents;
    if (this->eventCategory.empty()) this->classifyEvent();
    if (this->eventPassesEventSelection()){
    //Fill All Events in Category
	histos.fillMuonHist(strAfterDimuonSelection, this->muons.at(0).getTLorentz(), this->muons.at(1).getTLorentz(), pu_weight);
	histos.fillHiggsHist(this->eventCategory + strAllEvents, Dimuon, pu_weight);
    // Fill Event in RelIso/NonRelIsoFolder
	if (this->fsrPhotons.size() ==0){
		histos.fillHiggsHist(this->eventCategory  + relIsoEventTag + strNOFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
	} else if (this->fsrPhotons.size() ==1){
	    histos.fillHiggsHist(this->eventCategory  + relIsoEventTag + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
	    histos.fillHiggsHist(this->eventCategory  + relIsoEventTag + strFSREvents + tagDimuonPlusPH, Dimuon, fsrPhotons.at(0).getTLorentz(), pu_weight);
	} else if (this->fsrPhotons.size() ==2){
	    histos.fillHiggsHist(this->eventCategory  + relIsoEventTag + strFSREvents + tagOnlyDimuon, Dimuon, pu_weight);
	    histos.fillHiggsHist(this->eventCategory  + relIsoEventTag + strFSREvents + tagDimuonPlusPH, Dimuon, fsrPhotons.at(0).getTLorentz()+fsrPhotons.at(1).getTLorentz(), pu_weight);

	}

    }
    }
