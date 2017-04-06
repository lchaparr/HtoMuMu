//=============================================================
// H2Mu Analysis code, The code follow the cuts implementation
// provide by H2Mu CMS group.
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToMuMu
// created by: Luisa Fda Chaparro, Uniandes.
//             Gaston Lyons, FNAL.
//==============================================================


#include "MuonManager.h"
#include "JetManager.h"
#include "EventManager.h"
#include "ElectronManager.h"
#include "cmdparser.hpp"
#include <fstream>
#include <iostream>


using namespace std;

std::string configure_parser(cli::Parser& program_options);
double const MASS_MUON = 0.105658367;

int main(int argc, char **argv)
    {
    ios_base::sync_with_stdio(false);
    cli::Parser program_options(argc, argv);
    std::string running_configuration = configure_parser(program_options);
    cout << running_configuration << endl;
    std::string input_file = program_options.getAlternative<std::string>("input-file");
    std::string output_file_txt = "H2Mu_"+input_file.substr(input_file.find_last_of("/")+1, input_file.find_last_of(".") - input_file.find_last_of("/") -1) + "__" +  running_configuration + ".txt";
    cout << "Trying to open Root file " << program_options.getAlternative<std::string>("input-file") << endl;
    TFile* f = TFile::Open((program_options.getAlternative<std::string>("input-file")).c_str());
    f->cd("TNT");
    TTree* tree = (TTree*) f->Get("TNT/BOOM");
    EventManager event(tree);
    ElectronManager electrons(tree);
    MuonManager muons(tree, event.getIsRealData(), program_options);
    JetManager jets(tree, program_options);
    int cVbfTight = 0;
    ofstream outputTxt;
    outputTxt.open("output.txt");
    ofstream outputSummary;
    outputSummary.open(output_file_txt);
    int cggFTight = 0;
    int cVbfLoose = 0;
    int czero1jetTight = 0;
    int czero1jetLoose = 0;
    int trigger = 0;
    int eventElectronVeto = 0;
    int eventBjetVeto = 0;
    int selectedEventCounter = 0;
    int nentries = tree->GetEntries();
    //nentries = 100;
    for (int i = 0; i < nentries; i++)
	{
	if (i%10000 == 0 ) cout<< "Event " << i << endl;
	tree->GetEntry(i);
	int eventNumber = event.getEventNumber();
	//cout << event.getEventNumber() << endl;
	// First we check if event passses the trigger
	//cout << event.isRealData << "\tEvent Number: "<< (event.eventNumber) <<"\tLumi:" <<  (event.lumiBlock)<<"\tRun:"<< (event.runNumber) << endl;
	bool eventPassesTrigger = false;
	if (event.isTriggerMatched())
	    {
	    eventPassesTrigger = true;
	    trigger++;
	    }
	// if the event passes the trigger we check if it passes the bjet veto.
	bool eventisNotBJetVetoed = false;
	if (eventPassesTrigger)
	    {
	    if (!jets.bJetVeto())
		{
		eventisNotBJetVetoed = true;
		}
	    else
		{
		eventBjetVeto++;
		}
	    }
	// If event passes the bjetVeto, we see if it passes Muon selection

	bool eventPassesMuonSelection = true;
	std::vector<Muon> selectedMuons;
	TLorentzVector dimuon;

	try
	    {
	    selectedMuons = muons.selectMuonsPFIso();
	    dimuon = selectedMuons.at(0).getTLorentz() + selectedMuons.at(1).getTLorentz();
	    if (dimuon.M() < 12.0)
		{
		eventPassesMuonSelection = false;
		}

	    }
	catch (int e)
	    {
	    eventPassesMuonSelection = false;
	    }

	// If there is no exception, there are 2 muons in selectedMuons. check electron veto
	bool eventisNotElectronVeto = false;
	if (eventPassesMuonSelection)
	    {
	    if (!electrons.electronVeto(selectedMuons))
		{
		eventisNotElectronVeto = true;
		}
	    else
		{
		//cout << event.getEventNumber() << endl;
		eventElectronVeto++;
		}
	    }
	// Event Selection
	bool eventPassesEventSelection = false;
	if (eventPassesTrigger && eventPassesMuonSelection && eventisNotElectronVeto && eventisNotBJetVetoed)
	    {
	    eventPassesEventSelection = true;
	    selectedEventCounter++;
	    }
	// At this point we know if the event is selected!
	std::vector<Jet> selectedJets;
	bool eventPassesJetPreselection = true;
	if (eventPassesEventSelection)
	    {
	    try
		{
		selectedJets = jets.selectJets(selectedMuons);
		}
	    catch (int e)
		{
		eventPassesJetPreselection = false;
		}
	    }
	// At this point we know that the event is selected, and if at least there are two jets in selectedJets
	bool VbfTight = false;
	bool ggFTight = false;
	bool VbfLoose = false;
	bool zero1jetTight = false;
	bool zero1jetLoose = false;
	if (eventPassesEventSelection && eventisNotBJetVetoed)
	    {
	    if (eventPassesJetPreselection)
		{
		VbfLoose = true;
		for (auto jet1 : selectedJets)
		    {
		    if (VbfTight == true)
			break;
		    for (auto jet2 : selectedJets)
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
		cVbfTight++;
		outputTxt << "VbfTight\t" << eventNumber << endl;
		}
	    else if (ggFTight)
		{
		cggFTight++;
		outputTxt << "ggFTight\t" << eventNumber <<  endl;

		}
	    else if (VbfLoose)
		{
		cVbfLoose++;
		outputTxt << "cVbfLoose\t" << eventNumber <<  endl;

		}
	    else if (zero1jetTight)
		{
		czero1jetTight++;
		outputTxt << "fzero1jetTight\t" << eventNumber <<  endl;

		}
	    else if (zero1jetLoose)
		{
		czero1jetLoose++;
		outputTxt << "fzero1jetLoose\t" << eventNumber << endl;

		}
	    VbfTight = false;
	    ggFTight = false;
	    VbfLoose = false;
	    zero1jetTight = false;
	    zero1jetLoose = false;
	    }
	}
    outputSummary << running_configuration << "| Trigger |" << trigger << endl;
    outputSummary << running_configuration << "| NOT pass ElectronVeto |" << eventElectronVeto << endl;
    outputSummary << running_configuration << "| NOT pass BjetVeto |" << eventBjetVeto << endl;
    outputSummary << running_configuration << "| VbfTight |" << cVbfTight << endl;
    outputSummary << running_configuration << "| ggFTight |" << cggFTight << endl;
    outputSummary << running_configuration << "| VbfLoose |" << cVbfLoose << endl;
    outputSummary << running_configuration << "| zero1jetTight |" << czero1jetTight << endl;
    outputSummary << running_configuration << "| zero1jetLoose |" << czero1jetLoose << endl;
    outputSummary << running_configuration << "|total select |" << selectedEventCounter << endl;
    outputSummary << running_configuration << "| Total events |" << cVbfTight + cggFTight + cVbfLoose + czero1jetTight + czero1jetLoose << endl;
    outputSummary << running_configuration << "| Total in file | " << tree->GetEntries() << endl;
    outputTxt.close();
    outputSummary.close();
    delete tree;
    delete f;
    cout << "Execution finished" << endl;
    return 0;
    }

std::string configure_parser(cli::Parser& program_options)
    {
    program_options.set_optional<std::string>("i", "input-file", "GF_Moriond_Ntuple.root", "Root file to analize");
    program_options.set_optional<std::string>("o", "output-file", "output.root", "Root file to write results to");
    program_options.set_optional<double>("r", "maxRelIso", 0.25, "Minimum pt Value for Muons");
    program_options.set_optional<double>("p", "minMuonPt", 10.0, "Minimum pt Value for Muons");
    program_options.set_optional<double>("j", "minJetPt", 30.0, "Maximum pt Value for jet selection");
    program_options.set_optional<bool>("g", "noGenMatching", false, "Use Gen object matching in Rochester corrections");
    program_options.set_optional<bool>("c", "noRochesterCorrections", false, "Apply Rochester corrections");
    program_options.run_and_exit_if_error();
    std::ostringstream oss;
    if (!(program_options.isNonDefault<double>("minMuonPt").empty()))
	oss << program_options.isNonDefault<double>("minMuonPt") << "-";
    if (!(program_options.isNonDefault<double>("maxRelIso").empty()))
	oss << program_options.isNonDefault<double>("maxRelIso") << "-";
    if (!(program_options.isNonDefault<double>("minJetPt").empty()))
	oss << program_options.isNonDefault<double>("minJetPt")<<"-";
    if (!(program_options.isNonDefault<bool>("noGenMatching").empty()))
	oss << program_options.isNonDefault<bool>("noGenMatching")<<"-";
    if (!(program_options.isNonDefault<bool>("noRochesterCorrections").empty()))
	oss << program_options.isNonDefault<bool>("noRochesterCorrections")<<"-";
    return oss.str();
    }
