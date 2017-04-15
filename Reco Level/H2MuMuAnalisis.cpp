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
#include "FSRManager.h"
#include <TLorentzVector.h>
#include <fstream>
#include <iostream>


using namespace std;

std::string configure_parser(cli::Parser& program_options);
double const MASS_MUON = 0.1056583670;

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
    EventManager event(tree, program_options, running_configuration);
    ElectronManager electrons(tree);
    MuonManager muons(tree, event.getIsRealData(), program_options);
    JetManager jets(tree, program_options);
    FSRManager fsr(tree);
    ofstream outputTxt;
    outputTxt.open("output.txt");
    ofstream outputSummary;
    outputSummary.open(output_file_txt);
    int trigger = 0;
    int eventElectronVeto = 0;
    int eventBjetVeto = 0;
    int nentries = tree->GetEntries();
    map<string, int> counters;
    for (int i = 0; i < nentries; i++)
	{
	if (i%10000 == 0 ) cout<< "Event " << i << endl;
	tree->GetEntry(i);
	int eventNumber = event.getEventNumber();
	// First we check if event passses the trigger
	event.eventPassesTrigger = false;
	if (event.isTriggerMatched())
	    {
	    event.eventPassesTrigger = true;
	    trigger++;
	    }
	// if the event passes the trigger we check if it passes the bjet veto.
	event.eventisNotBJetVetoed = false;
	if (event.eventPassesTrigger)
	    {
	    if (!jets.bJetVeto())
		{
		event.eventisNotBJetVetoed = true;
		}
	    else
		{
		eventBjetVeto++;
		}
	    }
	// If event passes the bjetVeto, we see if it passes Muon selection

	event.eventPassesMuonSelection = true;
	event.eventPassesRelIso = true;
	TLorentzVector dimuon;
	try
	    {
	    event.muons = muons.selectMuonsPFIso();
	    event.Dimuon = event.muons.at(0).getTLorentz() + event.muons.at(1).getTLorentz();
	    if (event.Dimuon.M() < 12.0)
		{
		event.eventPassesMuonSelection = false;
		event.eventPassesRelIso = false;

		}

	    }
	catch (int e)
	    {
	    event.eventPassesMuonSelection = false;
	    }

	// If the event passes muon selection, let's do FSR Recovery
	if (program_options.getAlternative<bool>("fsr")){
	if (event.eventPassesMuonSelection){
	 std::vector<FSRPhoton> fsrPhotons2;
	 std::vector<FSRPhoton> fsrPhotons = fsr.selectFSRPhotons(event.muons);
	 TLorentzVector DimuonPH = event.Dimuon;
	 for (auto fsrph :fsrPhotons){
	     DimuonPH += fsrph.getTLorentz();
	 }
	 if (!(abs(125.0 - event.Dimuon.M()) <= abs(125.0 - DimuonPH.M())))
	 {
	     event.fsrPhotons = fsrPhotons;
	 }
	 else{
	     event.fsrPhotons = fsrPhotons2;
	 }
	}
	}
	// If there is no exception, there are 2 muons in event.muons. check electron veto
	event.eventisNotElectronVeto = false;
	if (event.eventPassesMuonSelection)
	    {
	    if (!electrons.electronVeto(event.muons))
		{
		event.eventisNotElectronVeto = true;
		}
	    else
		{
		eventElectronVeto++;
		}
	    }

	// At this point we know if the event is selected!
	event.eventPassesJetPreselection = true;
	if (event.eventPassesEventSelection())
	    {
	    try
		{
		event.jets = jets.selectJets(event.muons);
		}
	    catch (int e)
		{
		event.eventPassesJetPreselection = false;
		}
	    }
	// At this point we know that the event is selected, and if at least there are two jets in event.jets, let's classify it.
	    string classification = event.classifyEvent();
	    event.saveEvent();
	    counters[classification] = counters[classification] + 1;
	}
    for (auto counter : counters){
	outputSummary << running_configuration << "|" << counter.first << "|" << counter.second << endl;
    }

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
    program_options.set_optional<bool>("f", "fsr", false, "Use FSR photons");
    program_options.set_optional<bool>("c", "noRochesterCorrections", false,  "Apply Rochester corrections");
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
    if (!(program_options.isNonDefault<bool>("fsr").empty()))
    	oss << program_options.isNonDefault<bool>("fsr")<<"-";
    return oss.str();
    }
