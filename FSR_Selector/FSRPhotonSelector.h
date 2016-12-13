// 
// Authors:  Luisa Fda. Chaparro: Universidad de los Andes, Colombia.  
//
#ifndef __FSR_PHOTON_H_                                                                                                                                  
#define __FSR_PHOTON_H_

#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>                                                                    
#include <TClonesArray.h>

#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <iomanip>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
// user include files

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/VectorUtil.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "baseTree.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"


using namespace std;
using namespace pat;
using namespace edm;
//
// class declaration
//

class FSRPhotonSelector : public  baseTree{

public:

	FSRPhotonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
	~FSRPhotonSelector();
	void Fill(const edm::Event& iEvent);
	void SetBranches();
	void Clear();
private:
	FSRPhotonSelector(){};
  
	// ----------member data ---------------------------
  
	vector<float> FSRPhoton_pt , FSRPhoton_eta, FSRPhoton_phi, FSRPhoton_energy;
	vector<float> FSRPhoton_pX, FSRPhoton_pY, FSRPhoton_pZ;
	vector<float>  FSRPhoton_isoNH, FSRPhoton_isoCH, FSRPhoton_isoCHPU, FSRPhoton_isoPhot;
	vector<float> FSRPhoton_isoNHPhot;
      	vector<float> FSRPhoton_et;//, FSRPhoton_HoverE, FSRPhoton_phoR9, FSRPhoton_SigmaIEtaIEta;
	//vector<float> FSRPhoton_SigmaIPhiIPhi, FSRPhoton_PFChIso, FSRPhoton_PFPhoIso, FSRPhoton_PFNeuIso;
	//vector<int>   FSRPhoton_EleVeto, FSRPhoton_hasPixelSeed;
 
	// confit variables
	//  edm::InputTag _PhotonToken;
	edm::EDGetTokenT<edm::View<pat::PFParticle> > _FSRPhotonToken;
	//double _Photon_pt_min;
	//double _Photon_eta_max;
};

#endif
