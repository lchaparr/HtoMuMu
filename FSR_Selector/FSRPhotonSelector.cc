#include "NtupleMaker/BSM3G_TNT_Maker/interface/FSRPhotonSelector.h"

FSRPhotonSelector::FSRPhotonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC): baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the PhotonSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
//  _PhotonToken      = iConfig.getParameter<edm::InputTag>("photons");
  _FSRPhotonToken      = iCC.consumes<edm::View<pat::PFParticle> >(iConfig.getParameter<edm::InputTag>("fsrPhotonsSrc"));
  //_Photon_pt_min    = iConfig.getParameter<double>("Photon_pt_min");
  //_Photon_eta_max   = iConfig.getParameter<double>("Photon_eta_max");
  SetBranches();
}

FSRPhotonSelector::~FSRPhotonSelector(){
  delete tree_;
}

void FSRPhotonSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  edm::Handle<edm::View<pat::PFParticle> > FSRphoton_handle;
//  iEvent.getByLabel(_PhotonToken, photon_handle);
  iEvent.getByToken(_FSRPhotonToken, FSRphoton_handle);

  if(debug_) std::cout << "     FSRPhotonSelector: Cleared the vectors, grabbed the photon collection handle, and looping over photons." << std::endl;
 
  if (FSRphoton_handle.isValid()) { 
    for(edm::View<pat::PFParticle>::const_iterator ph = FSRphoton_handle->begin(); ph != FSRphoton_handle->end(); ph++){
      //if (ph->pt() < _Photon_pt_min) continue;
      //if (fabs(ph->eta()) > _Photon_eta_max) continue;  
      FSRPhoton_pt.push_back(ph->pt());
      FSRPhoton_eta.push_back(ph->eta());
      FSRPhoton_phi.push_back(ph->phi());
      FSRPhoton_energy.push_back(ph->energy());
      FSRPhoton_et.push_back(ph->et());
      /*FSRPhoton_HoverE.push_back(ph->hadTowOverEm());
      FSRPhoton_phoR9.push_back(ph->r9());
      FSRPhoton_SigmaIEtaIEta.push_back(ph->see());
      FSRPhoton_SigmaIPhiIPhi.push_back(ph->spp());
      FSRPhoton_PFChIso.push_back(ph->chargedHadronIso());
      FSRPhoton_PFPhoIso.push_back(ph->photonIso());
      FSRPhoton_PFNeuIso.push_back(ph->neutralHadronIso());
      FSRPhoton_EleVeto.push_back((int)ph->passElectronVeto());
      FSRPhoton_hasPixelSeed.push_back((int)ph->hasPixelSeed());*/
	  
	  FSRPhoton_pX.push_back(ph->px());
	  FSRPhoton_pY.push_back(ph->py());
	  FSRPhoton_pZ.push_back(ph->pz());
	  FSRPhoton_isoNH.push_back(ph->userFloat("fsrPhotonPFIsoNHad03"));
	  FSRPhoton_isoCH.push_back(ph->userFloat("fsrPhotonPFIsoChHad03pt02"));
	  FSRPhoton_isoCHPU.push_back(ph->userFloat("fsrPhotonPFIsoChHadPU03pt02"));
	  FSRPhoton_isoPhot.push_back(ph->userFloat("fsrPhotonPFIsoPhoton03"));
  	  FSRPhoton_isoNHPhot.push_back(ph->userFloat("fsrPhotonPFIsoNHadPhoton03"));
     }
  }
}

void FSRPhotonSelector::SetBranches(){
  if(debug_) std::cout << "     PhotonSelector: Setting branches by calling AddBranch of baseTree." << std::endl;
  
  AddBranch(&FSRPhoton_pt             ,"FSRPhoton_pt");
  AddBranch(&FSRPhoton_eta            ,"FSRPhoton_eta");
  AddBranch(&FSRPhoton_phi            ,"FSRPhoton_phi");
  AddBranch(&FSRPhoton_energy         ,"FSRPhoton_energy");
  AddBranch(&FSRPhoton_et             ,"FSRPhoton_et");
  /*AddBranch(&FSRPhoton_HoverE         ,"FSRPhoton_HoverE");
  AddBranch(&FSRPhoton_phoR9          ,"FSRPhoton_phoR9");
  AddBranch(&FSRPhoton_SigmaIEtaIEta  ,"FSRPhoton_SigmaIEtaIEta");
  AddBranch(&FSRPhoton_SigmaIPhiIPhi  ,"FSRPhoton_SigmaIPhiIPhi");
  AddBranch(&FSRPhoton_PFChIso        ,"FSRPhoton_PFChIso");
  AddBranch(&FSRPhoton_PFPhoIso       ,"FSRPhoton_PFPhoIso");
  AddBranch(&FSRPhoton_PFNeuIso       ,"FSRPhoton_PFNeuIso");
  AddBranch(&FSRPhoton_EleVeto        ,"FSRPhoton_EleVeto");
  AddBranch(&FSRPhoton_hasPixelSeed   ,"FSRPhoton_hasPixelSeed");*/
  
  AddBranch(&FSRPhoton_pX			  ,"FSRPhoton_pX");
  AddBranch(&FSRPhoton_pY			  ,"FSRPhoton_pY");
  AddBranch(&FSRPhoton_pZ			  ,"FSRPhoton_pZ");
  AddBranch(&FSRPhoton_isoNH		  ,"FSRPhoton_isoNH");
  AddBranch(&FSRPhoton_isoCH		  ,"FSRPhoton_isoCH");
  AddBranch(&FSRPhoton_isoCHPU		  ,"FSRPhoton_isoCHPU");
  AddBranch(&FSRPhoton_isoPhot		  ,"FSRPhoton_isoPhot");
  AddBranch(&FSRPhoton_isoNHPhot	  ,"FSRPhoton_isoNHPhot");

  if(debug_) std::cout << "     FSRPhotonSelector: Finished setting branches." << std::endl;
}

void FSRPhotonSelector::Clear(){
  
  FSRPhoton_pt.clear();
  FSRPhoton_eta.clear();
  FSRPhoton_phi.clear();
  FSRPhoton_energy.clear();
  FSRPhoton_et.clear();
  /*FSRPhoton_HoverE.clear();
  FSRPhoton_phoR9.clear();
  FSRPhoton_SigmaIEtaIEta.clear();
  FSRPhoton_SigmaIPhiIPhi.clear();
  FSRPhoton_PFChIso.clear();
  FSRPhoton_PFPhoIso.clear();
  FSRPhoton_PFNeuIso.clear();
  FSRPhoton_EleVeto.clear();
  FSRPhoton_hasPixelSeed.clear();*/
  
  FSRPhoton_pX.clear();
  FSRPhoton_pY.clear();
  FSRPhoton_pZ.clear();
  FSRPhoton_isoNH.clear();
  FSRPhoton_isoCH.clear();
  FSRPhoton_isoCHPU.clear();
  FSRPhoton_isoPhot.clear();
  FSRPhoton_isoNHPhot.clear();
}

