// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __MUON_MU_H_                                                                                                                                  
#define __MUON_MU_H_

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
#include <TClonesArray.h>
#include <boost/regex.hpp>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "baseTree.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
using namespace std;
using namespace pat;
using namespace edm;

//
// class declaration
//

class MuonSelector : public  baseTree{

public:

  MuonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~MuonSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtx);
  bool IsHltMatched(edm::View<pat::Muon>::const_iterator mu, const edm::Event& iEvent, const edm::EventSetup& iSetup,
		   const edm::Handle<pat::TriggerObjectStandAloneCollection>& _trigObjsHandle,
		   const edm::Handle<edm::TriggerResults>& _trigResultsHandle,
		   const std::string _desiredTrigName, const double _muon_trig_dR );
private:
  MuonSelector(){};
  
  // ----------member data ---------------------------

  vector <double> Muon_isoSum , Muon_isoCharParPt ;
  vector <double> Muon_pt ,Muon_eta,Muon_phi, Muon_dz_pv, Muon_energy, Muon_iso;
  vector <double> Muon_isoCharged, Muon_isoNeutralHadron , Muon_isoPhoton, Muon_isoPU, Muon_combinedIso, Muon_trackRe_iso;
  vector <double> Muon_charge, Muon_chi2, Muon_matchedStat, Muon_dxy_pv, Muon_dxy_bs; 
  vector <double> Muon_dz_bs, Muon_dzError, Muon_dxyError, Muon_ndof, Muon_vtx, Muon_vty, Muon_vtz; 
  vector <double> Muon_track_pt, Muon_track_ptError, Muon_validHits, Muon_validHitsInner, Muon_TLayers; 
  vector <bool> Muon_isTriggerMatched, Muon_loose, Muon_medium, Muon_tight, Muon_soft, Muon_isHighPt, Muon_pf, Muon_isGlobal;   
  vector<double> Muon_dB, Muon_besttrack_pt, Muon_besttrack_ptError, Muon_tunePBestTrack_pt;
  vector<double> Muon_isTrackerMuon, Muon_POGisGood;
  vector<double> Muon_chi2LocalPosition, Muon_trkKink, Muon_segmentCompatibility, Muon_validFraction, Muon_pixelLayersWithMeasurement, Muon_qualityhighPurity;
  vector<double> Muon_tunePBestTrackType;
  vector<double> Muon_track_PCAx_bs, Muon_track_PCAy_bs, Muon_track_PCAz_bs;
  vector<double> Muon_track_PCAx_pv, Muon_track_PCAy_pv, Muon_track_PCAz_pv;
  vector<double> Muon_trackFitErrorMatrix_00, Muon_trackFitErrorMatrix_01, Muon_trackFitErrorMatrix_02, Muon_trackFitErrorMatrix_11, Muon_trackFitErrorMatrix_12, Muon_trackFitErrorMatrix_22;

  // super tiny ntuple?
  bool _super_TNT;

  // confit variables
//  edm::InputTag _muonToken;
  edm::EDGetTokenT<edm::View<pat::Muon> > _muonToken;
  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _trigObjsToken;
//  edm::InputTag _vertexInputTag;
//  edm::InputTag _beamSpot; 
  edm::EDGetTokenT<reco::VertexCollection> _vertexInputTag;
  edm::EDGetTokenT<reco::BeamSpot> _beamSpot;

  double _Muon_pt_min;
  double _Muon_eta_max;
  int _Muon_vtx_ndof_min;
  int _Muon_vtx_rho_max;
  double _Muon_vtx_position_z_max;

};

#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
