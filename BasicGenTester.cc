#include <iostream>

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// essentials !!!
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepPDT/ParticleDataTable.hh"

class BasicGenTester : public edm::EDAnalyzer
{
  
public:
  
	//
	explicit BasicGenTester( const edm::ParameterSet& ) ;
	virtual ~BasicGenTester() {} // no need to delete ROOT stuff
	// as it'll be deleted upon closing TFile
  
	virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
	virtual void beginJob() ;
	virtual void beginRun( const edm::Run &, const edm::EventSetup& );
	virtual void endRun( const edm::Run&, const edm::EventSetup& ) ;
	virtual void endJob() ;
  
private:
  
	TH1D*       fPtHiggsProduced;
	TH1D*       fMassHiggsProduced;
	TH1D*       fEtaHiggsProduced;
	TH1D*       fPhiHiggsProduced ;
	TH1D*       fMassDimuon;
	TH1D*       fNumberPhotons;
	TH1D*       fEnergyofPhotons;
	TH1D*       fMassDimuonOver10;
	TH1D*       fEnergyPhotonsOver10OverHmass;
	TH1D*	    fDeltaRPhotonMuon1OverHmass;
	TH1D*	    fDeltaRPhotonMuon2OverHmass;
	TH1D*	    fDeltaRPhotonDimuonOverHmass;
	TH1D*       fEnergyPhotonsOver10UnderHmass;
	TH1D*	    fDeltaRPhotonMuon1UnderHmass;
	TH1D*	    fDeltaRPhotonMuon2UnderHmass;
	TH1D*	    fDeltaRPhotonDimuonUnderHmass;	
	TH1D*       fMassDimuonOnephoton;
	TH1D*       fEtaPositiveMuon;
	TH1D*	    fEtaNegativeMuon;
	TH1D*       fPtPositiveMuon;
	TH1D*       fPtNegativeMuon;


	int         mflag =0;
	int         photonFSR =0;
	int         fNPart;
	double      fPtMin;
	double      fPtMax;
	edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;
  
}; 

using namespace edm;
using namespace std;
using namespace reco;

BasicGenTester::BasicGenTester( const ParameterSet& pset )
	: fPtHiggsProduced(0),  fMassHiggsProduced(0), fEtaHiggsProduced(0), fPhiHiggsProduced(0),fMassDimuon(0),fNumberPhotons(0),fEnergyofPhotons(0),fMassDimuonOver10(0),fEnergyPhotonsOver10OverHmass(0),fDeltaRPhotonMuon1OverHmass(0),fDeltaRPhotonMuon2OverHmass(0),fDeltaRPhotonDimuonOverHmass(0),fEnergyPhotonsOver10UnderHmass(0), fDeltaRPhotonMuon1UnderHmass(0), fDeltaRPhotonMuon2UnderHmass(0),fDeltaRPhotonDimuonUnderHmass(0),fMassDimuonOnephoton(0), fEtaPositiveMuon(0), fEtaNegativeMuon(0), fPtPositiveMuon(0),fPtNegativeMuon(0)

	      		 
{
  
	fNPart = pset.getUntrackedParameter<int>( "NPartForHisto", 500 );
	fPtMin = pset.getUntrackedParameter<double>( "PtMinForHisto",  0. );
	fPtMax = pset.getUntrackedParameter<double>( "PtMaxForHisto", 25. );   
  
}

void BasicGenTester::beginJob()
{
  
	Service<TFileService> fs;
  
	fPtHiggsProduced = fs->make<TH1D>(  "PtHiggsProduced","Higgs Pt;P_T(GeV);Events", 150, 0, 150);
	fMassHiggsProduced = fs->make<TH1D>("MassHiggsProduced", "Higgs Mass;Mass(GeV);Events",250, 0, 250);
	fEtaHiggsProduced = fs->make<TH1D>("EtaHiggsProduced", ";#eta;Events", 80, -4, 4);
	fPhiHiggsProduced = fs->make<TH1D>("PhiHiggsProduced", ";#phi;Events", 70,-3.5,3.5);
	fMassDimuon = fs->make<TH1D>("MassDimuon","Dimuon Mass ;Mass (GeV);Events", 250, 0, 250);
	fNumberPhotons = fs->make<TH1D>("NumberofPhotons","Number of Photons; Photons per event; Events", 500,0,500);
	fEnergyofPhotons = fs->make<TH1D>("EnergyOfPhotons","Energy of Photons;Energy(GeV);Photons",150,0,150);
	fMassDimuonOver10 = fs->make<TH1D>("MassDimuon (Photons over 10GEV)","Dimuon mass ;Mass (GeV);Events", 250, 0, 250);
	fEnergyPhotonsOver10OverHmass = fs->make<TH1D>("EnergyPhotonsOver10OverHmass","Energy of Photons;Energy(GeV);Photons",150,0,150);
	fDeltaRPhotonMuon1OverHmass = fs->make<TH1D>("DeltaRPhotonMuon1OverHmass", "#Delta R (#gamma #mu_1); #Delta R; Events", 100, 0, 10);
	fDeltaRPhotonMuon2OverHmass = fs->make<TH1D>("DeltaRPhotonMuon2OverHmass", "#Delta R (#gamma #mu_2); #Delta R; Events", 100, 0, 10);
	fDeltaRPhotonDimuonOverHmass = fs->make<TH1D>("DeltaRPhotonDimuonOverHmass", "#Delta R (#gamma #mu#mu); #Delta R; Events", 100, 0, 10);
	fEnergyPhotonsOver10UnderHmass = fs->make<TH1D>("EnergyPhotonsOver10UnderHmass","Energy of Photons;Energy(GeV);Photons",150,0,150);
	fDeltaRPhotonMuon1UnderHmass = fs->make<TH1D>("DeltaRPhotonMuon1UnderHmass", "#Delta R (#gamma #mu_1); #Delta R; Events", 100, 0, 10);
	fDeltaRPhotonMuon2UnderHmass = fs->make<TH1D>("DeltaRPhotonMuon2UnderHmass", "#Delta R (#gamma #mu_2); #Delta R; Events", 100, 0, 10);
	fDeltaRPhotonDimuonUnderHmass = fs->make<TH1D>("DeltaRPhotonDimuonUnderHmass", "#Delta R (#gamma #mu#mu); #Delta R; Events", 100, 0, 10);
	fMassDimuonOnephoton = fs->make<TH1D>("MassDimuon (one photon FSR)","Dimuon mass + FSR Photon ;Mass (GeV);Events", 200, 0, 200);	
	fEtaPositiveMuon = fs->make<TH1D>("EtaPositiveMuon", ";#eta;Events", 80, -4, 4); 
	fEtaNegativeMuon = fs->make<TH1D>("EtaNegativeMuon", ";#eta;Events", 80, -4, 4);
	fPtPositiveMuon = fs->make<TH1D>(  "PtPositiveMuon","#mu+ pT;P_T(GeV);Events", 150, 0, 150);
	fPtNegativeMuon = fs->make<TH1D>(  "PtNegativeMuon","#mu- pT;P_T(GeV);Events", 150, 0, 150);

 

	return ;
  
}

void BasicGenTester::beginRun( const edm::Run& r, const edm::EventSetup& es )
{
  
	es.getData( fPDGTable ) ;
  
	return ;
  
}

void BasicGenTester::analyze( const Event& e, const EventSetup& )
{
  
	edm::Handle<edm::HepMCProduct > EvtHandle ;
	e.getByLabel( "generator", EvtHandle ) ;
  
	const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

	int iflag = 0, jflag = 0, phflag = 0; 
	int photon1 = 0, photon2 = 0;
	double px1 = 0, px2 = 0, py1 = 0, py2 = 0, pz1 = 0, pz2 = 0, en1 = 0, en2 = 0;
	double px3 = 0, py3 = 0, pz3 = 0, en3 = 0, pxp = 0, pyp = 0, pzp = 0;
	double ep = 0,pxp2 = 0, pyp2 = 0, pzp2 = 0, ep2 = 0,px4 = 0, py4 = 0, pz4 = 0, en4 = 0;
	double etamu1 = 0, etamu2 = 0, etap1 = 0; 
	double phimu1 = 0, phimu2 = 0, phip1 = 0; 
	double phi3 = 0, eta3 = 0;
	double massm = 0, massmp = 0, Egamma = 0;
	double deltar1 = 0, deltar2 = 0, deltar12 = 0;
	double Hmass = 198, Hsample=200.05;
	int phfsr = 0, NPhotons = 0,NPhotonOver10 = 0;


	for ( HepMC::GenEvent::particle_const_iterator part = Evt->particles_begin();
	part != Evt->particles_end(); ++part ) 
	{
      
		int pid = (*part)->pdg_id();
		int stat = (*part)->status();

		double px = (*part)->momentum().px();
		double py = (*part)->momentum().py();
		double pz = (*part)->momentum().pz();
		double e = (*part)->momentum().e();
		
		//double pt = ((*part)->momentum()).perp();
		
		double eta = (*part)->momentum().eta();
		double phi = (*part)->momentum().phi();
		double mass = sqrt(e*e -px*px - py*py - pz*pz);
		
		//----------Higgs production
		if ( abs(pid) == 25 && ((*part)->parent_event()) && (stat == 22) )
		{      
			fPtHiggsProduced->Fill( ((*part)->momentum()).perp() );
			fMassHiggsProduced->Fill(mass);
			fEtaHiggsProduced->Fill( ((*part)->momentum()).eta() );
			fPhiHiggsProduced->Fill( ((*part)->momentum()).phi() );
		}    

		//----------------Photons
		if(pid==22){

			NPhotons++;
			Egamma=e;
			fEnergyofPhotons->Fill(Egamma);
			
			if( ((*part)->momentum()).perp()>10){
				NPhotonOver10++;

				if(phflag==0){
					pxp = px;
					pyp = py;
					pzp = py;
					//ptp = pt;
					ep = e;
					etap1 = eta;
					phip1 = phi;
					photon1=1;
					
					
				}
				else{
					pxp2 = px;
					pyp2 = py;
					pzp2 = py;
					ep2 = e;
					//etap2 = eta;
					//phip2 = phi;
					photon2=1;
				}
				phfsr++;
				phflag++;
			}
		}
     
     
		//---------Muons
		if((pid==13||pid==-13) && (stat == 1) && ((*part)->production_vertex()) )
		{  
			//const HepPDT::ParticleData*  PData = fPDGTable->particle(HepPDT::ParticleID(pid));
			//double charge = PData->charge();
			if (pid==13 && iflag==0)
			{
				px1 = px;
				py1 = py;
				pz1 = pz;
				en1 = e;
				
				etamu1 = eta;
				phimu1 = phi;
				iflag++;
				
				fEtaPositiveMuon->Fill(etamu1); 
				fPtPositiveMuon->Fill(((*part)->momentum()).perp());
			}
	  
			if(pid==-13 && jflag==0)
			{
				px2 = px;
				py2 = py;
				pz2 = pz;
				en2 = e;
				
				etamu2 = eta;
				phimu2 = phi;
				jflag++;
				
				fEtaNegativeMuon->Fill(etamu2); 
				fPtNegativeMuon->Fill(((*part)->momentum()).perp());
			}
		}
	}//close loop part.


	//---------Dimuon Mass
	if(iflag==1 && jflag==1)
	{
		massm = sqrt(((en1 + en2)*(en1 + en2))-((px1 + px2)*(px1 + px2))-((py1 + py2)*(py1 + py2))-((pz1 +pz2)*(pz1+pz2)) );
		phi3 = phimu1+phimu2;
		eta3 = etamu1+etamu2;
	}
	fMassDimuon->Fill(massm);
	std::cout<<"masa:   "<<massm<<std::endl;
  
	//--------Photons Plots
 
	fNumberPhotons->Fill(NPhotons);
   
	if(photon1==1)
	{
		deltar1 = sqrt((etamu1-etap1)*(etamu1-etap1)+(phimu1-phip1)*(phimu1-phip1));
		deltar2 = sqrt((etamu2-etap1)*(etamu2-etap1)+(phimu2-phip1)*(phimu2-phip1));
		deltar12 = sqrt((eta3-etap1)*(eta3-etap1)+(phi3-phip1)*(phi3-phip1));
			
		fMassDimuonOver10->Fill(massm);
		std::cout<<"masa 10:   "<<massm<<std::endl;

		if(massm>Hmass){
			fEnergyPhotonsOver10OverHmass->Fill(en1);
			fDeltaRPhotonMuon1OverHmass->Fill(deltar1);
			fDeltaRPhotonMuon2OverHmass->Fill(deltar2);
			fDeltaRPhotonDimuonOverHmass->Fill(deltar12);
			
		}

		
		if(massm<Hmass && phfsr!=0)
		{
			mflag++;
			photonFSR = photonFSR+phfsr;
			
			fEnergyPhotonsOver10UnderHmass->Fill(en1);
			fDeltaRPhotonMuon1UnderHmass->Fill(deltar1);
			fDeltaRPhotonMuon2UnderHmass->Fill(deltar2);
			fDeltaRPhotonDimuonUnderHmass->Fill(deltar12);
			
			px3 = px1+px2;
			py3 = px1+py2;
			pz3 = pz1+pz2;
			en3 = en1+en2;
			
			massmp = sqrt( ((en3+ep)*(en3+ep))-((px3+pxp)*(px3+pxp))-((py3+pyp)*(py3+pyp))-((pz3+pzp)*(pz3+pzp)) ); 
			if(massmp<Hsample){	
				fMassDimuonOnephoton->Fill(massmp);
				std::cout<<"masa ph:   "<<massmp<<std::endl;
			}
		}
	}
	
  
	
  
	if((massm<118.5)&&(photon2==1))
	{
		px4 = px3+pxp;
		py4 = px3+pyp;
		pz4 = pz3+pzp;
		en4 = en3+ep;
		massm = sqrt( ((en4+ep2)*(en4+ep2))-((px4+pxp2)*(px4+pxp2))-((py4+pyp2)*(py4+pyp2))-((pz4+pzp2)*(pz4+pzp2)) ); 

  
	}
  
 
	
	std:: cout<<"No higgs mass events: "<<mflag<< "  FSR Photons per event: "<<phfsr<<"FSR Photons in total: "<<photonFSR<<std::endl;
  
  
	return ;
  
}

void BasicGenTester::endRun( const edm::Run& r, const edm::EventSetup& )
{
  
	return;
  
}


void BasicGenTester::endJob()
{
  
	return ;
}

DEFINE_FWK_MODULE(BasicGenTester);
