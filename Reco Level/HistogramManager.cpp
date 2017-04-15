/*
 * HistogramManager.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: gastonlp
 */

#include "HistogramManager.h"

using namespace std;
#include <iostream>
#include <libgen.h>
#include "TSystem.h"
#include <sstream>


HistogramManager::HistogramManager()
{
}

HistogramManager::HistogramManager(std::string filename)
{
	this->init(filename);
}

void HistogramManager::init(std::string filename)
{
	//We create the output root file
	rootFile = new TFile(filename.c_str(), "RECREATE");
}


void HistogramManager::addDirectory(std::string directoryName)
	{
	if (directories.count(directoryName) == 0)
		{
		std::size_t found = directoryName.find_last_of("/");
		if (found != std::string::npos)
			{
			this->addDirectory(directoryName.substr(0, found));
			}
		const char* directoryNameC = directoryName.c_str();
		TDirectory* createdDirectory = rootFile->mkdir(directoryNameC);
		createdDirectory->cd();
		directories[directoryName] = createdDirectory;
		cout << "Created directory: " << directoryName << endl;
		}
	}

void HistogramManager::addHistogram(std::string directoryName, std::string histogramName, std::string histogramTitle,
		Int_t nbinsx, Double_t xlow, Double_t xup)
{
	const char* directoryNameC = directoryName.c_str();
	const char* histogramNameC = histogramName.c_str();
	const char* histogramTitleC = histogramTitle.c_str();

	cout << "Received addHIst " << directoryName << endl;
	if (directories.count(directoryName) == 0)
	{
		this->addDirectory(directoryName);
	}
	rootFile->cd(directoryNameC);
	TH1* histogram = new TH1F(histogramNameC, histogramTitleC, nbinsx, xlow, xup);
	histograms[directoryName][histogramName] = histogram;
}

void HistogramManager::fillHistogram(std::string directoryName, std::string histogramName, Double_t x, Double_t w)
{
	histograms[directoryName][histogramName]->Fill(x, w);
}

void HistogramManager::fillHistogram(std::string directoryName, std::string histogramName, Double_t x)
{
	histograms[directoryName][histogramName]->Fill(x);
}

void HistogramManager::createMuonHistograms(std::string directoryName)
{
	this->addHistogram(directoryName, (char*) "twoLepton", (char*) "m_{#mu, #mu}", 600., 0., 300.);
	this->addHistogram(directoryName, (char*) "lead_muon_pT", (char*) "#mu p_{T}", 600, 0., 300.);
	this->addHistogram(directoryName, (char*) "slead_muon_pT", (char*) "#mu p_{T}", 600, 0., 300.);
	this->addHistogram(directoryName, (char*) "lead_muon_eta", (char*) "#mu #eta", 100, -2.5, 2.5);
	this->addHistogram(directoryName, (char*) "slead_muon_eta", (char*) "#mu #eta", 100, -2.5, 2.5);
	this->addHistogram(directoryName, (char*) "lead_muon_phi", (char*) "#mu #phi", 140, -3.5, 3.5);
	this->addHistogram(directoryName, (char*) "slead_muon_phi", (char*) "#mu #phi", 140, -3.5, 3.5);
}
void HistogramManager::createPhotonHistograms(std::string directoryName)
{
	cout << "Received name photon " << directoryName << endl;
	this->addHistogram(directoryName, (char*) "PhotonEnergy_Lep1", (char*) "E_{#gamma}", 200., 0., 100.);
	this->addHistogram(directoryName, (char*) "PhotonpT_Lep1", (char*) "pT_{#gamma}", 200., 0., 100.);
	this->addHistogram(directoryName, (char*) "PhotonEta_Lep1", (char*) "#eta_{#gamma}", 100, -2.5, 2.5);
	this->addHistogram(directoryName, (char*) "PhotonPhi_Lep1", (char*) "#phi_{#gamma}", 140, -3.5, 3.5);
	this->addHistogram(directoryName, (char*) "PhotonEnergy_Lep2", (char*) "E_{#gamma}", 200., 0., 100.);
	this->addHistogram(directoryName, (char*) "PhotonpT_Lep2", (char*) "pT_{#gamma}", 200., 0., 100.);
	this->addHistogram(directoryName, (char*) "PhotonEta_Lep2", (char*) "#eta_{#gamma}", 100, -2.5, 2.5);
	this->addHistogram(directoryName, (char*) "PhotonPhi_Lep2", (char*) "#phi_{#gamma}", 140, -3.5, 3.5);
}
void HistogramManager::createDimuonHistograms(std::string directoryName)
{
	this->addHistogram(directoryName, (char*) "diMuonMass", (char*) "m_{#mu, #mu}", 600., 0., 300.);
	this->addHistogram(directoryName, (char*) "diMuon_pT", (char*) "pT_{#mu,#mu}", 600, 0., 300.);
	this->addHistogram(directoryName, (char*) "diMuon_eta", (char*) "#eta_{#mu,#mu}", 100, -2.5, 2.5);
	this->addHistogram(directoryName, (char*) "diMuon_phi", (char*) "#phi_{#mu,#mu}", 140, -3.5, 3.5);
}
void HistogramManager::fillMuonHist(std::string directoryName, const TLorentzVector& Lepton1, const TLorentzVector& Lepton2,
		Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "twoLepton", (Lepton1 + Lepton2).M(), weight);
	this->fillHistogram(directoryName, (char*) "lead_muon_pT", Lepton1.Pt(), weight);
	this->fillHistogram(directoryName, (char*) "slead_muon_pT", Lepton2.Pt(), weight);
	this->fillHistogram(directoryName, (char*) "lead_muon_eta", Lepton1.Eta(), weight);
	this->fillHistogram(directoryName, (char*) "slead_muon_eta", Lepton2.Eta(), weight);
	this->fillHistogram(directoryName, (char*) "lead_muon_phi", Lepton1.Phi(), weight);
	this->fillHistogram(directoryName, (char*) "slead_muon_phi", Lepton2.Phi(), weight);
}

void HistogramManager::fillHiggsHist(std::string directoryName, const TLorentzVector& Dimuon, Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "diMuonMass", (Dimuon).M(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_pT", (Dimuon).Pt(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_eta", (Dimuon).Eta(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_phi", (Dimuon).Phi(), weight);
}

void HistogramManager::fillHiggsHist(std::string directoryName,const  TLorentzVector& Dimuon, const TLorentzVector& PhotonFSR,
		Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "diMuonMass", (Dimuon + PhotonFSR).M(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_pT", (Dimuon + PhotonFSR).Pt(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_eta", (Dimuon + PhotonFSR).Eta(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_phi", (Dimuon + PhotonFSR).Phi(), weight);
}

void HistogramManager::fillHiggsHist(std::string directoryName,const  TLorentzVector& Dimuon, const TLorentzVector& PhotonFSR1,
	const TLorentzVector& PhotonFSR2, Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "diMuonMass", (Dimuon + PhotonFSR1 + PhotonFSR2).M(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_pT", (Dimuon + PhotonFSR1 + PhotonFSR2).Pt(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_eta", (Dimuon + PhotonFSR1 + PhotonFSR2).Eta(), weight);
	this->fillHistogram(directoryName, (char*) "diMuon_phi", (Dimuon + PhotonFSR1 + PhotonFSR2).Phi(), weight);
}

void HistogramManager::fillPhotonHist(std::string directoryName, const TLorentzVector& FSR_Photon1_TL,
	const TLorentzVector& FSR_Photon2_TL, Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "PhotonEnergy_Lep1", FSR_Photon1_TL.E(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonpT_Lep1", FSR_Photon1_TL.Pt(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonEta_Lep1", FSR_Photon1_TL.Eta(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonPhi_Lep1", FSR_Photon1_TL.Phi(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonEnergy_Lep2", FSR_Photon2_TL.E(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonpT_Lep2", FSR_Photon2_TL.Pt(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonEta_Lep2", FSR_Photon2_TL.Eta(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonPhi_Lep2", FSR_Photon2_TL.Phi(), weight);
}

void HistogramManager::fillPhotonHist(std::string directoryName,const  TLorentzVector& FSR_Photon1_TL, Double_t weight)
{
	this->fillHistogram(directoryName, (char*) "PhotonEnergy_Lep1", FSR_Photon1_TL.E(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonpT_Lep1", FSR_Photon1_TL.Pt(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonEta_Lep1", FSR_Photon1_TL.Eta(), weight);
	this->fillHistogram(directoryName, (char*) "PhotonPhi_Lep1", FSR_Photon1_TL.Phi(), weight);
}

void HistogramManager::createPhotonHistograms(std::string directoryName, std::string tags)
{
	std::stringstream ss(tags);
	string token;
	while (ss >> token)
	{
		std::string newString = directoryName + "/" + token;
		cout << directoryName + "/" + token << endl;
		this->createPhotonHistograms(newString);
	}
}

void HistogramManager::createMuonHistograms(std::string directoryName, std::string tags)
{
	std::stringstream ss(tags);
	string token;
	while (ss >> token)
	{
		std::string newString = directoryName + "/" + token;
		cout << directoryName + "/" + token << endl;
		this->createMuonHistograms(newString);
	}
}

void HistogramManager::createDimuonHistograms(std::string directoryName, std::string tags)
{
	std::stringstream ss(tags);
	string token;
	while (ss >> token)
	{
		std::string newString = directoryName + "/" + token;
		cout << directoryName + "/" + token << endl;
		this->createDimuonHistograms(newString);
	}
}

void HistogramManager::writeTFile()
{
	//Loop over all the histograms, the data structure is histogram["directoryName"]["histogramName"] = TH1* histogram!
	for (auto dir : histograms)
	{
		for (auto hist : dir.second)
		{
			//cout << "dir.first: " << dir.first << endl;
			rootFile->cd(dir.first.c_str());
			//directories[dir.first]->cd();
			hist.second->Write();
		}
	}
	rootFile->Close();
}

void HistogramManager::printDirectories()
{
	for (auto dir : histograms)
	{
		cout << dir.first << endl;
	}
}

HistogramManager::~HistogramManager()
{
}

