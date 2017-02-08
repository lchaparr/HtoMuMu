/*
 * HistogramManager.h
 *
 *  Created on: Feb 6, 2017
 *      Author: gastonlp
 */

#ifndef HISTOGRAMMANAGER_H_
#define HISTOGRAMMANAGER_H_


#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TDirectory.h>


class HistogramManager
    {
public:
    // attributes
    TFile* rootFile;
    std::map<std::string, TDirectory* > directories;
    std::map<std::string, std::map<std::string, TH1*>> histograms;

    //methods
    HistogramManager(std::string filename);
    void addDirectory(std::string);
    void writeTFile();
    void printDirectories();
    void addHistogram(std::string, std::string histogramName, std::string histogramTitle, Int_t nbinsx, Double_t xlow, Double_t xup);
    void fillHistogram(std::string, std::string histogramName,Double_t x);
    void fillHistogram(std::string, std::string histogramName,Double_t x, Double_t w);
    void createMuonHistograms(std::string directoryName);
    void createPhotonHistograms(std::string directoryName);
    void createDimuonHistograms(std::string directoryName);
    void fillMuonHist(std::string directoryName, TLorentzVector& Lepton1, TLorentzVector& Lepton2, Double_t weight);
    void fillHiggsHist(std::string directoryName, TLorentzVector& Dimuom, Double_t weight);
    void fillHiggsHist(std::string directoryName, TLorentzVector& Dimuom, TLorentzVector& PhotonFSR, Double_t weight);
    void fillHiggsHist(std::string directoryName, TLorentzVector& Dimuon, TLorentzVector& PhotonFSR1, TLorentzVector& PhotonFSR2, Double_t weight);
    void fillPhotonHist(std::string directoryName,TLorentzVector& FSR_Photon1_TL, TLorentzVector& FSR_Photon2_TL, Double_t weight);
    void fillPhotonHist(std::string directoryName, TLorentzVector& FSR_Photon1_TL, Double_t weight);
    void createMuonHistograms(std::string directoryName, std::string tags);
    void createPhotonHistograms(std::string directoryName, std::string tags);
    void createDimuonHistograms(std::string directoryName, std::string tags);

    virtual ~HistogramManager();
private:

    };

#endif /* HISTOGRAMMANAGER_H_ */
