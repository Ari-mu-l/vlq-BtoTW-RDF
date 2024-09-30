#define rdf_cxx
#include "../lumiMask.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <TFile.h>
#include <TChain.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>
#include <algorithm> // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TVector2.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock
//#include "../correctionlib/include/correction.h"

using namespace std;
using namespace ROOT::VecOps;

void rdf::analyzer_check2016preVFP(TString testNum)
{
  TStopwatch time;
  time.Start();
  string sample = this->sample;
  string year = this->year;
  //bool isMC = this->isMC;
  //bool isSM = this->isSM;
  //bool isSE = this->isSE;
  //bool debug = false;

  cout << "Sample in cc: " << sample << endl;
  cout << "Year in cc: " << year << endl;

  std::string jsonfile = "../Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
  const auto myLumiMask = lumiMask::fromJSON(jsonfile);
  auto goldenjson = [myLumiMask](unsigned int &run, unsigned int &luminosityBlock){return myLumiMask.accept(run, luminosityBlock);};

  auto rdf_input = ROOT::RDataFrame("Events", files);
  
  auto METgeneralFilters = rdf_input.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters");
  
  auto truth = METgeneralFilters.Define("passesJSON", goldenjson, {"run","luminosityBlock"})
    .Filter("passesJSON == 1", "Data passes Golden JSON");
  
  auto LepDefs = truth.Define("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap, Electron_cutBased, Electron_pfRelIso03_all,Electron_eta,Electron_pt,Electron_sieie,Electron_eInvMinusPInv)")
    .Define("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
    .Define("TPassEl", Form("(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1%s",elHEMcut.c_str()))
    .Define("VetoMu", "TPassMu && (Muon_pt>25)")
    .Define("VetoEl", "TPassEl && (Electron_pt>25)")
    .Define("SignalIsoMu", "TPassMu && (Muon_pt>=55)")
    .Define("SignalIsoEl", "TPassEl && (Electron_pt>=55)")
    .Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
    .Define("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
    .Define("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
    .Define("VetoIsoMu", "(VetoMu == true && Muon_pt < 55)")
    .Define("VetoIsoEl", "(VetoEl == true && Electron_pt < 55)")
    .Define("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))");

  string tkmutrig = " || HLT_OldMu100 || HLT_TkMu100";
  string eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200";
  if(year == "2016APV" and (era == "A" || era == "B")){ tkmutrig = ""; eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175";}
  
  auto LepSelect = LepDefs.Define("isMu", Form("(nMuon>0) && (HLT_Mu50%s) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)",tkmutrig.c_str()))
    .Define("isEl", Form("(nElectron>0) && (%s) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)",eltrig.c_str()));

  cout << "-------------------------------------------------" << endl
       << ">>> Saving " << sample << " Snapshot..." << endl;
  //TString finalFile = "RDF_" + sample + "_" + year + "_" + testNum.Data() + ".root";
  TString finalFile = "RDF_" + sample + era + "_" + year + "_" + testNum.Data() + ".root";
  const char *stdfinalFile = finalFile;
  
  ROOT::RDF::RSnapshotOptions opts;
  Reconstruction.Snapshot("Events", stdfinalFile, {"Electron_r9", });
  cout << "Output File: " << finalFile << endl
       << "-------------------------------------------------" << endl;
  
  time.Stop();
  time.Print();
  
  cout << "Cut statistics:" << endl;
  Reconstruction.Report()->Print();

  std::cout << "Got past the report" << std::endl;
}
