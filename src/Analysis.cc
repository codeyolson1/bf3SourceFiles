// Source file for Analysis().
// Created by Codey Olson on June 1, 2021.

/// \file Analysis.cc
/// \brief Source code for Analysis class.

#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"
#include "Analysis.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SDManager.hh"
#include "g4root.hh"
#include "G4RootAnalysisManager.hh"
#include "G4ConvergenceTester.hh"

G4ThreadLocal Analysis* theAnalysis = 0;

namespace {
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
  G4ConvergenceTester* fConvTest1 = new G4ConvergenceTester("ConvTest1");
}

Analysis::Analysis()
{
  eDepHist1 = 0;
  eDepHist2 = 0;
  eDepHistTot = 0;
  primEneHist = 0;
  primPosHist = 0;
  convergenceName = "";
}

//
//

Analysis::~Analysis() 
{
}

//
//

Analysis* Analysis::GetAnalysis()
{
  if (!theAnalysis) {
    theAnalysis = new Analysis();
    G4AutoDelete::Register(theAnalysis);
  }
  return theAnalysis;
}

//
//

void Analysis::Book(G4String runName)
{
  convergenceName = runName;
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(2);
  #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
  #endif
  man->SetFirstNtupleId(0);
  man->SetFirstNtupleColumnId(0);

  eDepHist1 = man->CreateH1("BF3EnergyDep1", "BF3EnergyDep1", 512, 0., 5.);
  eDepHist2 = man->CreateH1("BF3EnergyDep2", "BF3EnergyDep2", 512, 0., 5.);
  eDepHistTot = man->CreateH1("BF3EnergyDepTot", "BF3EnergyDepTot", 512, 0., 5.);

  const std::vector<G4double> binEdges = {2.00E-07,7.28E-07,1.34E-06,2.84E-06,8.22E-06,1.64E-05,6.39E-05,1.65E-04,3.27E-04,6.94E-04,2.71E-03,5.12E-03,1.09E-02,2.29E-02,2.88E-02,2.65E-02,2.80E-02,3.75E-02,5.00E-02,9.39E-02,1.55E-01,2.55E-01,4.11E-01,4.41E-01,6.27E-01,7.18E-01,8.78E-01,9.03E-01,1.18E+00,1.70E+00,1.95E+00,2.19E+00,2.54E+00,2.47E+00,2.39E+00,2.57E+00,2.99E+00,3.29E+00,4.35E+00,6.26E+00,7.17E+00,8.75E+00,9.81E+00,1.14E+01,1.44E+01,1.62E+01,2.04E+01};

  primEneHist = man->CreateH1("PrimEnergy", "PrimaryEnergy", 50, 0., 20.);
  primPosHist = man->CreateH2("PrimaryPosition", "PrimaryPosition", 180, -9., 9., 110, -5.5, 5.5);
  
  return; 
}

//
//

void Analysis::OpenFile(const G4String& fileName)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->OpenFile(fileName.c_str());

  return;
}

//
//

void Analysis::Save()
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();

  return;
}

//
//

void Analysis::Close(G4bool reset)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->CloseFile(reset);

  return;
}

//
//

void Analysis::FillEDep1(G4double eDep)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(eDepHist1, eDep);
  G4AutoLock l(&aMutex);
  fConvTest1->AddScore(eDep);
  return;
}

//
//

void Analysis::FillEDep2(G4double eDep)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(eDepHist2, eDep);
  return;
}

//
//

void Analysis::FillEDepTot(G4double eDep)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(eDepHistTot, eDep);
  return;
}

//
//

void Analysis::FillPrimaryEne(G4double energy)
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(primEneHist, energy);
  return;
}

//
//

void Analysis::FillPrimaryPos(G4double xPos, G4double yPos)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH2(primPosHist, xPos, yPos);
  return;
}

//
//

void Analysis::CheckConvergence()
{
  std::ofstream convOutput;
  convOutput.open(convergenceName+"-conv.txt");
  fConvTest1->ShowResult(convOutput);
  fConvTest1->ShowHistory(convOutput);
  convOutput.close();

  return;
}