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
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  man->SetVerboseLevel(2);
  #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
  #endif
  man->SetFirstNtupleId(0);
  man->SetFirstNtupleColumnId(0);

  eDepHist1 = man->CreateH1("BF3EnergyDep1", "BF3EnergyDep1", 1000, 0., 5.);
  eDepHist2 = man->CreateH1("BF3EnergyDep2", "BF3EnergyDep2", 1000, 0., 5.);
  eDepHistTot = man->CreateH1("BF3EnergyDepTot", "BF3EnergyDepTot", 1000, 0., 5.);

  primEneHist = man->CreateH1("PrimaryEnergy", "PrimaryEnergy", 500, 0., 50.);
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