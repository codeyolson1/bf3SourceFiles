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

G4ThreadLocal Analysis* theAnalysis = 0;

Analysis::Analysis()
{
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

void Analysis::Book(G4bool isHe3)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  man->SetVerboseLevel(2);
  #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
  #endif
  man->SetFirstNtupleId(0);
  man->SetFirstNtupleColumnId(0);
  if (isHe3) {
    man->CreateNtuple("Helium3", "Helium3");
    man->CreateNtupleDColumn("EDep");
    man->FinishNtuple();
  } else {
    man->CreateNtuple("Boron Trifluoride", "Boron Trifluoride");
    man->CreateNtupleDColumn("EDep");
    man->FinishNtuple();
  }
  man->CreateNtuple("Primary", "Primary");
  man->CreateNtupleDColumn("Energy");
  man->CreateNtupleDColumn("XPos");
  man->CreateNtupleDColumn("YPos");
  man->FinishNtuple();
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

void Analysis::FillEDep(G4double eDep, G4int col)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(0, col, eDep);
  man->AddNtupleRow(0);
  return;
}

//
//

void Analysis::FillPrimary(G4double energy)
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(1, 0, energy);
}

//
//

void Analysis::FillPrimaryPos(G4double xPos, G4double yPos)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(1, 1, xPos);
  man->FillNtupleDColumn(1, 2, yPos);
  man->AddNtupleRow(1);
}