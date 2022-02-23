// Source file for RunAction() class.
// Created by Codey Olson on May 28, 2021.
/// \file RunAction.cc
/// \brief Source code for class RunAction.
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "Analysis.hh"
#include <G4WorkerThread.hh>
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4StatAnalysis.hh"
#include "G4RunManager.hh"
#include "G4ConvergenceTester.hh"
#include "tools/mpi/wrmpi"
#include <iostream>
#include <fstream>
RunAction::RunAction(G4bool He3) : G4UserRunAction()
{
  isHe3 = He3;
G4String primEnergy = "grp28";
  if (He3) {
    outFileName = "EDep_He3_" + primEnergy;
  } else {
    outFileName = "EDep_BF3_" + primEnergy;
  }
  eDepTest = new G4ConvergenceTester("EnergyDeposition");
}
//
//
RunAction::~RunAction()
{}
//
//
G4Run* RunAction::GenerateRun()
{
  return new Run(isHe3);
}
//
//
void RunAction::BeginOfRunAction(const G4Run*)
{
  Analysis* myAnalysis = Analysis::GetAnalysis();
  myAnalysis->Book(isHe3);
  myAnalysis->OpenFile(outFileName);
  if (IsMaster()) {
    eDepTest = new G4ConvergenceTester("EnergyDeposition");
  }
} 
//
//
void RunAction::EndOfRunAction(const G4Run* aRun)
{
  Analysis* myana = Analysis::GetAnalysis();
  if (IsMaster()) {
    G4cout << "End of Global Run" << G4endl;
    myana->Save();
    myana->Close();
    const Run* myRun = static_cast<const Run*>(aRun);
    for (G4int i = 0; i != myRun->GetNumEntries(); i++) {
      G4double eDep = myRun->GetEDepAtEvent(i);
      eDepTest->AddScore(eDep);
    }
    std::ofstream convergence;
    convergence.open("convergence.txt");
    eDepTest->ShowResult(convergence);
    eDepTest->ShowHistory(convergence);
    convergence.close();
  } else {
    myana->Save();
  }
}