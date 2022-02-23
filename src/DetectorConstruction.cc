// Definition of world geometry and detectors.
// Created by Codey Olson on May 10, 2021.

/// \file DetectorConstruction.cc
/// \brief Definition of world geometry and detectors.

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSCellFlux.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#define _USE_MATH_DEFINES 
#include <math.h>
#include <iomanip>
#include <iostream>
#include <string>


DetectorConstruction::DetectorConstruction(G4bool he3)
: G4VUserDetectorConstruction()
{
  fmats = {};
  isHe3 = he3;
  ConstructMaterials();
}

//
//

DetectorConstruction::~DetectorConstruction()
{}

//
//

void DetectorConstruction::ConstructMaterials()
{
  // Get instance of nist material manager:
  G4NistManager* nist = G4NistManager::Instance();

  // Create materials and input into dictionary.
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  fmats["air"] = air;
  G4Material* poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  fmats["poly"] = poly;
  // Fill gases:
  G4Material* bf3 = new G4Material("Boron Trifluoride", 2.73e-3*g/cm3, 2, kStateGas, 297.039*kelvin, 1.*atmosphere); // From Walker Dissertation
  G4Element* boron = nist->FindOrBuildElement(5, true);
  G4Element* fluorine = nist->FindOrBuildElement(9, true);
  G4Element* hydrogen = nist->FindOrBuildElement(1, true);
  bf3->AddElement(boron, 1);
  bf3->AddElement(fluorine, 3);
  fmats["bf3"] = bf3;

  G4Material* bf3En = new G4Material("Boron Trifluoride (Enriched)", 2.73e-3*g/cm3, 2, kStateGas, 297.039*kelvin, 1.*atmosphere); // From Walker Dissertai
  G4Element* enrBoron = new G4Element("Enriched Boron", "B", 2);
  G4Isotope* boron10 = new G4Isotope("Boron10", 5, 10, 10.012936862*g/mole); // 
  G4Isotope* boron11 = new G4Isotope("Boron11", 5, 11, 11.009305167*g/mole); //
  enrBoron->AddIsotope(boron10, 96.*perCent);
  enrBoron->AddIsotope(boron11, 4.*perCent);
  bf3En->AddElement(enrBoron, 1);
  bf3En->AddElement(fluorine, 3);
  fmats["enrBF3"] = bf3En;


  G4Material* he3 = new G4Material("Helium 3", 5.39e-4*g/cm3, 1, kStateGas, 297.039*kelvin, 4.*atmosphere); // From Walker Dissertai
  G4Element* helium = new G4Element("Helium", "He", 1);
  G4Isotope* helium3 = new G4Isotope("Helium3", 2, 3, 3.01602932197*g/mole); // from IAEA
  helium->AddIsotope(helium3, 100.*perCent);
  he3->AddElement(helium, 1);
  fmats["he3"] = he3;

  G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fmats["steel"] = steel;

  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  fmats["aluminum"] = aluminum;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;

  //
  // World:
  // Params:
  G4double worldX, worldY, worldZ;
  if (isHe3) {
    worldX = 11.54*cm; 
    worldY = 9.54*cm; 
    worldZ = 11.*cm;
  } else {
    worldX = 18.3*cm; 
    worldY = 11.4*cm; 
    worldZ = 11.*cm;
  }
  // Construction:
  G4Box* solidWorld = new G4Box("World", 0.5*worldX, 0.5*worldY,0.5*worldZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fmats["air"], "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

  // 
  // Specify here if the detector will be he3
  G4double tubeDiam;
  G4double tubeHeight;
  G4double modx, mody, modz;

  // Tube and moderator dimensions:
  if (isHe3) {
    tubeDiam = 2.54*cm;
    tubeHeight = 10*cm;
    modx = 6.54*cm; mody = 4.54*cm; modz = tubeHeight;
  } else {
    tubeDiam = 4.4*cm;
    tubeHeight = 10*cm;
    modx = tubeDiam*2. + 4.5*cm; mody = tubeDiam + 2.*cm; modz = tubeHeight;
  }

  if (isHe3) {
    // Construct He3 Detector:
    // Detector:
    // SS Shell:
    G4Tubs* ssShellSolid = new G4Tubs("SS Shell", 0, 0.5*tubeDiam, 0.5*tubeHeight, 0, 360.*deg);
    G4LogicalVolume* ssShellLogic = new G4LogicalVolume(ssShellSolid, fmats["aluminum"], "SS Shell");
    new G4PVPlacement(0, G4ThreeVector(), ssShellLogic, "SS Shell", logicWorld, false, 0, checkOverlaps); 
    // helium3 fill gas:
    G4Tubs* he3GasSolid = new G4Tubs("He3 Gas", 0, 0.5*(tubeDiam - 2*mm), 0.5*(tubeHeight - 2*mm), 0, 360.*deg);
    G4LogicalVolume* he3GasLogic = new G4LogicalVolume(he3GasSolid, fmats["he3"], "He3 Gas");
    new G4PVPlacement(0, G4ThreeVector(), he3GasLogic, "He3 Tube", ssShellLogic, false, 0, checkOverlaps); 
    //Moderator:
    // Dummies for subtraction solid:
    G4Box* moderatorDummy = new G4Box("He3 Moderator Dummy", 0.5*modx, 0.5*mody, 0.5*modz);
    G4Tubs* moderatorVoidDummy = new G4Tubs("He3 Void Dummy", 0, 0.5*tubeDiam, 0.5*(tubeHeight + 1*cm), 0, 360.*deg);
    // Final solid:
    G4VSolid* he3ModeratorSolid = new G4SubtractionSolid("He3 Moderator", moderatorDummy, moderatorVoidDummy, 0, G4ThreeVector());
    G4LogicalVolume* he3ModeratorLogic = new G4LogicalVolume(he3ModeratorSolid, fmats["poly"], "He3 Moderator");
    new G4PVPlacement(0, G4ThreeVector(), he3ModeratorLogic, "He3 Moderator", logicWorld, false, 0, checkOverlaps);
  } else {
    // Construct BF3 Detectors:
    // SS Shells
    G4Tubs* bf3ShellSolid1 = new G4Tubs("BF3 Shell1", 0, 0.5*tubeDiam, 0.5*tubeHeight, 0, 360.*deg);
    G4LogicalVolume* bf3ShellLogic1 = new G4LogicalVolume(bf3ShellSolid1, fmats["aluminum"], "BF3 Shell1");
    new G4PVPlacement(0, G4ThreeVector(tubeDiam*0.5 + 0.5*cm, 0, 0), bf3ShellLogic1, "BF3 Shell1", logicWorld, false, 0, checkOverlaps);
    G4Tubs* bf3ShellSolid2 = new G4Tubs("BF3 Shell2", 0, 0.5*tubeDiam, 0.5*tubeHeight, 0, 360.*deg);
    G4LogicalVolume* bf3ShellLogic2 = new G4LogicalVolume(bf3ShellSolid2, fmats["aluminum"], "BF3 Shell2");
    new G4PVPlacement(0, G4ThreeVector(-tubeDiam*0.5 - 0.5*cm, 0, 0), bf3ShellLogic2, "BF3 Shell2", logicWorld, false, 0, checkOverlaps);
    // BF3 fill gas:
    G4Tubs* bf3GasSolid1 = new G4Tubs("BF3 Gas1", 0, 0.5*(tubeDiam - 2.*mm), 0.5*(tubeHeight - 2.*mm), 0, 360.*deg);
    G4LogicalVolume* bf3GasLogic1 = new G4LogicalVolume(bf3GasSolid1, fmats["enrBF3"], "BF3 Gas1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), bf3GasLogic1, "BF3 Gas1", bf3ShellLogic1, false, 0, checkOverlaps);
    G4Tubs* bf3GasSolid2 = new G4Tubs("BF3 Gas2", 0, 0.5*(tubeDiam - 2.*mm), 0.5*(tubeHeight - 2.*mm), 0, 360.*deg);
    G4LogicalVolume* bf3GasLogic2 = new G4LogicalVolume(bf3GasSolid2, fmats["enrBF3"], "BF3 Gas2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), bf3GasLogic2, "BF3 Gas2", bf3ShellLogic2, false, 0, checkOverlaps);
    // Moderator:
    G4Box* moderatorDummy1 = new G4Box("BF3 Moderator Dummy", 0.5*modx, 0.5*mody, 0.5*modz);
    G4Tubs* moderatorVoidDummy1 = new G4Tubs("BF3 Moderator Void Dummy", 0, 0.5*tubeDiam, 0.5*(tubeHeight + 1.*cm), 0, 360.*deg);
    G4VSolid* bf3ModeratorTemp = new G4SubtractionSolid("Mod Temp", moderatorDummy1, moderatorVoidDummy1, 0, G4ThreeVector(tubeDiam*0.5 + 0.5*cm, 0, 0));
    G4VSolid* bf3ModeratorSolid = new G4SubtractionSolid("BF3 Moderator", bf3ModeratorTemp, moderatorVoidDummy1, 0, G4ThreeVector(-tubeDiam*0.5 - 0.5*cm, 0, 0));
    G4LogicalVolume* moderatorBF3Logic = new G4LogicalVolume(bf3ModeratorSolid, fmats["poly"], "ModeratorBF3");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), moderatorBF3Logic, "ModeratorBF3", logicWorld, false, 0, checkOverlaps);

  }

  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  if (isHe3) {
    G4MultiFunctionalDetector* he3Detector = new G4MultiFunctionalDetector("Helium-3");
    G4SDManager::GetSDMpointer()->AddNewDetector(he3Detector);
    G4VPrimitiveScorer* energyDep = new G4PSEnergyDeposit("EnergyDep");
    he3Detector->RegisterPrimitive(energyDep);
    SetSensitiveDetector("He3 Gas", he3Detector);
  } else {
    G4MultiFunctionalDetector* bf3Detector1 = new G4MultiFunctionalDetector("BF31");
    G4SDManager::GetSDMpointer()->AddNewDetector(bf3Detector1);
    G4VPrimitiveScorer* energyDep1 = new G4PSEnergyDeposit("EnergyDep1");
    bf3Detector1->RegisterPrimitive(energyDep1);
    SetSensitiveDetector("BF3 Gas1", bf3Detector1);

    G4MultiFunctionalDetector* bf3Detector2 = new G4MultiFunctionalDetector("BF32");
    G4SDManager::GetSDMpointer()->AddNewDetector(bf3Detector2);
    G4VPrimitiveScorer* energyDep2 = new G4PSEnergyDeposit("EnergyDep2");
    bf3Detector2->RegisterPrimitive(energyDep2);
    SetSensitiveDetector("BF3 Gas2", bf3Detector2);
  }
}