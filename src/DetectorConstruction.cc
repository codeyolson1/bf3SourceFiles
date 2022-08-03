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


DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  fmats = {};
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
  G4Material* bf3 = new G4Material("Boron Trifluoride", 2.73e-3*g/cm3, 2, kStateGas, 293.*kelvin, 1.*atmosphere); // From Walker Dissertation
  G4Element* boron = nist->FindOrBuildElement(5, true);
  G4Element* fluorine = nist->FindOrBuildElement(9, true);
  bf3->AddElement(boron, 1);
  bf3->AddElement(fluorine, 3);
  fmats["bf3"] = bf3;

  G4Material* bf3En = new G4Material("Boron Trifluoride (Enriched)", 2.73e-3*g/cm3, 2, kStateGas, 293.*kelvin, 1.*atmosphere); // From Walker Dissertai
  G4Element* enrBoron = new G4Element("Enriched Boron", "B", 2);
  G4Isotope* boron10 = new G4Isotope("Boron10", 5, 10, 10.012936862*g/mole); // 
  G4Isotope* boron11 = new G4Isotope("Boron11", 5, 11, 11.009305167*g/mole); //
  enrBoron->AddIsotope(boron10, 99.*perCent);
  enrBoron->AddIsotope(boron11, 1.*perCent);
  bf3En->AddElement(enrBoron, 1);
  bf3En->AddElement(fluorine, 3);
  //fmats["enrBF3"] = bf3En;


// Material info from :
// https://gitlab.cern.ch/clemenci/Geant4-srcs/-/blob/92686251452762ac5947193b5f02ba43b77f546b/examples/extended/hadronic/FissionFragment/src/FFDetectorConstruction.cc
    G4double const B10Enrichment = 0.96;
    G4double const B11Enrichment = 0.04;
    G4Isotope* const iB10
        = new G4Isotope("iB10",                         // name
                        5,                              // ZZZ
                        10,                             // AAA
                        10.0129370 * (g / mole));       // molecular weight
    G4Isotope* const iB11
        = new G4Isotope("iB11",                         // name
                        5,                              // ZZZ
                        11,                             // AAA
                        11.0093054 * (g / mole));       // molecular weight
    // Now create the elements and add the isotopes
    G4Element* const B10
        = new G4Element("B10",                          // name
                        "B10",                          // symbol
                        1);                             // number of isotopes
    B10->AddIsotope(iB10,                               // isotope
                     1.0);                              // abundance
    G4Element* const B11
        = new G4Element("B11",                          // name
                        "B11",                          // symbol
                        1);                             // number of isotopes
    B11->AddIsotope(iB11,                               // isotope
                     1.0);                              // abundance
    G4Element* const flouride = nist->FindOrBuildElement("F");
    // Calculate the mass fractions
    const G4double BF3MolecularWeight = B10->GetA() * B10Enrichment
                                        + B11->GetA() * B11Enrichment
                                        + flouride->GetA() * 3;
    const G4double B10MassFraction = (B10->GetA() * B10Enrichment)
                                     / BF3MolecularWeight;
    const G4double B11MassFraction = (B11->GetA() * B11Enrichment)
                                     / BF3MolecularWeight;
    const G4double flourideMassFraction = (flouride->GetA() * 3)
                                          / BF3MolecularWeight;
    // create the material and add the elements
    fmats["enrBF3"] = new G4Material("BF3_96E",                // name
                              2.73E-3 * (g/cm3),          // density
                              3, kStateGas, 293.*kelvin, 1.0*atmosphere);                       // number of components
    fmats["enrBF3"]->AddElement(B10,                           // element
                         B10MassFraction);              // mass fraction
    fmats["enrBF3"]->AddElement(B11,                           // element
                         B11MassFraction);              // mass fraction
    fmats["enrBF3"]->AddElement(flouride,                      // element
                         flourideMassFraction);         // mass fraction


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

  worldX = 18.3*cm; 
  worldY = 11.4*cm; 
  worldZ = 11.*cm;
  
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
  tubeDiam = 4.4*cm;
  tubeHeight = 10.0*cm;
  modx = tubeDiam*2. + 4.5*cm; mody = tubeDiam + 2.*cm; modz = tubeHeight;

  // Construct BF3 Detectors:
  // SS Shells
  G4Tubs* bf3ShellSolid1 = new G4Tubs("BF3 Shell1", 0, 0.5*(tubeDiam + 0.2*cm), 0.5*(tubeHeight + 0.2*cm), 0., 360.*deg);
  G4LogicalVolume* bf3ShellLogic1 = new G4LogicalVolume(bf3ShellSolid1, fmats["steel"], "BF3 Shell1");
  //new G4PVPlacement(0, G4ThreeVector(tubeDiam*0.5 + 0.5*cm, 0, 0), bf3ShellLogic1, "BF3 Shell1", logicWorld, false, 0, checkOverlaps);
  G4Tubs* bf3ShellSolid2 = new G4Tubs("BF3 Shell2", 0, 0.5*(tubeDiam + 0.2*cm), 0.5*(tubeHeight + 0.2*cm), 0., 360.*deg);
  G4LogicalVolume* bf3ShellLogic2 = new G4LogicalVolume(bf3ShellSolid2, fmats["steel"], "BF3 Shell2");
  //new G4PVPlacement(0, G4ThreeVector(-tubeDiam*0.5 - 0.5*cm, 0, 0), bf3ShellLogic2, "BF3 Shell2", logicWorld, false, 0, checkOverlaps);
  // Visual Stuff for shells
  G4VisAttributes* shellAttr = new G4VisAttributes(G4Colour(192., 192., 192.)); // silver
  shellAttr->SetForceWireframe(true);
  bf3ShellLogic1->SetVisAttributes(shellAttr);
  bf3ShellLogic2->SetVisAttributes(shellAttr);
  // BF3 fill gas:
  G4Tubs* bf3GasSolid1 = new G4Tubs("BF3 Gas1", 0, 0.5*(tubeDiam), 0.5*(tubeHeight), 0, 360.*deg);
  G4LogicalVolume* bf3GasLogic1 = new G4LogicalVolume(bf3GasSolid1, fmats["enrBF3"], "BF3 Gas1");
  new G4PVPlacement(0, G4ThreeVector(tubeDiam*0.5 + 0.25*cm, 0, 0), bf3GasLogic1, "BF3 Gas1", logicWorld, false, 0, checkOverlaps);
  G4Tubs* bf3GasSolid2 = new G4Tubs("BF3 Gas2", 0, 0.5*(tubeDiam), 0.5*(tubeHeight), 0, 360.*deg);
  G4LogicalVolume* bf3GasLogic2 = new G4LogicalVolume(bf3GasSolid2, fmats["enrBF3"], "BF3 Gas2");
  new G4PVPlacement(0, G4ThreeVector(-(tubeDiam)*0.5 - 0.25*cm, 0, 0), bf3GasLogic2, "BF3 Gas2", logicWorld, false, 0, checkOverlaps);
  G4cout << "BF3 gas volume: " << bf3GasSolid1->GetCubicVolume()/cm3 + bf3GasSolid2->GetCubicVolume()/cm3 << G4endl;
  // Visual Stuff for gas
  G4VisAttributes* gasAttr = new G4VisAttributes(G4Colour(255., 0., 0.)); // red
  gasAttr->SetForceSolid(true);
  bf3GasLogic1->SetVisAttributes(gasAttr);
  bf3GasLogic2->SetVisAttributes(gasAttr);

  // Moderator:
  G4Box* moderatorDummy1 = new G4Box("BF3 Moderator Dummy", 0.5*modx, 0.5*mody, 0.5*modz);
  G4Tubs* moderatorVoidDummy1 = new G4Tubs("BF3 Moderator Void Dummy", 0, 0.5*(tubeDiam), 0.5*(tubeHeight + 1.*cm), 0, 360.*deg);
  G4VSolid* bf3ModeratorTemp = new G4SubtractionSolid("Mod Temp", moderatorDummy1, moderatorVoidDummy1, 0, G4ThreeVector((tubeDiam)*0.5 + 0.25*cm, 0, 0));
  G4VSolid* bf3ModeratorSolid = new G4SubtractionSolid("BF3 Moderator", bf3ModeratorTemp, moderatorVoidDummy1, 0, G4ThreeVector(-(tubeDiam)*0.5 - 0.25*cm, 0, 0));
  G4LogicalVolume* moderatorBF3Logic = new G4LogicalVolume(bf3ModeratorSolid, fmats["poly"], "ModeratorBF3");
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), moderatorBF3Logic, "ModeratorBF3", logicWorld, false, 0, checkOverlaps);
  G4cout << "Moderator volume: " << bf3ModeratorSolid->GetCubicVolume()/cm3 << G4endl;
  // Visual Stuff for moderator
  G4VisAttributes* moderatorAttr = new G4VisAttributes(G4Colour()); // white
  moderatorAttr->SetForceSolid(true);
  moderatorBF3Logic->SetVisAttributes(moderatorAttr);

  // Air Source
  G4Box* airSourceDummy = new G4Box("AirSourceDummy", (modx + 4.*cm)*0.5, (mody + 4.*cm)*0.5, 0.5*(modz));
  G4Box* moderatorDummy2 = new G4Box("ModeratorDummy2", 0.5*modx, 0.5*mody, 0.5*(modz + 0.5*cm));
  G4VSolid* airSource = new G4SubtractionSolid("AirSource", airSourceDummy, moderatorDummy2, 0, G4ThreeVector(0, 0, 0));
  G4LogicalVolume* airLogic = new G4LogicalVolume(airSource, fmats["air"], "AirSource");
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), airLogic, "AirSource", logicWorld, false, 0, checkOverlaps);
  G4cout << "Air source volume: " << airSource->GetCubicVolume()/cm3 << G4endl;
  // visual Stuff for Air source
  G4VisAttributes* airAttr = new G4VisAttributes(G4Colour(0., 255., 0.)); // green
  airAttr->SetForceSolid(true);
  airLogic->SetVisAttributes(airAttr);
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDParticleFilter* nFilter = new G4SDParticleFilter("NeutronFilter");
  nFilter->add("alpha");
  nFilter->addIon(3,7); // Li7
  nFilter->addIon(3,6); // Li6
  nFilter->add("proton");
  nFilter->addIon(4,10); // Be10
  nFilter->add("deuteron");
  nFilter->addIon(4,9); // Be9
  nFilter->add("triton");
  nFilter->addIon(4, 8); // Be8
  nFilter->addIon(5,10); // B-10
  nFilter->addIon(5,11); // B-11
  nFilter->add("neutron");

  
  G4MultiFunctionalDetector* bf3Detector1 = new G4MultiFunctionalDetector("BF31");
  G4SDManager::GetSDMpointer()->AddNewDetector(bf3Detector1);
  G4VPrimitiveScorer* energyDep1 = new G4PSEnergyDeposit("EnergyDep1");
  bf3Detector1->RegisterPrimitive(energyDep1);
  energyDep1->SetFilter(nFilter);
  SetSensitiveDetector("BF3 Gas1", bf3Detector1);

  G4MultiFunctionalDetector* bf3Detector2 = new G4MultiFunctionalDetector("BF32");
  G4SDManager::GetSDMpointer()->AddNewDetector(bf3Detector2);
  G4VPrimitiveScorer* energyDep2 = new G4PSEnergyDeposit("EnergyDep2");
  bf3Detector2->RegisterPrimitive(energyDep2);
  energyDep2->SetFilter(nFilter);
  SetSensitiveDetector("BF3 Gas2", bf3Detector2);
  
}