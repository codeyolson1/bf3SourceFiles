// Created by Codey Olson - 7/28/2021
// Copied from example Hadr03

/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4HadronElasticPhysicsPHP.hh"
#include "G4HadronElasticPhysicsXS.hh"

#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsShieldingLEND.hh"
#include "G4HadronElasticPhysicsLEND.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"

#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
//#include "GammaNuclearPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPCapture.hh"
#include "G4ThermalNeutrons.hh"
// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "NeutronHPphysics.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{
  G4int verb = 1;  
  SetVerboseLevel(verb);
  SetCutValue(0, "proton");
  
  //add new units for cross sections
  // 
  new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);  
  
  // Hadron Elastic scattering
  //
  //
  RegisterPhysics( new G4HadronElasticPhysicsLEND(verb));
  //RegisterPhysics( new NeutronHPphysics("neutronHP"));
  
  // Hadron Inelastic physics
  //
  ////RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
  //RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
  //RegisterPhysics( new G4HadronPhysicsQGSP_BIC_AllHP(verb));
  //RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC(verb));  
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));
  ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
  //RegisterPhysics( new G4HadronPhysicsShieldingLEND(verb));
  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_AllHP());
  RegisterPhysics( new G4ThermalNeutrons(verb));  
  // Electromagnetic physics:
  RegisterPhysics( new G4EmStandardPhysics_option4(verb));
  RegisterPhysics( new G4EmExtraPhysics(verb));
  
  // Ion Elastic scattering
  //
  RegisterPhysics( new G4IonElasticPhysics(verb));
  
  // Ion Inelastic physics
  //
  ///RegisterPhysics( new G4IonPhysicsXS(verb));
  RegisterPhysics( new G4IonPhysicsPHP(verb));
  RegisterPhysics( new G4StoppingPhysics(verb));
  ////RegisterPhysics( new G4IonQMDPhysics(verb));
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));

  // Gamma physics
  //
  //RegisterPhysics( new GammaNuclearPhysics("gamma"));
  
  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());
  RegisterPhysics(new G4DecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
   //SetCutValue(0.1*mm, "proton");
   //SetCutValue(0.1*mm, "alpha");
   //SetCutValue(0.1*mm, "triton");
   //SetCutValue(0.1*mm, "deuteron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
