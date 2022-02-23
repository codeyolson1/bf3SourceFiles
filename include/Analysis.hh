// Header file for analysis class.
// Created by Codey Olson on May 28, 2021.

/// \file Analysis.hh
/// \brief Selection of analysis classes

#ifndef Analysis_h
#define Analysis_h 1

#include <tools/histo/h1d>
#include <tools/histo/h2d>


#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class Analysis {
  public:
    ~Analysis();

    static Analysis* GetAnalysis();

    void Book(G4bool);
    void EndOfRun();

    void OpenFile(const G4String& fname);
    void Save();
    void Close(G4bool reset = true);

    void FillDose(const G4String& detName, const G4String& collName, G4double dose);
    void FillEDep(G4double eDep, G4int col);
    G4int FindNtupleID(const G4String& detName);
    G4int FindNtupleColID(const G4String& collName);
    void FillPrimary(G4double);
    void FillPrimaryPos(G4double, G4double);

  private:
    Analysis();
    DISALLOW_COPY_AND_ASSIGN(Analysis);

    std::vector<G4String> fdetNames;
    std::map<G4String, G4int> fntupleIDs;
    std::map<G4String, G4int> fcolIDs;
    G4int nPrims;
};

#endif