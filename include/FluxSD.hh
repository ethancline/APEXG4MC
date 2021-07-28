#ifndef FluxSD_h
#define FluxSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4VHitsCollection.hh"
#include "FluxHit.hh"

//---------------------------------------------------------------------------

class FluxSD : public G4VSensitiveDetector
{
public:
  
  FluxSD( G4String, G4int );
  ~FluxSD();
  
  void   Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  G4bool ProcessHits_constStep(const G4Step* astep,G4TouchableHistory* ROHist);
  void   EndOfEvent(G4HCofThisEvent*);
  
private:
  
  FluxHitsCollection*  fCollection;  
  G4int                    fDetID;
  G4int                    fNelements;
  G4int                    fNhits;
  G4int*                   fhitID;
  G4int*                   fHits;

};
#endif

//---------------------------------------------------------------------------






