#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4THitsCollection.hh"

#include "FluxSD.hh"
#include "FluxHit.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

FluxSD::FluxSD(G4String name, G4int)
  :G4VSensitiveDetector(name)
{
  collectionName.insert(G4String("SDHits")+name);
  fCollection = NULL;
  fNelements  = 10000;
  fNhits      = 0;
  fhitID      = new G4int[fNelements];
  fHits       = new G4int[fNelements];
  for(G4int i=0; i<fNelements; i++) fhitID[i] = -1;
  for(G4int i=0; i<fNelements; i++) fHits[i]  = 0;
}

//---------------------------------------------------------------------------

FluxSD::~FluxSD()
{
}

//---------------------------------------------------------------------------

void FluxSD::Initialize(G4HCofThisEvent* HCE)
{

  fCollection = new FluxHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;

  if( HCID < 0 )  {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  }  
  HCE->AddHitsCollection( HCID, (G4VHitsCollection*)fCollection ); 
}

//---------------------------------------------------------------------------

G4bool FluxSD::ProcessHits( G4Step* aStep,G4TouchableHistory* )
{ 
  return false;
}

//---------------------------------------------------------------------------

G4bool FluxSD::ProcessHits_constStep(const G4Step* aStep,G4TouchableHistory*)
{

  G4Track*              aTrack       = aStep->GetTrack();
  G4TouchableHistory*   theTouchable = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
  G4VPhysicalVolume*    volume       = theTouchable->GetVolume();
  G4int                 id           = volume->GetCopyNo();
  G4String              ParticleName = aTrack->GetDefinition()->GetParticleName();
  G4double              ene          = aTrack->GetKineticEnergy();
  G4int                 AtomicMass   = aTrack->GetDefinition()->GetAtomicMass(); 
  G4int                 AtomicNumber = aTrack->GetDefinition()->GetAtomicNumber();  
  G4int                 pid          = aTrack->GetParentID(); 
  G4int                 tid          = aTrack->GetTrackID(); 

  FluxHit* Hit = new FluxHit;
  Hit->SetEnergy(ene);
  Hit->SetMomentum(aTrack->GetMomentum());
  Hit->SetPrePosition(aStep->GetPreStepPoint()->GetPosition());
  Hit->SetPostPosition(aStep->GetPostStepPoint()->GetPosition());
  Hit->SetPDef(aTrack->GetDefinition());
  Hit->SetID(id);
  Hit->SetTrackID(tid);
  Hit->SetParentID(pid);
  Hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
  Hit->SetAtomicMass(AtomicMass);
  Hit->SetAtomicNumber(AtomicNumber);
  fhitID[id] = fCollection->insert(Hit) -1;
  fHits[fNhits++]=id;
  return true;

}

//---------------------------------------------------------------------------

void FluxSD::EndOfEvent(G4HCofThisEvent*)
{

  for (G4int i=0;i<fNhits;i++) {
    fhitID[fHits[i]] = -1;
    fHits[i]         = 0;
  }
  fNhits = 0;
}

//---------------------------------------------------------------------------














