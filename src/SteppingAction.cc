#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "FluxSD.hh"

#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

SteppingAction::SteppingAction(DetectorConstruction* det)
  :detector(det)
{;}

//---------------------------------------------------------------------------

SteppingAction::~SteppingAction()
{;}

//---------------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* endPoint = aStep->GetPostStepPoint();

  FluxSD* fluxSD = (FluxSD*)detector->GetFluxSD();
  G4Track* track = aStep->GetTrack();
  
  if( prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detector->GetBlockerVol() )
    track->SetTrackStatus(fStopAndKill);
  
  else {
    for(int i=0; i<=detector->GetNoSD(); i++) {
      
      if ( prePoint->GetTouchableHandle()->GetVolume() != detector->GetDetVol(i) &&
	   endPoint->GetTouchableHandle()->GetVolume() == detector->GetDetVol(i) ) {
	if(fluxSD) 
	  fluxSD->ProcessHits_constStep(aStep,NULL);
      }
    }
  }
}

//---------------------------------------------------------------------------


