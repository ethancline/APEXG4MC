#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

#include "EventAction.hh"
#include "FluxHit.hh"
#include "EnergyDepositHit.hh"
#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"

//---------------------------------------------------------------------------

EventAction::EventAction( AnalysisManager* ana, PrimaryGeneratorAction* pga )
  :fAnaManager(ana), fPGA(pga)
{;}

//---------------------------------------------------------------------------

EventAction::~EventAction()
{;}

//---------------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event* evt)
{ 
  if( evt->GetEventID() == 0 )
    fAnaManager->InitOutput();
  
}

//---------------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event* evt)
{

  G4int event_id   = evt->GetEventID();  
  if ( event_id%1 == 0 ) 
    G4cout <<"Event " << event_id << G4endl;

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  FluxHit* hit;
  G4int nfluxhits = 0;
  EnergyDepositHit* edhit;
  G4int nedephits = 0;

  if(HCE) {
    G4int CollSize = HCE->GetNumberOfCollections();
    G4int hci      = 0;
    
    fAnaManager->ZeroArray(); 

    for(G4int i = 0; i < CollSize; i++) {
      FluxHitsCollection* hc;
      while(!(hc = static_cast<FluxHitsCollection*>(HCE->GetHC(hci++))));
      G4int hc_nhits = hc->entries();
      
      if( hc_nhits == 0 ) continue;
      
      // Fill output hit arrays for particles fluxes
      else if( hc->GetName().contains("Flux")  ) 
	{
	  for(G4int j = 0; j < hc_nhits; j++) {
	    hit = static_cast<FluxHit*>( hc->GetHit(j) );
	    fAnaManager->SetFluxPDef( (G4ParticleDefinition*) hit->GetPDef() );
	    fAnaManager->SetFluxP3( (G4ThreeVector) hit->GetMom() );
	    fAnaManager->SetFluxPosPre( (G4ThreeVector) hit->GetPosPre() );
	    fAnaManager->SetFluxPosPost( (G4ThreeVector) hit->GetPosPost() );
	    fAnaManager->SetFluxTime( (G4double) hit->GetTime() );
	    fAnaManager->SetFluxEnergy( (G4double) hit->GetEdep() );
	    fAnaManager->SetFluxID( (G4int) hit->GetID() );
	    fAnaManager->SetFluxPID( (G4int) hit->GetParentID() );
	    fAnaManager->SetFluxTID( (G4int) hit->GetTrackID() );
	    fAnaManager->SetFluxANum( (G4int) hit->GetAtomicNumber() );
	    fAnaManager->SetFluxAMass( (G4int) hit->GetAtomicMass() );

	    fAnaManager->FillFluxArray( nfluxhits ); 
	    nfluxhits++;
	  }
	}

    }
    
    if( nfluxhits != 0 || nedephits != 0 ) {
      // Fill output data from PrimaryGeneratorAction
      fAnaManager->SetPrimaryDirection ( (G4ThreeVector)fPGA->GetDirection() );
      fAnaManager->SetPrimaryEnergy    ( (G4double)fPGA->GetEnergy() );
      fAnaManager->SetPrimaryTime      ( (G4double)fPGA->GetTime() );
      fAnaManager->SetPrimaryPDef      ( (G4ParticleDefinition*)fPGA->GetPrimPDef() );
      
      fAnaManager->FillTree(); 
    }
  }

}

//---------------------------------------------------------------------------
