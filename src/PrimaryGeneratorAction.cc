#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4String.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
  :fDetector(det)
{
  fMode              = 1;
  fBeamE             = 2.2;
  fRasterX           = 2.0;
  fRasterY           = 2.0;
  fThMin             = 2.0;
  fThMax             = 8.0;
  fPhMin             = -5.0;
  fPhMax             = 5.0;
  fDeltaRange        = 0.2;
  
  fGunMessenger      = new PrimaryGeneratorMessenger(this);

  fParticleSource    = new G4GeneralParticleSource();
  fParticleGun       = new G4ParticleGun(1);
}

//---------------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fParticleSource;
  delete fGunMessenger;
}

//---------------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4int itarg;
  G4double th, ph, p;
  G4ThreeVector p3_lab;
  
  switch(fMode) {

  case EPGA_GPS:

    fParticleSource->GeneratePrimaryVertex(anEvent);

    // Get the electron primary variables for the output tree (for now)
    fVx          = fParticleSource->GetParticlePosition().getX();
    fVy          = fParticleSource->GetParticlePosition().getY();
    fVz          = fParticleSource->GetParticlePosition().getZ();
    fPxp         = fParticleSource->GetParticleMomentumDirection().getX();
    fPyp         = fParticleSource->GetParticleMomentumDirection().getY();
    fPzp         = fParticleSource->GetParticleMomentumDirection().getZ();
    fEp          = fParticleSource->GetParticleEnergy()/GeV;
    fTp          = fParticleSource->GetParticleTime();
    fPDefinition = fParticleSource->GetParticleDefinition();

    break;

  case EPGA_APEX:

    // first check that raster is not wider than target foil
    if( fRasterX*mm > fDetector->GetTargetWidth() )
      fRasterX = fDetector->GetTargetWidth();
    if( fRasterY*mm > fDetector->GetTargetWidth() )
      fRasterY = fDetector->GetTargetWidth();

    // restrict raster for wire targets
    if( fDetector->GetTargetType() == APEX::kVWires )
      fRasterX = fDetector->GetTargetThick();
    if( fDetector->GetTargetType() == APEX::kHWires )
      fRasterY = fDetector->GetTargetThick();

    // now generate vertex
    itarg = CLHEP::RandFlat::shootInt( fDetector->GetNTargets() );
    fVx   = fDetector->GetTargetXPos()[itarg] + CLHEP::RandFlat::shoot(-fRasterX/2.0, fRasterX/2.0 );
    fVy   = fDetector->GetTargetYPos()[itarg] + CLHEP::RandFlat::shoot(-fRasterY/2.0, fRasterY/2.0 );
    fVz   = -fDetector->GetDistTarPivot()*cm + fDetector->GetTargetZPos()[itarg] +
      CLHEP::RandFlat::shoot(-fDetector->GetTargetThick()/2.0, fDetector->GetTargetThick()/2.0 );

    // Generate electron in the lab
    th = acos( CLHEP::RandFlat::shoot(cos(fThMax*deg), cos(fThMin *deg)) );
    ph = CLHEP::RandFlat::shoot(fPhMin *deg, fPhMax *deg);
    p  = CLHEP::RandFlat::shoot((1-fDeltaRange/2.)*fDetector->GetHRSMomentum()*GeV,
				(1+fDeltaRange/2.)*fDetector->GetHRSMomentum()*GeV);

    p3_lab.setX( p*sin(th)*cos(ph) );
    p3_lab.setY( p*sin(th)*sin(ph) );
    p3_lab.setZ( p*cos(th) );
    
    fPDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");      
    
    fParticleGun->SetParticlePosition( G4ThreeVector( fVx, fVy, fVz) );
    fParticleGun->SetParticleMomentumDirection( p3_lab.unit() );
    fParticleGun->SetParticleEnergy( p );
    fParticleGun->SetParticleDefinition( fPDefinition );
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Get the electron primary variables for the output tree (for now)
    fVx          = fParticleGun->GetParticlePosition().getX();
    fVy          = fParticleGun->GetParticlePosition().getY();
    fVz          = fParticleGun->GetParticlePosition().getZ();
    fPxp         = fParticleGun->GetParticleMomentumDirection().getX();
    fPyp         = fParticleGun->GetParticleMomentumDirection().getY();
    fPzp         = fParticleGun->GetParticleMomentumDirection().getZ();
    fEp          = fParticleGun->GetParticleEnergy()/GeV;
    fTp          = fParticleGun->GetParticleTime();
    fPDefinition = fParticleGun->GetParticleDefinition();

    // now generate the partner positron (keeping it simple for now)
    fPDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e+");      

    p3_lab.setX( (fBeamE*GeV - p)*sin(-th)*cos(-ph) );
    p3_lab.setY( (fBeamE*GeV - p)*sin(-th)*sin(-ph) );
    p3_lab.setZ( (fBeamE*GeV - p)*cos(-th) );
    
    fParticleGun->SetParticleMomentumDirection( p3_lab.unit() );
    fParticleGun->SetParticleEnergy( fBeamE*GeV - p );
    fParticleGun->SetParticleDefinition( fPDefinition );
    fParticleGun->GeneratePrimaryVertex(anEvent);

    break;

  default:
    G4cout << "Unknown mode given to PrimiaryGeneratorAction (0 for gps, 1 for apex)" << G4endl;
  }			       
}

//---------------------------------------------------------------------------
