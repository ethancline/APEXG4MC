#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh" 
#include "G4ThreeVector.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class DetectorConstruction;
class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

enum { EPGA_GPS, EPGA_APEX};

//---------------------------------------------------------------------------

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction( DetectorConstruction* );    
  ~PrimaryGeneratorAction();

  void GeneratePrimaries (G4Event*);

  void SetMode(G4int mo)          { fMode       = mo; }
  void SetBeamE(G4double eb)      { fBeamE      = eb; }
  void SetRasterX(G4double rx)    { fRasterX    = rx; }
  void SetRasterY(G4double ry)    { fRasterY    = ry; }
  void SetThetaMin(G4double tm)   { fThMin      = tm; }
  void SetThetaMax(G4double tx)   { fThMax      = tx; }
  void SetPhiMin(G4double pm)     { fPhMin      = pm; }
  void SetPhiMax(G4double px)     { fPhMax      = px; }
  void SetDeltaRange(G4double dr) { fDeltaRange = dr; }

  G4int         GetMode()             { return fMode; }
  G4ThreeVector GetVertex()           { return G4ThreeVector(fVx,  fVy,  fVz);  }
  G4ThreeVector GetDirection()        { return G4ThreeVector(fPxp, fPyp, fPzp); }
  G4double      GetEnergy()           { return fEp; }
  G4double      GetTime()             { return fTp; }
  G4ParticleDefinition* GetPrimPDef() { return fPDefinition; }

private:
  PrimaryGeneratorMessenger* fGunMessenger;

  DetectorConstruction*      fDetector;
  G4ParticleGun*             fParticleGun;    
  G4GeneralParticleSource*   fParticleSource;
  G4ParticleDefinition*      fPDefinition;
  G4int                      fMode;   
  
  G4double                   fVx;
  G4double                   fVy;
  G4double                   fVz;
  G4double                   fPxp;
  G4double                   fPyp;
  G4double                   fPzp;
  G4double                   fEp;
  G4double                   fTp;

  G4double                   fBeamE;
  G4double                   fRasterX;
  G4double                   fRasterY;
  G4double                   fThMin;
  G4double                   fThMax;
  G4double                   fPhMin;
  G4double                   fPhMax;
  G4double                   fDeltaRange;
};
#endif

//---------------------------------------------------------------------------


