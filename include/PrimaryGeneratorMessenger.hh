#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

//---------------------------------------------------------------------------

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
   ~PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  PrimaryGeneratorAction* fAction;
  G4UIdirectory*          fGunDir;

  G4UIcmdWithAnInteger*   fSetModeCmd;
  G4UIcmdWithADouble*     fBeamECmd;
  G4UIcmdWithADouble*     fRasterXCmd;
  G4UIcmdWithADouble*     fRasterYCmd;
  G4UIcmdWithADouble*     fThMinCmd;
  G4UIcmdWithADouble*     fThMaxCmd;
  G4UIcmdWithADouble*     fPhMinCmd;
  G4UIcmdWithADouble*     fPhMaxCmd;
  G4UIcmdWithADouble*     fDeltaCmd;
};

#endif

//---------------------------------------------------------------------------
