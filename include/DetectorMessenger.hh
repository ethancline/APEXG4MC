#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

//---------------------------------------------------------------------------

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  DetectorConstruction*        fDetector;
  G4UIdirectory*               fDetectorDir;

  G4UIcmdWithADouble*          fHRSAngleCmd;
  G4UIcmdWithADouble*          fHRSMomCmd;
  G4UIcmdWithADouble*          fDistTPCmd;
  G4UIcmdWithADouble*          fDistPQ1Cmd;
  G4UIcmdWithADouble*          fSepScaleCmd;
  G4UIcmdWithAString*          fFieldMapCmd;  
  G4UIcmdWithABool*            fSieveOnCmd;
  G4UIcmdWithADouble*          fSieveAngleCmd;
  G4UIcmdWithAString*          fTargetCmd;
  
  G4UIcommand*                 fUpdateCmd;
};

#endif

//---------------------------------------------------------------------------

