#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;

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

  G4UIcmdWithADouble*          fTargetLengthCmd;
  G4UIcmdWithADouble*          fBigBiteDistanceCmd;
  G4UIcmdWithADouble*          fBigBiteAngleCmd;
  G4UIcmdWithADouble*          fBigBiteFieldCmd;
  G4UIcmdWithADouble*          fBigCalDistanceCmd;
  G4UIcmdWithADouble*          fBigCalAngleCmd;
  G4UIcmdWithADouble*          fHMSAngleCmd;

  G4UIcommand*                 fUpdateCmd;
};

#endif

//---------------------------------------------------------------------------

