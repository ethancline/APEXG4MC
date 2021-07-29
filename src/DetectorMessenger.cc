#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

//---------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger(DetectorConstruction* Detect)
  :fDetector(Detect)
{
  fDetectorDir        = new G4UIdirectory("/APEXG4MC/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fTargetLengthCmd     = new G4UIcmdWithADouble("/APEXG4MC/detector/TargetLength",this);
  fTargetLengthCmd->SetGuidance("Set the LH2 target length in cm.");

  fHRSAngleCmd        = new G4UIcmdWithADouble("/APEXG4MC/detector/HRSAngle",this);
  fHRSAngleCmd->SetGuidance("Set the HRS angle in degrees.");

  fUpdateCmd          = new G4UIcommand("/APEXG4MC/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");  
}

//---------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete fDetectorDir;

  delete fTargetLengthCmd;
  delete fHRSAngleCmd;

  delete fUpdateCmd;
}

//---------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == fUpdateCmd)
    fDetector->UpdateGeometry();
  else if (command == fTargetLengthCmd)
    fDetector->SetTargetLength(fTargetLengthCmd->GetNewDoubleValue(newValue));
  else if (command == fHRSAngleCmd)
   fDetector->SetHRSAngle(fHRSAngleCmd->GetNewDoubleValue(newValue));
}

//---------------------------------------------------------------------------
