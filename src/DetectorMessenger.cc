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
  fDetectorDir        = new G4UIdirectory("/ComptonMC/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fTargetLengthCmd     = new G4UIcmdWithADouble("/ComptonMC/detector/TargetLength",this);
  fTargetLengthCmd->SetGuidance("Set the LH2 target length in cm.");

  fBigBiteDistanceCmd = new G4UIcmdWithADouble("/ComptonMC/detector/BigBiteDistance",this);
  fBigBiteDistanceCmd->SetGuidance("Set the BigBite distance between EFB and the hall pivot in m.");
  fBigBiteAngleCmd    = new G4UIcmdWithADouble("/ComptonMC/detector/BigBiteAngle",this);
  fBigBiteAngleCmd->SetGuidance("Set the BigBite angle in degrees.");
  fBigBiteFieldCmd    = new G4UIcmdWithADouble("/ComptonMC/detector/BigBiteField",this);
  fBigBiteFieldCmd->SetGuidance("Set the magnetic flux density in tesla (Bx).");

  fBigCalDistanceCmd  = new G4UIcmdWithADouble("/ComptonMC/detector/BigCalDistance",this);
  fBigCalDistanceCmd->SetGuidance("Set the BigCal distance between front face and the hall pivot in m.");
  fBigCalAngleCmd     = new G4UIcmdWithADouble("/ComptonMC/detector/BigCalAngle",this);
  fBigCalAngleCmd->SetGuidance("Set the BigCal angle in degrees.");

  fHMSAngleCmd        = new G4UIcmdWithADouble("/ComptonMC/detector/HMSAngle",this);
  fHMSAngleCmd->SetGuidance("Set the HMS angle in degrees.");

  fUpdateCmd          = new G4UIcommand("/ComptonMC/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");  
}

//---------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete fDetectorDir;

  delete fTargetLengthCmd;
  delete fBigBiteDistanceCmd;
  delete fBigBiteAngleCmd;
  delete fBigBiteFieldCmd;
  delete fBigCalDistanceCmd;
  delete fBigCalAngleCmd;
  delete fHMSAngleCmd;

  delete fUpdateCmd;
}

//---------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == fUpdateCmd)
    fDetector->UpdateGeometry();
  else if (command == fTargetLengthCmd)
    fDetector->SetTargetLength(fTargetLengthCmd->GetNewDoubleValue(newValue));
  else if (command == fBigBiteDistanceCmd)
    fDetector->SetBigBiteDistance(fBigBiteDistanceCmd->GetNewDoubleValue(newValue));
  else if (command == fBigBiteAngleCmd)
   fDetector->SetBigBiteAngle(fBigBiteAngleCmd->GetNewDoubleValue(newValue));
  else if (command == fBigBiteFieldCmd)
   fDetector->SetBigBiteField(fBigBiteFieldCmd->GetNewDoubleValue(newValue));
  else if (command == fBigCalDistanceCmd)
    fDetector->SetBigCalDistance(fBigCalDistanceCmd->GetNewDoubleValue(newValue));
  else if (command == fBigCalAngleCmd)
   fDetector->SetBigCalAngle(fBigCalAngleCmd->GetNewDoubleValue(newValue));
  else if (command == fHMSAngleCmd)
   fDetector->SetHMSAngle(fHMSAngleCmd->GetNewDoubleValue(newValue));
}

//---------------------------------------------------------------------------
