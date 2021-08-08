#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"

//---------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger(DetectorConstruction* Detect)
  :fDetector(Detect)
{
  fDetectorDir        = new G4UIdirectory("/APEXG4MC/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fHRSAngleCmd        = new G4UIcmdWithADouble("/APEXG4MC/detector/HRSAngle",this);
  fHRSAngleCmd->SetGuidance("Set the HRS angle in degrees.");

  fHRSMomCmd        = new G4UIcmdWithADouble("/APEXG4MC/detector/HRSMomentum",this);
  fHRSMomCmd->SetGuidance("Set the HRS central momentum in GeV/c.");
  
  fDistTPCmd          = new G4UIcmdWithADouble("/APEXG4MC/detector/DistTarPivot",this);
  fDistTPCmd->SetGuidance("Set the distance between centre of scattering chamber and hall pivot in cm.");

  fDistPQ1Cmd         = new G4UIcmdWithADouble("/APEXG4MC/detector/DistPivotQ1",this);
  fDistPQ1Cmd->SetGuidance("Set the distance between hall pivot and Q1 entrance in cm.");

  fSepScaleCmd        = new G4UIcmdWithADouble("/APEXG4MC/detector/SeptumFieldScale",this);
  fSepScaleCmd->SetGuidance("Set the septum field scale (from nominal central momentum).");

  fFieldMapCmd = new G4UIcmdWithAString("/APEXG4MC/detector/SeptumFieldMap",this);
  fFieldMapCmd->SetGuidance("Set the full name and path of the septum field map");

  fUpdateCmd          = new G4UIcommand("/APEXG4MC/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");  
}

//---------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete fDetectorDir;

  delete fHRSAngleCmd;
  delete fHRSMomCmd;
  delete fDistTPCmd;
  delete fDistPQ1Cmd;
  delete fSepScaleCmd;
  delete fFieldMapCmd;
  
  delete fUpdateCmd;
}

//---------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == fUpdateCmd)
    fDetector->UpdateGeometry();
  else if (command == fHRSAngleCmd)
   fDetector->SetHRSAngle(fHRSAngleCmd->GetNewDoubleValue(newValue));
    else if (command == fHRSMomCmd)
   fDetector->SetHRSMomentum(fHRSMomCmd->GetNewDoubleValue(newValue));
  else if (command == fDistTPCmd)
    fDetector->SetDistTarPivot(fDistTPCmd->GetNewDoubleValue(newValue));
  else if (command == fDistPQ1Cmd)
    fDetector->SetDistPivotQ1(fDistPQ1Cmd->GetNewDoubleValue(newValue));
  else if (command == fSepScaleCmd)
    fDetector->SetSeptumScale(fSepScaleCmd->GetNewDoubleValue(newValue));
  else if(command == fFieldMapCmd)
    fDetector->SetSeptFieldMap(newValue.data());

}

//---------------------------------------------------------------------------
