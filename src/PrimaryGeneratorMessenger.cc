#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "Randomize.hh"

//---------------------------------------------------------------------------

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger( PrimaryGeneratorAction* Gun )
:fAction(Gun)
{
  fGunDir = new G4UIdirectory("/APEXG4MC/generator/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fSetModeCmd = new G4UIcmdWithAnInteger("/APEXG4MC/generator/Mode",this);
  fSetModeCmd->SetGuidance("Set the mode of the generator (0 for GPS or 1 for APEX).");

  fBeamECmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/BeamEnergy",this);
  fBeamECmd->SetGuidance("Set the beam energy GeV.");

  fRasterXCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/RasterX",this);
  fRasterXCmd->SetGuidance("Set the beam x raster in mm.");

  fRasterYCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/RasterY",this);
  fRasterYCmd->SetGuidance("Set the beam y raster in mm.");

  fThMinCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/ThMin",this);
  fThMinCmd->SetGuidance("Set the minimum polar angle for the electron in degrees.");

  fThMaxCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/ThMax",this);
  fThMaxCmd->SetGuidance("Set the maximum polar angle for the electron in degrees.");

  fPhMinCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/PhMin",this);
  fPhMinCmd->SetGuidance("Set the minimum azimuthal angle for the electron in degrees.");

  fPhMaxCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/PhMax",this);
  fPhMaxCmd->SetGuidance("Set the maximum azimuthal angle for the electron in degrees.");

  fDeltaCmd  = new G4UIcmdWithADouble("/APEXG4MC/generator/DeltaRange",this);
  fDeltaCmd->SetGuidance("Set the generated range in delta (absolute).");
  
}

//---------------------------------------------------------------------------


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fGunDir;
  delete fSetModeCmd;
}

//---------------------------------------------------------------------------


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fSetModeCmd )
     { fAction->SetMode(fSetModeCmd->GetNewIntValue(newValue));}
  else if (command == fBeamECmd)
   fAction->SetBeamE(fBeamECmd->GetNewDoubleValue(newValue));
  else if (command == fRasterXCmd)
    fAction->SetRasterX(fRasterXCmd->GetNewDoubleValue(newValue));
  else if (command == fRasterYCmd)
    fAction->SetRasterY(fRasterYCmd->GetNewDoubleValue(newValue));
  else if (command == fThMinCmd)
    fAction->SetThetaMin(fThMinCmd->GetNewDoubleValue(newValue));
  else if (command == fThMaxCmd)
    fAction->SetThetaMax(fThMaxCmd->GetNewDoubleValue(newValue));
  else if (command == fPhMinCmd)
    fAction->SetPhiMin(fPhMinCmd->GetNewDoubleValue(newValue));
  else if (command == fPhMaxCmd)
    fAction->SetPhiMax(fPhMaxCmd->GetNewDoubleValue(newValue));
  else if (command == fDeltaCmd)
    fAction->SetDeltaRange(fDeltaCmd->GetNewDoubleValue(newValue));
  
}

//---------------------------------------------------------------------------


