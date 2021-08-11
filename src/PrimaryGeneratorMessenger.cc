#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "Randomize.hh"

//---------------------------------------------------------------------------

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger( PrimaryGeneratorAction* Gun )
:fAction(Gun)
{
  fGunDir = new G4UIdirectory("/APEXG4MC/generator/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fSetInputCmd = new G4UIcmdWithAString("/APEXG4MC/generator/InputFile",this);
  fSetInputCmd->SetGuidance("Set the full name and path of the input ROOT tree.");
  fSetInputCmd->SetParameterName("inputfile",false);
  fSetInputCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSetModeCmd = new G4UIcmdWithAnInteger("/APEXG4MC/generator/Mode",this);
  fSetModeCmd->SetGuidance("Set the mode of the generator (0 for GPS or 1 for ROOT tree).");
  fSetModeCmd->SetParameterName("Mode",false);
  fSetModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSetSeedCmd = new G4UIcmdWithAnInteger("/APEXG4MC/generator/Seed",this);
  fSetSeedCmd->SetGuidance("Set the random seed for the generator");
  fSetSeedCmd->SetParameterName("Seed",false);
  fSetSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//---------------------------------------------------------------------------


PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fGunDir;
  delete fSetInputCmd;
  delete fSetModeCmd;
  delete fSetSeedCmd;
}

//---------------------------------------------------------------------------


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fSetInputCmd )
    { fAction->SetUpROOTInput(static_cast<TString>(newValue));}
   
  if( command == fSetModeCmd )
     { fAction->SetMode(fSetModeCmd->GetNewIntValue(newValue));}

 //  if( command == fSetSeedCmd )
//     { CLHEP::HepRandom::setTheSeed(fSetSeedCmd->GetNewIntValue(newValue));}

}

//---------------------------------------------------------------------------


