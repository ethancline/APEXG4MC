#include "AnalysisMessenger.hh"
#include "AnalysisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

AnalysisMessenger::AnalysisMessenger(AnalysisManager* anMana)
:fAnalysisManager(anMana)
{
  fAnalysisDir = new G4UIdirectory("/APEXG4MC/analysis/");
  fAnalysisDir->SetGuidance("Analysis output control");

  fOutFileCmd = new G4UIcmdWithAString("/APEXG4MC/analysis/setOutputFile",this);
  fOutFileCmd->SetGuidance("Set the full name and path of the output file");
  fOutFileCmd->SetParameterName("choice",true);
  fOutFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
} 

//---------------------------------------------------------------------------

AnalysisMessenger::~AnalysisMessenger()
{
}

//---------------------------------------------------------------------------

void AnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == fOutFileCmd)
    {fAnalysisManager->SetOutFileName(newValue.data());}
}

//---------------------------------------------------------------------------
