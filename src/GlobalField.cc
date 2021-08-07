#include "GlobalField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "BField_Septum_New.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//---------------------------------------------------------------------------

GlobalField::GlobalField(G4double mom, G4double scale, const char* mapname) {

  fHRSMomentum = mom;
  fSeptumFieldScale = scale; 
  
  fMapField = new BField_Septum_New( mapname );
}

//---------------------------------------------------------------------------

GlobalField::~GlobalField() {
  delete fMapField;
}

//---------------------------------------------------------------------------

void GlobalField::GetFieldValue(const double Point[3], double *Bfield) const {

  unsigned int i;
  for( i = 0; i < 3; i++ ){ Bfield[i] = 0.0; }

  if (Point[2]>(-100.)*cm && Point[2]<(250)*cm)
    {
	  G4double pos_sept[3]={Point[0], Point[1], Point[2]};
	  G4double B_sept[3]={0,0,0};
	  fMapField->GetBField(pos_sept, B_sept);
	  B_sept[0]*=2.2/2.140045;
	  B_sept[1]*=2.2/2.140045;
	  B_sept[2]*=2.2/2.140045;

 	  Bfield[0]+=fSeptumFieldScale*B_sept[0]*(fHRSMomentum/GeV)/2.2;
	  Bfield[1]+=fSeptumFieldScale*B_sept[1]*(fHRSMomentum/GeV)/2.2;
 	  Bfield[2]+=fSeptumFieldScale*B_sept[2]*(fHRSMomentum/GeV)/2.2;
    }
  return;
}
