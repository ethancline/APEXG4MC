#ifndef GlobalField_hh
#define GlobalField_hh

#include "globals.hh"
#include <vector>

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"

class BField_Septum_New;

class GlobalField : public G4MagneticField {
public:
  GlobalField(G4double, G4double, const char*);
  ~GlobalField();

  void GetFieldValue( const  double Point[3], double *Bfield ) const;

  BField_Septum_New *fMapField;
  G4double fSeptumFieldScale;
  G4double fHRSMomentum;


};

#endif
