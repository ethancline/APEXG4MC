// BField_Septum_New.h: interface for the BField_Septum_New class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BFIELD_Septum_New_H)
#define BFIELD_Septum_New_H
#include <math.h>
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>
#include "G4MagneticField.hh"
using namespace CLHEP;

class BField_Septum_New 

{
public:
  BField_Septum_New( const char *mapfile="Septa-JB_map.table" );
  virtual ~BField_Septum_New();
  void GetBField(double Pos[3],double B[3]);

private:
  void ReadMap(const char *filename);
  
private:
  double ****mBField;
  float ****Btxt;
  float ****BCoord;
  int    nx;
  int    ny;
  int    nz;
  
};

#endif // !defined(BFIELD_Septum_H)
