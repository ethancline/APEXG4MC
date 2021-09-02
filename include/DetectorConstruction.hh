#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "DetectorMessenger.hh"

//---------------------------------------------------------------------------

class G4VPhysicalVolume;
class FluxSD;
class EnergyDepositSD;

//---------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  inline G4VPhysicalVolume* GetExpHall()       { return fExpHall;   };
  inline G4VPhysicalVolume* GetDetVol(G4int i) { return fDetVol[i]; };
  inline G4VPhysicalVolume* GetSurfVol()       { return fSurfVol;   };
  inline FluxSD*            GetFluxSD()        { return fFluxSD;    };
  inline G4int              GetNoSD()          { return fNSD;       };

  void SetHRSAngle      ( G4double a  ) { fHRSAngle     = a;  }
  void SetHRSMomentum   ( G4double p  ) { fHRSMomentum  = p;  }
  void SetDistTarPivot  ( G4double tp ) { fDistTarPivot = tp; }
  void SetDistPivotQ1   ( G4double pq ) { fDistPivotQ1  = pq; }
  void SetSeptumScale   ( G4double ss ) { fScaleSeptum  = ss; }
  void SetSeptFieldMap  ( G4String mf ) { fFieldMapFile = mf; }
  void SetSieveOn       ( G4bool o    ) { fSieveOn      = o;  }
  void SetSieveAngle    ( G4double sa ) { fSieveAngle   = sa; }
  void SetTarget        ( G4String t  ) { fTarget       = t;  }
  
  private:

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  static const G4int fMaxNSD = 50;

  G4int fNSD;
  
  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDetVol[fMaxNSD];
  G4VPhysicalVolume* fSurfVol;

  FluxSD*            fFluxSD;
  EnergyDepositSD*   fEdepSD;

  G4double           fHRSAngle;
  G4double           fDistTarPivot;
  G4double           fDistPivotQ1;
  G4double           fHRSMomentum;
  G4double           fScaleSeptum;
  G4String           fFieldMapFile;
  G4bool             fSieveOn;
  G4double           fSieveAngle;
  G4String           fTarget;

  G4int              fTargetType;
  enum Targ_t { kProdW, kProdC, kOptics1, kOptics2, kOptics3, kVWires, kHWires };

};
#endif

//---------------------------------------------------------------------------

