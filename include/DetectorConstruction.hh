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

namespace APEX {
  enum Targ_t { kProdW, kProdC, kOptics1, kOptics2, kOptics3, kVWires, kHWires };
}

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  
  inline G4VPhysicalVolume* GetExpHall()       { return fExpHall;     };
  inline G4VPhysicalVolume* GetDetVol(G4int i) { return fDetVol[i];   };
  inline G4LogicalVolume*   GetBlockerVol()    { return fBlockerLog;  };
  inline FluxSD*            GetFluxSD()        { return fFluxSD;      };
  inline G4int              GetNoSD()          { return fNSD;         };
  inline G4int              GetTargetType()    { return fTargetType;  };
  inline G4int              GetNTargets()      { return fNtargs;      };
  inline G4double           GetTargetWidth()   { return fTargwidth;   };
  inline G4double           GetTargetThick()   { return fTargthick;   };
  inline G4double*          GetTargetXPos()    { return fTargxpos;    };
  inline G4double*          GetTargetYPos()    { return fTargypos;    };
  inline G4double*          GetTargetZPos()    { return fTargzpos;    };
  inline G4double           GetHRSMomentum()   { return fHRSMomentum; };
  inline G4double           GetDistTarPivot()  { return fDistTarPivot; };

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
  G4LogicalVolume*   fBlockerLog;

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

  static const G4int fMaxTargs = 10;
  G4int       fNtargs;
  G4double    fTargwidth;
  G4double    fTargthick;
  G4Material* fTargMaterial;
  G4double    fTargxpos[fMaxTargs];
  G4double    fTargypos[fMaxTargs];
  G4double    fTargzpos[fMaxTargs];
  
};
#endif

//---------------------------------------------------------------------------

