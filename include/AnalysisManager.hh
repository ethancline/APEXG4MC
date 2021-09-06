#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "globals.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMessenger.hh"

#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "Rtypes.h"
#include "TVector3.h"
#include "TString.h"

class TTree;
class TFile;

//---------------------------------------------------------------------------

class AnalysisManager {

public:

  AnalysisManager();
  ~AnalysisManager();

  void InitOutput();

  void ZeroArray();
  void FillFluxArray( Int_t );
  void FillTree();

  void SetOutFileName  ( TString fname )           { fOutFileName  = fname; }

  void SetPrimaryEnergy   (G4double       ene  )       { fPEne  = ene;  }
  void SetPrimaryTime     (G4double       time )       { fPTime = time; }
  void SetPrimaryPDef     (G4ParticleDefinition* pdef) { fPPDef = pdef; }
  void SetPrimaryDirection(G4ThreeVector  dir  )       { fPdir  = dir;  }
  void SetPrimaryVertex   (G4ThreeVector  vtx  )       { fPvtx  = vtx;  }

  void SetFluxPDef      ( G4ParticleDefinition* sp ) { fFluxpdef = sp;    }
  void SetFluxPosPre    ( G4ThreeVector  spos )      { fFluxpospre  = spos;  }
  void SetFluxPosPost   ( G4ThreeVector  spos )      { fFluxpospost  = spos;  }
  void SetFluxP3        ( G4ThreeVector  smom )      { fFluxp3   = smom;  }
  void SetFluxPolarization ( G4ThreeVector  pol )    { fFluxpol  = pol;  }
  void SetFluxTime      ( G4double       stime )     { fFluxtime = stime; }
  void SetFluxID        ( G4int          sid )       { fFluxid   = sid;   }
  void SetFluxEnergy    ( G4double       sedep)      { fFluxedep = sedep; }
  void SetFluxANum      ( G4int          ann )       { fFluxANum  = ann;   }
  void SetFluxAMass     ( G4int          amm )       { fFluxAMass = amm;   }
  void SetFluxPID       ( G4int          pid )       { fFluxPid   = pid;   }
  void SetFluxTID       ( G4int          tid )       { fFluxTid   = tid;   }

private:
  
  AnalysisMessenger*    fAnaMessenger;
  TString               fOutFileName;
  TFile*                fROOTfile;
  TTree*                fROOTtree;
  
  // Primary
  Float_t               fPEne;
  Float_t               fPth;
  Float_t               fPph;
  Float_t               fPTime;
  Float_t               fPVx;
  Float_t               fPVy;
  Float_t               fPVz;
  G4ParticleDefinition* fPPDef;
  Int_t                 fPpdg;
  G4ThreeVector         fPdir;
  G4ThreeVector         fPvtx;

  // Flux raw
  G4ParticleDefinition* fFluxpdef;
  G4ThreeVector         fFluxp3;
  G4ThreeVector         fFluxpol;
  G4ThreeVector         fFluxpospre;
  G4ThreeVector         fFluxpospost;
  G4double              fFluxtime;
  G4int                 fFluxid;
  G4int                 fFluxPid;
  G4int                 fFluxTid;
  G4int                 fFluxANum;
  G4int                 fFluxAMass;
  G4double              fFluxedep;

  static const Int_t    fMaxhits = 50000;

  Int_t                 fFlux_Nhits;
  Int_t                 fFlux_id[fMaxhits];
  Float_t               fFlux_time[fMaxhits];
  Float_t               fFlux_Edep[fMaxhits];
  Int_t                 fFlux_pdg[fMaxhits];
  Int_t                 fFlux_tid[fMaxhits];//Track ID 
  Int_t                 fFlux_pid[fMaxhits];// Process ID 
  Int_t                 fFlux_anum[fMaxhits];//Atomic Number
  Int_t                 fFlux_amass[fMaxhits];//Atomic Mass
  Float_t               fFlux_mass[fMaxhits];
  Float_t               fFlux_mom[fMaxhits];
  Float_t               fFlux_px[fMaxhits];
  Float_t               fFlux_py[fMaxhits];
  Float_t               fFlux_pz[fMaxhits];
  Float_t               fFlux_polx[fMaxhits];
  Float_t               fFlux_poly[fMaxhits];
  Float_t               fFlux_polz[fMaxhits];
  Float_t               fFlux_th[fMaxhits];
  Float_t               fFlux_ph[fMaxhits];
  Float_t               fFlux_xpre[fMaxhits];
  Float_t               fFlux_ypre[fMaxhits];
  Float_t               fFlux_zpre[fMaxhits];
  Float_t               fFlux_xpost[fMaxhits];
  Float_t               fFlux_ypost[fMaxhits];
  Float_t               fFlux_zpost[fMaxhits];
  Float_t               fFlux_Energy[fMaxhits];

};

#endif

//---------------------------------------------------------------------------
