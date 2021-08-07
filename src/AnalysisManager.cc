#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

//---------------------------------------------------------------------------

AnalysisManager::AnalysisManager()
{
  ZeroArray();

  fOutFileName = TString("output/out_default.root");

  fAnaMessenger = new AnalysisMessenger(this);

}

//---------------------------------------------------------------------------

AnalysisManager::~AnalysisManager()
{
  if(fROOTfile) {
    fROOTtree->Write();
    fROOTfile->Close();
  }
}

//---------------------------------------------------------------------------

void AnalysisManager::InitOutput()
{ 

  fROOTfile = new TFile(fOutFileName,"RECREATE","fROOTfile",1);
  fROOTtree = new TTree("T","Output Tree");
  fROOTtree->SetAutoSave();

  // Set Primary Branches
  fROOTtree->Branch("Prim_E",      &fPEne,   "Prim_E/F"      );
  fROOTtree->Branch("Prim_Th",     &fPth,    "Prim_Th/F"     );
  fROOTtree->Branch("Prim_Ph",     &fPph,    "Prim_Ph/F"     );
  fROOTtree->Branch("Prim_pdg",    &fPpdg,   "Prim_pdg/I"    );

  // Set Flux Hit Branches
  fROOTtree->Branch("Flux_Nhits", &fFlux_Nhits, "Flux_Nhits/I");  
  fROOTtree->Branch("Flux_pdg",   fFlux_pdg,    "Flux_pdg[Flux_Nhits]/I");
  fROOTtree->Branch("Flux_id",    fFlux_id,     "Flux_id[Flux_Nhits]/I");
  fROOTtree->Branch("Flux_x",     fFlux_xpre,   "Flux_x[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_y",     fFlux_ypre,   "Flux_y[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_z",     fFlux_zpre,   "Flux_z[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_polx",     fFlux_polx,   "Flux_polx[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_poly",     fFlux_poly,   "Flux_poly[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_polz",     fFlux_polz,   "Flux_polz[Flux_Nhits]/F"  );
  fROOTtree->Branch("Flux_E",     fFlux_Edep,   "Flux_E[Flux_Nhits]/F" );
  fROOTtree->Branch("Flux_A",     fFlux_amass,  "Flux_A[Flux_Nhits]/I");
  fROOTtree->Branch("Flux_Z",     fFlux_anum,   "Flux_Z[Flux_Nhits]/I");  
  fROOTtree->Branch("Flux_tid",   fFlux_tid,    "Flux_tid[Flux_Nhits]/I");
  fROOTtree->Branch("Flux_pid",   fFlux_pid,    "Flux_pid[Flux_Nhits]/I");
  fROOTtree->Branch("Flux_th",    &fFlux_th,    "Flux_th[Flux_Nhits]/F"     );
  fROOTtree->Branch("Flux_ph",    &fFlux_ph,    "Flux_ph[Flux_Nhits]/F"     );
}

//---------------------------------------------------------------------------

void AnalysisManager::ZeroArray()
{
  // Primary
  G4ThreeVector zero(0.,0.,0.);
  fPEne   = 9999;
  fPdir   = (zero);
  fPth    = 9999;
  fPph    = 9999;
  fPTime  = 9999;
  fPPDef  = NULL;
  fPpdg   = 9999;

  // Raw Hits
  fFlux_Nhits  = 0;
  fFluxpdef    = NULL;
  fFluxp3      = (zero);
  fFluxpol     = (zero);
  fFluxpospre  = (zero);
  fFluxpospost = (zero);
  fFluxid      = 0;
  fFluxtime    = 0;
  fFluxedep    = 0;
  fFluxTid     = 0;
  fFluxANum    = 0;
  fFluxAMass   = 0;

  for ( Int_t i = 0; i < fMaxhits; i++ ) {
    fFlux_id[i]      = 9999;  
    fFlux_time[i]    = 9999;
    fFlux_Edep[i]    = 9999;
    fFlux_pdg[i]     = 9999;
    fFlux_mass[i]    = 9999;
    fFlux_mom[i]     = 9999;
    fFlux_px[i]      = 9999;
    fFlux_py[i]      = 9999;
    fFlux_pz[i]      = 9999;
    fFlux_polx[i]      = 9999;
    fFlux_poly[i]      = 9999;
    fFlux_polz[i]      = 9999;
    fFlux_th[i]      = 9999;
    fFlux_ph[i]      = 9999;

    fFlux_xpre[i]    = 9999;
    fFlux_ypre[i]    = 9999;
    fFlux_zpre[i]    = 9999;
    fFlux_xpost[i]   = 9999;
    fFlux_ypost[i]   = 9999;
    fFlux_zpost[i]   = 9999;
    fFlux_Energy[i]  = 9999;
    fFlux_pid[i]     = 9999;  
    fFlux_tid[i]     = 9999;  
    fFlux_amass[i]   = 9999;  
    fFlux_anum[i]    = 9999; 

  }

}

//---------------------------------------------------------------------------

void AnalysisManager::FillFluxArray( Int_t hitn ) 
{
    fFlux_Nhits++;
    fFlux_id[hitn]     = (Int_t)fFluxid;
    fFlux_tid[hitn]    = (Int_t)fFluxTid;
    fFlux_pid[hitn]    = (Int_t)fFluxPid;
    fFlux_anum[hitn]   = (Int_t)fFluxANum;
    fFlux_amass[hitn]  = (Int_t)fFluxAMass;
    fFlux_pdg[hitn]    = (Int_t)fFluxpdef->GetPDGEncoding();
    fFlux_mass[hitn]   = (Float_t)fFluxpdef->GetPDGMass();
    fFlux_time[hitn]   = (Float_t)fFluxtime;                                   
    fFlux_mom[hitn]    = (Float_t)fFluxp3.mag();                             
    fFlux_px[hitn]     = (Float_t)fFluxp3.getX();                             
    fFlux_py[hitn]     = (Float_t)fFluxp3.getY();                             
    fFlux_pz[hitn]     = (Float_t)fFluxp3.getZ();
    fFlux_polx[hitn]     = (Float_t)fFluxpol.getX();                             
    fFlux_poly[hitn]     = (Float_t)fFluxpol.getY();                             
    fFlux_polz[hitn]     = (Float_t)fFluxpol.getZ();
    fFlux_th[hitn]     = (Float_t)fFluxp3.getTheta();
    fFlux_ph[hitn]     = (Float_t)fFluxp3.getPhi();                             
    fFlux_xpre[hitn]   = (Float_t)fFluxpospre.getX();                             
    fFlux_ypre[hitn]   = (Float_t)fFluxpospre.getY();                             
    fFlux_zpre[hitn]   = (Float_t)fFluxpospre.getZ();                             
    fFlux_xpost[hitn]  = (Float_t)fFluxpospost.getX();                             
    fFlux_ypost[hitn]  = (Float_t)fFluxpospost.getY();                             
    fFlux_zpost[hitn]  = (Float_t)fFluxpospost.getZ();                             
    fFlux_Edep[hitn]   = (Float_t)fFluxedep;
    fFlux_Energy[hitn] = TMath::Sqrt( fFlux_mom[hitn]*fFlux_mom[hitn] 
				     + fFlux_mass[hitn]*fFlux_mass[hitn] );

    fFlux_Energy[hitn] = fFlux_mom[hitn] - fFlux_mass[hitn];

    fFlux_xpre[hitn]   = (fFlux_xpre[hitn] + fFlux_xpost[hitn])/2.;
    fFlux_ypre[hitn]   = (fFlux_ypre[hitn] + fFlux_ypost[hitn])/2.;
    fFlux_zpre[hitn]   = (fFlux_zpre[hitn] + fFlux_zpost[hitn])/2.;

}

//---------------------------------------------------------------------------

void AnalysisManager::FillTree()
{
  // Primary Variables
  fPTime  = (Float_t)fPTime;
  fPth    = (Float_t)fPdir.getTheta();                         
  fPph    = (Float_t)fPdir.getPhi();                                                      
  fPEne   = (Float_t)fPEne;                         
  fPpdg   = (Int_t)  fPPDef->GetPDGEncoding();
  
  fROOTtree->Fill();
}

//---------------------------------------------------------------------------
