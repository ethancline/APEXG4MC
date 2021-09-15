#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

#include <fstream>

const Float_t phrs   = 1.063;  // GeV/c 
const Float_t thhrs  = 12.5;   // deg
const Float_t q1_ent = 1720.5; // mm

void MergeG4SIMC( TString fname="merged.root" )
{
  // ---------------------------------------------------------------------------
  // Get G4 tree
  // ---------------------------------------------------------------------------
  
  TFile *fG4     = new TFile("batch_apex.root");
  TTree* G4_tree = (TTree*)fG4->Get("T");

  const Int_t MaxHits = 500;
  
  Float_t         Prim_E;
  Float_t         Prim_Th;
  Float_t         Prim_Ph;
  Float_t         Prim_x;
  Float_t         Prim_y;
  Float_t         Prim_z;
  Int_t           Prim_pdg;
  Int_t           Flux_Nhits;
  Int_t           Flux_pdg[500];
  Int_t           Flux_id[500];
  Float_t         Flux_x[500];
  Float_t         Flux_y[500];
  Float_t         Flux_z[500];
  Float_t         Flux_E[500];
  Int_t           Flux_tid[500];
  Int_t           Flux_pid[500];
  Float_t         Flux_th[500];
  Float_t         Flux_ph[500];

  G4_tree->SetBranchAddress("Prim_E",&Prim_E);
  G4_tree->SetBranchAddress("Prim_Th",&Prim_Th);
  G4_tree->SetBranchAddress("Prim_Ph",&Prim_Ph);
  G4_tree->SetBranchAddress("Prim_x",&Prim_x);
  G4_tree->SetBranchAddress("Prim_y",&Prim_y);
  G4_tree->SetBranchAddress("Prim_z",&Prim_z);
  G4_tree->SetBranchAddress("Prim_pdg",&Prim_pdg);
  G4_tree->SetBranchAddress("Flux_Nhits",&Flux_Nhits);
  G4_tree->SetBranchAddress("Flux_pdg",Flux_pdg);
  G4_tree->SetBranchAddress("Flux_id",Flux_id);
  G4_tree->SetBranchAddress("Flux_x",Flux_x);
  G4_tree->SetBranchAddress("Flux_y",Flux_y);
  G4_tree->SetBranchAddress("Flux_z",Flux_z);
  G4_tree->SetBranchAddress("Flux_E",Flux_E);
  G4_tree->SetBranchAddress("Flux_tid",Flux_tid);
  G4_tree->SetBranchAddress("Flux_pid",Flux_pid);
  G4_tree->SetBranchAddress("Flux_th",Flux_th);
  G4_tree->SetBranchAddress("Flux_ph",Flux_ph);

  // ---------------------------------------------------------------------------
  // Get Left SIMC tree
  // ---------------------------------------------------------------------------
  
  TFile *fLSC     = new TFile("extendedl.root");
  TTree* LSC_tree = (TTree*)fLSC->Get("h1");

  Float_t         L_xfp;
  Float_t         L_yfp;
  Float_t         L_xpfp;
  Float_t         L_ypfp;
  Float_t         L_ytari;
  Float_t         L_deltai;
  Float_t         L_yptari;
  Float_t         L_xptari;
  Float_t         L_ytar;
  Float_t         L_delta;
  Float_t         L_yptar;
  Float_t         L_xptar;
  Float_t         L_good;
  
  LSC_tree->SetBranchAddress("hsxfp",&L_xfp);
  LSC_tree->SetBranchAddress("hsyfp",&L_yfp);
  LSC_tree->SetBranchAddress("hsxpfp",&L_xpfp);
  LSC_tree->SetBranchAddress("hsypfp",&L_ypfp);
  LSC_tree->SetBranchAddress("hsytari",&L_ytari);
  LSC_tree->SetBranchAddress("hsdeltai",&L_deltai);
  LSC_tree->SetBranchAddress("hsyptari",&L_yptari);
  LSC_tree->SetBranchAddress("hsxptari",&L_xptari);
  LSC_tree->SetBranchAddress("hsytar",&L_ytar);
  LSC_tree->SetBranchAddress("hsdelta",&L_delta);
  LSC_tree->SetBranchAddress("hsyptar",&L_yptar);
  LSC_tree->SetBranchAddress("hsxptar",&L_xptar);
  LSC_tree->SetBranchAddress("ok_spec",&L_good);

  // ---------------------------------------------------------------------------
  // Get Right SIMC tree
  // ---------------------------------------------------------------------------

  TFile *fRSC     = new TFile("extendedr.root");
  TTree* RSC_tree = (TTree*)fRSC->Get("h1");

  Float_t         R_xfp;
  Float_t         R_yfp;
  Float_t         R_xpfp;
  Float_t         R_ypfp;
  Float_t         R_ytari;
  Float_t         R_deltai;
  Float_t         R_yptari;
  Float_t         R_xptari;
  Float_t         R_ytar;
  Float_t         R_delta;
  Float_t         R_yptar;
  Float_t         R_xptar;
  Float_t         R_good;
  
  RSC_tree->SetBranchAddress("hsxfp",&R_xfp);
  RSC_tree->SetBranchAddress("hsyfp",&R_yfp);
  RSC_tree->SetBranchAddress("hsxpfp",&R_xpfp);
  RSC_tree->SetBranchAddress("hsypfp",&R_ypfp);
  RSC_tree->SetBranchAddress("hsytari",&R_ytari);
  RSC_tree->SetBranchAddress("hsdeltai",&R_deltai);
  RSC_tree->SetBranchAddress("hsyptari",&R_yptari);
  RSC_tree->SetBranchAddress("hsxptari",&R_xptari);
  RSC_tree->SetBranchAddress("hsytar",&R_ytar);
  RSC_tree->SetBranchAddress("hsdelta",&R_delta);
  RSC_tree->SetBranchAddress("hsyptar",&R_yptar);
  RSC_tree->SetBranchAddress("hsxptar",&R_xptar);
  RSC_tree->SetBranchAddress("ok_spec",&R_good);

  if( G4_tree->GetEntries() != LSC_tree->GetEntries() && LSC_tree->GetEntries() != RSC_tree->GetEntries() )
    std::cerr << "ERROR: number of entries in the tree do not match" << std::endl;

  Long64_t nentries = G4_tree->GetEntries();
  
  // ---------------------------------------------------------------------------

  for(Long64_t i=0; i<nentries;i++) {

    G4_tree->GetEntry(i);
    LSC_tree->GetEntry(i);
    RSC_tree->GetEntry(i);
    
    for( Int_t j=0; j<Flux_Nhits; j++) {
    }
    
  }

}
