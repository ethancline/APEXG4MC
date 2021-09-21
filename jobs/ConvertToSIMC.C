#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

#include <fstream>

const Float_t phrs   = 1.063;  // GeV/c 
const Float_t thhrs  = 12.5;   // deg
const Float_t q1_ent = 1720.5; // mm

void ConvertToSIMC( TString fname="batch_apex.root" )
{
  TFile *f = new TFile(fname);
  TTree* T = (TTree*)gDirectory->Get("T");

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

  T->SetBranchAddress("Prim_E",&Prim_E);
  T->SetBranchAddress("Prim_Th",&Prim_Th);
  T->SetBranchAddress("Prim_Ph",&Prim_Ph);
  T->SetBranchAddress("Prim_x",&Prim_x);
  T->SetBranchAddress("Prim_y",&Prim_y);
  T->SetBranchAddress("Prim_z",&Prim_z);
  T->SetBranchAddress("Prim_pdg",&Prim_pdg);
  T->SetBranchAddress("Flux_Nhits",&Flux_Nhits);
  T->SetBranchAddress("Flux_pdg",Flux_pdg);
  T->SetBranchAddress("Flux_id",Flux_id);
  T->SetBranchAddress("Flux_x",Flux_x);
  T->SetBranchAddress("Flux_y",Flux_y);
  T->SetBranchAddress("Flux_z",Flux_z);
  T->SetBranchAddress("Flux_E",Flux_E);
  T->SetBranchAddress("Flux_tid",Flux_tid);
  T->SetBranchAddress("Flux_pid",Flux_pid);
  T->SetBranchAddress("Flux_th",Flux_th);
  T->SetBranchAddress("Flux_ph",Flux_ph);
  
  Long64_t nentries = T->GetEntries();
  
  ofstream outfileL("LQ1_TCS.dat");
  ofstream outfileR("RQ1_TCS.dat");   

  TH2F* hLTp = new TH2F("hLTp","",100,-100, 100, 100, -100, 100);
  TH2F* hLTa = new TH2F("hLTa","",100,-100, 100, 100, -100, 100);
  TH2F* hRTp = new TH2F("hRTp","",100,-100, 100, 100, -100, 100);
  TH2F* hRTa = new TH2F("hRTa","",100,-100, 100, 100, -100, 100);

  TH1F* phi = new TH1F("phi","",100,-10,10);

  for(Long64_t i=0; i<nentries;i++) {
    T->GetEntry(i);

    Float_t delta     = 9999.;
    Float_t Ldydz_tcs = 9999.;
    Float_t Ldxdz_tcs = 9999.;
    Float_t Rdydz_tcs = 9999.;
    Float_t Rdxdz_tcs = 9999.;
    Float_t LxQ1_tcs  = 9999.;
    Float_t LyQ1_tcs  = 9999.;
    Float_t RxQ1_tcs  = 9999.;
    Float_t RyQ1_tcs  = 9999.;;
    
    for( Int_t j=0; j<Flux_Nhits; j++) {

      // calculate delta
      Float_t m     = 0.511/1000.;
      Float_t E     = Flux_E[j]/1000. + m;
      Float_t p     = TMath::Sqrt(E*E - m*m);
      delta         = (p-phrs)/phrs*100.;
      
      Float_t px    = p*TMath::Sin(Flux_th[j])*TMath::Cos(Flux_ph[j]);
      Float_t py    = p*TMath::Sin(Flux_th[j])*TMath::Sin(Flux_ph[j]);
      Float_t pz    = p*TMath::Cos(Flux_th[j]);
      
      TVector3 plab;
      TVector3 poslab;
      if(Flux_id[j]==3)
        phi->Fill(Flux_th[j]*180/3.14);
      if( Flux_id[j]  == 3 && Flux_pdg[j] == 11 ) { // look only at electrons in LHRS Q1
	
	// convert lab angles to TCS angles
	plab.SetXYZ( px, py, pz );
	plab.RotateY( -12.5 *TMath::DegToRad() );
	Ldydz_tcs = TMath::ATan( plab.X()/plab.Z() );
	Ldxdz_tcs = TMath::ATan( plab.Y()/plab.Z() );
	
	// convert x and y positions in lab at Q1 to TCS (x is dispersive)
	poslab.SetXYZ( Flux_x[j], Flux_y[j], Flux_z[j] );
	poslab.RotateY( -12.5 *TMath::DegToRad() );
	LxQ1_tcs = 0.1*(poslab.Y() + (q1_ent - poslab.Z())*Ldxdz_tcs);
	LyQ1_tcs = 0.1*(poslab.X() + (q1_ent - poslab.Z())*Ldydz_tcs);
	
	hLTp->Fill(LyQ1_tcs, LxQ1_tcs);
	hLTa->Fill(Ldydz_tcs*1000., Ldxdz_tcs*1000.);

      }
      else if( Flux_id[j]  == 4 && Flux_pdg[j] == -11 ) { // look only at positrons in RHRS Q1

	// convert lab angles to TCS angles
	plab.SetXYZ( px, py, pz );
	plab.RotateY( 12.5 *TMath::DegToRad() );
	Rdydz_tcs = TMath::ATan( plab.X()/plab.Z() );
	Rdxdz_tcs = TMath::ATan( plab.Y()/plab.Z() );
	
	// convert x and y positions in lab at Q1 to TCS (x is dispersive)
	poslab.SetXYZ( Flux_x[j], Flux_y[j], Flux_z[j] );
	poslab.RotateY( 12.5 *TMath::DegToRad() );
	RxQ1_tcs = 0.1*(poslab.Y() + (q1_ent - poslab.Z())*Rdxdz_tcs);
	RyQ1_tcs = 0.1*(poslab.X() + (q1_ent - poslab.Z())*Rdydz_tcs);
	
	hRTp->Fill(RyQ1_tcs, RxQ1_tcs);
	hRTa->Fill(Rdydz_tcs*1000., Rdxdz_tcs*1000.);
      
      }
    }
    
    outfileL <<" 172.05\t"<< delta << "\t" << 1000.*Ldydz_tcs << "\t" << 1000.*Ldxdz_tcs << "\t" << LyQ1_tcs << "\t" << LxQ1_tcs << endl;
    outfileR <<" 172.05\t"<< delta << "\t" << 1000.*Rdydz_tcs << "\t" << 1000.*Rdxdz_tcs << "\t" << RyQ1_tcs << "\t" << RxQ1_tcs << endl;

  }

  TCanvas* c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,2);
  c1->cd(1);
  hLTp->Draw("colz");
  c1->cd(2);
  hLTa->Draw("colz");
  c1->cd(3);
  hRTp->Draw("colz");
  c1->cd(4);
  hRTa->Draw("colz");



  outfileL.close();
  outfileR.close();
}
