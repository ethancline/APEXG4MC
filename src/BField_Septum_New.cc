// ********************************************************************
// $Id: BField_Septum_new.cc,v 3.0, 2011/1/19  G2P Exp $
// Implementation of the BField_Septum_new class.
//
//////////////////////////////////////////////////////////////////////

#include "BField_Septum_New.hh"
#include "G4UImanager.hh"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include  <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"


//#define CREATE_MAP_NTUPLE 1
//#define BFIELD_SEPTUM_NEW_DEBUG 1

#ifdef BFIELD_SEPTUM_NEW_DEBUG
#include "GlobalDebuger.hh"
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BField_Septum_New::BField_Septum_New(double pMomentumL,double pMomentumR, const char *mapfile)
{

  //-------------------------------------------------------------
  // DJH Mod
  //-------------------------------------------------------------

  fLHRSMomentum = pMomentumL;
  fSeptumFieldScale = 1.05;
  
  //-------------------------------------------------------------

  ReadMap(mapfile);

}

//////////////////////////////////////////////////////////////////////

BField_Septum_New::~BField_Septum_New()
{
//	delete mRotL2F;
//	delete mRotF2L;
  delete mBField;
  delete Btxt;
  delete BCoord;
}

//////////////////////////////////////////////////////////////////////

void BField_Septum_New::ReadMap(const char *filename)
{
//    Btxt= new float[101][51][201][3];
//    cout<<"kas stex1?"<<endl;
    nx=101;
    ny=51;
    nz=201;
    Btxt   =new float ***[nx];
    BCoord =new float ***[nx];
    for(int i=0;i<101;i++)
    {
        Btxt[i]   =new float **[ny];
        BCoord[i] =new float **[ny];
        for (int j=0;j<51;j++)
        {
            Btxt[i][j]   =new float *[nz];
            BCoord[i][j] =new float *[nz];
            for (int k=0;k<201;k++)
            {
              Btxt[i][j][k]   =new float [0];
              BCoord[i][j][k] =new float [0];
              //initial the array
              for (int l=0;l<3;l++)
              {
                    Btxt[i][j][k][l] = 0.0;
                  BCoord[i][j][k][l] = 0.0;
              }
            }
        }
    }
    string line;
    ifstream myfile(filename,ios::in);
    if (myfile.is_open())
    {
      clog<<"reading septum new map"<<endl;
      while (! myfile.eof() )
      {
        getline(myfile,line);
        if (line.size() > 90)
        {
          float x, y, z;
          float bx, by, bz;
          istringstream vars(line);
          vars>>x>>y>>z>>bx>>by>>bz;
          int ix=int(x);
          int iy=int(y);
          int iz=int((z+200.)/2.);
//          cout<<setw(15)<<ix<<setw(15)<<iy<<setw(15)<<iz<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<setw(15)<<bx<<setw(15)<<by<<setw(15)<<bz<<endl;
          Btxt[ix][iy][iz][0]=bx;
          Btxt[ix][iy][iz][1]=by;
          Btxt[ix][iy][iz][2]=bz;
          BCoord[ix][iy][iz][0]=x;
          BCoord[ix][iy][iz][1]=y;
          BCoord[ix][iy][iz][2]=z;
        }
      }
      clog<<"septum new map has been read"<<endl;
    }
    else 
    {
      clog<<"Septum field is missing. we skip this field"<<endl;
    }

//    cout<<"kas stex3?"<<endl;
/*
    for (int i=0; i<101; i++)
    for (int j=0; j<51; j++)
    for (int k=0; k<201; k++)
    for (int l=0; l<3; l++)
    {
       Btxt[i][j][k][l]=i+j+k*l;
    }
    */
    cout<<"arrrray       "<<Btxt[29][9][99][0]<<endl;
    cout<<"arrrray       "<<Btxt[19][9][99][1]<<endl;
    cout<<"arrrray       "<<Btxt[29][9][99][2]<<endl;
    cout<<"arrrray       "<<BCoord[29][9][99][0]<<endl;
    cout<<"arrrray       "<<BCoord[19][9][99][1]<<endl;
    cout<<"arrrray       "<<BCoord[29][9][99][2]<<endl;
}

//////////////////////////////////////////////////////////////////////

void  BField_Septum_New::GetFieldValue(const double point[4],
				 double *Bfield ) const
{
  
  //-------------------------------------------------------------
  // DJH Mod
  //-------------------------------------------------------------

  double fPos[3];
  double fB[3];

  fPos[0] = point[0];
  fPos[1] = point[1];
  fPos[2] = point[2];
  
  //-------------------------------------------------------------

//    double sep_cent= 69.612752*cm;  //distance from target to septum is 175 cm 
    double sep_cent= 68.6457*cm;  //distance from target to septum is 173.939 cm
    double x_loc=fabs(fPos[0]);
    double y_loc=fabs(fPos[1]);
    double z_loc=sep_cent-fPos[2];

    {
//        double q1_center = 207.0649*cm;
//        double q1_entry = 160.0*cm;
//        double q1_exit = 254.13*cm;
//        double ztmp_quad1=zn-q1_center;
        int ix=int(x_loc/cm);
        int iy=int(y_loc/cm);
        int iz=int((z_loc/cm + 200)/2.);

//        cout<<"locxyz = "<<ix<<setw(15)<<iy<<setw(15)<<iz<<setw(15)<<x_loc/cm<<setw(15)<<y_loc/cm<<setw(15)<<z_loc/cm<<endl;

        if ((ix>=0) && (ix<nx-1) && (iy>=0) && (iy<ny-1) && (iz>=0) && (iz<nz-1))
        {
            double x1=BCoord[ix][iy][iz][0]*cm;
            double x2=BCoord[ix+1][iy][iz][0]*cm;
            double y1=BCoord[ix][iy][iz][1]*cm;
            double y2=BCoord[ix][iy+1][iz][1]*cm;
            double z1=BCoord[ix][iy][iz][2]*cm;
            double z2=BCoord[ix][iy][iz+1][2]*cm;


            double Bx_y1z1 = Btxt[ix][iy][iz][0] + ((Btxt[ix+1][iy][iz][0]-Btxt[ix][iy][iz][0])/(x2-x1))*(x_loc - x1);
//            cout<<"Bx_y1z1="<<Bx_y1z1<<endl;
            double Bx_y2z1 = Btxt[ix][iy+1][iz][0] + ((Btxt[ix+1][iy+1][iz][0]-Btxt[ix][iy+1][iz][0])/(x2-x1))*(x_loc - x1);
//            cout<<"Bx_y2z1="<<Bx_y2z1<<endl;
            double Bx_z1   = Bx_y1z1 + (y_loc-y1)*(Bx_y2z1-Bx_y1z1)/(y2-y1);
//            cout<<"Bx_z1="<<Bx_z1<<endl;

            double Bx_y1z2 = Btxt[ix][iy][iz+1][0] + ((Btxt[ix+1][iy][iz+1][0]-Btxt[ix][iy][iz+1][0])/(x2-x1))*(x_loc - x1);
//            cout<<"Bx_y1z2="<<Bx_y1z2<<endl;
            double Bx_y2z2 = Btxt[ix][iy+1][iz+1][0] + ((Btxt[ix+1][iy+1][iz+1][0]-Btxt[ix][iy+1][iz+1][0])/(x2-x1))*(x_loc - x1);
//            cout<<"Bx_y2z2="<<Bx_y2z2<<endl;
            double Bx_z2   = Bx_y1z2 + (y_loc-y1)*(Bx_y2z2-Bx_y1z2)/(y2-y1);
//            cout<<"Bx_z2="<<Bx_z2<<endl;

            double Bx_z = Bx_z1 + (z_loc - z1)*(Bx_z2-Bx_z1)/(z2-z1);
//            cout<<"Bx_z="<<Bx_z<<endl;



//             cout<<"coord x  "<<x_q1[i111]<<"     "<<x_q1[i121]<<"      "<<x_q1[i211]<<"      "<<x_q1[i221]<<endl;
//             cout<<"coord y  "<<y_q1[i111]<<"     "<<y_q1[i121]<<"      "<<y_q1[i211]<<"      "<<y_q1[i221]<<endl;
//             cout<<"field   "<<Bxq1[i111]<<"     "<<Bxq1[i121]<<"      "<<Bxq1[i211]<<"      "<<Bxq1[i221]<<endl;
//             cout<<" temp 1 "<<Bx_y1z1<<"        "<<Bx_y2z1<<"    "<<(xn/cm - x_q1[i111])<<endl;
//             cout<<"field Bx= "<<Bx_z1<<endl;



            double By_y1z1 = Btxt[ix][iy][iz][1] + ((Btxt[ix+1][iy][iz][1]-Btxt[ix][iy][iz][1])/(x2-x1))*(x_loc - x1);
            double By_y2z1 = Btxt[ix][iy+1][iz][1] + ((Btxt[ix+1][iy+1][iz][1]-Btxt[ix][iy+1][iz][1])/(x2-x1))*(x_loc - x1);
            double By_z1   = By_y1z1 + (y_loc-y1)*(By_y2z1-By_y1z1)/(y2-y1);

            double By_y1z2 = Btxt[ix][iy][iz+1][1] + ((Btxt[ix+1][iy][iz+1][1]-Btxt[ix][iy][iz+1][1])/(x2-x1))*(x_loc - x1);
            double By_y2z2 = Btxt[ix][iy+1][iz+1][1] + ((Btxt[ix+1][iy+1][iz+1][1]-Btxt[ix][iy+1][iz+1][1])/(x2-x1))*(x_loc - x1);
            double By_z2   = By_y1z2 + (y_loc-y1)*(By_y2z2-By_y1z2)/(y2-y1);

            double By_z = By_z1 + (z_loc - z1)*(By_z2-By_z1)/(z2-z1);



            double Bz_y1z1 = Btxt[ix][iy][iz][2] + ((Btxt[ix+1][iy][iz][2]-Btxt[ix][iy][iz][2])/(x2-x1))*(x_loc - x1);
            double Bz_y2z1 = Btxt[ix][iy+1][iz][2] + ((Btxt[ix+1][iy+1][iz][2]-Btxt[ix][iy+1][iz][2])/(x2-x1))*(x_loc - x1);
            double Bz_z1   = Bz_y1z1 + (y_loc-y1)*(Bz_y2z1-Bz_y1z1)/(y2-y1);

            double Bz_y1z2 = Btxt[ix][iy][iz+1][2] + ((Btxt[ix+1][iy][iz+1][2]-Btxt[ix][iy][iz+1][2])/(x2-x1))*(x_loc - x1);
            double Bz_y2z2 = Btxt[ix][iy+1][iz+1][2] + ((Btxt[ix+1][iy+1][iz+1][2]-Btxt[ix][iy+1][iz+1][2])/(x2-x1))*(x_loc - x1);
            double Bz_z2   = Bz_y1z2 + (y_loc-y1)*(Bz_y2z2-Bz_y1z2)/(y2-y1);

            double Bz_z = Bz_z1 + (z_loc - z1)*(Bz_z2-Bz_z1)/(z2-z1);

            fB[0]=Bx_z/10000.*tesla;
            fB[1]=By_z/10000.*tesla;
            fB[2]=Bz_z/10000.*tesla;
            if (fPos[0] < 0 ) {fB[0]*= -1.;};
            if (fPos[1] < 0 ) {fB[0]*= -1.; fB[2]*= -1.;};
	    

	    //-------------------------------------------------------------
	    // DJH Mod
	    //-------------------------------------------------------------

	    if(fB[0] > 0.001 ) 
	      G4cout << "-----" << fPos[0] << "\t" <<  fPos[1] << "\t" <<  fPos[2] << "\t" << G4endl;

	    fB[0] *=2.2/2.140045;
	    fB[1] *=2.2/2.140045;
	    fB[2] *=2.2/2.140045;

	    Bfield[0] =  fSeptumFieldScale *fB[0] * (fLHRSMomentum/GeV)/2.2;
	    Bfield[1] =  fSeptumFieldScale *fB[1] * (fLHRSMomentum/GeV)/2.2;	    
	    Bfield[2] =  fSeptumFieldScale *fB[2] * (fLHRSMomentum/GeV)/2.2;


        }

    }
//    return true;
}
