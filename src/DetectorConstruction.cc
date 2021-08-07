#include "DetectorConstruction.hh"
#include "FluxSD.hh"
#include "EnergyDepositSD.hh"
#include "DetectorMessenger.hh"
#include "GlobalField.hh"

#include "G4Material.hh"
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"

#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4ClassicalRK4.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4GenericTrap.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"
#include "G4TwoVector.hh"

using namespace CLHEP;
using namespace std;

//---------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction()
{
  fNistManager  = G4NistManager::Instance();
  fDetMessenger = new DetectorMessenger(this);

  fHRSAngle     = 12.5;
  fDistTarPivot = 105.29;
  fDistPivotQ1  = 159.03;

  fHRSMomentum   = 1.1;
  fScaleSeptum  = 1.0;
  fFieldMapFile = "Septa-JB_map.table";
    
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/control/execute macros/DetectorSetup.mac";
  UI->ApplyCommand(command);

}

//---------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() 
{
}

//---------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  double LarmStepLimit=2.000 * mm; //why
  G4UserLimits* LarmStepLimits = new G4UserLimits(LarmStepLimit);
  G4int nonSDcounter = 100;
  G4int SDcounter    = 0;

  //---------------------------------------------------------------------------
  // Set up magnetic field
  //---------------------------------------------------------------------------

  G4MagneticField *fMagField = new GlobalField(fHRSMomentum, fScaleSeptum, fFieldMapFile);
  G4FieldManager *fm = new G4FieldManager(fMagField);
  G4Mag_SpinEqRhs* fBMTequation = new G4Mag_SpinEqRhs(fMagField);
  G4MagIntegratorStepper *pStepper = new G4ClassicalRK4(fBMTequation,12);
  G4ChordFinder *cftemp = new G4ChordFinder(fMagField, 0.01*mm, pStepper);
  fm->SetChordFinder(cftemp);

  //---------------------------------------------------------------------------
  // Define Materials
  //---------------------------------------------------------------------------
  
  G4Element*  N   = fNistManager->FindOrBuildElement(7);
  G4Element*  O   = fNistManager->FindOrBuildElement(8);
  
  G4Material* Beamline = new G4Material("Beam", 1.e-5*g/cm3, 2, kStateGas, STP_Temperature, 2.e-2*bar );
  Beamline->AddElement(N, 0.7);
  Beamline->AddElement(O, 0.3);

  //---------------------------------------------------------------------------
  // Create Experimental Hall
  //---------------------------------------------------------------------------

  G4Box* expHall_box           = new G4Box("expHall_box",
					   2.5 *m, 1.0 *m, 4.0 *m );

  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     fNistManager->FindOrBuildMaterial("G4_AIR"),
						     "expHall_log", 0, 0, 0);

  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, nonSDcounter);

  //---------------------------------------------------------------------------
  // Create Scattering Chamber
  //---------------------------------------------------------------------------

  G4double scat_inrad    = 602.25 *mm;
  G4double scat_outrad   = 603.758 *mm;
  G4double scat_height   = 190.5   *mm;

  G4RotationMatrix* scat_rm  = new G4RotationMatrix();
  scat_rm->rotateX(90. *deg);
  
  //---------------------------------------------------------------------------

  G4Tubs* scatv_tubs         = new G4Tubs("scatv_tubs",
					  0.0, scat_inrad, scat_height,
					  0.0, 360.0 *deg );
  
  G4LogicalVolume* scatv_log = new G4LogicalVolume(scatv_tubs,
						   Beamline,
						   "scatv_log", 0, 0, 0);
  
  new G4PVPlacement(scat_rm, G4ThreeVector(0,0,-fDistTarPivot*cm),scatv_log, "scatv", expHall_log, false, nonSDcounter++);
  
  //---------------------------------------------------------------------------
  
  G4Tubs* scat_tubs         = new G4Tubs("scat_tubs",
					 scat_inrad, scat_outrad, scat_height,
					 0.0, 360.0 *deg );
  
  G4LogicalVolume* scat_log = new G4LogicalVolume(scat_tubs,
						  fNistManager->FindOrBuildMaterial("G4_Al"),
						  "scat_log", 0, 0, 0);
  
  new G4PVPlacement(scat_rm, G4ThreeVector(0,0,-fDistTarPivot*cm),scat_log, "scat", expHall_log, false, nonSDcounter++);

  //--------------------------------------------------------------------------- 
  // Create Septum 
  //--------------------------------------------------------------------------- 
  double inch           = 2.54*cm;
  double y_en = 2.44*inch;
  double y_ex = 4.7 *inch;
  double zlength        = 173.939*cm;
  double z_sep_cen      = zlength*tan(5.*deg)/tan(12.5*deg);
  double z_real_tar     = z_sep_cen-zlength;
  double z_sept_en_min1 = z_real_tar + (17.31*inch+22.50*inch);
  double z_sept_en_max1 = z_sept_en_min1 - 3.15*inch*tan(5.*deg);
  double x_cen_sep_en   = 3.966*inch;
  double x_width_sep_en = 3.15*inch;
  double xmin_sep_en1   = x_cen_sep_en -x_width_sep_en*0.5*cos(5.*deg);
  double xmax_sep_en1   = x_cen_sep_en +x_width_sep_en*0.5*cos(5.*deg);
  double ang_en_min_1   = 5.*deg;
  double ang_en_max_1   = 5.*deg;
  double ang_en_min_2   = 6.6*deg;
  double ang_en_max_2   = 10.8*deg;
  double length_max_1   = 50.19* cm; // 19.76*inch
  double length_min_1   = 52.1 * cm; // 19.76*inch
  double length_min_2   = 59.82* cm; // 23.55*inch
  double length_max_2   = 60.3 * cm; // 23.74*inch
  
  double xmin_sep_ex1   = xmin_sep_en1 + length_min_1 * sin(ang_en_min_1);
  double xmax_sep_ex1   = xmax_sep_en1 + length_max_1 * sin(ang_en_max_1);
  double z_sept_ex_min1 = z_sept_en_min1 + length_min_1 * cos(ang_en_min_1);
  double z_sept_ex_max1 = z_sept_en_max1 + length_max_1 * cos(ang_en_max_1);

  double xmin_sep_en2   = xmin_sep_ex1;
  double xmax_sep_en2   = xmax_sep_ex1;
  double z_sept_en_min2 = z_sept_ex_min1;
  double z_sept_en_max2 = z_sept_ex_max1;

  double xmin_sep_ex2   = xmin_sep_en2 + length_min_2 * sin(ang_en_min_2);
  double xmax_sep_ex2   = xmax_sep_en2 + length_max_2 * sin(ang_en_max_2);
  double z_sept_ex_min2 = z_sept_en_min2 + length_min_2 * cos(ang_en_min_2);
  double z_sept_ex_max2 = z_sept_en_max2 + length_max_2 * cos(ang_en_max_2);

  double ymin_sep_ex1 = y_en + (z_sept_ex_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_ex1 = y_en + (z_sept_ex_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymin_sep_en1 = y_en + (z_sept_en_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_en1 = y_en + (z_sept_en_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);


  double ymin_sep_ex2 = y_en + (z_sept_ex_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_ex2 = y_en + (z_sept_ex_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymin_sep_en2 = y_en + (z_sept_en_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_en2 = y_en + (z_sept_en_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);

  G4double  pDz_1=0.31*inch/2.;
  G4double  pDz_2=0.25*inch/2.;
  
  //---------------------------------LEFT SEPTUM HALF----------------------------------------------
  //LHS Back Right Vertices
  vector<G4TwoVector> vertices_en1_min;
  vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
  vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
  vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
  vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );
  vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
  vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
  vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
  vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );
 

  //LHS Back Right Trap
  G4VSolid* testTrap_en1_min = new G4GenericTrap("testTrap_en1_min", pDz_1, vertices_en1_min);
  G4LogicalVolume* trrap_en1_min = new G4LogicalVolume(testTrap_en1_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_min",0,0,LarmStepLimits);
  trrap_en1_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en1_min=new G4RotationMatrix();
  pRotY90deg_en1_min->rotateY(90.*deg-ang_en_min_1);
  new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(pDz_1+xmin_sep_en1-0.02*mm,0,z_sept_en_min1), 
		    trrap_en1_min,"trrap_en1_min",expHall_log,0,nonSDcounter++);


  //LHS Back Left Vertices
  vector<G4TwoVector> vertices_en1_max;
  vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
  vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
  vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
  vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );
  vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
  vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
  vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
  vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );

  //LHS Back Left Trap
  G4VSolid* testTrap_en1_max = new G4GenericTrap("testTrap_en1_max", pDz_1, vertices_en1_max);
  G4LogicalVolume* trrap_en1_max = new G4LogicalVolume(testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_max",0,0,LarmStepLimits);
  trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en1_max=new G4RotationMatrix();
  pRotY90deg_en1_max->rotateY(90.*deg-ang_en_min_1);
  new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1+xmax_sep_en1+0.02*mm,0,z_sept_en_max1), 
		    trrap_en1_max,"trrap_en1_max",expHall_log,0,nonSDcounter++);


  //LHS Front Left Vertices
  vector<G4TwoVector> vertices_en2_max;
  vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
  vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
  vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
  vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );
  vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
  vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
  vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
  vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );


  //LHS Front Left Trap
  G4VSolid* testTrap_en2_max = new G4GenericTrap("testTrap_en2_max", pDz_1, vertices_en2_max);
  G4LogicalVolume* trrap_en2_max = new G4LogicalVolume(testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_max",0,0,LarmStepLimits);
  trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en2_max=new G4RotationMatrix();
  pRotY90deg_en2_max->rotateY(90.*deg-ang_en_max_2);
  new G4PVPlacement(pRotY90deg_en2_max,G4ThreeVector(pDz_1+xmax_sep_en2+0.02*mm,0,z_sept_en_max2), 
		    trrap_en2_max,"trrap_en2_max",expHall_log,0,nonSDcounter++);


  //LHS Front Right Vertices
  vector<G4TwoVector> vertices_en2_min;
  vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
  vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
  vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
  vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );
  vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
  vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
  vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
  vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );


  //LHS Front Right Trap
  G4VSolid* testTrap_en2_min = new G4GenericTrap("testTrap_en2_min", pDz_1, vertices_en2_min);
  G4LogicalVolume* trrap_en2_min = new G4LogicalVolume(testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_min",0,0,LarmStepLimits);
  trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en2_min=new G4RotationMatrix();
  pRotY90deg_en2_min->rotateY(90.*deg-ang_en_min_2);
  new G4PVPlacement(pRotY90deg_en2_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en2-0.02*mm,0,z_sept_en_min2), 
		    trrap_en2_min,"trrap_en2_min",expHall_log,0,nonSDcounter++);



  //LHS Back Top and Bottom Vertices
  vector<G4TwoVector> vertices_cov1_up;
  vertices_cov1_up.push_back( G4TwoVector(xmin_sep_en1-1.95*pDz_1, z_sept_en_min1) );
  vertices_cov1_up.push_back( G4TwoVector(xmax_sep_en1+1.95*pDz_1, z_sept_en_max1) );
  vertices_cov1_up.push_back( G4TwoVector(xmax_sep_ex1+1.95*pDz_1, z_sept_ex_max1) );
  vertices_cov1_up.push_back( G4TwoVector(xmin_sep_ex1-1.95*pDz_1, z_sept_ex_min1) );
  vertices_cov1_up.push_back( G4TwoVector(xmin_sep_en1-1.95*pDz_1, z_sept_en_min1) );
  vertices_cov1_up.push_back( G4TwoVector(xmax_sep_en1+1.95*pDz_1, z_sept_en_max1) );
  vertices_cov1_up.push_back( G4TwoVector(xmax_sep_ex1+1.95*pDz_1, z_sept_ex_max1) );
  vertices_cov1_up.push_back( G4TwoVector(xmin_sep_ex1-1.95*pDz_1, z_sept_ex_min1) );


  //LHS Back Top Trap
  G4VSolid* testTrap_cov_1_up = new G4GenericTrap("testTrap_cov_1_up", pDz_2, vertices_cov1_up);
  G4LogicalVolume* trrap_cov_1_up = new G4LogicalVolume(testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_up",0,0,LarmStepLimits);
  trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_1_up=new G4RotationMatrix();
  pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(pRotX90deg_cov_1_up,G4ThreeVector(0,0.02*mm+ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.), 
		    trrap_cov_1_up,"cov1_up",expHall_log,0,nonSDcounter++);



  //LHS Back Bottom Trap
  G4VSolid* testTrap_cov_1_bot = testTrap_cov_1_up;
  G4LogicalVolume* trrap_cov_1_bot = new G4LogicalVolume(testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_bot",0,0,LarmStepLimits);
  trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_1_bot=new G4RotationMatrix();
  pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(pRotX90deg_cov_1_bot,G4ThreeVector(0,-0.02*mm-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.), 
		    trrap_cov_1_bot,"cov1_bot",expHall_log,0,nonSDcounter++);
  trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));


  //LHS Front Top and Bottom Vertices
  vector<G4TwoVector> vertices_cov2_up;
  vertices_cov2_up.push_back( G4TwoVector(xmin_sep_en2-1.95*pDz_1, z_sept_en_min2) );
  vertices_cov2_up.push_back( G4TwoVector(xmax_sep_en2+1.95*pDz_1, z_sept_en_max2) );
  vertices_cov2_up.push_back( G4TwoVector(xmax_sep_ex2+1.95*pDz_1, z_sept_ex_max2) );
  vertices_cov2_up.push_back( G4TwoVector(xmin_sep_ex2-1.95*pDz_1, z_sept_ex_min2) );
  vertices_cov2_up.push_back( G4TwoVector(xmin_sep_en2-1.95*pDz_1, z_sept_en_min2) );
  vertices_cov2_up.push_back( G4TwoVector(xmax_sep_en2+1.95*pDz_1, z_sept_en_max2) );
  vertices_cov2_up.push_back( G4TwoVector(xmax_sep_ex2+1.95*pDz_1, z_sept_ex_max2) );
  vertices_cov2_up.push_back( G4TwoVector(xmin_sep_ex2-1.95*pDz_1, z_sept_ex_min2) );


  //LHS Front Top Trap
  G4VSolid* testTrap_cov_2_up = new G4GenericTrap("testTrap_cov_2_up", pDz_2, vertices_cov2_up);
  G4LogicalVolume* trrap_cov_2_up = new G4LogicalVolume(testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_up",0,0,LarmStepLimits);
  trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_2_up=new G4RotationMatrix();
  pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(pRotX90deg_cov_2_up,G4ThreeVector(0, 0.02*mm+ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.), 
		    trrap_cov_2_up,"cov2_up",expHall_log,0,nonSDcounter++);

  
  //RHS Front Bottom Trap
  G4VSolid* testTrap_cov_2_bot = testTrap_cov_2_up;
  G4LogicalVolume* trrap_cov_2_bot = new G4LogicalVolume(testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_bot",0,0,LarmStepLimits);
  trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_2_bot=new G4RotationMatrix();
  pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(pRotX90deg_cov_2_bot,G4ThreeVector(0,-0.02*mm-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.), 
		    trrap_cov_2_bot,"cov2_bot",expHall_log,0,nonSDcounter++);

  //----------------------------------RIGHT SEPTUM HALF----------------------------------------
  //RHS Back Left Vertices
  vector<G4TwoVector> r_vertices_en1_min;
  r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
  r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
  r_vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
  r_vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );
  r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1, -1.0*ymin_sep_en1) );
  r_vertices_en1_min.push_back( G4TwoVector(0.*length_min_1,  1.0*ymin_sep_en1) );
  r_vertices_en1_min.push_back( G4TwoVector(   length_min_1,  1.0*ymin_sep_ex1) );
  r_vertices_en1_min.push_back( G4TwoVector(   length_min_1, -1.0*ymin_sep_ex1) );


  //RHS Back Left Trap
  G4VSolid* r_testTrap = new G4GenericTrap("vac Box 1", pDz_1, r_vertices_en1_min);
  G4LogicalVolume* r_trrap = new G4LogicalVolume(r_testTrap,fNistManager->FindOrBuildMaterial("G4_Al"),"vac Box 1",0,0,LarmStepLimits);
  r_trrap->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en1_min=new G4RotationMatrix();
  r_pRotY90deg_en1_min->rotateY(90*deg+ang_en_max_1);
  new G4PVPlacement(r_pRotY90deg_en1_min,G4ThreeVector(1.*pDz_1-xmin_sep_en1,0,z_sept_en_min1), 
		    r_trrap,"vac Box en1 min",expHall_log,0,nonSDcounter++);


  //RHS Back Right Vertices
  vector<G4TwoVector> r_vertices_en1_max;
  r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
  r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
  r_vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
  r_vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );
  r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1, -1.0*ymax_sep_en1) );
  r_vertices_en1_max.push_back( G4TwoVector(0.*length_max_1,  1.0*ymax_sep_en1) );
  r_vertices_en1_max.push_back( G4TwoVector(   length_max_1,  1.0*ymax_sep_ex1) );
  r_vertices_en1_max.push_back( G4TwoVector(   length_max_1, -1.0*ymax_sep_ex1) );


  //RHS Back Right Trap
  G4VSolid* r_testTrap_en1_max = new G4GenericTrap("testTrap_en1_max", pDz_1, r_vertices_en1_max);
  G4LogicalVolume* r_trrap_en1_max = new G4LogicalVolume(r_testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_max",0,0,LarmStepLimits);
  r_trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en1_max=new G4RotationMatrix();
  r_pRotY90deg_en1_max->rotateY(90.*deg+ang_en_min_1);
  new G4PVPlacement(r_pRotY90deg_en1_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en1),0,z_sept_en_max1),
		    r_trrap_en1_max,"trrap_en1_max",expHall_log,0,nonSDcounter++);


  //RHS Front Right Vertices
  vector<G4TwoVector> r_vertices_en2_max;
  r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
  r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
  r_vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
  r_vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );
  r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2, -1.0*ymax_sep_en2) );
  r_vertices_en2_max.push_back( G4TwoVector(0.*length_max_2,  1.0*ymax_sep_en2) );
  r_vertices_en2_max.push_back( G4TwoVector(   length_max_2,  1.0*ymax_sep_ex2) );
  r_vertices_en2_max.push_back( G4TwoVector(   length_max_2, -1.0*ymax_sep_ex2) );


  //RHS Font Right Trap
  G4VSolid* r_testTrap_en2_max = new G4GenericTrap("testTrap_en2_max", pDz_1, vertices_en2_max);
  G4LogicalVolume* r_trrap_en2_max = new G4LogicalVolume(r_testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_max",0,0,LarmStepLimits);
  r_trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en2_max=new G4RotationMatrix();
  r_pRotY90deg_en2_max->rotateY(90.*deg+ang_en_max_2);
  new G4PVPlacement(r_pRotY90deg_en2_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en2),0,z_sept_en_max2),
		    r_trrap_en2_max,"trrap_en2_max",expHall_log,0,nonSDcounter++);


  //RHS Front Left Trap Vertices
  vector<G4TwoVector> r_vertices_en2_min;
  r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
  r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
  r_vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
  r_vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );
  r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2, -1.0*ymin_sep_en2) );
  r_vertices_en2_min.push_back( G4TwoVector(0.*length_min_2,  1.0*ymin_sep_en2) );
  r_vertices_en2_min.push_back( G4TwoVector(   length_min_2,  1.0*ymin_sep_ex2) );
  r_vertices_en2_min.push_back( G4TwoVector(   length_min_2, -1.0*ymin_sep_ex2) );


  //RHS Front Left Trap
  G4VSolid* r_testTrap_en2_min = new G4GenericTrap("testTrap_en2_min", pDz_1, vertices_en2_min);
  G4LogicalVolume* r_trrap_en2_min = new G4LogicalVolume(r_testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_min",0,0,LarmStepLimits);
  r_trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en2_min=new G4RotationMatrix();
  r_pRotY90deg_en2_min->rotateY(90.*deg+ang_en_min_2);
  new G4PVPlacement(r_pRotY90deg_en2_min,G4ThreeVector(pDz_1-xmin_sep_en2,0,z_sept_en_min2),
		    r_trrap_en2_min,"trrap_en2_min",expHall_log,0,nonSDcounter++);


  //RHS Back Top and Bottom Vertices
  vector<G4TwoVector> r_vertices_cov1_up;
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_en1-1.95*pDz_1), z_sept_en_min1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_en1+1.95*pDz_1), z_sept_en_max1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_ex1+1.95*pDz_1), z_sept_ex_max1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_ex1-1.95*pDz_1), z_sept_ex_min1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_en1-1.95*pDz_1), z_sept_en_min1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_en1+1.95*pDz_1), z_sept_en_max1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmax_sep_ex1+1.95*pDz_1), z_sept_ex_max1) );
  r_vertices_cov1_up.push_back( G4TwoVector(-1.*(xmin_sep_ex1-1.95*pDz_1), z_sept_ex_min1) );


  //RHS Back Top Trap
  G4VSolid* r_testTrap_cov_1_up = new G4GenericTrap("testTrap_cov_1_up", pDz_2, r_vertices_cov1_up);
  G4LogicalVolume* r_trrap_cov_1_up = new G4LogicalVolume(r_testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_up",0,0,LarmStepLimits);
  r_trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_1_up=new G4RotationMatrix();
  r_pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(r_pRotX90deg_cov_1_up,G4ThreeVector(0,ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.),
		    r_trrap_cov_1_up,"cov1_up",expHall_log,0,nonSDcounter++);

  
  //RHS Back Bottom Trap
  G4VSolid* r_testTrap_cov_1_bot = r_testTrap_cov_1_up;
  G4LogicalVolume* r_trrap_cov_1_bot = new G4LogicalVolume(r_testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_bot",0,0,LarmStepLimits);
  r_trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_1_bot=new G4RotationMatrix();
  r_pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(r_pRotX90deg_cov_1_bot,G4ThreeVector(0,-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.),
		    r_trrap_cov_1_bot,"cov1_bot",expHall_log,0,nonSDcounter++);
  r_trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));


  //RHS Front Top and Bottom Vertices
  vector<G4TwoVector> r_vertices_cov2_up;
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_en2-1.95*pDz_1), z_sept_en_min2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_en2+1.95*pDz_1), z_sept_en_max2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_ex2+1.95*pDz_1), z_sept_ex_max2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_ex2-1.95*pDz_1), z_sept_ex_min2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_en2-1.95*pDz_1), z_sept_en_min2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_en2+1.95*pDz_1), z_sept_en_max2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmax_sep_ex2+1.95*pDz_1), z_sept_ex_max2) );
  r_vertices_cov2_up.push_back( G4TwoVector(-1.*(xmin_sep_ex2-1.95*pDz_1), z_sept_ex_min2) );


  //RHS Front Top Trap
  G4VSolid* r_testTrap_cov_2_up = new G4GenericTrap("testTrap_cov_2_up", pDz_2, r_vertices_cov2_up);
  G4LogicalVolume* r_trrap_cov_2_up = new G4LogicalVolume(r_testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_up",0,0,LarmStepLimits);
  r_trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_up=new G4RotationMatrix();
  r_pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_up,G4ThreeVector(0,     ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_up,"cov2_up",expHall_log,0,nonSDcounter++);


  //RHS Front Bottom Trap
  G4VSolid* r_testTrap_cov_2_bot = r_testTrap_cov_2_up;
  G4LogicalVolume* r_trrap_cov_2_bot = new G4LogicalVolume(r_testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_bot",0,0,LarmStepLimits);
  r_trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_bot=new G4RotationMatrix();
  r_pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_bot,G4ThreeVector(0,-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_bot,"cov2_bot",expHall_log,0,nonSDcounter++);
  
  //------------------------------End of Septum----------------------------------------
  
  //--------------------------------------------------------------------------- 
  // Create Q1 Virtual Detectors
  //--------------------------------------------------------------------------- 
  
 
  G4Box* Q1_box           = new G4Box("Q1_box", 
				       0.20 *m, 0.20 *m, 0.05 *m ); 
   
  G4LogicalVolume* Q1_log = new G4LogicalVolume(Q1_box, 
						fNistManager->FindOrBuildMaterial("G4_AIR"), 
						"Q1_log", 0, 0, 0); 

  //--------------------------------------------------------------------------- 
    
  G4double LQ1_th           = fHRSAngle *deg; 
  G4double LQ1_d            = fDistPivotQ1 *cm; 
  G4double LQ1_xprime       = -LQ1_d * std::sin(LQ1_th); 
  G4double LQ1_zprime       = LQ1_d * std::cos(LQ1_th); 
  G4RotationMatrix* LQ1_rm  = new G4RotationMatrix(); 
  LQ1_rm->rotateY(LQ1_th); 

  fDetVol[0]                    = new G4PVPlacement(LQ1_rm, G4ThreeVector(LQ1_xprime,0.,LQ1_zprime), 
						    Q1_log, "LQ1", expHall_log, false, SDcounter);

  //--------------------------------------------------------------------------- 
  
  G4double RQ1_th           = -fHRSAngle *deg; 
  G4double RQ1_d            = fDistPivotQ1 *cm;
  G4double RQ1_xprime       = -RQ1_d * std::sin(RQ1_th); 
  G4double RQ1_zprime       = RQ1_d * std::cos(RQ1_th); 
  G4RotationMatrix* RQ1_rm  = new G4RotationMatrix(); 
  RQ1_rm->rotateY(RQ1_th); 
     
  fDetVol[1]                    = new G4PVPlacement(RQ1_rm, G4ThreeVector(RQ1_xprime,0.,RQ1_zprime), 
						    Q1_log, "RQ1", expHall_log, false, SDcounter++); 

  //---------------------------------------------------------------------------
  // Set Logical Attributes
  //---------------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fFluxSD = new FluxSD("FluxSD", fNSD);
  SDman->AddNewDetector( fFluxSD );
  Q1_log->SetSensitiveDetector( fFluxSD );

  //---------------------------------------------------------------------------

  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0)   );
  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  Q1_log->SetVisAttributes(blue);

  //---------------------------------------------------------------------------
  
  expHall_log->SetFieldManager(fm, true);

  //---------------------------------------------------------------------------

  return fExpHall;
}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------
