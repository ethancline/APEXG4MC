#include "DetectorConstruction.hh"
#include "FluxSD.hh"
#include "EnergyDepositSD.hh"
#include "DetectorMessenger.hh"
#include "BField_Septum_New.hh"

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
#include "G4ChordFinder.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction()
{
  fNistManager  = G4NistManager::Instance();
  fDetMessenger = new DetectorMessenger(this);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/control/execute macros/DetectorSetup.mac";
  UI->ApplyCommand(command);


  G4FieldManager   *pFieldMgr;      
  fSeptumField = new BField_Septum_New( 2.2, 2.2, "Septa-JB_map.table" );
  pFieldMgr=G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4ChordFinder *pChordFinder = new G4ChordFinder(fSeptumField);
  pFieldMgr->SetChordFinder( pChordFinder );
  pFieldMgr->SetDetectorField(fSeptumField);

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

  //---------------------------------------------------------------------------
  // Define Materials
  //---------------------------------------------------------------------------
  
  G4Element*  N   = fNistManager->FindOrBuildElement(7);
  G4Element*  O   = fNistManager->FindOrBuildElement(8);
  
  G4Material* Air      = new G4Material("Air", 1.290*mg/cm3, 2 );
  Air->AddElement(N, 0.7);
  Air->AddElement(O, 0.3);

  G4Material* Beamline = new G4Material("Beam", 1.e-5*g/cm3, 2, kStateGas, STP_Temperature, 2.e-2*bar );
  Beamline->AddElement(N, 0.7);
  Beamline->AddElement(O, 0.3);

  //---------------------------------------------------------------------------
  // Create Experimental Hall
  //---------------------------------------------------------------------------

  G4Box* expHall_box           = new G4Box("expHall_box",
					   15. *m, 12.0 *m, 15. *m );

  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     fNistManager->FindOrBuildMaterial("G4_AIR"),
						     "expHall_log", 0, 0, 0);

  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, 0);

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
  
  new G4PVPlacement(scat_rm, G4ThreeVector(0,0,0),scatv_log, "scatv", expHall_log, false, 0);
  
  //---------------------------------------------------------------------------
  
  G4Tubs* scat_tubs         = new G4Tubs("scat_tubs",
					 scat_inrad, scat_outrad, scat_height,
					 0.0, 360.0 *deg );
  
  G4LogicalVolume* scat_log = new G4LogicalVolume(scat_tubs,
						  fNistManager->FindOrBuildMaterial("G4_Al"),
						  "scat_log", 0, 0, 0);
  
  new G4PVPlacement(scat_rm, G4ThreeVector(0,0,0),scat_log, "scat", expHall_log, false, 0);

  //---------------------------------------------------------------------------
  // Create Vertical Wire Targets
  //---------------------------------------------------------------------------
  
  G4double vwire_outrad   = 0.5 *mm;
  G4double vwire_height   = 10.0 *mm;
  
  G4Tubs* vwire_tubs         = new G4Tubs("vwire_tubs",
					  0.0, vwire_outrad, vwire_height,
					  0.0, 360.0 *deg );
  
  G4LogicalVolume* vwire_log = new G4LogicalVolume(vwire_tubs,
						   fNistManager->FindOrBuildMaterial("G4_Al"),
						   "vwire_log", 0, 0, 0);
  
  new G4PVPlacement(0, G4ThreeVector(0,-20.0*mm,0.0), vwire_log, "v1", scatv_log, false, 0); 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), vwire_log, "v2", scatv_log, false, 0);  
  new G4PVPlacement(0, G4ThreeVector(0,20.0*mm,0.0), vwire_log, "v3", scatv_log, false, 0); 

  //--------------------------------------------------------------------------- 
  // Create Left Septum
  //--------------------------------------------------------------------------- 
  
  G4double LSeptum_th           = 5.0 *deg; 
  G4double LSeptum_d            = 1.5 *m; 
  G4double LSeptum_xprime       = -LSeptum_d * std::sin(LSeptum_th); 
  G4double LSeptum_zprime       = LSeptum_d * std::cos(LSeptum_th); 
  G4RotationMatrix* LSeptum_rm  = new G4RotationMatrix(); 
  LSeptum_rm->rotateY(LSeptum_th); 
 
  //--------------------------------------------------------------------------- 
  
  G4Box* LSeptum_box           = new G4Box("LSeptum_box", 
					   0.05 *m, 0.20 *m, 0.3 *m ); 
  
  G4LogicalVolume* LSeptum_log = new G4LogicalVolume(LSeptum_box, 
						     fNistManager->FindOrBuildMaterial("G4_AIR"), 
						     "LSeptum_log", 0, 0, 0); 
  
  new G4PVPlacement(LSeptum_rm, G4ThreeVector(LSeptum_xprime,0.,LSeptum_zprime), LSeptum_log, "LSeptum", expHall_log, false, 0); 
 
  
  //--------------------------------------------------------------------------- 
  // Create Left Q1 Virtual Detector
  //--------------------------------------------------------------------------- 
  
  G4double LQ1_th           = 12.5 *deg; 
  G4double LQ1_d            = 2.5 *m; 
  G4double LQ1_xprime       = -LQ1_d * std::sin(LQ1_th); 
  G4double LQ1_zprime       = LQ1_d * std::cos(LQ1_th); 
  G4RotationMatrix* LQ1_rm  = new G4RotationMatrix(); 
  LQ1_rm->rotateY(LQ1_th); 
 
  //--------------------------------------------------------------------------- 
  
  G4Box* LQ1_box           = new G4Box("LQ1_box", 
				       0.20 *m, 0.20 *m, 0.05 *m ); 
   
  G4LogicalVolume* LQ1_log = new G4LogicalVolume(LQ1_box, 
						 fNistManager->FindOrBuildMaterial("G4_AIR"), 
						 "LQ1_log", 0, 0, 0); 
     
  fDetVol[0]                    = new G4PVPlacement(LQ1_rm, G4ThreeVector(LQ1_xprime,0.,LQ1_zprime), 
						    LQ1_log, "LQ1", expHall_log, false, 0); 

  //--------------------------------------------------------------------------- 
  // Create Right Septum
  //--------------------------------------------------------------------------- 
  
  G4double RSeptum_th           = -5.0 *deg; 
  G4double RSeptum_d            = 1.5 *m; 
  G4double RSeptum_xprime       = -RSeptum_d * std::sin(RSeptum_th); 
  G4double RSeptum_zprime       = RSeptum_d * std::cos(RSeptum_th); 
  G4RotationMatrix* RSeptum_rm  = new G4RotationMatrix(); 
  RSeptum_rm->rotateY(RSeptum_th); 
 
  //--------------------------------------------------------------------------- 
  
  G4Box* RSeptum_box           = new G4Box("RSeptum_box", 
					   0.05 *m, 0.20 *m, 0.3 *m ); 
  
  G4LogicalVolume* RSeptum_log = new G4LogicalVolume(RSeptum_box, 
						     fNistManager->FindOrBuildMaterial("G4_AIR"), 
						     "RSeptum_log", 0, 0, 0); 
  
  new G4PVPlacement(RSeptum_rm, G4ThreeVector(RSeptum_xprime,0.,RSeptum_zprime), RSeptum_log, "RSeptum", expHall_log, false, 0); 
  
  //--------------------------------------------------------------------------- 
  // Create Right Q1 Virtual Detector
  //--------------------------------------------------------------------------- 
  
  G4double RQ1_th           = -12.5 *deg; 
  G4double RQ1_d            = 2.5 *m; 
  G4double RQ1_xprime       = -RQ1_d * std::sin(RQ1_th); 
  G4double RQ1_zprime       = RQ1_d * std::cos(RQ1_th); 
  G4RotationMatrix* RQ1_rm  = new G4RotationMatrix(); 
  RQ1_rm->rotateY(RQ1_th); 
 
  //--------------------------------------------------------------------------- 
  
  G4Box* RQ1_box           = new G4Box("RQ1_box", 
				       0.20 *m, 0.20 *m, 0.05 *m ); 
   
  G4LogicalVolume* RQ1_log = new G4LogicalVolume(RQ1_box, 
						 fNistManager->FindOrBuildMaterial("G4_AIR"), 
 						 "RQ1_log", 0, 0, 0); 
     
  fDetVol[1]                    = new G4PVPlacement(RQ1_rm, G4ThreeVector(RQ1_xprime,0.,RQ1_zprime), 
						    RQ1_log, "RQ1", expHall_log, false, 0); 

  //---------------------------------------------------------------------------
  // Set Step Limits, Sensitive Detector and Visualisation
  //---------------------------------------------------------------------------

  G4double maxStep = 0.5 *mm;;
  G4UserLimits* stepLimit = new G4UserLimits(maxStep);

  //---------------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fFluxSD = new FluxSD("FluxSD", fNSD);
  SDman->AddNewDetector( fFluxSD );
  LQ1_log->SetSensitiveDetector( fFluxSD );
  RQ1_log->SetSensitiveDetector( fFluxSD );

  //---------------------------------------------------------------------------

  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0)   );
  G4VisAttributes* red     = new G4VisAttributes( G4Colour(1.0,0.0,0.0)   );

  LSeptum_log->SetVisAttributes(red);
  RSeptum_log->SetVisAttributes(red);
  LQ1_log->SetVisAttributes(blue);
  RQ1_log->SetVisAttributes(blue);

  //---------------------------------------------------------------------------

  return fExpHall;
}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------
