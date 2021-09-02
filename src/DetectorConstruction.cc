#include "DetectorConstruction.hh"
#include "FluxSD.hh"
#include "DetectorMessenger.hh"
#include "GlobalField.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh" 
#include "G4GenericTrap.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4ClassicalRK4.hh"

#include "G4UserLimits.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

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

  // Set detector parameters defaults
  fHRSAngle     = 12.5;
  fDistTarPivot = 105.29;  // from original APEX G4
  fDistPivotQ1  = 171.095; // Q1 SOS (from original APEX G4)
  fHRSMomentum  = 1.063;
  fScaleSeptum  = 1.05;
  fFieldMapFile = "Septa-JB_map.table"; // new septum geometry but not latest map
  fSieveOn      = false;
  fSieveAngle   = 5 *deg; // 5.81 degrees?
  fTarget       = "ProdW";
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/control/execute macros/DetectorSetup.mac";
  UI->ApplyCommand(command);
}

//---------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() 
{
  delete fDetMessenger;
}

//---------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  double LarmStepLimit=2.000 * mm; // why? -- from original APEX_G4
  G4UserLimits* LarmStepLimits = new G4UserLimits(LarmStepLimit);

  G4int nonSDcounter = 100;       // detector copy id counters
  G4int SDcounter    = 0;

  //---------------------------------------------------------------------------
  // Set up magnetic field
  //---------------------------------------------------------------------------

  G4MagneticField*        magField     = new GlobalField(fHRSMomentum, fScaleSeptum, fFieldMapFile);
  G4FieldManager*         fm           = new G4FieldManager(magField);
  G4Mag_SpinEqRhs*        BMTequation  = new G4Mag_SpinEqRhs(magField);
  G4MagIntegratorStepper* pStepper     = new G4ClassicalRK4(BMTequation,12);
  G4ChordFinder*          cftemp       = new G4ChordFinder(magField, 0.01*mm, pStepper);
  fm->SetChordFinder(cftemp);

  //---------------------------------------------------------------------------
  // Define non-NIST materials
  //---------------------------------------------------------------------------
  
  G4Element*  N   = fNistManager->FindOrBuildElement(7);
  G4Element*  O   = fNistManager->FindOrBuildElement(8);
  
  G4Material* Beamline = new G4Material("Beam", 1.e-5*g/cm3, 2, kStateGas, STP_Temperature, 2.e-2*bar );
  Beamline->AddElement(N, 0.7);
  Beamline->AddElement(O, 0.3);

  //---------------------------------------------------------------------------
  // Create experimental hall -- default to vacuum for now
  //---------------------------------------------------------------------------

  G4Box* expHall_box           = new G4Box("expHall_box",
					   2.5 *m, 1.0 *m, 4.0 *m );
  
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     Beamline,                    // DJH: default to vacuum for now
						     "expHall_log", 0, 0, 0);
  
  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, nonSDcounter);

  //--------------------------------------------------------------------------- 
  // Create upstream beamline (copy from "standard" (GMn) config in g4sbs)
  //--------------------------------------------------------------------------- 

  G4bool ChkOverlaps = false;
  G4double inch      = 2.54*cm;
  
  double sc_entbeampipeflange_dist = 25.375*2.54*cm; // entrance pipe flange distance from hall center 
  
  G4double ent_len = 2*m; 
  G4double ent_rin = 31.75*mm; 
  G4double ent_rou = ent_rin+0.120*mm; 
  
  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg ); 
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg ); 
  
  G4LogicalVolume *entLog    = new G4LogicalVolume(ent_tube,
						   fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), "ent_log", 0, 0, 0); 
  G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac,
						   Beamline, "entvac_log", 0, 0, 0); 
  
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist-fDistTarPivot*cm),
		    entLog, "ent_phys", expHall_log, false, ++nonSDcounter , ChkOverlaps); 
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist-fDistTarPivot*cm),
		    entvacLog, "entvac_phys", expHall_log,false, ++nonSDcounter , ChkOverlaps); 

  //---------------------------------------------------------------------------
  // Create scattering chamber (copy from "standard" (GMn) config in g4sbs, just the tank and vacuum for now)
  //---------------------------------------------------------------------------

  G4LogicalVolume *logicScatChamberTank = 0;
  G4LogicalVolume *logicScatChamber     = 0;

  // Scattering chamber tank:
  // basic volume:
  G4double SCHeight        = 44.75*inch;
  G4double SCRadius        = 20.0*inch;
  G4double SCTankThickness = 2.5*inch;
  G4double SCTankRadius    = SCRadius+SCTankThickness;
  G4double SCTankHeight    = SCHeight;
  G4double SCOffset        = 3.75*inch;
  
  G4Tubs* solidSCTank_0 = new G4Tubs("SCTank_0", SCRadius, SCTankRadius, 0.5*SCTankHeight, 0.0*deg, 360.0*deg);
  
  // exit flange:
  G4double SCExitFlangePlateHLength = 22.5*sin(25.5*atan(1)/45.0)*inch;
  G4double SCExitFlangePlateHeight  = 11.0*inch;
  G4double SCExitFlangePlateThick   = 1.25*inch;
  G4double SCExitFlangeHAngleApert  = atan(SCExitFlangePlateHLength/(SCTankRadius+SCExitFlangePlateThick));
  G4double SCExitFlangeMaxRad       = SCExitFlangePlateHLength/sin(SCExitFlangeHAngleApert);
  
  G4Tubs* solidSCExitFlangetubs =  new G4Tubs("SCExFlange_tubs", SCRadius, SCExitFlangeMaxRad, 
					      0.5*SCExitFlangePlateHeight, 0.0, 2.0*SCExitFlangeHAngleApert);
  
  G4RotationMatrix* rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHAngleApert);
  
  G4UnionSolid* solidSCTank_0_exft =  new G4UnionSolid("solidSCTank_0_exft", solidSCTank_0, solidSCExitFlangetubs,
						       rot_temp, G4ThreeVector(0,0,SCOffset));
  
  G4Box* ExitFlangeHeadCut = new G4Box("ExitFlangeHeadCut", 0.5*m, 0.5*m, 0.5*m); 
    
  G4SubtractionSolid* solidSCTank_0_exf = new G4SubtractionSolid("solidSCTank_0_exf", solidSCTank_0_exft, ExitFlangeHeadCut,
								 0, G4ThreeVector(-SCTankRadius-0.5*m,0,0));
  
  // exit flange hole:
  G4double SCExitFlangeHoleHeight    = 7.85*inch;
  G4double SCExitFlangeHoleAngleApert = 38.25*deg;
  
  G4Tubs* solidSCExFH = new G4Tubs("SCExFH", SCRadius-1.0*cm, SCTankRadius+1.5*inch,
				   0.5*SCExitFlangeHoleHeight, 0.0, SCExitFlangeHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCTank_0_exfh = new G4SubtractionSolid("solidSCTank_0_exfh", solidSCTank_0_exf, solidSCExFH,
								  rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // windows holes: 
  G4double SCWindowHeight      = 18.0*inch;
  G4double SCWindowAngleApert  = 149.0*deg;
  G4double SCWindowAngleOffset = 11.0*deg;
  
  G4Tubs* solidSCWindow = new G4Tubs("SCWindow", SCRadius-1.0*cm, SCTankRadius+1.0*cm,
				     0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wf = new G4SubtractionSolid("solidSCTank_0_wf", solidSCTank_0_exfh, solidSCWindow,
								rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wb = new G4SubtractionSolid("solidSCTank_0_wb", solidSCTank_0_wf, solidSCWindow,
								rot_temp, G4ThreeVector(0,0,SCOffset));

  //---------------------------------------------------------------------------  
  // Scattering chamber tank (aluminum)
  
  G4VSolid* solidScatChamberTank = solidSCTank_0_wb;
  logicScatChamberTank =   new G4LogicalVolume(solidScatChamberTank,
					       fNistManager->FindOrBuildMaterial("G4_Al"), "ScatChamberTank_log");
  
  G4RotationMatrix* rotSC = new G4RotationMatrix();
  rotSC->rotateX(-90.0*deg);
  
  G4ThreeVector* SCPlacement = new G4ThreeVector(0,-fDistTarPivot*cm,-SCOffset);
  SCPlacement->rotateX(90*deg);
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamberTank, "ScatChamberTankPhys", expHall_log, false, ++nonSDcounter, ChkOverlaps);

  
  //---------------------------------------------------------------------------  
  // Scattering chamber volume (vacuum)

  G4Tubs* solidExitBeamPipeHole     = new G4Tubs("solidBackViewPipeHole", 
						 0.0, 175.0*mm, 7.903*inch, 0.0, 360.0*deg);  // DJH: artificially increased radius for now
  
  G4Tubs* solidEntranceBeamPipeHole = new G4Tubs("solidEntranceBeamPipeHole",
						 0.0, 1.0*inch, 5.375*inch, 0.0, 360.0*deg);
  
  G4Tubs* solidScatChamber_0        = new G4Tubs("SC", 0.0, SCRadius, 0.5* SCHeight, 0.0*deg, 360.0*deg);
  
  G4Tubs* solidSCWindowVacuum       = new G4Tubs("SCWindowVacuumFront", SCRadius-1.0*cm, SCTankRadius,
						 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4UnionSolid* solidScatChamber_0_wbv = new G4UnionSolid("solidScatChamber_0_wbv", solidScatChamber_0, solidSCWindowVacuum,
							  rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);

  G4UnionSolid* solidScatChamber_0_entbp = new G4UnionSolid("solidScatChamber_0_entbp", solidScatChamber_0_wbv,
							     solidEntranceBeamPipeHole, 
							     rot_temp, G4ThreeVector(0, -SCRadius, SCOffset));

  G4UnionSolid* solidScatChamber_0_exbp  = new G4UnionSolid("solidScatChamber_0_exbp", solidScatChamber_0_entbp,
							    solidExitBeamPipeHole, 
							    rot_temp, G4ThreeVector(0, +SCRadius, SCOffset));
  
  G4VSolid* solidScatChamber = solidScatChamber_0_exbp;
  logicScatChamber = new G4LogicalVolume(solidScatChamber, Beamline, "ScatChamber_log");
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamber, "ScatChamberPhys",
		    expHall_log, false, ++nonSDcounter, ChkOverlaps);

  //--------------------------------------------------------------------------- 
  // Create APEX targets (from Silviu's CAD model -- need to add offsets, etc)
  //--------------------------------------------------------------------------- 

  G4int       ntargs;
  G4double    targwidth, targthick;
  G4Material* targMaterial;
  G4double    targzpos[10] = { 0.0 };
  G4double    targxpos[10] = { 0.0 };
  G4double    targypos[10] = { 0.0 };
  
  if( fTarget.compareTo("ProdW") == 0 )
    fTargetType = kProdW;
  else if( fTarget.compareTo("ProdC") == 0 )
    fTargetType = kProdC;
  else if( fTarget.compareTo("Optics1") == 0 )
    fTargetType = kOptics1;
  else if( fTarget.compareTo("Optics2") == 0 )
    fTargetType = kOptics2;
  else if( fTarget.compareTo("Optics3") == 0 )
    fTargetType = kOptics3;
  else if( fTarget.compareTo("VWires") == 0 )
    fTargetType = kVWires;
  else if( fTarget.compareTo("HWires") == 0 )
    fTargetType = kHWires;
  else {
    cout << "Unknown target type " << fTarget << " setting to W production" << endl;
    fTargetType = kProdW;
  }
  
  switch(fTargetType){
  case kProdW:
    ntargs       = 10;
    targwidth    = 2.5*mm;
    targthick    = 0.01*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_W");
    targzpos[0]  = -247.4 *mm;
    targzpos[1]  = -192.4 *mm;
    targzpos[2]  = -137.4 *mm;
    targzpos[3]  = -82.4 *mm;
    targzpos[4]  = -27.4 *mm;
    targzpos[5]  = 27.6 *mm;
    targzpos[6]  = 82.6 *mm;
    targzpos[7]  = 137.6 *mm;
    targzpos[8]  = 192.6 *mm;
    targzpos[9]  = 247.6 *mm;
    break;
  case kProdC:
    ntargs       = 10;
    targwidth    = 2.5*mm;
    targthick    = 0.125*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_C");
    targzpos[0]  = -247.4 *mm;
    targzpos[1]  = -192.4 *mm;
    targzpos[2]  = -137.4 *mm;
    targzpos[3]  = -82.4 *mm;
    targzpos[4]  = -27.4 *mm;
    targzpos[5]  = 27.6 *mm;
    targzpos[6]  = 82.6 *mm;
    targzpos[7]  = 137.6 *mm;
    targzpos[8]  = 192.6 *mm;
    targzpos[9]  = 247.6 *mm;
    break;
  case kOptics1:
    ntargs       = 4;
    targwidth    = 5.0*mm;
    targthick    = 0.2*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_C");
    targzpos[0]  = -300.0 *mm;
    targzpos[1]  = -150.0 *mm;
    targzpos[2]  = 75.0 *mm;
    targzpos[3]  = 219.0 *mm;
    break;
  case kOptics2:
    ntargs       = 4;
    targwidth    = 5.0*mm;
    targthick    = 0.2*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_C");
    targzpos[0]  = -300.0 *mm;
    targzpos[1]  = -219.0 *mm;
    targzpos[2]  = -150.0 *mm;
    targzpos[3]  = -75.0 *mm;
    targzpos[4]  = 75.0 *mm;
    targzpos[5]  = 150.0 *mm;
    targzpos[6]  = 219.0 *mm;
    targzpos[7]  = 300.0 *mm;
    break;
  case kOptics3:
    ntargs       = 4;
    targwidth    = 5.0*mm;
    targthick    = 0.2*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_C");
    targzpos[0]  = -219.0 *mm;
    targzpos[1]  = -75.0 *mm;
    targzpos[2]  = 150 *mm;
    targzpos[3]  = 300 *mm;
    break;
  case kVWires:
    ntargs       = 3;
    targwidth    = 30.0*mm;
    targthick    = 0.1*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_W");
    targzpos[0]  = -200.0 *mm;
    targzpos[1]  = 0.0 *mm;
    targzpos[2]  = 200. *mm;
    targxpos[0]  = -2.5 *mm;
    targxpos[1]  = 0.0 *mm;
    targxpos[2]  = 2.5 *mm;
    break;
  case kHWires:
    ntargs       = 4;
    targwidth    = 30.0*mm;
    targthick    = 0.1*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_W");
    targzpos[0]  = -250.0 *mm;
    targzpos[1]  = -100.0 *mm;
    targzpos[2]  = 100. *mm;
    targzpos[3]  = 250. *mm;
    targypos[0]  = 0.0 *mm;
    targypos[1]  = 5.0 *mm;
    targypos[2]  = 10. *mm;
    targypos[3]  = 15. *mm;
    break;
  default: // kProdW
    ntargs       = 10;
    targwidth    = 2.5*mm;
    targthick    = 0.01*mm;
    targMaterial = fNistManager->FindOrBuildMaterial("G4_W");
    targzpos[0]  = -247.4 *mm;
    targzpos[1]  = -192.4 *mm;
    targzpos[2]  = -137.4 *mm;
    targzpos[3]  = -82.4 *mm;
    targzpos[4]  = -27.4 *mm;
    targzpos[5]  = 27.6 *mm;
    targzpos[6]  = 82.6 *mm;
    targzpos[7]  = 137.6 *mm;
    targzpos[8]  = 192.6 *mm;
    targzpos[9]  = 247.6 *mm;
    break;
  }

  // Target mother volume
  G4Box *TargMother_box = new G4Box("solidTarg", targwidth, targwidth, 350. *mm );
  G4LogicalVolume *TargMother_log = new G4LogicalVolume( TargMother_box, Beamline, "TargMother_log" );
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);  
  new G4PVPlacement( rot_temp, G4ThreeVector( 0., 0., 0.), TargMother_log, "TargMother", logicScatChamber, false, ++nonSDcounter );

  // Target solids
  char solidname[100];
  G4VSolid *Target_solid;
  for( G4int itarg=0; itarg<ntargs; itarg++ ) {

    sprintf(solidname, "Target_solid%d", itarg );
    if( fTargetType == kVWires || fTargetType == kHWires )
      Target_solid = new G4Tubs(G4String(solidname), 0., targthick/2., targwidth/2., 0.0, 360.0 *deg );
    else
      Target_solid = new G4Box(G4String(solidname), targwidth/2., targwidth/2., targthick/2. );
    
    G4String logname = solidname;
    logname += "_log";
    G4LogicalVolume *Target_log = new G4LogicalVolume( Target_solid, targMaterial, logname );

    G4String physname = solidname;
    physname += "_phys";
    if( fTargetType == kVWires ) {
      rot_temp = new G4RotationMatrix();
      rot_temp->rotateX(90.0*deg);  
      new G4PVPlacement( rot_temp, G4ThreeVector( targxpos[itarg], targypos[itarg], targzpos[itarg]), Target_log, physname, TargMother_log, false, ++nonSDcounter );
    }
    else if( fTargetType == kHWires ) {
      rot_temp = new G4RotationMatrix();
      rot_temp->rotateY(90.0*deg);  
      new G4PVPlacement( rot_temp, G4ThreeVector( targxpos[itarg], targypos[itarg], targzpos[itarg]), Target_log, physname, TargMother_log, false, ++nonSDcounter );
    }
    else
      new G4PVPlacement( 0, G4ThreeVector( targxpos[itarg], targypos[itarg], targzpos[itarg]), Target_log, physname, TargMother_log, false, ++nonSDcounter );
  }

  //--------------------------------------------------------------------------- 
  // Create Septum (from original APEX G4)
  //--------------------------------------------------------------------------- 

  double y_en           = 2.44*inch;
  double y_ex           = 4.7 *inch;
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

  double ymin_sep_ex1   = y_en + (z_sept_ex_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_ex1   = y_en + (z_sept_ex_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymin_sep_en1   = y_en + (z_sept_en_min1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_en1   = y_en + (z_sept_en_max1-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);

  double ymin_sep_ex2   = y_en + (z_sept_ex_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_ex2   = y_en + (z_sept_ex_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymin_sep_en2   = y_en + (z_sept_en_min2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
  double ymax_sep_en2   = y_en + (z_sept_en_max2-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);

  G4double  pDz_1       = 0.31*inch/2.;
  G4double  pDz_2       = 0.25*inch/2.;
  
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
  G4LogicalVolume* trrap_en1_min = new G4LogicalVolume(testTrap_en1_min, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en1_min",0,0,LarmStepLimits);
  trrap_en1_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en1_min=new G4RotationMatrix();
  pRotY90deg_en1_min->rotateY(90.*deg-ang_en_min_1);
  new G4PVPlacement(pRotY90deg_en1_min,G4ThreeVector(pDz_1+xmin_sep_en1-0.02*mm,0,z_sept_en_min1), 
		    trrap_en1_min,"trrap_en1_min",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* trrap_en1_max = new G4LogicalVolume(testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en1_max",0,0,LarmStepLimits);
  trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en1_max=new G4RotationMatrix();
  pRotY90deg_en1_max->rotateY(90.*deg-ang_en_min_1);
  new G4PVPlacement(pRotY90deg_en1_max,G4ThreeVector(pDz_1+xmax_sep_en1+0.02*mm,0,z_sept_en_max1), 
		    trrap_en1_max,"trrap_en1_max",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* trrap_en2_max = new G4LogicalVolume(testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en2_max",0,0,LarmStepLimits);
  trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en2_max=new G4RotationMatrix();
  pRotY90deg_en2_max->rotateY(90.*deg-ang_en_max_2);
  new G4PVPlacement(pRotY90deg_en2_max,G4ThreeVector(pDz_1+xmax_sep_en2+0.02*mm,0,z_sept_en_max2), 
		    trrap_en2_max,"trrap_en2_max",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* trrap_en2_min = new G4LogicalVolume(testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en2_min",0,0,LarmStepLimits);
  trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotY90deg_en2_min=new G4RotationMatrix();
  pRotY90deg_en2_min->rotateY(90.*deg-ang_en_min_2);
  new G4PVPlacement(pRotY90deg_en2_min,G4ThreeVector(-1.*pDz_1+xmin_sep_en2-0.02*mm,0,z_sept_en_min2), 
		    trrap_en2_min,"trrap_en2_min",expHall_log,0,++nonSDcounter);



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
  G4LogicalVolume* trrap_cov_1_up = new G4LogicalVolume(testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_1_up",0,0,LarmStepLimits);
  trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_1_up=new G4RotationMatrix();
  pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(pRotX90deg_cov_1_up,G4ThreeVector(0,0.02*mm+ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.), 
		    trrap_cov_1_up,"cov1_up",expHall_log,0,++nonSDcounter);



  //LHS Back Bottom Trap
  G4VSolid* testTrap_cov_1_bot = testTrap_cov_1_up;
  G4LogicalVolume* trrap_cov_1_bot = new G4LogicalVolume(testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_1_bot",0,0,LarmStepLimits);
  trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_1_bot=new G4RotationMatrix();
  pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(pRotX90deg_cov_1_bot,G4ThreeVector(0,-0.02*mm-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.), 
		    trrap_cov_1_bot,"cov1_bot",expHall_log,0,++nonSDcounter);
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
  G4LogicalVolume* trrap_cov_2_up = new G4LogicalVolume(testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_2_up",0,0,LarmStepLimits);
  trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_2_up=new G4RotationMatrix();
  pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(pRotX90deg_cov_2_up,G4ThreeVector(0, 0.02*mm+ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.), 
		    trrap_cov_2_up,"cov2_up",expHall_log,0,++nonSDcounter);

  
  //RHS Front Bottom Trap
  G4VSolid* testTrap_cov_2_bot = testTrap_cov_2_up;
  G4LogicalVolume* trrap_cov_2_bot = new G4LogicalVolume(testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_2_bot",0,0,LarmStepLimits);
  trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_2_bot=new G4RotationMatrix();
  pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(pRotX90deg_cov_2_bot,G4ThreeVector(0,-0.02*mm-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.), 
		    trrap_cov_2_bot,"cov2_bot",expHall_log,0,++nonSDcounter);

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
  G4LogicalVolume* r_trrap = new G4LogicalVolume(r_testTrap,fNistManager->FindOrBuildMaterial("G4_W"),"vac Box 1",0,0,LarmStepLimits);
  r_trrap->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en1_min=new G4RotationMatrix();
  r_pRotY90deg_en1_min->rotateY(90*deg+ang_en_max_1);
  new G4PVPlacement(r_pRotY90deg_en1_min,G4ThreeVector(1.*pDz_1-xmin_sep_en1,0,z_sept_en_min1), 
		    r_trrap,"vac Box en1 min",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* r_trrap_en1_max = new G4LogicalVolume(r_testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en1_max",0,0,LarmStepLimits);
  r_trrap_en1_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en1_max=new G4RotationMatrix();
  r_pRotY90deg_en1_max->rotateY(90.*deg+ang_en_min_1);
  new G4PVPlacement(r_pRotY90deg_en1_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en1),0,z_sept_en_max1),
		    r_trrap_en1_max,"trrap_en1_max",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* r_trrap_en2_max = new G4LogicalVolume(r_testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en2_max",0,0,LarmStepLimits);
  r_trrap_en2_max->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en2_max=new G4RotationMatrix();
  r_pRotY90deg_en2_max->rotateY(90.*deg+ang_en_max_2);
  new G4PVPlacement(r_pRotY90deg_en2_max,G4ThreeVector(-1.*(pDz_1+xmax_sep_en2),0,z_sept_en_max2),
		    r_trrap_en2_max,"trrap_en2_max",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* r_trrap_en2_min = new G4LogicalVolume(r_testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_en2_min",0,0,LarmStepLimits);
  r_trrap_en2_min->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotY90deg_en2_min=new G4RotationMatrix();
  r_pRotY90deg_en2_min->rotateY(90.*deg+ang_en_min_2);
  new G4PVPlacement(r_pRotY90deg_en2_min,G4ThreeVector(pDz_1-xmin_sep_en2,0,z_sept_en_min2),
		    r_trrap_en2_min,"trrap_en2_min",expHall_log,0,++nonSDcounter);


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
  G4LogicalVolume* r_trrap_cov_1_up = new G4LogicalVolume(r_testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_1_up",0,0,LarmStepLimits);
  r_trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_1_up=new G4RotationMatrix();
  r_pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(r_pRotX90deg_cov_1_up,G4ThreeVector(0,ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.),
		    r_trrap_cov_1_up,"cov1_up",expHall_log,0,++nonSDcounter);

  
  //RHS Back Bottom Trap
  G4VSolid* r_testTrap_cov_1_bot = r_testTrap_cov_1_up;
  G4LogicalVolume* r_trrap_cov_1_bot = new G4LogicalVolume(r_testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_1_bot",0,0,LarmStepLimits);
  r_trrap_cov_1_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_1_bot=new G4RotationMatrix();
  r_pRotX90deg_cov_1_bot->rotateX(-90.*deg-atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(r_pRotX90deg_cov_1_bot,G4ThreeVector(0,-1.*ymax_sep_en1-pDz_2-fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.),
		    r_trrap_cov_1_bot,"cov1_bot",expHall_log,0,++nonSDcounter);
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
  G4LogicalVolume* r_trrap_cov_2_up = new G4LogicalVolume(r_testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_2_up",0,0,LarmStepLimits);
  r_trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_up=new G4RotationMatrix();
  r_pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_up,G4ThreeVector(0,     ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_up,"cov2_up",expHall_log,0,++nonSDcounter);


  //RHS Front Bottom Trap
  G4VSolid* r_testTrap_cov_2_bot = r_testTrap_cov_2_up;
  G4LogicalVolume* r_trrap_cov_2_bot = new G4LogicalVolume(r_testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_W"),"trrap_cov_2_bot",0,0,LarmStepLimits);
  r_trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_bot=new G4RotationMatrix();
  r_pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_bot,G4ThreeVector(0,-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_bot,"cov2_bot",expHall_log,0,++nonSDcounter);
  
  //--------------------------------------------------------------------------- 
  // Create Sieve (from original APEX G4)
  //---------------------------------------------------------------------------
 
  //sieve geometry is 8.38 x 8.38 x 1 inches, gap distance is 0.492 inch (horizontal) and 0.984 inch (vertical), gap diameter is 0.157e-3 inch, two gaps with 0.236e-3 inch

  double pSieveSlitX     = 3.25*inch; 
  double pSieveSlitY     = 4.25*inch; 
  double sieve_thickness = 0.5;
  double pSieveSlitZ     = sieve_thickness*inch;

  double sieve_distance  = 31.23*inch;

  double sieve_pos_z     = -fDistTarPivot*cm+(sieve_distance + pSieveSlitZ/2.)*cos(fSieveAngle*deg); // DJH: was 105.30001365 cm and 5.81 degrees???
  double sieve_pos_x     = (sieve_distance + pSieveSlitZ/2.)*sin(fSieveAngle*deg);                   // DJH: was 5.81 degrees???

  double startphi        = 0.0*deg;
  double deltaphi        = 360.0*deg;

  double pSieveSlitHoleR       = 0.5*0.055*inch;      //radius of small hole 0.055/2 inch
  double pSieveSlitLargeHoleR  = 0.5*0.106*inch;      //radius of large hole 0.106/2 inch
  double pSieveSlitMediumHoleR = 1.*0.5*0.075*inch;   //radius of small hole 0.055/2 inch
  double pSieveSlitHoldL       = pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

  //the hole positions relative to the slit center 
  double pSieveSlitHolePosH[15]={0.295, 0.485, 0.675, 0.865, 1.055, 1.245, 1.435, 1.625, 1.815, 2.005, 2.195, 2.385, 2.575, 2.765, 2.955};
  double pSieveSlitHolePosV[9] = {0.285, 0.745, 1.205, 1.665, 2.125, 2.585, 3.045, 3.505, 3.965};
  for(int ii=0;ii<15;ii++)
    {
      pSieveSlitHolePosH[ii] = (pSieveSlitHolePosH[ii]-1.625)*inch;
    }
  for(int ii=0;ii<9;ii++)
    {
      pSieveSlitHolePosV[ii] = (pSieveSlitHolePosV[ii]-2.125)*inch;
    }

  //now start to build box then subtract 63 holes 
  G4VSolid* sieveSlitWholeSolid      = new G4Box("sieveSlitWholeBox",pSieveSlitX/2.0,
						 pSieveSlitY/2.0,pSieveSlitZ/2.0);

  G4VSolid* sieveSlitHoleSolid       = new G4Tubs("sieveSlitHoleTubs",0,pSieveSlitHoleR,
						  pSieveSlitHoldL/2.0,startphi,deltaphi); 
  
  G4VSolid* sieveSlitLargeHoleSolid  = new G4Tubs("sieveSlitLargeHoleTubs",0,
						  pSieveSlitLargeHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 
  
  G4VSolid* sieveSlitMediumHoleSolid = new G4Tubs("sieveSlitMediumHoleTubs",0,
						  pSieveSlitMediumHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 

  // DJH: added else if because medium holes were overlapping small holes in original
  G4SubtractionSolid* sieveSlitSolid = (G4SubtractionSolid*)sieveSlitWholeSolid;
  char strName[100];
  for(int ih=0;ih<15;ih++) {
      for(int iv=0;iv<9;iv++) {
	sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
	if((ih==7 && iv==4) || (ih==3 && iv==2)) {
	  //now dig large holes in the block
	  sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
						sieveSlitLargeHoleSolid,0,
						G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
	}
	else if ((iv>=7) || (iv<2)) {
	  sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
						sieveSlitMediumHoleSolid,0,
						G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv], 0.3*inch));
	}
	else {
	  //now dig small holes in the block
	  sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
						sieveSlitHoleSolid,0,
						G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv], 0));
	}
      }
  }

  // DJH: commented out for now -- not sure what this next block is doing ?????
  // for(int ih=0;ih<14-2;ih++)
  //   {
  //     for(int iv=0;iv<8;iv++)
  // 	{
  // 	  sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
  // 	  if ( !( ( ih == 6 ) && (iv>0) && (iv<7) ) )
  // 	    {
  // 	      //now dig small holes in the block
  // 	      sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
  // 						    sieveSlitHoleSolid,0,
  // 						    G4ThreeVector(0.5*(pSieveSlitHolePosH[ih]+pSieveSlitHolePosH[ih+1]), 0.5*(pSieveSlitHolePosV[iv]+pSieveSlitHolePosV[iv+1]), 0));
  // 	    }
  // 	  if ((iv==7) || (iv==0))
  // 	    {
  // 	      sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
  // 						    sieveSlitMediumHoleSolid,0,
  // 						    G4ThreeVector(0.5*(pSieveSlitHolePosH[ih]+pSieveSlitHolePosH[ih+1]), 0.5*(pSieveSlitHolePosV[iv]+pSieveSlitHolePosV[iv+1]), 0.3*inch));
  // 	    }
  // 	}
  // }

  sieveSlitSolid->SetName("sieveSlitSolid");
  
  G4LogicalVolume* sieveSlitLogical = new G4LogicalVolume(sieveSlitSolid,
							  fNistManager->FindOrBuildMaterial("G4_AIR"),
							  "sieveSlitLogical",0,0,LarmStepLimits);

  G4RotationMatrix *R_RotY90deg_sieve_slit=new G4RotationMatrix();
  R_RotY90deg_sieve_slit->rotateY(-fSieveAngle*deg);
  G4RotationMatrix *L_RotY90deg_sieve_slit=new G4RotationMatrix();
  L_RotY90deg_sieve_slit->rotateY(fSieveAngle*deg);

  if( fSieveOn ) {
    new G4PVPlacement(R_RotY90deg_sieve_slit,G4ThreeVector(sieve_pos_x, 0, sieve_pos_z),
		      sieveSlitLogical,"RSievePhys", expHall_log, 0, 0, 0);
    
    
    new G4PVPlacement(L_RotY90deg_sieve_slit,G4ThreeVector(-sieve_pos_x, 0, sieve_pos_z),
		      sieveSlitLogical,"LSievePhys", expHall_log, 0, 0, 0);
  }
  
  //--------------------------------------------------------------------------- 
  // Create Q1 SOS (from original APEX G4)
  //--------------------------------------------------------------------------- 

  double pQ1Rin          = 12.827 *cm;
  double pQ1Length       = 1.0 *cm;
  
  G4VSolid* Q1_tubs      = new G4Tubs("Q1Front",
				      0, 2*(pQ1Rin+1.*cm), pQ1Length/2., // make virtual detector radius 1cm larger
				      0.0, 360.0*deg);
  
  G4LogicalVolume* Q1_log = new G4LogicalVolume(Q1_tubs,
						fNistManager->FindOrBuildMaterial("G4_AIR"), 
						"Q1_log", 0, 0, 0); 

  //--------------------------------------------------------------------------- 
    
  G4double LQ1_th           = fHRSAngle *deg; 
  G4double LQ1_d            = fDistPivotQ1 *cm - pQ1Length; 
  G4double LQ1_xprime       = -LQ1_d * std::sin(LQ1_th); 
  G4double LQ1_zprime       = LQ1_d * std::cos(LQ1_th); 
  G4RotationMatrix* LQ1_rm  = new G4RotationMatrix(); 
  LQ1_rm->rotateY(LQ1_th); 

  fDetVol[0]                    = new G4PVPlacement(LQ1_rm, G4ThreeVector(LQ1_xprime,0.,LQ1_zprime), 
						    Q1_log, "LQ1", expHall_log, false, SDcounter);

  //--------------------------------------------------------------------------- 
  
  G4double RQ1_th           = -fHRSAngle *deg; 
  G4double RQ1_d            = fDistPivotQ1 *cm - pQ1Length; 
  G4double RQ1_xprime       = -RQ1_d * std::sin(RQ1_th); 
  G4double RQ1_zprime       = RQ1_d * std::cos(RQ1_th); 
  G4RotationMatrix* RQ1_rm  = new G4RotationMatrix(); 
  RQ1_rm->rotateY(RQ1_th); 
     
  fDetVol[1]                    = new G4PVPlacement(RQ1_rm, G4ThreeVector(RQ1_xprime,0.,RQ1_zprime), 
						    Q1_log, "RQ1", expHall_log, false, ++SDcounter); 

  //---------------------------------------------------------------------------
  // Set Logical Attributes
  //---------------------------------------------------------------------------

  // Senstive detector
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  fFluxSD = new FluxSD("FluxSD", fNSD);
  SDman->AddNewDetector( fFluxSD );
  Q1_log->SetSensitiveDetector( fFluxSD );

  // Magnetic field
  expHall_log->SetFieldManager(fm, true);

  // Visualisation
  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0) );
  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  TargMother_log->SetVisAttributes(G4VisAttributes::Invisible);
  Q1_log->SetVisAttributes(blue);

  //---------------------------------------------------------------------------

  return fExpHall;
}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------
