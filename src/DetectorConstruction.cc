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
#include "G4UnionSolid.hh" 
#include "G4IntersectionSolid.hh" 

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
  fDistTarPivot = 105.29;  // from original APEX G4
  fDistPivotQ1  = 171.095; // Q1 SOS (from original APEX G4)

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
  // Create Upstream Beamline (copy from "standard" (GMn) config in g4sbs)
  //--------------------------------------------------------------------------- 

  G4bool ChkOverlaps = false;
  G4double inch = 2.54*cm;
  
  double sc_entbeampipeflange_dist = 25.375*2.54*cm; // entrance pipe flange distance from hall center 
  
  G4double ent_len = 2*m; 
  G4double ent_rin = 31.75*mm; 
  G4double ent_rou = ent_rin+0.120*mm; 
  
  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg ); 
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg ); 
  
  G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), "ent_log", 0, 0, 0); 
  G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, Beamline, "entvac_log", 0, 0, 0); 
   
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist-fDistTarPivot*cm), entLog, "ent_phys", expHall_log, false, ++nonSDcounter , ChkOverlaps); 
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist-fDistTarPivot*cm), entvacLog, "entvac_phys", expHall_log,false, ++nonSDcounter , ChkOverlaps); 

  //---------------------------------------------------------------------------
  // Create Scattering Chamber (copy from "standard" (GMn) config in g4sbs)
  //---------------------------------------------------------------------------

  G4LogicalVolume *logicScatChamberTank =0;
  G4LogicalVolume *logicScatChamberExitFlangePlate =0;
  G4LogicalVolume *logicScatChamberFrontClamshell =0;
  G4LogicalVolume *logicScatChamberBackClamshell =0;
  G4LogicalVolume *logicScatChamberLeftSnoutWindow =0;
  G4LogicalVolume *logicScatChamberLeftSnoutWindowFrame =0;
  G4LogicalVolume *logicScatChamberRightSnoutWindow =0;
  G4LogicalVolume *logicScatChamberRightSnoutWindowFrame =0;
  G4LogicalVolume *logicScatChamber =0;

  // Scattering chamber tank:
  // basic volume:
  G4double SCHeight = 44.75*inch;
  G4double SCRadius = 20.0*inch;
  G4double SCTankThickness = 2.5*inch;
  G4double SCTankRadius = SCRadius+SCTankThickness;
  G4double SCTankHeight = SCHeight;
  G4double SCOffset = 3.75*inch;
  
  G4Tubs* solidSCTank_0 = 
    new G4Tubs("SCTank_0", SCRadius, SCTankRadius, 0.5*SCTankHeight, 0.0*deg, 360.0*deg);
  
  // exit flange:
  G4double SCExitFlangePlateHLength = 22.5*sin(25.5*atan(1)/45.0)*inch;
  G4double SCExitFlangePlateHeight = 11.0*inch;
  G4double SCExitFlangePlateThick = 1.25*inch;
  G4double SCExitFlangeHAngleApert = atan(SCExitFlangePlateHLength/(SCTankRadius+SCExitFlangePlateThick));
  G4double SCExitFlangeMaxRad = SCExitFlangePlateHLength/sin(SCExitFlangeHAngleApert);
  
  G4Tubs* solidSCExitFlangetubs = 
    new G4Tubs("SCExFlange_tubs", SCRadius, SCExitFlangeMaxRad, 
  	       0.5*SCExitFlangePlateHeight, 0.0, 2.0*SCExitFlangeHAngleApert);
  
  G4RotationMatrix* rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHAngleApert);
  
  G4UnionSolid* solidSCTank_0_exft = 
    new G4UnionSolid("solidSCTank_0_exft", solidSCTank_0, solidSCExitFlangetubs,
  		     rot_temp, G4ThreeVector(0,0,SCOffset));
  
  G4Box* ExitFlangeHeadCut = new G4Box("ExitFlangeHeadCut", 0.5*m, 0.5*m, 0.5*m); 
  
  
  G4SubtractionSolid* solidSCTank_0_exf = 
    new G4SubtractionSolid("solidSCTank_0_exf", solidSCTank_0_exft, ExitFlangeHeadCut,
  			   0, G4ThreeVector(-SCTankRadius-0.5*m,0,0));
  
  // exit flange hole:
  G4double SCExitFlangeHoleHeight = 7.85*inch;
  G4double SCExitFlangeHoleAngleApert = 38.25*deg;
  
  G4Tubs* solidSCExFH = 
    new G4Tubs("SCExFH", SCRadius-1.0*cm, SCTankRadius+1.5*inch,
  	       0.5*SCExitFlangeHoleHeight, 0.0, SCExitFlangeHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCTank_0_exfh = 
    new G4SubtractionSolid("solidSCTank_0_exfh", solidSCTank_0_exf, solidSCExFH,
  			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // windows holes: 
  G4double SCWindowHeight = 18.0*inch;
  G4double SCWindowAngleApert = 149.0*deg;
  G4double SCWindowAngleOffset = 11.0*deg;
  
  G4Tubs* solidSCWindow = 
    new G4Tubs("SCWindow", SCRadius-1.0*cm, SCTankRadius+1.0*cm, 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wf = 
    new G4SubtractionSolid("solidSCTank_0_wf", solidSCTank_0_exfh, solidSCWindow,
  			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wb = 
    new G4SubtractionSolid("solidSCTank_0_wb", solidSCTank_0_wf, solidSCWindow,
  			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // Solid Scattering chamber tank
  G4VSolid* solidScatChamberTank = solidSCTank_0_wb;
  
  // Logic scat chamber tank
  logicScatChamberTank = 
    new G4LogicalVolume(solidScatChamberTank, fNistManager->FindOrBuildMaterial("G4_Al"), "ScatChamberTank_log");
  
  // Scattering chamber tank placement:
  G4RotationMatrix* rotSC = new G4RotationMatrix();
  rotSC->rotateX(-90.0*deg);
  
  G4ThreeVector* SCPlacement = new G4ThreeVector(0,-fDistTarPivot*cm,-SCOffset);
  SCPlacement->rotateX(90*deg);
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamberTank, "ScatChamberTankPhys", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  //Exit Flange Plate
  G4Box* solidSCExitFlangePlate = 
    new G4Box("SCExitFlangePlate_sol", SCExitFlangePlateThick/2.0, 
  	      SCExitFlangePlateHeight/2.0, SCExitFlangePlateHLength); 
  
  logicScatChamberExitFlangePlate = 
    new G4LogicalVolume(solidSCExitFlangePlate, fNistManager->FindOrBuildMaterial("G4_Al"), "SCExitFlangePlate_log");

  new G4PVPlacement(0, G4ThreeVector(-SCTankRadius-SCExitFlangePlateThick/2.0,0,-fDistTarPivot*cm), 
  		    logicScatChamberExitFlangePlate, "SCExitFlangePlate", expHall_log, false, ++nonSDcounter, ChkOverlaps); 
  
  // Front and back Clamshells...
  // Basic solid: 
  G4double SCClamHeight = 20.0*inch;
  G4double SCClamThick = 1.25*inch;
  G4double SCClamAngleApert = 151.0*deg;
  
  G4Tubs* solidSCClamshell_0 = 
    new G4Tubs("solidSCClamshell_0", SCTankRadius, SCTankRadius+SCClamThick, 
  	       0.5*SCClamHeight, 0.0, SCClamAngleApert);
  
  // Front Clamshell:
  G4double SCFrontClamOuterRadius = SCTankRadius+SCClamThick;
  G4double SCBeamExitAngleOffset = 64.5*deg;
  G4double SCLeftSnoutAngle = -24.2*deg;
  G4double SCLeftSnoutAngleOffset = SCBeamExitAngleOffset+SCLeftSnoutAngle;
  //G4double SCRightSnoutAngle = 50.1*deg;
  G4double SCRightSnoutAngle = 47.5*deg;
  G4double SCRightSnoutAngleOffset = SCBeamExitAngleOffset+SCRightSnoutAngle;
  
  // Snouts: NB: similarly to the GEp scattering chamber, 
  // the "left" and "right" is defined as viewed from downstream.
  // In other words, the "left" snout is actually on the right side of the beam 
  // (looking on the right direction) and vice-versa.
  
  // Right snout opening:
  G4double SCRightSnoutDepth = 15.0*inch;// x
  G4double SCRightSnoutWidth = 26.0*inch;// y
  G4double SCRightSnoutHeight = 18.0*inch; //z
  G4double SCRightSnoutHoleHeight = 14.0*inch;
  G4double SCRightSnoutHoleAngleApert = 55.0*deg;
  G4double SCRightSnoutBoxAngleOffset = 9.40*deg;
  
  // Basic "box"
  G4Box* RightSnoutBox = 
    new G4Box("RightSnoutBox", SCRightSnoutDepth*0.5, SCRightSnoutWidth*0.5, SCRightSnoutHeight*0.5);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCRightSnoutAngleOffset-SCRightSnoutBoxAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_rsb = 
    new G4UnionSolid("solidSCClamshell_0_rsb", solidSCClamshell_0, 
  		     RightSnoutBox, rot_temp, 
  		     G4ThreeVector(SCFrontClamOuterRadius*sin(90.0*deg-SCRightSnoutAngleOffset),
  		      		   SCFrontClamOuterRadius*cos(90.0*deg-SCRightSnoutAngleOffset), 
  		     		   0)
  		     );
  
  // Basic box cut: remove the surplus outside of the clamshell
  // NB: the surplus inside is removed along with the one of the left snout
  G4Box* RightSnoutBox_cut = 
    new G4Box("RightSnoutBox_cut", 
  	      SCRightSnoutDepth+0.01*inch, SCRightSnoutWidth, SCRightSnoutHeight*0.5+0.01*inch);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCRightSnoutAngleOffset);
  
  G4SubtractionSolid* solidSCFrontClam_0_rsbc = 
    new G4SubtractionSolid("solidSCFrontClam_0_rsbc", solidSCFrontClam_0_rsb,
  			   RightSnoutBox_cut, rot_temp, 
  			   G4ThreeVector((SCFrontClamOuterRadius+SCRightSnoutDepth)
  					 *sin(90.0*deg-SCRightSnoutAngleOffset), 
  					 (SCFrontClamOuterRadius+SCRightSnoutDepth)
  					 *cos(90.0*deg-SCRightSnoutAngleOffset),
  					 0)
  			   );
  
  // Cutting the hole...
  G4Tubs* RightSnoutApertCut = 
    new G4Tubs("RightSnoutApertCut", 0.0, 30.0*inch, SCRightSnoutHoleHeight/2.0, 
  	       0.0, SCRightSnoutHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-SCRightSnoutAngleOffset+SCRightSnoutHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCFrontClam_0_rs =
    new G4SubtractionSolid("solidSCClamshell_0_rs", solidSCFrontClam_0_rsbc, 
  			   RightSnoutApertCut, rot_temp, G4ThreeVector(0, 0, 0));
  
  
  // Right snout window+frame:
  G4double SCRightSnoutWindowWidth = 26.351*inch;
  G4double SCSnoutWindowThick = 0.02*inch;
  G4double SCSnoutWindowFrameThick = 0.75*inch;
  
  G4double SCRightSnoutWindowDist = 23.74*inch+SCSnoutWindowThick*0.5;
  //G4double SCRightSnoutHoleWidth = 21.855*inch;
  G4double SCRightSnoutHoleWidth = 24.47*inch;
  G4double SCRightSnoutHoleCurvRad = 2.1*inch;
  G4double SCRightSnoutWindowFrameDist = SCRightSnoutWindowDist+SCSnoutWindowThick*0.5+SCSnoutWindowFrameThick*0.5;
  
  // window
  G4Box* solidRightSnoutWindow = 
    new G4Box("RightSnoutWindow_sol", SCRightSnoutWindowWidth*0.5, 
  	      SCRightSnoutHeight*0.5, SCSnoutWindowThick*0.5);

  logicScatChamberRightSnoutWindow = 
    new G4LogicalVolume(solidRightSnoutWindow, fNistManager->FindOrBuildMaterial("G4_Al"), "SCRightSnoutWindow_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCRightSnoutAngle);
  
  new G4PVPlacement(rot_temp, 
  		    G4ThreeVector(SCRightSnoutWindowDist*sin(SCRightSnoutAngle),
  				  0,
  				  SCRightSnoutWindowDist*cos(SCRightSnoutAngle)-fDistTarPivot*cm), 
  		    logicScatChamberRightSnoutWindow, "SCRightSnoutWindow", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  // basic window frame
  G4Box* solidRightSnoutWindowFrame_0 = 
    new G4Box("RightSnoutWindowFrame_0", SCRightSnoutWidth*0.5, 
  	      SCRightSnoutHeight*0.5, SCSnoutWindowFrameThick*0.5);
  
  // + lots of cut out solids...
  G4Tubs* solidRightSnoutWindowFrame_cut_0 = 
    new G4Tubs("RightSnoutWindowFrame_cut_0", 0.0, SCRightSnoutHoleCurvRad, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_htr = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_htr", solidRightSnoutWindowFrame_0,
  			   solidRightSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
  					 SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad,
  					 0)
  			   );

  G4SubtractionSolid* solidRightSnoutWindowFrame_0_hbr = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_hbr", solidRightSnoutWindowFrame_0_htr,
  			   solidRightSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
  					 -SCRightSnoutHoleHeight*0.5+SCRightSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_hbl = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_hbl", solidRightSnoutWindowFrame_0_hbr,
  			   solidRightSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(-SCRightSnoutHoleWidth*0.5+SCRightSnoutHoleCurvRad, 
  					 -SCRightSnoutHoleHeight*0.5+SCRightSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_htl = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_htl", solidRightSnoutWindowFrame_0_hbl,
  			   solidRightSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(-SCRightSnoutHoleWidth*0.5+SCRightSnoutHoleCurvRad, 
  					 SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4Box* solidRightSnoutWindowFrame_cut_1 = 
    new G4Box("RightSnoutWindowFrame_cut_1", SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
  	      SCRightSnoutHoleHeight*0.5, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_1 = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame", solidRightSnoutWindowFrame_0_htl,
  			   solidRightSnoutWindowFrame_cut_1, 0, G4ThreeVector());
  
  G4Box* solidRightSnoutWindowFrame_cut_2 = 
    new G4Box("RightSnoutWindowFrame_cut_2", SCRightSnoutHoleWidth*0.5,
  	      SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame", solidRightSnoutWindowFrame_1,
  			   solidRightSnoutWindowFrame_cut_2, 0, G4ThreeVector());

  logicScatChamberRightSnoutWindowFrame = 
    new G4LogicalVolume(solidRightSnoutWindowFrame, fNistManager->FindOrBuildMaterial("G4_Al"), "SCRightSnoutWindowFrame_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCRightSnoutAngle);
  
  new G4PVPlacement(rot_temp, 
  		    G4ThreeVector(SCRightSnoutWindowFrameDist*sin(SCRightSnoutAngle),
  				  0,
  				  SCRightSnoutWindowFrameDist*cos(SCRightSnoutAngle)-fDistTarPivot*cm), 
  		    logicScatChamberRightSnoutWindowFrame, "SCRightSnoutWindowFrame", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
		    
  // Left snout opening:
  G4double SCLeftSnoutDepth = 4.0*inch;// x
  G4double SCLeftSnoutWidth = 16.338*inch;// y
  G4double SCLeftSnoutHeight = 11.0*inch; //z
  G4double SCLeftSnoutHoleHeight = 7.0*inch;
  G4double SCLeftSnoutHoleAngleApert = 30.0*deg;
  G4double SCLeftSnoutYOffset = 0.50*inch;
  
  // Basic box
  G4Box* LeftSnoutBox = 
    new G4Box("LeftSnoutBox", SCLeftSnoutDepth*0.5, SCLeftSnoutWidth*0.5, SCLeftSnoutHeight*0.5);

  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCLeftSnoutAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_lsb = 
    new G4UnionSolid("solidSCClamshell_0_lsb", solidSCFrontClam_0_rs, 
  		     LeftSnoutBox, rot_temp, 
  		     G4ThreeVector((SCFrontClamOuterRadius-1.85*inch)*sin(90.0*deg-SCLeftSnoutAngleOffset),
  		     		   (SCFrontClamOuterRadius-1.85*inch)*cos(90.0*deg-SCLeftSnoutAngleOffset), 
  		     		   SCLeftSnoutYOffset)
  		     );
  
  // remove all surplus material inside the scat chamber
  G4Tubs* SnoutsInnerCut = 
    new G4Tubs("SnoutsInnerCut", 0.0, SCTankRadius, SCClamHeight/2.0, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_lsbc =
    new G4SubtractionSolid("solidSCClamshell_0_lsbc", solidSCFrontClam_0_lsb, 
  			   SnoutsInnerCut, 0, G4ThreeVector());
  
  // Cut the hole
  G4Tubs* LeftSnoutApertCut = 
    new G4Tubs("LeftSnoutApertCut", 0.0, 30.0*inch, SCLeftSnoutHoleHeight/2.0, 
  	       0.0, SCLeftSnoutHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-SCLeftSnoutAngleOffset+SCLeftSnoutHoleAngleApert*0.5);
 
  G4SubtractionSolid* solidSCFrontClam_0_ls =
    new G4SubtractionSolid("solidSCClamshell_0_ls", solidSCFrontClam_0_lsbc, 
  			   LeftSnoutApertCut, rot_temp, G4ThreeVector(0, 0, SCLeftSnoutYOffset));
  
  // Left snout window+frame:
  G4double SCLeftSnoutWindowDist = 23.9*inch+SCSnoutWindowThick*0.5;
  G4double SCLeftSnoutWindowFrameDist = SCLeftSnoutWindowDist+SCSnoutWindowThick*0.5+SCSnoutWindowFrameThick*0.5;
  G4double SCLeftSnoutHoleWidth = 12.673*inch;
  G4double SCLeftSnoutHoleCurvRad = 1.05*inch;
  
  // window
  G4Box* solidLeftSnoutWindow_0 = 
    new G4Box("LeftSnoutWindow_sol", SCLeftSnoutWidth*0.5, SCLeftSnoutHeight*0.5, SCSnoutWindowThick*0.5);
  
  G4Tubs* solidLeftSnoutWindow_cut_0 = 
    new G4Tubs("LeftSnoutWindow_cut_0", 0.0, 2.579*inch, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidLeftSnoutWindow = 
    new G4SubtractionSolid("solidLeftSnoutWindow", solidLeftSnoutWindow_0, solidLeftSnoutWindow_cut_0,
  			   0, G4ThreeVector(10.287*inch, -SCLeftSnoutYOffset, 0));
  
  logicScatChamberLeftSnoutWindow = 
    new G4LogicalVolume(solidLeftSnoutWindow, fNistManager->FindOrBuildMaterial("G4_Al"), "SCLeftSnoutWindow_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCLeftSnoutAngle);
  
  new G4PVPlacement(rot_temp, 
  		    G4ThreeVector(SCLeftSnoutWindowDist*sin(SCLeftSnoutAngle),
  				  SCLeftSnoutYOffset,
  				  SCLeftSnoutWindowDist*cos(SCLeftSnoutAngle)-fDistTarPivot*cm), 
  		    logicScatChamberLeftSnoutWindow, "SCLeftSnoutWindow", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  // window frame
  G4Box* solidLeftSnoutWindowFrame_0 = 
    new G4Box("LeftSnoutWindowFrame_0", SCLeftSnoutWidth*0.5, 
  	      SCLeftSnoutHeight*0.5, SCSnoutWindowFrameThick*0.5);
  
  // + lots of cut out solids...
  G4Tubs* solidLeftSnoutWindowFrame_cut_0 = 
    new G4Tubs("LeftSnoutWindowFrame_cut_0", 0.0, SCLeftSnoutHoleCurvRad, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_htr = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_htr", solidLeftSnoutWindowFrame_0,
  			   solidLeftSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
  					 SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad,
  					 0)
  			   );

  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_hbr = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_hbr", solidLeftSnoutWindowFrame_0_htr,
  			   solidLeftSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
  					 -SCLeftSnoutHoleHeight*0.5+SCLeftSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_hbl = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_hbl", solidLeftSnoutWindowFrame_0_hbr,
  			   solidLeftSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(-SCLeftSnoutHoleWidth*0.5+SCLeftSnoutHoleCurvRad, 
  					 -SCLeftSnoutHoleHeight*0.5+SCLeftSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_htl = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_htl", solidLeftSnoutWindowFrame_0_hbl,
  			   solidLeftSnoutWindowFrame_cut_0, 0, 
  			   G4ThreeVector(-SCLeftSnoutHoleWidth*0.5+SCLeftSnoutHoleCurvRad, 
  					 SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad,
  					 0)
  			   );
  
  G4Box* solidLeftSnoutWindowFrame_cut_1 = 
    new G4Box("LeftSnoutWindowFrame_cut_1", SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
  	      SCLeftSnoutHoleHeight*0.5, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_1 = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame", solidLeftSnoutWindowFrame_0_htl,
  			   solidLeftSnoutWindowFrame_cut_1, 0, G4ThreeVector());
  
  G4Box* solidLeftSnoutWindowFrame_cut_2 = 
    new G4Box("LeftSnoutWindowFrame_cut_2", SCLeftSnoutHoleWidth*0.5,
  	      SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_2 = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_2", solidLeftSnoutWindowFrame_1,
  			   solidLeftSnoutWindowFrame_cut_2, 0, G4ThreeVector());

  G4SubtractionSolid* solidLeftSnoutWindowFrame = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame", solidLeftSnoutWindowFrame_2, solidLeftSnoutWindow_cut_0,
  			   0, G4ThreeVector(10.287*inch, -SCLeftSnoutYOffset, 0));
  
  logicScatChamberLeftSnoutWindowFrame = 
    new G4LogicalVolume(solidLeftSnoutWindowFrame, fNistManager->FindOrBuildMaterial("G4_Al"), "SCLeftSnoutWindowFrame_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCLeftSnoutAngle);
  
  new G4PVPlacement(rot_temp, 
  		    G4ThreeVector(SCLeftSnoutWindowFrameDist*sin(SCLeftSnoutAngle),
  				  SCLeftSnoutYOffset,
  				  SCLeftSnoutWindowFrameDist*cos(SCLeftSnoutAngle)-fDistTarPivot*cm), 
  		    logicScatChamberLeftSnoutWindowFrame, "SCLeftSnoutWindowFrame", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  // Exit Beam Pipe:
  // Should come after left snout
  G4Tubs* solidExitBeamPipe = new G4Tubs("solidExitBeamPipe", 
  					 0.0, 55.0*mm, 2.150*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg+SCBeamExitAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCFrontClam_0_ls, 
  		     solidExitBeamPipe, rot_temp, 
  		     G4ThreeVector((SCFrontClamOuterRadius+2.003*inch)*sin(90.0*deg-SCBeamExitAngleOffset),
  				   (SCFrontClamOuterRadius+2.003*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
  				   0.0)
  		     );
  
  // remove extra material from the left snout around the chamber
  G4Tubs* solidExitBeamPipeSurroundCut = new G4Tubs("solidExitBeamPipeSurroundCut", 
  						    55.0*mm, 2.803*inch, 1.684*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_ebps =
    new G4SubtractionSolid("solidSCClamshell_0_ebps", solidSCFrontClam_0_ebp, 
  			   solidExitBeamPipeSurroundCut, rot_temp, 
  			   G4ThreeVector((SCFrontClamOuterRadius+1.684*inch)*sin(90.0*deg-SCBeamExitAngleOffset),
  					 (SCFrontClamOuterRadius+1.684*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
  					 0.0)
  			   );
   
  G4Tubs* solidExitBeamPipeFlange = new G4Tubs("solidExitBeamPipeFlange", 
  					       0.0, 2.985*inch, 0.3925*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCFrontClam_0_ebpf =
    new G4UnionSolid("solidSCClamshell_0_ebpf", solidSCFrontClam_0_ebps, 
  		     solidExitBeamPipeFlange, rot_temp, 
  		     G4ThreeVector((SCFrontClamOuterRadius+3.7605*inch)*sin(90.0*deg-SCBeamExitAngleOffset), 
  				   (SCFrontClamOuterRadius+3.7605*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
  				   0.0));
  
  G4Tubs* solidExitBeamPipeHole = new G4Tubs("solidBackViewPipeHole", 
  					     0.0, 50.0*mm, 7.903*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCFrontClam_0_ebpf, 
  			   solidExitBeamPipeHole, rot_temp, 
  			   G4ThreeVector(SCFrontClamOuterRadius*sin(90.0*deg-SCBeamExitAngleOffset),
  					 SCFrontClamOuterRadius*cos(90.0*deg-SCBeamExitAngleOffset), 
  					 0.0)
  			   );
  
  // Placing Front ClamShell (at last)
  G4VSolid* solidSCFrontClamShell = solidSCFrontClam_0_ebph;
  
  logicScatChamberFrontClamshell = 
    new G4LogicalVolume(solidSCFrontClamShell, fNistManager->FindOrBuildMaterial("G4_Al"), "SCFrontClamshell_log");
  //new G4LogicalVolume(solidSCFrontClam_0_ebps, fNistManager->FindOrBuildMaterial("G4_Al"), "SCFrontClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(90.0*deg+SCClamAngleApert*0.5-SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp, G4ThreeVector(0,0,-fDistTarPivot*cm), logicScatChamberFrontClamshell, 
  		    "SCFrontClamshell", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  // Back Clamshell:
  G4double SCBackClamThick = 0.80*inch;
  G4double SCBackClamHeight = 12.0*inch;
  G4double SCBackClamAngleApert = 145.7*deg;
  G4double SCBackClamAngleOffset = -2.5*deg;
  G4double SCBackClamOuterRadius = SCTankRadius+SCClamThick+SCBackClamThick;
  G4double SCBeamEntranceAngleOffset = SCBackClamAngleOffset-84*deg;
  G4double SCBackViewPipeAngleOffset = SCBackClamAngleOffset-39.9*deg;
  
  G4Tubs* solidSCAddBackClam = 
    new G4Tubs("solidSCAddBackClam", SCTankRadius+SCClamThick-0.5*cm, SCBackClamOuterRadius, 
  	       0.5*SCBackClamHeight, 0.0, SCBackClamAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(SCBackClamAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_abc = 
    new G4UnionSolid("solidSCClamshell_0_abc", solidSCClamshell_0, solidSCAddBackClam,
  		     rot_temp, G4ThreeVector(0,0,0));
  
  // Entrance beam pipe
  G4Tubs* solidEntranceBeamPipe = 
    new G4Tubs("solidEntranceBeamPipe", 0.0, 2.25*inch, 0.825*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBeamEntranceAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCBackClam_0_abc, solidEntranceBeamPipe, rot_temp, 
  		     G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBeamEntranceAngleOffset),
  				   SCBackClamOuterRadius*cos(90.0*deg+SCBeamEntranceAngleOffset),
  				   0.0)
  		     );
  
  G4Tubs* solidEntranceBeamPipeHole = 
    new G4Tubs("solidEntranceBeamPipeHole", 0.0, 1.0*inch, 5.375*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCBackClam_0_ebp, 
  			   solidEntranceBeamPipeHole, rot_temp, 
  			   G4ThreeVector(SCBackClamOuterRadius*
  					 sin(90.0*deg+SCBeamEntranceAngleOffset),
  					 SCBackClamOuterRadius*
  					 cos(90.0*deg+SCBeamEntranceAngleOffset),
  					 0.0)
  			   );
  
  // Back view pipe
  G4Tubs* solidBackViewPipe = new G4Tubs("solidBackViewPipe", 
  					 0.0, 2.12*inch, 1.555*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBackViewPipeAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_bvp =
    new G4UnionSolid("solidSCClamshell_0_bvp", solidSCBackClam_0_ebph, 
  		     solidBackViewPipe, rot_temp, 
  		     G4ThreeVector((SCBackClamOuterRadius+1.501*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset),
  				   (SCBackClamOuterRadius+1.501*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
  				   0.0)
  		     );
  
  G4Tubs* solidBackViewPipeFlange = new G4Tubs("solidBackViewPipeFlange", 
  					       0.0, 3.0*inch, 0.5*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCBackClam_0_bvpf =
    new G4UnionSolid("solidSCClamshell_0_bvpf", solidSCBackClam_0_bvp, 
  		     solidBackViewPipeFlange, rot_temp, 
  		     G4ThreeVector((SCBackClamOuterRadius+2.556*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset), 
  				   (SCBackClamOuterRadius+2.556*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
  				   0.0));
  
  G4Tubs* solidBackViewPipeHole = new G4Tubs("solidBackViewPipeHole", 
  					     0.0, 2.0*inch, 4.0*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_bvph =
    new G4SubtractionSolid("solidSCClamshell_0_bvph", solidSCBackClam_0_bvpf, 
  			   solidBackViewPipeHole, rot_temp, 
  			   G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBackViewPipeAngleOffset),
  					 SCBackClamOuterRadius*cos(90.0*deg+SCBackViewPipeAngleOffset), 
  					 0.0)
  			   );
  
  
  // Placing back clamshell
  G4VSolid* solidSCBackClamShell = solidSCBackClam_0_bvph;
    
  logicScatChamberBackClamshell = new G4LogicalVolume(solidSCBackClamShell, 
  						      fNistManager->FindOrBuildMaterial("G4_Al"), 
  						      "SCBackClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(-90.0*deg+SCClamAngleApert*0.5+SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp, G4ThreeVector(0,0,-fDistTarPivot*cm), logicScatChamberBackClamshell, 
  		    "SCBackClamshell", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  // Scattering chamber volume
  //
  G4Tubs* solidScatChamber_0 = new G4Tubs("SC", 0.0, SCRadius, 0.5* SCHeight, 0.0*deg, 360.0*deg);
  
  G4Tubs* solidSCWindowVacuum = 
    new G4Tubs("SCWindowVacuumFront", SCRadius-1.0*cm, SCTankRadius, 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4UnionSolid* solidScatChamber_0_wbv = 
    new G4UnionSolid("solidScatChamber_0_wbv", solidScatChamber_0, solidSCWindowVacuum,
  		     rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  // G4UnionSolid* solidScatChamber_0_wfv = 
  //   new G4UnionSolid("solidScatChamber_0_wfv", solidScatChamber_0_wbv, solidSCWindowVacuum,
  // 		     rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  G4UnionSolid* solidScatChamber_0_entbp = 
    new G4UnionSolid("solidScatChamber_0_entbp", solidScatChamber_0_wbv, solidEntranceBeamPipeHole, 
  		     rot_temp, G4ThreeVector(0, -SCRadius, SCOffset));
  //solidExitBeamPipeHole
  

  G4UnionSolid* solidScatChamber_0_exbp = 
    new G4UnionSolid("solidScatChamber_0_exbp", solidScatChamber_0_entbp, solidExitBeamPipeHole, 
  		     rot_temp, G4ThreeVector(0, +SCRadius, SCOffset));
  
  G4VSolid* solidScatChamber = solidScatChamber_0_exbp;

  logicScatChamber = new G4LogicalVolume(solidScatChamber, Beamline, "ScatChamber_log");
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamber, "ScatChamberPhys", expHall_log, false, ++nonSDcounter, ChkOverlaps);
  
  G4VisAttributes* Invisible  = new G4VisAttributes(G4Colour(0.,0.,0.)); 
  Invisible->SetVisibility(false);
  G4VisAttributes* colourDarkGrey = new G4VisAttributes(G4Colour(0.3,0.3,0.3)); 
  colourDarkGrey->SetForceWireframe(true);
  G4VisAttributes* colourGrey = new G4VisAttributes(G4Colour(0.7,0.7,0.7)); 
  colourGrey->SetForceWireframe(true);
  G4VisAttributes* colourCyan = new G4VisAttributes(G4Colour(0.,1.,1.)); 
  
  logicScatChamberTank->SetVisAttributes(colourDarkGrey);
  logicScatChamberFrontClamshell->SetVisAttributes(colourGrey);
  logicScatChamberLeftSnoutWindow->SetVisAttributes(Invisible);
  logicScatChamberLeftSnoutWindowFrame->SetVisAttributes(colourCyan);
  logicScatChamberRightSnoutWindow->SetVisAttributes(Invisible);
  logicScatChamberRightSnoutWindowFrame->SetVisAttributes(colourCyan);
  logicScatChamberBackClamshell->SetVisAttributes(colourGrey);
  logicScatChamberExitFlangePlate->SetVisAttributes(colourGrey);
  logicScatChamber->SetVisAttributes(Invisible);

  //--------------------------------------------------------------------------- 
  // Create Septum 
  //--------------------------------------------------------------------------- 

  //  double inch           = 2.54*cm;
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
  G4LogicalVolume* trrap_en1_min = new G4LogicalVolume(testTrap_en1_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_min",0,0,LarmStepLimits);
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
  G4LogicalVolume* trrap_en1_max = new G4LogicalVolume(testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_max",0,0,LarmStepLimits);
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
  G4LogicalVolume* trrap_en2_max = new G4LogicalVolume(testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_max",0,0,LarmStepLimits);
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
  G4LogicalVolume* trrap_en2_min = new G4LogicalVolume(testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_min",0,0,LarmStepLimits);
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
  G4LogicalVolume* trrap_cov_1_up = new G4LogicalVolume(testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_up",0,0,LarmStepLimits);
  trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_1_up=new G4RotationMatrix();
  pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(pRotX90deg_cov_1_up,G4ThreeVector(0,0.02*mm+ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.), 
		    trrap_cov_1_up,"cov1_up",expHall_log,0,++nonSDcounter);



  //LHS Back Bottom Trap
  G4VSolid* testTrap_cov_1_bot = testTrap_cov_1_up;
  G4LogicalVolume* trrap_cov_1_bot = new G4LogicalVolume(testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_bot",0,0,LarmStepLimits);
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
  G4LogicalVolume* trrap_cov_2_up = new G4LogicalVolume(testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_up",0,0,LarmStepLimits);
  trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *pRotX90deg_cov_2_up=new G4RotationMatrix();
  pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(pRotX90deg_cov_2_up,G4ThreeVector(0, 0.02*mm+ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.), 
		    trrap_cov_2_up,"cov2_up",expHall_log,0,++nonSDcounter);

  
  //RHS Front Bottom Trap
  G4VSolid* testTrap_cov_2_bot = testTrap_cov_2_up;
  G4LogicalVolume* trrap_cov_2_bot = new G4LogicalVolume(testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_bot",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap = new G4LogicalVolume(r_testTrap,fNistManager->FindOrBuildMaterial("G4_Al"),"vac Box 1",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap_en1_max = new G4LogicalVolume(r_testTrap_en1_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en1_max",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap_en2_max = new G4LogicalVolume(r_testTrap_en2_max, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_max",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap_en2_min = new G4LogicalVolume(r_testTrap_en2_min, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_en2_min",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap_cov_1_up = new G4LogicalVolume(r_testTrap_cov_1_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_up",0,0,LarmStepLimits);
  r_trrap_cov_1_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_1_up=new G4RotationMatrix();
  r_pRotX90deg_cov_1_up->rotateX(-90.*deg+atan((ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1)));
  new G4PVPlacement(r_pRotX90deg_cov_1_up,G4ThreeVector(0,ymax_sep_en1+pDz_2+fabs(z_sept_en_min1)*(ymax_sep_ex1-ymax_sep_en1)/(z_sept_ex_max1-z_sept_en_max1),0.),
		    r_trrap_cov_1_up,"cov1_up",expHall_log,0,++nonSDcounter);

  
  //RHS Back Bottom Trap
  G4VSolid* r_testTrap_cov_1_bot = r_testTrap_cov_1_up;
  G4LogicalVolume* r_trrap_cov_1_bot = new G4LogicalVolume(r_testTrap_cov_1_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_1_bot",0,0,LarmStepLimits);
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
  G4LogicalVolume* r_trrap_cov_2_up = new G4LogicalVolume(r_testTrap_cov_2_up, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_up",0,0,LarmStepLimits);
  r_trrap_cov_2_up->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_up=new G4RotationMatrix();
  r_pRotX90deg_cov_2_up->rotateX(-90.*deg+atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_up,G4ThreeVector(0,     ymin_sep_en2+pDz_2-fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_up,"cov2_up",expHall_log,0,++nonSDcounter);


  //RHS Front Bottom Trap
  G4VSolid* r_testTrap_cov_2_bot = r_testTrap_cov_2_up;
  G4LogicalVolume* r_trrap_cov_2_bot = new G4LogicalVolume(r_testTrap_cov_2_bot, fNistManager->FindOrBuildMaterial("G4_Al"),"trrap_cov_2_bot",0,0,LarmStepLimits);
  r_trrap_cov_2_bot->SetVisAttributes((G4Colour(0.8, 0.8, 1.0)));
  G4RotationMatrix *r_pRotX90deg_cov_2_bot=new G4RotationMatrix();
  r_pRotX90deg_cov_2_bot->rotateX(-90.*deg-atan((ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2)));
  new G4PVPlacement(r_pRotX90deg_cov_2_bot,G4ThreeVector(0,-1.*ymin_sep_en2-pDz_2+fabs(z_sept_en_min2)*(ymax_sep_ex2-ymax_sep_en2)/(z_sept_ex_max2-z_sept_en_max2),0.),
		    r_trrap_cov_2_bot,"cov2_bot",expHall_log,0,++nonSDcounter);
  
  //------------------------------End of Septum----------------------------------------
  
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
