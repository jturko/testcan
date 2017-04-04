// G.A.P. Cirrone, May 2012

/**
 * @file
 * @brief Implements mandatory user class DetectorConstruction.
 */

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

DetectorConstruction::DetectorConstruction() :
scoringVol(0)
{
  //--------- Material definition ---------
  DefineMaterials();

  //--------- Sizes of the principal geometrical components (solids)  ---------
  ComputeParameters();

  detMessenger = new DetectorMessenger(this);
}
 
DetectorConstruction::~DetectorConstruction()
{delete detMessenger;}

void DetectorConstruction::DefineMaterials() 
{
  //Get Materials from NIST database
    
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(0);

  // define NIST materials
    
  vacuum  = man->FindOrBuildMaterial("G4_Galactic");
  air     = man->FindOrBuildMaterial("G4_AIR");
  Germanium = man->FindOrBuildMaterial("G4_Ge");

    
    G4double a, z, density , fractionmass;
     G4String name , symbol;
     G4int nComp;
    
    //hydrogen
    a = 1.01*CLHEP::g/CLHEP::mole;
    G4Element * elH = new G4Element(name="Hydrogen",symbol="H",z=1.,a);
    
    //carbon
    a = 12.01*CLHEP::g/CLHEP::mole;
    G4Element * elC = new G4Element(name="Carbon",symbol="C",z=6.,a);
    
    //deuterium
   a = 2.02*CLHEP::g/CLHEP::mole;
G4Isotope* De = new G4Isotope("De", 1, 2, a);
   G4Element * D = new G4Element(name="Deuterium",symbol="D",nComp=1);
  D->AddIsotope(De,fractionmass=100.*CLHEP::perCent);
    
    //aluminum
    a = 26.981538*CLHEP::g/CLHEP::mole;
    G4Element * elAl = new G4Element(name="Aluminum",symbol="Al",z=13.,a);
    
    density = 2.70*CLHEP::g/CLHEP::cm3;
    Alum = new G4Material("Aluminum",density,1);
    Alum->AddElement(elAl,fractionmass=100.*CLHEP::perCent);
    G4MaterialPropertiesTable * alumMPT = new G4MaterialPropertiesTable();
    const G4int NUMB = 3;
    G4double phoEnergy[NUMB] = {1240/399.9*CLHEP::eV,1240/450.86*CLHEP::eV,1240/506.06*CLHEP::eV};
    G4double refract[NUMB] = { 0.49,0.6203,0.789 };

    alumMPT->AddProperty("RINDEX",phoEnergy,refract,NUMB)->SetSpline(true);
    Alum->SetMaterialPropertiesTable(alumMPT);
    
    
    //Scin
    
    
    density = 0.874*CLHEP::g/CLHEP::cm3;
    Scin = new G4Material("non-deuterated scintillator",density,2);
    Scin->AddElement(elH,fractionmass=9.2497*CLHEP::perCent);
    Scin->AddElement(elC,fractionmass=90.7503*CLHEP::perCent);
    
    //DeuScin
    density = 0.954*CLHEP::g/CLHEP::cm3;
    DeuScin = new G4Material("Deuterated Scintillator",density,3);
    DeuScin->AddElement(elH,fractionmass=0.0625*CLHEP::perCent);
    DeuScin->AddElement(elC,fractionmass=85.7326*CLHEP::perCent);
    DeuScin->AddElement(D,fractionmass=14.2049*CLHEP::perCent);
    
    G4MaterialPropertiesTable * MPT = new G4MaterialPropertiesTable();
    const G4int NUM = 6;
    G4double photon_energies[NUM] = {3.1*CLHEP::eV,2.88*CLHEP::eV,2.82*CLHEP::eV,2.695*CLHEP::eV,2.58*CLHEP::eV,2.48*CLHEP::eV};
    G4double emission_spectra[NUM] = {0.05,1.,0.7,0.37,0.2,0.1};
    MPT->AddConstProperty("RINDEX",1.498);
    MPT->AddConstProperty("SCINTILLATIONYIELD",9200./CLHEP::MeV); 

    G4double electron_energy[4] = {1.*CLHEP::keV,1.*CLHEP::MeV,10.*CLHEP::MeV,100*CLHEP::MeV};
    G4double electron[4] = {1.,9200.,92000.,920000.};
     MPT->AddProperty("ELECTRONSCINTILLATIONYIELD",electron_energy,electron,4);
     // MPT->AddProperty("DEUTERONSCINTILLATIONYIELD",1.*CLHEP::MeV,2000./CLHEP::MeV,1);
    /*
    MPT->AddConstProperty("IONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("DEUTERONSCINTILLATIONYIELD",9200./CLHEP::MeV);
    MPT->AddConstProperty("TRITONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("ALPHASCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("PROTONSCINTILLATIONYIELD",9200./CLHEP::MeV);
    MPT->AddConstProperty("ELECTRONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    */
    MPT->AddConstProperty("RESOLUTIONSCALE",25.);
    

    MPT->AddConstProperty("ABSLENGTH",3.*CLHEP::m);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // values for fast and slow time component taken from :                                              //
    //     http://research.physics.lsa.umich.edu/twinsol/Publications/Ojaruegathesis.pdf                 //
    //                                                                                                   //
    //     Pulse shape analysis of liquid scintillators for neutron studies,  S. Marrone et al.          //
    MPT->AddConstProperty("FASTTIMECONSTANT",3.5*CLHEP::ns);
    MPT->AddProperty("FASTCOMPONENT",photon_energies,emission_spectra,NUM)->SetSpline(true);
    // MPT->AddConstProperty("SLOWTIMECONSTANT",50.7*CLHEP::ns);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    MPT->AddConstProperty("YIELDRATIO",1.);
    DeuScin->SetMaterialPropertiesTable(MPT);
    
    //CARBON
    density = 12.01*CLHEP::g/CLHEP::cm3;
    CARBON = new G4Material("CARBON",density,1);
    CARBON->AddElement(elC,fractionmass=100.*CLHEP::perCent);
    
    //HYDROGEN
    density = 1.0079*CLHEP::g/CLHEP::cm3;
    HYDROGEN = new G4Material("HYDROGEN",density,1);
    HYDROGEN->AddElement(elH,fractionmass=100.*CLHEP::perCent);
    
    //DEUTERIUM
    density = 2.014*CLHEP::g/CLHEP::cm3;
    DEUTERIUM = new G4Material("DEUTERIUM",density,1);
    DEUTERIUM->AddElement(D,fractionmass=100.*CLHEP::perCent);
    
    

}
 
void DetectorConstruction::ComputeParameters() 
{
  //This function defines the defaults
  //of the geometry construction

  // ** world **
  halfWorldLength = 8* CLHEP::m;
    
    //test can
    
    posTestCan = G4ThreeVector(0.0,0.0,50.0*CLHEP::cm);
    radTestCan = 5.55*CLHEP::cm;
    zTestCan = (2.22*CLHEP::cm)/2.;
    
    caseThickness = (0.16*CLHEP::cm)/2;

    caseShift = 48.78*CLHEP::cm;
    caseUp = 0;
    posCaseing = G4ThreeVector(0,caseUp,caseShift);
    
    shellRadius = 5.71*CLHEP::cm;
    
    
    

}


G4VPhysicalVolume * DetectorConstruction::Construct(){return ConstructTestCan();} 

G4VPhysicalVolume* DetectorConstruction::ConstructTestCan()
{
  //This function is called by G4 when the detector has to be created
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  
  //------------------------------
  // World
  //------------------------------

 G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
 
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);

  G4Box * solidWorld= new G4Box("world",halfWorldLength,halfWorldLength,halfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, vacuum, "World", 0, 0, 0);
 
  G4VPhysicalVolume * physiWorld = new G4PVPlacement(0,               // no rotation
						     G4ThreeVector(), // at (0,0,0)
						     logicWorld,      // its logical volume
						     "World",         // its name
						     0,               // its mother  volume
						     false,           // no boolean operations
						     0);              // copy number
				
    
    //ConstructTestCan();
    
    G4double ang1, ang2;
    ang1 = 0;
    ang2 = 360*CLHEP::deg;
    
    
    G4Tubs * physTestCan = new G4Tubs("PhysicalCan",0.,radTestCan,zTestCan,ang1,ang2);
    
    YourChoice = DeuScin;
    // ...000000000000000000000000000000000000000000000000000000000000000000000000...
    logicTestCan = new G4LogicalVolume(physTestCan,YourChoice,"LogicalTestCan");
    // ...000000000000000000000000000000000000000000000000000000000000000000000000...
    
    G4VPhysicalVolume * testcan = new G4PVPlacement(0,posTestCan,logicTestCan,"Test Can",logicWorld,false,0);
    
    
    
    
    
    
   G4Tubs * caseTestCan = new G4Tubs("CaseCan",0.,shellRadius,caseThickness,ang1,ang2);
   logicCase = new G4LogicalVolume(caseTestCan,Alum,"LogicalCase");
   //  G4VPhysicalVolume * caseing = new G4PVPlacement(0,posCaseing,logicCase,"Caseing",logicWorld,false,0);
     caseing = new G4PVPlacement(0,posCaseing,logicCase,"Caseing",logicWorld,false,0);
    
    
    G4Tubs * shellTestCan = new G4Tubs("ShellCan",radTestCan,shellRadius,zTestCan,ang1,ang2);
    G4LogicalVolume * logicShell = new G4LogicalVolume(shellTestCan,Alum,"LogicalShell");
    
    G4VPhysicalVolume * shell = new G4PVPlacement(0,posTestCan,logicShell,"Shell",logicWorld,false,0);
    
    
    // G4Tubs * physTestCannon = new G4Tubs("PhysicalCannon",0.,radTestCan/5,2.5*CLHEP::cm,ang1,ang2);
    
    // G4LogicalVolume * logicCannon = new G4LogicalVolume(physTestCannon,air,"LogicalTestCan");
    
    //logicTestCan->SetVisAttributes(new G4VisAttributes(green));
    //G4VPhysicalVolume * launcher = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicCannon,"Cannon",logicWorld,false,0);
    
    
    G4MultiFunctionalDetector * sensTestCan = new G4MultiFunctionalDetector("SensitiveTestCan");
    G4SDManager * sdman = G4SDManager::GetSDMpointer();
    sdman->AddNewDetector(sensTestCan);
    logicTestCan->SetSensitiveDetector(sensTestCan);
    
    G4VPrimitiveScorer * eDeposit = new G4PSEnergyDeposit("Energy Deposit",0);
    //eDeposit->SetMultiFunctionalDetector(sensTestCan);
    sensTestCan->RegisterPrimitive(eDeposit);
    
    scoringVol = logicTestCan;

  //--------- Visualization attributes -------------------------------

  G4Color
    green(0.0,1.0,0.0),
    blue(0.0,0.0,1.0),
    brown(1.0,1.0,0.1),
    white(1.0,1.0,1.0),
    magenta(1.0, 0.0, 1.0),
    red(1.0,0.0,0.0);

  logicWorld -> SetVisAttributes(new G4VisAttributes(white));
  logicWorld -> SetVisAttributes(G4VisAttributes::Invisible);
    logicTestCan->SetVisAttributes(new G4VisAttributes(blue));
    // logicCannon->SetVisAttributes(new G4VisAttributes(brown));
    logicCase->SetVisAttributes(new G4VisAttributes(red));
    //logicShell->SetVisAttributes(new G4VisAttributes(brown));
    
  //always return the physical World
  //
  return physiWorld;
}


// Setting the Detector Material
void DetectorConstruction::SetTestCanMaterial(G4String materialChoice){

  if(materialChoice == "Deuterium"){YourChoice = DEUTERIUM;}
  else if(materialChoice == "Hydrogen"){YourChoice = HYDROGEN;}
  else if(materialChoice == "Carbon"){YourChoice = CARBON;}
  else if(materialChoice == "BC501A"){YourChoice = Scin;}
  else if(materialChoice == "BC537"){YourChoice = DeuScin;}
  else if(materialChoice == "Germanium"){YourChoice = Germanium;}
  else{ G4cout << " that wasnt right..." << G4endl; return; }

  if(logicTestCan){
    logicTestCan->SetMaterial(YourChoice);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  }
}

//Setting if the aluminum shielding is present or not
void DetectorConstruction::SetAlumCase(G4bool choice)
{

  if(choice == true){
    //if(caseing == NULL) caseing = new G4PVPlacement(0,posCaseing,logicCase,"Caseing",logicWorld,false,0);
    // G4cout << "\n\n the casing was just created!!! \n";
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
    caseUp = 0;
    posCaseing = G4ThreeVector(0,0,caseShift); 
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
	//this->UpdateGeometry();
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  }
  else if(choice == false){
    caseUp = 20*CLHEP::cm;
    posCaseing = G4ThreeVector(0,caseUp,caseShift);
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
	//this->UpdateGeometry();
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  }

  else G4cout << "this is the final else... something is weird..." << G4endl;
}

// Setting the thickness of the detector
void DetectorConstruction::SetzTestCan(G4double z)
{
  zTestCan = z/2;
  caseShift = (50*CLHEP::cm - (zTestCan) - 0.1*CLHEP::cm);
  posCaseing = G4ThreeVector(0.0,0.0,caseShift);
  //G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4RunManager::GetRunManager()->GeometryHasBeenModified();
	//this->UpdateGeometry();
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}



