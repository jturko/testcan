// $Id: DetectorConstruction.hh 33 2010-01-14 17:08:18Z adotti $
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

/**
 * @file
 * @brief Defines mandatory user class DetectorConstruction.
 */

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*!
\brief This mandatory user class defines the geometry.

It is responsible for
 - Definition of material, and
 - Construction of geometry

\sa Construct()
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //! Constructor
  DetectorConstruction();
  //! Destructor
  ~DetectorConstruction();
public:
  //! Construct geometry of the setup
  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume * ConstructTestCan();
  G4VPhysicalVolume * GetAlum() { return caseing;}

    
    //Set and Get functions for the testcan position
    G4ThreeVector TestCanPosition() const { return posTestCan; }
    G4ThreeVector SetCanPosition(const G4ThreeVector & pos) {
    return posTestCan = pos; }
    
    //Set and get functions for the testcan aluminum casing position
    G4ThreeVector CaseingPosition() const { return posCaseing; }
    G4ThreeVector SetCaseingPosition(const G4ThreeVector & pos) {
    return posCaseing = pos; }
    
    //Set and get functions for the z length of the detector ( the function for SetzTestCan is located in DetectorConstruction.cc
    G4double GetzTestCan(){return zTestCan;}
    void SetzTestCan(G4double z);
    
    // Unused scoring volume get function
    G4LogicalVolume * GetScoringVolume() const { return scoringVol;}

    // functions for setting the material type and getting its name to verify that we are using the correct material
    G4String GetDetectorMaterialName(){G4String name = YourChoice->GetName();
	return name;}
    void SetTestCanMaterial(G4String);

  void SetAlumCase(G4bool);
    

  //@}
private:
  //! define needed materials
  void DefineMaterials();
  //! initialize geometry parameters
  void ComputeParameters();

private:

    G4Material * vacuum;
    G4Material* air;
    G4Material * Scin;
    G4Material * DeuScin;
    G4Material * Alum;
    
    G4Material * CARBON;
    G4Material * DEUTERIUM;
    G4Material * HYDROGEN;
    G4Material * Germanium;

  G4Material * YourChoice;
  DetectorMessenger * detMessenger;
    
    G4double radTestCan;
    G4double zTestCan;
    G4ThreeVector posTestCan;
  
    G4ThreeVector posCaseing;
    G4double caseThickness;
    G4double caseShift;
  G4double caseUp;
    
    G4double shellRadius;
    

    G4LogicalVolume * logicWorld;
    G4double halfWorldLength;
    
    //for stepping action class

    G4LogicalVolume * scoringVol;
  G4LogicalVolume * logicTestCan;
  G4LogicalVolume * logicCase;
  G4VPhysicalVolume * caseing;

  };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
