


#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class DetectorMessenger: public G4UImessenger
{
 
public:
  DetectorMessenger(DetectorConstruction*);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*,G4String);


private:
  DetectorConstruction * detectorConstruct;

  G4UIdirectory * testcanDir;
  G4UIdirectory * detectorDir;
  G4UIcmdWithAString * detectorMaterialCmd;
   G4UIcmdWithABool * detectorAlumCmd;
  G4UIcmdWithADoubleAndUnit * canThicknessCmd;

};

#endif
