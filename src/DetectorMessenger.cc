

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
//#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"


DetectorMessenger::DetectorMessenger( DetectorConstruction * detector )
: G4UImessenger(), 
  detectorConstruct(detector),
  testcanDir(0),
  detectorDir(0), 
  detectorMaterialCmd(0),
  detectorAlumCmd(0),
  canThicknessCmd(0)
{

  testcanDir = new G4UIdirectory("/testcan/");
  testcanDir->SetGuidance("UI commands for the testcan setup");

  G4bool broadcast = false;
  detectorDir = new G4UIdirectory("/testcan/det/",broadcast);
  detectorDir->SetGuidance("detector control is yours commander...");

  detectorMaterialCmd = new G4UIcmdWithAString("/testcan/det/setDetectorMaterial",this);
  detectorMaterialCmd->SetGuidance("set the material of the testcan:   Carbon, Deuterium, Hydrogen, BC501A, BC537, Germanium");
  detectorMaterialCmd->SetParameterName("choice",false);
  detectorMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  canThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testcan/det/setDetectorThickness",this);
  canThicknessCmd->SetGuidance("set the length of the testcan");
  canThicknessCmd->SetParameterName("Size",false);
  canThicknessCmd->SetRange("Size>0.");
  canThicknessCmd->SetUnitCategory("Length");
  canThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  detectorAlumCmd = new G4UIcmdWithABool("/testcan/det/setDetectorAlum",this);
  detectorAlumCmd->SetGuidance("set the alum casing on or off");
  detectorAlumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

DetectorMessenger::~DetectorMessenger()
{

  delete canThicknessCmd;
  delete detectorMaterialCmd;
  delete detectorAlumCmd;
  delete detectorDir;
  delete testcanDir;

}



void DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{

  if(command == canThicknessCmd){
    detectorConstruct->SetzTestCan(canThicknessCmd->GetNewDoubleValue(newValue));}
    
  if(command == detectorMaterialCmd){
    detectorConstruct->SetTestCanMaterial(newValue);}

    if(command == detectorAlumCmd){
      detectorConstruct->SetAlumCase(detectorAlumCmd->GetNewBoolValue(newValue));}
}
