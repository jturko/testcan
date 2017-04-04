#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsListMessenger: public G4UImessenger
{
public:
    
    PhysicsListMessenger(PhysicsList* );
    virtual ~PhysicsListMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
private:
    
    PhysicsList* fPPhysicsList;
    
    G4UIdirectory*             fPhysDir;
    G4UIcmdWithADoubleAndUnit* fGammaCutCmd;
    G4UIcmdWithADoubleAndUnit* fElectCutCmd;
    G4UIcmdWithADoubleAndUnit* fProtoCutCmd;
    G4UIcmdWithADoubleAndUnit* fAllCutCmd;
    G4UIcmdWithADoubleAndUnit* fMCutCmd;
    G4UIcmdWithADoubleAndUnit* fECutCmd;
    G4UIcmdWithAString*        fPListCmd;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

