

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicsConstructor;
class PhysicsListMessenger;
class G4ProductionCuts;
class G4Scintillation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:
    PhysicsList();
    virtual ~PhysicsList();
    
    virtual void ConstructParticle();
    
    virtual void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    
    void SelectPhysicsList(const G4String& name);
    virtual void ConstructProcess();
   
    void ConstructOp();

    void SetTargetCut(G4double val);
    void SetDetectorCut(G4double val);
    
private:
    
    void AddExtraBuilders(G4bool flagHP);
    
    // hide assignment operator
    PhysicsList & operator=(const PhysicsList &right);
    PhysicsList(const PhysicsList&);
    
    G4double fCutForGamma;
    G4double fCutForElectron;
    G4double fCutForPositron;
    
    G4VPhysicsConstructor*  fEmPhysicsList;
    G4VPhysicsConstructor*  fRaddecayList;
    G4VPhysicsConstructor*  fParticleList;
    G4VPhysicsConstructor*  fHadPhysicsList;
    
    std::vector<G4VPhysicsConstructor*>  fHadronPhys;
    G4int fNhadcomp;
    
    PhysicsListMessenger* fPMessenger;
    G4ProductionCuts* fDetectorCuts;
    G4ProductionCuts* fTargetCuts;

  G4Scintillation * scintProcess;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

