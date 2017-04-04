
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4LogicalVolume;


class SteppingAction : public G4UserSteppingAction
{
    
public:
    //SteppingAction(){}
    SteppingAction(EventAction*);
    virtual ~SteppingAction();
    
    virtual void UserSteppingAction(const G4Step*);
    
private:
    EventAction * evtAction;
    G4LogicalVolume * logicVol;
    G4double Pi;
    
};

#endif
