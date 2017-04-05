

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//if you change this you also must change HistoManager.cc
#define MAX_PHOTONS 100000

class HistoManager;




class EventAction : public G4UserEventAction
{
public:
    EventAction(HistoManager *);
    virtual ~EventAction();
    
    virtual void BeginOfEventAction(const G4Event* );
    virtual void EndOfEventAction(const G4Event* );
    
    static G4int GetnEvents(){return nEvents;}
    static G4int GetnScatteringEvents(){return nScatteringEvents;}
    static G4int GetnEdepEvents(){return nEdepEvents;}
  static void OneMoreScatteringEvent(){nScatteringEvents += 1;}
    
  static void ResetEvents(){nEvents = 0; nScatteringEvents = 0; nEdepEvents = 0;}
    
    
  void AddEdep(G4double ener) { Edeposit += ener; }
  void AddStep(G4double stepl) {
    if(!stepLength) stepLength += stepl; }
  void AddTime(G4double time) { totalTime += time;}
  G4double GetTime(){return totalTime;}

   G4double labScatteringAngle;
  G4double getLabScatAngle(){return labScatteringAngle;}
  void setLabScatAngle(G4double angle){ labScatteringAngle = angle; }

  G4double cmScatteringAngle;
  G4double getcmScatAngle(){return cmScatteringAngle;}
  void setcmScatAngle(G4double angle){ cmScatteringAngle = angle;}

  G4double NumScintPhotons() {return numScintPhotons;}
  void OneMoreScintPhoton(G4double time){
    // if(numScintPhotons<MAX_PHOTONS)
    // photonTime[numScintPhotons] = time;
    numScintPhotons++;
  }
void OneLessScintPhoton(){
   numScintPhotons = numScintPhotons - .8;
}
 
    void SetPrimaryEnergy(G4double energy) { primaryEnergy = energy; }
    G4double GetPrimaryEnergy() { return primaryEnergy; }
    
private:
   static G4int nEvents;
   static G4int nScatteringEvents;
   static G4int nEdepEvents;
    
  G4double Edeposit;
  G4double stepLength;
  G4double totalTime;
  G4double numScintPhotons;
  G4double photonTime[MAX_PHOTONS];

    G4double primaryEnergy;

  HistoManager * histoMan;
    
};



#endif
