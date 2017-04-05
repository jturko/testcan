
#include "G4Step.hh"
                      #include "G4StepPoint.hh"
                      #include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "PhysicsUtilities.hh"
#include "G4TrackVector.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include <math.h>  
#include <stdlib.h>
#include <time.h>


SteppingAction::SteppingAction(EventAction * evtact) :
G4UserSteppingAction(),
evtAction(evtact),
logicVol(0)
{}

SteppingAction::~SteppingAction()
{G4cout << "im dead!" << G4endl;}

void SteppingAction::UserSteppingAction(const G4Step * step)
{
  if (!logicVol)
    {
      const DetectorConstruction * detectorConstruct =
	static_cast<const DetectorConstruction*>
	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      logicVol = detectorConstruct->GetScoringVolume();
    }
    
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
  if (volume != logicVol) return;

  G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  // G4cout << "process: " << processName << G4endl;
  G4String HadronicElastic = "hadElastic";

  G4int numSecondaries = step->GetSecondary()->size();
    
  G4double runAngle = evtAction->getcmScatAngle();
  //  G4cout << "whats my angle? should be -1...... but it actually is: " << runAngle << G4endl;

    G4double pEnergy = evtAction->GetPrimaryEnergy();

  if(processName == HadronicElastic && numSecondaries == 1 && runAngle == -1 && step->GetPreStepPoint()->GetKineticEnergy()/MeV == pEnergy)
  //if(processName == HadronicElastic && numSecondaries == 1 && runAngle == -1)
    { 

      G4ThreeVector PreScatMomentum = step->GetPreStepPoint()->GetMomentum();
      // G4cout << "the track of the step is a : " << step->GetTrack()->GetDefinition()->GetParticleName()<< G4endl;
      G4ThreeVector PostScatMomentum = step->GetPostStepPoint()->GetMomentum();
      G4double rad = PostScatMomentum.angle(PreScatMomentum);
      G4double deg = (rad*180)/3.14159265359;

      G4double T_pre = step->GetPreStepPoint()->GetKineticEnergy();
      G4double T_post = step->GetPostStepPoint()->GetKineticEnergy();
      //  G4cout << "kinetic energy after step: " << G4BestUnit(T_post,"Energy") << G4endl;
      // G4cout << "lab angle: " << deg << G4endl;


      const G4TrackVector * trackVector = step->GetSecondary();
  
      // G4cout << "number of secondaries generated: " << numSecondaries << G4endl;

      G4double massProjectile = step->GetPreStepPoint()->GetMass();
      G4double massRecoil = (*trackVector)[0]->GetDefinition()->GetPDGMass();
      // G4cout << "recoil mass: " << massRecoil << G4endl;
      // G4cout << "mass: " << mass << G4endl;
      G4double cosineCM_Angle = Physics::CosineThetaEjectileCm(T_pre,massProjectile,massRecoil,massProjectile,massRecoil,0.,rad);
      G4double cmAngle = 180*acos(cosineCM_Angle)/3.14159265359;
     
      // G4cout << "----------------->>>>>>>>>>>>>> cosine of CM scattering angle: " << cosineCM_Angle << G4endl;
      // G4cout << "----------------->>>>>>>>>>>>>> the CM scattering angle      : " << cmAngle << G4endl;
     	

      if(T_post > 0){
	evtAction->setLabScatAngle(deg);
	evtAction->setcmScatAngle(cmAngle);
	//evtAction->OneMoreScatteringEvent();
      
    
      }
   
    
      G4double energy = -(step->GetDeltaEnergy());
      // G4cout << "energy of the step: " << step->GetDeltaEnergy() << G4endl;
      evtAction->AddEdep(energy);

      G4double stepl = step->GetStepLength();
      evtAction->AddStep(stepl);

      G4double deltaTime = step->GetDeltaTime();
      evtAction->AddTime(deltaTime);
  
  
    }


  G4Track * track = step->GetTrack();
  G4String ParticleName = track->GetDynamicParticle()->
    GetParticleDefinition()->GetParticleName();
  G4String C12 = "C12";

  //srand(time(NULL));
  //G4int rand = rand() % 100;

  // G4cout << rand << G4endl;

  // killing the track for c12
 
      //if(ParticleName!="opticalphoton") G4cout << "the particles name is: " << ParticleName << G4endl;

  const std::vector<const G4Track*>* secondaries =
    step->GetSecondaryInCurrentStep();

  if (secondaries->size()>0) {//G4cout << "step 1 \n";
    for(unsigned int i=0; i<secondaries->size(); ++i) {//G4cout << "step 2 \n"; 
      if (secondaries->at(i)->GetParentID()>0) {//G4cout << "step 3 \n";
	if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
	   == G4OpticalPhoton::OpticalPhotonDefinition()){
	  if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
	      == "Scintillation")evtAction->OneMoreScintPhoton(track->GetLocalTime());
	 if(ParticleName == C12)
	   {
	     evtAction->OneLessScintPhoton();
	   }
	  track->SetTrackStatus(fKillTrackAndSecondaries);
	   
	}}}}
}

