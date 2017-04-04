

#include "EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
//#include "Hit.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4HCofThisEvent.hh"
#include "HistoManager.hh"

#include <iostream>


EventAction::EventAction(HistoManager * histo) :
G4UserEventAction(),
Edeposit(0.),
stepLength(0.),
totalTime(0.),
labScatteringAngle(0.),
cmScatteringAngle(0.),
histoMan(histo),
numScintPhotons(0)
{}

G4int EventAction::nEvents = 0;
G4int EventAction::nScatteringEvents = 0;
G4int EventAction::nEdepEvents = 0;

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*evt)
{Edeposit = 0;
  stepLength = 0;
  totalTime = 0.;
   cmScatteringAngle = -1 ;
   numScintPhotons = 0.;

  G4int evtNb = evt->GetEventID();
  if((evtNb%1000) == 0) {
	G4cout << "        Processing Event: " << evtNb/1000 << " 000\r" << std::flush;
	//printf("---> Event Number %5d\r",evtNb);
	}
	 

}

void EventAction::EndOfEventAction(const G4Event* evt)
{
    //nEvents += 1;
    //G4SDManager * sdman = G4SDManager::GetSDMpointer();
    
   G4VHitsCollection * HC = evt->GetHCofThisEvent()->GetHC(0);
   G4int hitcount = HC->GetSize();
    
   
    
    nEvents += 1;
    if(hitcount != 0) nScatteringEvents += 1;
    
  //  if(Edeposit > 0) nEdepEvents += 1;

    // G4cout << "distance before interaction: " << G4BestUnit(stepLength,"Length") << G4endl;
    
    //  G4cout << "energy deposit : " << G4BestUnit(Edeposit,"Energy") << G4endl;


    

 		if(Edeposit > 0*CLHEP::MeV)
   		 histoMan->FillHisto(1,Edeposit);

  		  if(stepLength > 0 && stepLength < 15*CLHEP::cm)
   		   histoMan->FillHisto(2,stepLength);

   		 histoMan->FillHisto(3,totalTime);

  		  if(labScatteringAngle > 0)
  		  histoMan->FillHisto(4,labScatteringAngle);

 		   if(cmScatteringAngle != -1)
  		    histoMan->FillHisto(5,cmScatteringAngle);
   		 // G4cout << "this should be the bin content of 1 degree: " << (*(histoMan->GetHisto(5))).GetBinContent(90) << G4endl;
    

   		 //   G4cout << "Total time of this event was: " << G4BestUnit(GetTime(),"Time") << G4endl;
    
    		//  G4cout << " numScintPhotons: " << numScintPhotons << G4endl;
   		 if(numScintPhotons > 0)histoMan->FillHisto(7,numScintPhotons);
  		  histoMan->FillNtuple(Edeposit,labScatteringAngle,totalTime,stepLength,numScintPhotons,photonTime);



}
