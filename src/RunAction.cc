
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"


#include "G4RunManager.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include <math.h>


G4int RunAction::j = 0;

//void RunAction::SetCrossArray(G4double crosSect, G4int i){
//    this->CrossArray[i]=crosSect;
//}

//G4double RunAction::GetCrossArray(G4int i){
//    return CrossArray[i];
//}


void RunAction::BeginOfRunAction(const G4Run * run)
{
    EventAction::ResetEvents();
    histoMan->book();
	G4int runID = run->GetRunID() + 1;
	G4cout << "------> Starting Run " << runID << G4endl;

	//  	Check if we are using the correct detector material and initial particle energy
	//			DOESN'T WORK HERE!!!!!!!!!!!!!!!!!!!!
	//     the BeginOfRunAction is called BEFORE the UI commands are sent, therefore this will just give the energy and material of the PREVIOUS run... I think?

	//G4double energy = primaryGen->GetGPSEnergy();
	//G4cout << "Initial Paricle Energy: " << G4BestUnit(energy,"Energy") << G4endl;
	//G4String detectorMaterial = detConstruction->GetDetectorMaterialName();
	//G4cout << "Detector Material: " << detectorMaterial << G4endl;

}

void RunAction::EndOfRunAction(const G4Run * run )
{
	G4int runID = run->GetRunID() + 1;
	G4cout << "\n\tRun " << runID << " Summary:			" << G4endl;

	G4double energy = primaryGen->GetGPSEnergy();
	G4cout << "Initial Paricle Energy: " << G4BestUnit(energy,"Energy") << G4endl;
	G4String detectorMaterial = detConstruction->GetDetectorMaterialName();
	G4cout << "Detector Material: " << detectorMaterial << G4endl;

	SetEnergyVector(energy);

	// FOR THE CALCULATION OF THE TOTAL CROSS SECTION
	// THIS WILL ONLY WORK FOR A 1 cm THICK TARGET WITH SPECIFICALLY DEFINED DENSITY: SEE JAMES WONG MSC THESIS FOR DETAILS
    G4double scat = EventAction::GetnScatteringEvents();
    G4double inc = EventAction::GetnEvents();
    
    G4double avo = 6.022e+23;
    G4double ratio;
    ratio = (scat/inc);
  
    G4double nlog = double (log(1-ratio));
    G4double crossSect = -nlog/avo;
    G4double bCross = crossSect*1e+24;
    G4double errCross = (1e+24/avo)*(1/(inc-scat))*sqrt(scat);

	//G4cout << "ratio: " << ratio << G4endl;
	G4cout << "Cross Section: " << bCross << G4endl;
    
    SetCrossArray(bCross,j);
    SetHitArray(scat,j);
    SeterrArray(errCross,j);
    j++;
    
    G4double num = 0;
  G4int i = 0;
    for( i = 1; i < 180; i++){
      
      num = ((*(histoMan->GetHisto(5))).GetBinContent(i))/sin(i*3.14159265359/180);
      //HYDROGEN ELASTIC CROSS SECTION 10 MeV
      num = (num/inc)*2.05101;
      (*(histoMan->GetHisto(6))).SetBinContent(i,num);

    }
    

    // histoMan->GetHisto(6)->Scale(1/(180*1.05276),"");
    histoMan->save();

	numberOfRuns += 1;
    
}
