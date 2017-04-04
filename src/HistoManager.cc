#include <algorithm>

 #include "TH1D.h"
 #include "TFile.h" 
 #include "TTree.h"
 #include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandGauss.h"

 #include "HistoManager.hh"
 #include "G4UnitsTable.hh"

HistoManager::HistoManager() : 
  rootFile(0),Ntuple(0),eDeposit(0),labAngle(0),totalTime(0),pathLength(0),scintLight(0)
{
  for(G4int i=0;i<MaxHisto;i++){histo[i]=0;}
  Ntuple = 0;
}

HistoManager::~HistoManager()
{
  if (rootFile) delete rootFile;
}

void HistoManager::book()
{
  G4String fileName = "testcan.root";
  rootFile = new TFile(fileName,"RECREATE");
  if(!rootFile){
 G4cout << "Problem booking the histogram!!" << G4endl;
    return;
               } 

 histo[1] = new TH1D("h1","Edep in testcan",400,0*CLHEP::MeV,10.*CLHEP::MeV);
 if(!histo[1]) G4cout << "\n can't create hist 1..." << G4endl;

 histo[2] = new TH1D("h2","mean free path",100,0.,15*CLHEP::cm);
 if(!histo[2]) G4cout << "\n can't create histo 2..." << G4endl;

 histo[3] = new TH1D("h3","total time of event",5000,0.,10*CLHEP::nanosecond);
 if(!histo[3]) G4cout << " \n can't create histo 3..." << G4endl;

 histo[4] = new TH1D("h4","lab scattering angle",180 ,0,181);
 if(!histo[4]) G4cout << " \n cant create histo 4..." << G4endl;

 histo[5] = new TH1D("h5","maybe CM scattering angle...?",180,0,181);
 if(!histo[5]) G4cout << " \n cant create histo 5 ..." << G4endl;

 histo[6] = new TH1D("h6","cm scattering angle / sin(theta)",180,0,181);
 if(!histo[6]) G4cout << " \n cant create histo 6 .. " << G4endl;

 histo[7] = new TH1D("h7","scint light",10000,0,100000);


 Ntuple = new TTree("tree","rootTree");
 Ntuple->Branch("eDep",&eDeposit,"eDeposit/D");
 Ntuple->Branch("lAng",&labAngle,"labAngle/D");
 Ntuple->Branch("tTime",&totalTime,"totalTime/D");
 Ntuple->Branch("pLength",&pathLength,"pathLength/D");
 Ntuple->Branch("scintLight",&scintLight,"scintLight/I");
 Ntuple->Branch("photonTime",photonTime,"photonTime[scintLight]/D");


//G4cout << "\n\n -------> histogram file is opened in " << fileName;
G4cout << G4endl;

}


void HistoManager::save()
{
  if(rootFile) {
    rootFile->Write();
    rootFile->Close();
    //G4cout << "\n------>>>>>> Histogram Tree is saved...? Lets hope so\n";
  }
}


void HistoManager::FillHisto(G4int ih,G4double xbin,G4double weight)
{
  if (ih >= MaxHisto){
    G4cout << "HistoManager::FillHisto -> histo" << ih << "doesnt exist";
    G4cout << G4endl;
      return;
  }
  if (histo[ih]) { histo[ih]->Fill(xbin,weight);}
}

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if(ih>=MaxHisto) {
    G4cout << "---> HistoManager::Normalize : histo " << ih
	   << "does not exist" << G4endl;
  }
  if(histo[ih]) histo[ih]->Scale(fac);
}

void HistoManager::FillNtuple(G4double energy, G4double angle, G4double time, G4double path, G4int sLight, G4double* pTime){

  if(energy > 0 && time > 0){
    eDeposit = energy;
    labAngle = angle;
    totalTime = time;
    scintLight = sLight;

    // scintLight = sLight +  CLHEP::RandGauss::shoot(0, sqrt(sLight));
    //std::copy(pTime,pTime+MAX_PHOTONS,photonTime);
    pathLength = path;}
  else if(time > 0){ labAngle = angle;
    totalTime = time;
    scintLight = sLight;

    //scintLight = sLight +  CLHEP::RandGauss::shoot(0, sqrt(sLight));
    //std::copy(pTime,pTime+MAX_PHOTONS,photonTime);
    pathLength = path;}
  else{ labAngle = angle;
    scintLight = sLight;

    //scintLight = sLight +  CLHEP::RandGauss::shoot(0, sqrt(sLight));
    //std::copy(pTime,pTime+MAX_PHOTONS,photonTime);
  }



  if(Ntuple) Ntuple->Fill();

}


void HistoManager::PrintStatistic()
{
  if(histo[1])
    {
      G4cout << "\n\n -----> printing histogram statistics" << G4endl;
      G4cout << "Edep: mean = " << G4BestUnit(histo[1]->GetMean(),"Energy");
      G4cout << "       rms = " << G4BestUnit(histo[1]->GetRMS(),"Energy") 
	     << G4endl;
   
      G4cout << "Mean Free Path : mean = " << G4BestUnit(histo[2]->GetMean(),"Length");
      G4cout << "                  rms = " << G4BestUnit(histo[2]->GetRMS(),"Length") 
	     << G4endl;



    }

}












