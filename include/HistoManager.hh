

#ifndef HistoManager_h
#define HistoManager_h 1


#include "globals.hh"
#include "TH1.h"

class TFile;
class TTree;
//class TH1D;

//if you change this you also must change EventAction.cc
#define MAX_PHOTONS 100000


const G4int MaxHisto = 8;



class HistoManager
{
public:
  HistoManager();
  ~HistoManager();

  void book();
  void save();

  void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
  void Normalize(G4int id, G4double fac);
  void FillNtuple(G4double energy, G4double angle, G4double time, G4double path, G4int scintPhotons, G4double* photonTime);

  void PrintStatistic();
  TH1D * GetHisto(G4int i){return histo[i];}

private:

  TFile * rootFile;
  TH1D* histo[MaxHisto];

  TTree * Ntuple;

  G4double eDeposit;
  G4double labAngle;
  G4double totalTime;
  G4double pathLength;
  G4int scintLight;
  G4double photonTime[MAX_PHOTONS];
};

#endif
