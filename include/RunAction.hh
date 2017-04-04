#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class G4Run;

class HistoManager;
class PrimaryGeneratorAction;
class DetectorConstruction;

class RunAction : public G4UserRunAction
{
public:
  RunAction(HistoManager* hist, PrimaryGeneratorAction* pga, DetectorConstruction* detConst) : G4UserRunAction(), histoMan(hist), primaryGen(pga), detConstruction(detConst), numberOfRuns(0)  {}
    virtual ~RunAction(){}
    
    //virtual G4Run * GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
     //000000000000000000000000oooooooo0000000000000000000000
    
     void SetCrossArray(G4double crosSect, G4int i) {
        CrossArray[i]=crosSect;
	CrossVector.push_back(crosSect);}

     G4double GetCrossArray(G4int i) {
        return CrossArray[i];}

    G4double GetCrossVector(G4int i){
	G4double crossSection = CrossVector.at(i);
	return crossSection;}
    //000000000000000000000000oooooooo0000000000000000000000
    
    void SetHitArray(G4int hits,G4int i) {
        HitArray[i]=hits;
	HitVector.push_back(hits);}

    G4int GetHitArray(G4int i) {
        return HitArray[i];}

    G4int GetHitVector(G4int i) {
	G4int hits = HitVector.at(i);	
	return hits;}
    
     //000000000000000000000000oooooooo0000000000000000000000
    
    void SeterrArray(G4double errCrosSect,G4int i) {
        errArray[i]=errCrosSect;
	errVector.push_back(errCrosSect);}

    G4double GeterrArray(G4int i) {
        return errArray[i];}

    G4double GeterrVector(G4int i) {
	G4double error = errVector.at(i);
	return error;}
     //000000000000000000000000oooooooo0000000000000000000000
    
    void SetMeanPathArray(G4double meanPath,G4int i) {
        MeanPathArray[i]=meanPath;
    }
    G4double GetMeanPathArray(G4int i) {
        return MeanPathArray[i];
    }
    
     //000000000000000000000000oooooooo0000000000000000000000

    void SetEnergyVector(G4double energy) {
	G4double GPSenergy = energy;
	energyVector.push_back(GPSenergy);}

   G4double GetEnergyVector(G4int i){
	G4double energy = energyVector.at(i);
	return energy;}

     //000000000000000000000000oooooooo0000000000000000000000

    G4int GetVectorSize(){
	G4int vSize = CrossVector.size();
	return vSize;}

    G4int GetNumberOfRuns(){return numberOfRuns;}

     //000000000000000000000000oooooooo0000000000000000000000
    
    
private:
    G4double CrossArray[20];
    G4int HitArray[20];
    G4double errArray[20];
    G4double MeanPathArray[20];

	// Vectors storing data needed to calculate cross sections
    std::vector<G4double> CrossVector;
    std::vector<G4int> HitVector;
    std::vector<G4double> errVector;
    std::vector<G4double> energyVector;
    
  HistoManager * histoMan;
  PrimaryGeneratorAction * primaryGen;
  DetectorConstruction * detConstruction;

   G4int numberOfRuns;

    static G4int j;
};


#endif

