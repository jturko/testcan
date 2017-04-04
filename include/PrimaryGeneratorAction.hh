// $Id: PrimaryGeneratorAction.hh 94 2010-01-26 13:18:30Z adotti $

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

/**
 * @file
 * @brief Defines mandatory user class PrimaryGeneratorAction.
 */

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include <fstream>


class G4VPrimaryGenerator;
class G4GeneralParticleSource;
 
/*!
\brief This mandatory user class provides the primary particle generator

Geant4 provides a number of predefined primary particle generator, to be utilised by the user.
 - G4ParticleGun
 - G4GeneralParticleSource

\sa GeneratePrimaries()
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  //! constructor
  PrimaryGeneratorAction();
  //! destructor
  ~PrimaryGeneratorAction();
  //! defines primary particles (mandatory)
  void GeneratePrimaries(G4Event*);

	//G4GeneralParticleSource * GetGPS(){
	//return gun;}
  G4double GetGPSEnergy(){
	G4double energy = gun->GetParticleEnergy();
	return energy;}

private:  
  G4GeneralParticleSource* InitializeGPS();
private:
  //G4VPrimaryGenerator* gun;
  G4GeneralParticleSource * gun;
  //std::ofstream * outfile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
