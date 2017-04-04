// $Id: PrimaryGeneratorAction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief implements mandatory user class PrimaryGeneratorAction
 */

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
  //: outfile(0)
{
  gun = InitializeGPS();
  
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  gun->GeneratePrimaryVertex(anEvent);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete gun;
}

G4GeneralParticleSource* PrimaryGeneratorAction::InitializeGPS()
{
  G4GeneralParticleSource * gps = new G4GeneralParticleSource();
  
  // setup details easier via UI commands see gps.mac
    //G4int neutNum = 2112;
  // particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  gps->GetCurrentSource()->SetParticleDefinition(neutron);

  // set energy distribution
  G4SPSEneDistribution *eneDist = gps->GetCurrentSource()->GetEneDist() ;
  eneDist->SetEnergyDisType("Mono"); // or gauss
  //eneDist->SetMonoEnergy(10*CLHEP::MeV);

  // set position distribution
  G4SPSPosDistribution *posDist = gps->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Beam");  // or Point,Plane,Volume,Beam
  posDist->SetCentreCoords(G4ThreeVector(0.0*CLHEP::cm,0.0*CLHEP::cm,0.0*CLHEP::cm));
  posDist->SetBeamSigmaInX(0.0*CLHEP::mm);
  posDist->SetBeamSigmaInY(0.0*CLHEP::mm);

  // set angular distribution
  G4SPSAngDistribution *angDist = gps->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection( G4ThreeVector(0., 0., 0.) );
  angDist->SetAngDistType("beam2d");
  angDist->SetBeamSigmaInAngX(0.0*CLHEP::mrad);
  angDist->SetBeamSigmaInAngY(0.0*CLHEP::mrad);
  angDist->DefineAngRefAxes("angref1",G4ThreeVector(-1.,0.,0.));

  return gps;
}

