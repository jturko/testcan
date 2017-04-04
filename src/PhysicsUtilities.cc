#include "PhysicsUtilities.hh"

using namespace std;
using namespace TMath;
using namespace Physics;

double Physics::bwe(int lam, int A, int outputflag) {
  /* output */
  if(outputflag==1) {
    printf("\nbwe:\n*****\n");
    printf("calculates the Weisskopf single particle value for given massnumber\n");
    printf("and multipolarity for electric multipole transitions in e^2  fm^2lam\n");
  }

  /* calculate Weisskopf single particle value */
  double bwe = -1;
  bwe = 1/(4*(TMath::Pi())) * pow((3/(double)(lam+3)),2) * pow(1.2,2*lam) * pow(A,2*lam/3.);

  /* output */
  if(outputflag==1) {
    printf("\nWeisskopf unit for A = %d is: B_w(E%d) = %.3f e^2 fm^%d\n\n", A, lam, bwe, 2*lam);
  }

  return bwe;
}

double Physics::SolidAngleConversionRecoil(double BetaCm, double RecoilBetaCm, double RecoilThetaCm) {
  //calculate conversion of solid angle from lab to cm (dOmega_lab/dOmega_cm)
  //from heiko's phd thesis:
  //(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
  //                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
  //d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
  //beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
  //but I want to calculate (dsigma/dOmega)_cm from (dsigma/dOmega)_lab of the recoil => switch ejectile with recoil
  //below everything w/o Recoil means the cm-frame
  RecoilThetaCm = Pi() - RecoilThetaCm;//different definitions of theta
  return SolidAngleConversionEjectile(BetaCm, RecoilBetaCm, RecoilThetaCm);
}

double Physics::SolidAngleConversionRecoil(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm) {
  recoilThetaCm = Pi() - recoilThetaCm;//different definitions of theta
  return SolidAngleConversionEjectile(beamEnergy, projectileMass, targetMass, ejectileMass, recoilMass, eEx, recoilThetaCm);
}

double Physics::SolidAngleConversionEjectile(double BetaCm, double ejectileBetaCm, double ejectileThetaCm) {
  //calculate conversion of solid angle from lab to cm (dOmega_lab/dOmega_cm)
  //from heiko's phd thesis:
  //(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
  //                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
  //d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
  //beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
  //below everything w/o ejectile means the cm-frame
  double GammaCmSquared = 1/(1 - Power(BetaCm,2));

  double Denominator = Power(GammaCmSquared*Power(Cos(ejectileThetaCm)+BetaCm/ejectileBetaCm,2)+Power(Sin(ejectileThetaCm),2),3./2.);
  //double Denominator = Power(Power(Cos(ejectileThetaCm)+BetaCm/ejectileBetaCm,2)+Power(Sin(ejectileThetaCm),2),3./2.);
  
  if(Denominator != 0)
    return Sqrt(GammaCmSquared) * (1 + Cos(ejectileThetaCm)*BetaCm/ejectileBetaCm) / Denominator;
  else
    return -1.;
}

double Physics::SolidAngleConversionEjectile(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double ejectileThetaCm) {
  double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
  double eCm = ECm(beamEnergy, projectileMass, targetMass);
  double ejectileBetaCm = BetaEjectileCm(eCm, eEx, ejectileMass, recoilMass);
  //calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
  //from heiko's phd thesis:
  //(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
  //                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
  //d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
  //beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
  //below everything w/o ejectile means the cm-frame
  double gammaCmSquared = 1/(1 - Power(betaCm,2));

  double denominator = Power(gammaCmSquared*Power(Cos(ejectileThetaCm)+betaCm/ejectileBetaCm,2)+Power(Sin(ejectileThetaCm),2),3./2.);
  
  if(denominator != 0)
    return Sqrt(gammaCmSquared) * (1 + Cos(ejectileThetaCm)*betaCm/ejectileBetaCm) / denominator;
  else
    return -1.;
}

double Physics::SolidAngleConversionRecoil(double* x, double* par) {
  //calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
  //parameters: 0 = beta of cm system, 1 = beta of recoil in cm system

  return SolidAngleConversionRecoil(par[0], par[1], x[0]);
}

double Physics::SolidAngleConversionEjectile(double* x, double* par) {
  //calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
  //parameters: 0 = beta of cm system, 1 = beta of recoil in cm system

  return SolidAngleConversionEjectile(par[0], par[1], x[0]);
}


//WARNING: theta_lab is 180 degree for zero degree cm scattering (=> recoil)
double Physics::ThetaLab(double BetaCm, double RecoilBetaCm, double RecoilThetaCm) {
  //calculate theta_lab
  double CosineThetaLab = 0.;

  //from heiko's phd thesis:
  //cos(theta_lab) = gamma*(cos(theta_cm)+beta/beta_cm)/sqrt(sin^2(theta_cm)+gamma^2*(cos(theta_cm)+beta/beta_cm)^2)
  //beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
  double GammaCmSquared = 1/(1 - Power(BetaCm,2));

  double CosinePlusBetaRatio = Cos(RecoilThetaCm)+BetaCm/RecoilBetaCm;

  double Denominator = Sqrt(Power(Sin(RecoilThetaCm),2)+GammaCmSquared*Power((CosinePlusBetaRatio),2));

  if(Denominator != 0)
    CosineThetaLab = Sqrt(GammaCmSquared)*(CosinePlusBetaRatio) / Denominator;
  else
    return -1.;

  return Pi() - ACos(CosineThetaLab);
}

double Physics::ThetaLab(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm) {
  double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
  double eCm = ECm(beamEnergy, projectileMass, targetMass);
  double recoilBetaCm = BetaRecoilCm(eCm, eEx, ejectileMass, recoilMass);
  //calculate theta_lab
  double cosineThetaLab = 0.;

  //from heiko's phd thesis:
  //cos(theta_lab) = gamma*(cos(theta_cm)+beta/beta_cm)/sqrt(sin^2(theta_cm)+gamma^2*(cos(theta_cm)+beta/beta_cm)^2)
  //beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
  double gammaCmSquared = 1/(1 - Power(betaCm,2));

  double cosinePlusBetaRatio = Cos(recoilThetaCm)+betaCm/recoilBetaCm;

  double denominator = Sqrt(Power(Sin(recoilThetaCm),2)+gammaCmSquared*Power((cosinePlusBetaRatio),2));

  if(denominator != 0)
    cosineThetaLab = Sqrt(gammaCmSquared)*(cosinePlusBetaRatio)/denominator;
  else
    return -1.;

  return Pi() - ACos(cosineThetaLab);
}

//what's the difference to ThetaRecoilCm??? this one gives unreasonable results???
double Physics::ThetaCm(double BetaCm, double RecoilBetaCm, double RecoilThetaLab) {
  //double Resolution = 1.0e-3;
  ////calculate theta_cm by recursion
  //double BetaRatio = BetaCm/RecoilBetaCm;
  //double GammaCmSquared = 1/(1 - Power(BetaCm,2));
  //
  //double ThetaLabMinimum = 0.;
  //
  ////check whether the given theta_lab can be reached with the given beta ratio
  //if(BetaRatio == 1)//maximum theta_lab is 90 degree
  //  {
  //    ThetaLabMinimum = Pi()/2.;
  //  }
  //else if(BetaRatio > 1)
  //  {
  //    ThetaLabMinimum = Pi() - ATan2(RecoilBetaCm,Sqrt(GammaCmSquared*(Power(BetaCm,2)-Power(RecoilBetaCm,2))));
  //
  //    cerr<<"Warning, theta_cm is not unambiguous, will only return one value!"<<endl;
  //  }
  //
  //if(RecoilThetaLab < ThetaLabMinimum)
  //  {
  //    cerr<<__PRETTY_FUNCTION__<<": this theta_lab can't be reached with the given beta's"<<endl;
  //
  //    return -1;
  //  }
  //
  //double ThetaCmLow = 0.;
  //double ThetaCmHigh = Pi();
  //double ThetaCm = Pi()/2.;
  //
  //while(Abs(ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab) > Resolution)
  //  {
  //    //check which two value pairs bracket the desired theta_lab, ThetaCm and ThetaCmLow or ThetaCm and ThetaCmHigh
  //    if((ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab)*(ThetaLab(BetaCm, RecoilBetaCm, ThetaCmLow) - RecoilThetaLab) < 0)
  //	{
  //	  ThetaCmHigh = ThetaCm;
  //	  ThetaCm = (ThetaCmLow + ThetaCmHigh)/2.;
  //	}
  //    else if((ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab)*(ThetaLab(BetaCm, RecoilBetaCm, ThetaCmHigh) - RecoilThetaLab) < 0)
  //	{
  //	  ThetaCmLow = ThetaCm;
  //	  ThetaCm = (ThetaCmLow + ThetaCmHigh)/2.;
  //	}
  //    else
  //	{
  //	  cerr<<"bracketing method failed, return average of high and low bracket"<<endl;
  //
  //	  return (ThetaLab(BetaCm, RecoilBetaCm, ThetaCmLow) + ThetaLab(BetaCm, RecoilBetaCm, ThetaCmHigh))/2.;
  //	}
  //  }
  //
  //return ThetaLab(BetaCm, RecoilBetaCm, ThetaCm);

  //x4 = betacm/betarecoil
  //thetacm = acos((-x4*gammacm^2*tan(thetalab_recoil)^2+-sqrt(1+gammacm^2*tan(thetalab_recoil)^2*(1-x4^2)))/(1+gammacm^2*tan(thetalab_recoil)^2));
  double BetaRatio = BetaCm/RecoilBetaCm;
  double GammaCmSquared = 1/(1 - Power(BetaCm,2));
  double GammaTangensSquared = GammaCmSquared*Power(RecoilThetaLab,2);

  //check whether the denominator will be zero
  if(GammaTangensSquared == -1)
    {
      cerr<<__PRETTY_FUNCTION__<<": GammaTangensSquared is -1"<<endl;
      return -1.;
    }

  double CosineThetaCm;

  //WARNING: this is for recoils i.e. theta_lab = 180 degree => theta_cm = 0 degree!
  if(RecoilThetaLab > Pi()/2.)
    {
      CosineThetaCm = (-BetaRatio*GammaTangensSquared+Sqrt(1+GammaTangensSquared*(1-Power(BetaRatio,2)))) / (1+GammaTangensSquared);
    }
  else
    {
      CosineThetaCm = (-BetaRatio*GammaTangensSquared-Sqrt(1+GammaTangensSquared*(1-Power(BetaRatio,2)))) / (1+GammaTangensSquared);
    }
  
  return ACos(CosineThetaCm);
}

double Physics::RelativisticGamma(double Beta) {
  if(Beta < 0 || Beta > 1)
    {
      cerr<<__PRETTY_FUNCTION__<<": Beta is not in range 0-1: "<<Beta<<endl;
      return -1;
    }

  return 1./Sqrt(1.-Power(Beta,2));
}

//total kinetic energy in the ingoing channel of the cm-system
double Physics::TiCm(double BeamEnergy, double ProjectileMass, double TargetMass) {
  return Sqrt(2.*(ProjectileMass+BeamEnergy)*TargetMass + Power(ProjectileMass, 2) + Power(TargetMass, 2)) - ProjectileMass - TargetMass;
}

//total kinetic energy in the outgoing channel of the cm-system
double Physics::TfCm(double ECm, double Eex, double EjectileMass, double RecoilMass) {
  return ECm - Eex - (EjectileMass+RecoilMass);
}

//kinetic energy of the projectile in the cm-system
double Physics::TProjectileCm(double BeamEnergy, double ProjectileMass, double TargetMass) {
  if(ECm(BeamEnergy, ProjectileMass, TargetMass) > 0.) {
    return TiCm(BeamEnergy, ProjectileMass, TargetMass)/2.*(TiCm(BeamEnergy, ProjectileMass, TargetMass) + 2.*TargetMass)/ECm(BeamEnergy, ProjectileMass, TargetMass);
  }
  return 0.;
}

//kinetic energy of the recoil in the cm-system
double Physics::TRecoilCm(double ECm, double Eex, double EjectileMass, double RecoilMass) {
  double tfCm = TfCm(ECm,Eex,EjectileMass,RecoilMass);

  if(ECm == 0.)
    {
      cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<endl;
      return -1.;
    }

  //cout<<"tfCm = "<<tfCm<<", EjectileMass = "<<EjectileMass<<", Eex = "<<Eex<<", ECm = "<<ECm<<" => TRecoilCm = "<<tfCm*(tfCm+2*(EjectileMass+Eex))/(2.*ECm)<<endl;

  return tfCm*(tfCm+2.*(EjectileMass+Eex))/(2.*ECm);
}

//kinetic energy of the ejectile in the cm-system
double Physics::TEjectileCm(double ECm, double Eex, double EjectileMass, double RecoilMass) {
  double tfCm = TfCm(ECm,Eex,EjectileMass,RecoilMass);

  if(ECm == 0.)
    {
      cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<endl;
      return -1.;
    }

  //cout<<"tfCm = "<<tfCm<<", EjectileMass = "<<EjectileMass<<", Eex = "<<Eex<<", ECm = "<<ECm<<" => TRecoilCm = "<<tfCm*(tfCm+2*(EjectileMass+Eex))/(2.*ECm)<<endl;

  return tfCm*(tfCm+2.*(RecoilMass+Eex))/(2.*ECm);
}

//beta of the recoil in the cm-system
double Physics::BetaRecoilCm(double ECm, double Eex, double EjectileMass, double RecoilMass) {
  double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);

  //check that neither denominator gets zero or sqrt gets imaginary
  //this is fulfilled if trecoilcm and RecoilMass are larger than zero (which they have to be anyway)
  if(tRecoilCm < 0 || RecoilMass < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0"<<endl;
      return -1;
    }

  return Sqrt(Power(tRecoilCm,2) + 2.*tRecoilCm*RecoilMass)/(tRecoilCm + RecoilMass);
}

//beta of the ejectile in the cm-system
double Physics::BetaEjectileCm(double ECm, double Eex, double EjectileMass, double RecoilMass) {
  double tEjectileCm = TEjectileCm(ECm,Eex,EjectileMass,RecoilMass);

  //check that neither denominator gets zero or sqrt gets imaginary
  //this is fulfilled if tEjectileCm and EjectileMass are larger than zero (which they have to be anyway)
  if(tEjectileCm < 0 || EjectileMass < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": either the ejectile mass ("<<EjectileMass<<") or the ejectile kinetic energy in the cm system ("<<tEjectileCm<<") are 0"<<endl;
      return -1;
    }

  return Sqrt(Power(tEjectileCm,2) + 2.*tEjectileCm*EjectileMass)/(tEjectileCm + EjectileMass);
}

//theta of recoil in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::ThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
  double GammaTangensSquared = Power(Tan(ThetaRecoilLab),2)/(1.-Power(Beta,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(GammaTangensSquared == -1 || (1.+GammaTangensSquared*(1.-Power(BetaRatio,2))) < 0)
    {
      //cerr<<__PRETTY_FUNCTION__<<": either the denominator gets zero or the sqrt will be imaginary: "<<GammaTangensSquared<<", "<<(1.+GammaTangensSquared*(1.-Power(Beta,2)))<<endl;
      return -1.;
    }

  //this seems to be the right combination, but why?
  if(ThetaRecoilLab < Pi()/2.)
    return  ACos((BetaRatio*GammaTangensSquared - Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared));
  
  return ACos((BetaRatio*GammaTangensSquared + Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared));
}

//cosine of theta of recoil in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
  double GammaTangensSquared = Power(Tan(ThetaRecoilLab),2)/(1.-Power(Beta,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(GammaTangensSquared == -1)
    {
      cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<endl;
      return -2.;
    }

  if((1.+GammaTangensSquared*(1.-Power(BetaRatio,2))) < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-Power(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-Power("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<endl;
      return -2.;
    }

  //this seems to be the right combination, but why?
  if(ThetaRecoilLab < Pi()/2.)
    {
      return (BetaRatio*GammaTangensSquared - Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
    }
  
  return (BetaRatio*GammaTangensSquared + Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
}

//cosine of theta of ejectile in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaEjectileCm(double Beta, double ECm, double Eex, double ThetaEjectileLab, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaEjectileCm(ECm,Eex,EjectileMass,RecoilMass);
  double GammaTangensSquared = Power(Tan(ThetaEjectileLab),2)/(1.-Power(Beta,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(GammaTangensSquared == -1)
    {
      cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<endl;
      return -2.;
    }

  if((1.+GammaTangensSquared*(1.-Power(BetaRatio,2))) < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-Power(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-Power("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<endl;
      return -2.;
    }

  //this seems to be the right combination, but why?
  if(ThetaEjectileLab < Pi()/2.)
    {
      return (BetaRatio*GammaTangensSquared - Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
    }
  
  return (BetaRatio*GammaTangensSquared + Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
}

//cosine of theta of ejectile in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaEjectileCm(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double thetaEjectileLab) {
  double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
  double eCm = ECm(beamEnergy, projectileMass, targetMass);
  double betaRatio = betaCm/BetaEjectileCm(eCm, eEx, ejectileMass, recoilMass);
  double gammaTangensSquared = Power(Tan(thetaEjectileLab),2)/(1.-Power(betaCm,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(gammaTangensSquared == -1) {
    cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<gammaTangensSquared<<endl;
    return -2.;
  }

  if((1.+gammaTangensSquared*(1.-Power(betaRatio,2))) < 0) {
    //cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+gammaTangensSquared*(1.-Power(betaRatio,2)))<<" (1+"<<gammaTangensSquared<<"*(1.-Power("<<betaRatio<<",2)), Beta = "<<betaCm<<")"<<endl;
    return -2.;
  }

  //this seems to be the right combination, but why?
  if(thetaEjectileLab < Pi()/2.) {
    return (betaRatio*gammaTangensSquared - Sqrt(1.+gammaTangensSquared*(1.-Power(betaRatio,2))))/(1+gammaTangensSquared);
  }
  
  return (betaRatio*gammaTangensSquared + Sqrt(1.+gammaTangensSquared*(1.-Power(betaRatio,2))))/(1+gammaTangensSquared);
}

//beta of recoil in cm-system for given lab and cm-angle and beta of cm-system
double Physics::BetaRecoilCm(double Beta, double ThetaRecoilLab, double ThetaRecoilCm) {
  //second solution: plus sign in numerator (does not give the right results)
  double Denominator = 1-Power(Beta,2)-Power(Cos(ThetaRecoilCm),2)*(1-Power(Beta,2) + Power(Tan(ThetaRecoilLab),2));

  //check whether the denominator will be zero
  if(Denominator == -1 || Denominator == 0 || Beta < 0 || Beta > 1)
    {
      cerr<<__PRETTY_FUNCTION__<<": Beta is not in range 0-1: "<<Beta<<" or denominator is wrong: "<<Denominator<<endl;
      return -1.;
    }

  return (Beta*cos(ThetaRecoilCm)*Power(Tan(ThetaRecoilLab),2) - Sin(ThetaRecoilCm)*Tan(ThetaRecoilLab)*Beta*Sqrt(1-Power(Beta,2))) / Denominator;
}

//excitation energy needed to reach certain lab and cm-angle combinations
double Physics::Eex(double ECm, double Beta, double EjectileMass, double RecoilMass, double ThetaRecoilLab, double ThetaRecoilCm) {
  //alternative solutions: minus and plus signs before the two Sqrt terms (i.e. four combinations)
  //and of course the second solutions from BetaRecoilCm => total eight different results
  //this one looks best in gnuplot (somewhat matches orbits output) and this and the other BetaRecoilCm solution give the only reasonable Eex-values (i.e. neither around 1e6 nor around -1e7)

  if((1-Power(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2)) < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": beta of recoil in cm is larger than 1: "<<BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm)<<endl;
      return -1.;
    }

  if((Power(ECm,2) + Power(RecoilMass,2) - 2*ECm*RecoilMass/Sqrt(1-Power(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2))) < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": something is wrong"<<endl;
      return -1.;
    }

  return (- EjectileMass + Sqrt(Power(ECm,2) + Power(RecoilMass,2) - 2*ECm*RecoilMass/Sqrt(1-Power(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2))));
}

//theta of recoil in lab-system, Beta = beta of cm-system, Eex = excitation energy
double Physics::ThetaRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);

  return ATan((Sqrt(1.-Power(Beta,2))*Sin(ThetaRecoilCm))/(Cos(ThetaRecoilCm)+BetaRatio));
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double gamma = RelativisticGamma(Beta);
  double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
  double cosThetaRecoilCm = CosineThetaRecoilCm(Beta, ECm, Eex, ThetaRecoilLab, EjectileMass, RecoilMass);

  //this happens only if the computation of ThetaRecoilCm failed
  if(cosThetaRecoilCm == -2.)
    {
      return 0.;
    }

  if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.)
    {
      cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<endl;
      return -1;
    }

  return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*Sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab2(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass) {
  double gamma = RelativisticGamma(Beta);
  double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
  double cosThetaRecoilCm = Cos(ThetaRecoilCm);

  if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.)
    {
      cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<endl;
      return -1;
    }

  return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*Sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

//kinetic energy of recoil in the lab system - doesn't work???
double Physics::TEjectileLab(double Beta, double ECm, double Eex, double ThetaEjectileLab, double EjectileMass, double RecoilMass) {
  double gamma = RelativisticGamma(Beta);
  double tEjectileCm = TEjectileCm(ECm,Eex,EjectileMass,RecoilMass);
  double cosThetaEjectileCm = CosineThetaEjectileCm(Beta, ECm, Eex, ThetaEjectileLab, EjectileMass, RecoilMass);

  //this happens only if the computation of ThetaEjectileCm failed
  if(cosThetaEjectileCm == -2.) {
    return 0.;
  }

  if(tEjectileCm < 0. || RecoilMass < 0. || cosThetaEjectileCm < -1. || cosThetaEjectileCm > 1.) {
    cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tEjectileCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaEjectileCm<<endl;
    return -1;
  }

  return (gamma-1)*RecoilMass + gamma*tEjectileCm - gamma*Beta*Sqrt(tEjectileCm*(tEjectileCm+2.*RecoilMass))*cosThetaEjectileCm;
}

//energy in the cm-system
double Physics::ECm(double BeamEnergy, double ProjectileMass, double TargetMass) {
  return Sqrt(2.*(ProjectileMass+BeamEnergy)*TargetMass + Power(ProjectileMass, 2) + Power(TargetMass, 2));
}

//beta of cm-system
double Physics::BetaCm(double BeamEnergy, double ProjectileMass, double TargetMass) {
  return Sqrt(2*ProjectileMass*BeamEnergy + Power(BeamEnergy,2))/(ProjectileMass + TargetMass + BeamEnergy);
}

double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double LabAngle, double LabEnergy) {
  //LabAngle has to be in rad, everything else (masses and energy) in the same unit as the result

  //from mathematica (fourth out of four solution, first two give imaginary results, third is always negative)
  //(-8*m1 - 4*m2 + t4Lab*Power(Sec(aLab),2) + Sec(aLab)*Sqrt(16*Power(m2,2)*Power(Cos(aLab),2) + t4Lab*(t4Lab*Power(Sec(aLab),2)*Power(Cos(aLab) - Sin(aLab),4) + 8*m2*(3 + Sin(2*aLab)))) - 2*t4Lab*Tan(aLab) + Sqrt((Power(Sec(aLab),4)*(24*Power(m2,2)*t4Lab*Power(Cos(aLab),2) + 32*Power(m1,2)*m2*Power(Cos(aLab),4) + 2*m2*Power(t4Lab,2)*Power(Cos(aLab) - Sin(aLab),4) + 16*Power(m2,2)*t4Lab*Power(Cos(aLab),3)*Sin(aLab) - 8*Power(m1,2)*t4Lab*Power(Cos(aLab),2)*(-1 + Sin(2*aLab)) + 2*Power(m1,2)*Cos(3*aLab)*Sqrt(2*(4*Power(m2,2) + 12*m2*t4Lab + Power(t4Lab,2)) + (8*Power(m2,2) - 2*Power(t4Lab,2))*Cos(2*aLab) + t4Lab*(t4Lab*Power(Sec(aLab),2) + 8*m2*Sin(2*aLab) - 4*t4Lab*Tan(aLab))) + 2*Cos(aLab)*(3*Power(m1,2) + m2*t4Lab - m2*t4Lab*Sin(2*aLab))*Sqrt(2*(4*Power(m2,2) + 12*m2*t4Lab + Power(t4Lab,2)) + (8*Power(m2,2) - 2*Power(t4Lab,2))*Cos(2*aLab) + t4Lab*(t4Lab*Power(Sec(aLab),2) + 8*m2*Sin(2*aLab) - 4*t4Lab*Tan(aLab)))))/m2))/8.

  double ts = Power(TargetMass,2);
  double ps = Power(ProjectileMass,2);
  double cs = Power(Cos(LabAngle),2);
  double es = Power(LabEnergy,2);
  double te = TargetMass*LabEnergy;

  return (-8.*ProjectileMass - 4.*TargetMass + LabEnergy/cs - 2.*LabEnergy*Tan(LabAngle) 
	  + Sqrt(16.*ts*cs + LabEnergy*(LabEnergy/cs*Power(Cos(LabAngle) - Sin(LabAngle),4) + 8.*TargetMass*(3. + Sin(2.*LabAngle))))/Cos(LabAngle) 
	  + Sqrt((24.*ts*LabEnergy*cs + 
		  32.*ps*TargetMass*Power(Cos(LabAngle),4) + 
		  2.*TargetMass*es*Power(Cos(LabAngle) - Sin(LabAngle),4) + 
		  16.*ts*LabEnergy*Power(Cos(LabAngle),3)*Sin(LabAngle) - 
		  8.*ps*LabEnergy*cs*(Sin(2.*LabAngle) - 1) + 
		  2.*ps*Cos(3.*LabAngle)*Sqrt(2.*(4.*ts + 12.*te + es) + 
					    (8.*ts - 2.*es)*Cos(2.*LabAngle) + 
					    LabEnergy*(LabEnergy/cs + 8.*TargetMass*Sin(2.*LabAngle) - 4.*LabEnergy*Tan(LabAngle))) + 
		  2.*Cos(LabAngle)*(3.*ps + te - te*Sin(2.*LabAngle))*Sqrt(2.*(4.*ts + 12.*te + es) + (8.*ts - 2.*es)*Cos(2.*LabAngle) + LabEnergy*(LabEnergy/cs + 8.*TargetMass*Sin(2.*LabAngle) - 4.*LabEnergy*Tan(LabAngle))))/(Power(Cos(LabAngle),4)*TargetMass)))/8.;
}

//-------------------- no excitation
//total kinetic energy in the outgoing channel of the cm-system
double Physics::TfCm(double ECm, double EjectileMass, double RecoilMass) {
  return ECm - (EjectileMass+RecoilMass);
}

//theta of recoil in cm-system, Beta = beta of cm-system, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::ThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaRecoilCm(ECm,0.,EjectileMass,RecoilMass);
  double GammaTangensSquared = Power(Tan(ThetaRecoilLab),2)/(1.-Power(Beta,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(GammaTangensSquared == -1 || (1.+GammaTangensSquared*(1.-Power(BetaRatio,2))) < 0)
    {
      //cerr<<__PRETTY_FUNCTION__<<": either the denominator gets zero or the sqrt will be imaginary: "<<GammaTangensSquared<<", "<<(1.+GammaTangensSquared*(1.-Power(Beta,2)))<<endl;
      return -1.;
    }

  //this seems to be the right combination, but why?
  if(ThetaRecoilLab < Pi()/2.)
    return  ACos((BetaRatio*GammaTangensSquared - Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared));
  
  return ACos((BetaRatio*GammaTangensSquared + Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared));
}

//cosine of theta of recoil in cm-system, Beta = beta of cm-system, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double BetaRatio = Beta/BetaRecoilCm(ECm,0.,EjectileMass,RecoilMass);
  double GammaTangensSquared = Power(Tan(ThetaRecoilLab),2)/(1.-Power(Beta,2));

  //check whether the denominator will be zero or sqrt gets imaginary
  if(GammaTangensSquared == -1)
    {
      cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<endl;
      return -2.;
    }

  if((1.+GammaTangensSquared*(1.-Power(BetaRatio,2))) < 0)
    {
      cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-Power(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-Power("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<endl;
      return -2.;
    }

  //this seems to be the right combination, but why?
  if(ThetaRecoilLab < Pi()/2.)
    {
      return (BetaRatio*GammaTangensSquared - Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
    }
  
  return (BetaRatio*GammaTangensSquared + Sqrt(1.+GammaTangensSquared*(1.-Power(BetaRatio,2))))/(1+GammaTangensSquared);
}

//kinetic energy of the recoil in the cm-system
double Physics::TRecoilCm(double ECm, double EjectileMass, double RecoilMass) {
  double tfCm = TfCm(ECm,EjectileMass,RecoilMass);

  if(ECm == 0.)
    {
      cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<endl;
      return -1.;
    }

  return tfCm*(tfCm+2*(EjectileMass))/(2.*ECm);
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass) {
  double gamma = RelativisticGamma(Beta);
  double tRecoilCm = TRecoilCm(ECm,EjectileMass,RecoilMass);
  double cosThetaRecoilCm = CosineThetaRecoilCm(Beta, ECm, ThetaRecoilLab, EjectileMass, RecoilMass);

  //this happens only if the computation of ThetaRecoilCm failed
  if(cosThetaRecoilCm == -2.)
    {
      return 0.;
    }

  if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.)
    {
      cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are less than 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<endl;
      return -1;
    }

  return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*Sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double LabEnergy) {
  //everything (masses and energy) in the same unit as the result

  //from mathematica (fourth out of four solution, first two give imaginary results, third is always negative)
  //(-8*m1 - 4*m2 + t4Lab*Power(Sec(aLab),2) + Sec(aLab)*Sqrt(16*Power(m2,2)*Power(Cos(aLab),2) + t4Lab*(t4Lab*Power(Sec(aLab),2)*Power(Cos(aLab) - Sin(aLab),4) + 8*m2*(3 + Sin(2*aLab)))) - 2*t4Lab*Tan(aLab) + Sqrt((Power(Sec(aLab),4)*(24*Power(m2,2)*t4Lab*Power(Cos(aLab),2) + 32*Power(m1,2)*m2*Power(Cos(aLab),4) + 2*m2*Power(t4Lab,2)*Power(Cos(aLab) - Sin(aLab),4) + 16*Power(m2,2)*t4Lab*Power(Cos(aLab),3)*Sin(aLab) - 8*Power(m1,2)*t4Lab*Power(Cos(aLab),2)*(-1 + Sin(2*aLab)) + 2*Power(m1,2)*Cos(3*aLab)*Sqrt(2*(4*Power(m2,2) + 12*m2*t4Lab + Power(t4Lab,2)) + (8*Power(m2,2) - 2*Power(t4Lab,2))*Cos(2*aLab) + t4Lab*(t4Lab*Power(Sec(aLab),2) + 8*m2*Sin(2*aLab) - 4*t4Lab*Tan(aLab))) + 2*Cos(aLab)*(3*Power(m1,2) + m2*t4Lab - m2*t4Lab*Sin(2*aLab))*Sqrt(2*(4*Power(m2,2) + 12*m2*t4Lab + Power(t4Lab,2)) + (8*Power(m2,2) - 2*Power(t4Lab,2))*Cos(2*aLab) + t4Lab*(t4Lab*Power(Sec(aLab),2) + 8*m2*Sin(2*aLab) - 4*t4Lab*Tan(aLab)))))/m2))/8.

  double ts = Power(TargetMass,2);
  double ps = Power(ProjectileMass,2);
  double es = Power(LabEnergy,2);
  double te = TargetMass*LabEnergy;

  return (-8.*ProjectileMass - 4.*TargetMass + LabEnergy
	  + Sqrt(16.*ts + LabEnergy*(LabEnergy + 8.*TargetMass*3.))
	  + Sqrt((24.*ts*LabEnergy + 
		  32.*ps*TargetMass + 
		  2.*TargetMass*es - 
		  8.*ps*LabEnergy + 
		  2.*ps*Sqrt(2.*(4.*ts + 12.*te + es) + 
			     (8.*ts - 2.*es) + 
			     Power(LabEnergy,2)) + 
		  2.*(3.*ps + te)*Sqrt(2.*(4.*ts + 12.*te + es) + (8.*ts - 2.*es) + Power(LabEnergy,2)))/TargetMass))/8.;
}

//no excitation and theta_lab of recoil = 0 degree; two solutions: +- Sqrt, +Sqrt yields a much higher beam energty (too high)
double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double EjectileMass, double RecoilMass, double LabEnergy) {
  double nominator = (-Power(TargetMass,3) + TargetMass*Power(EjectileMass,2) + 3.*Power(TargetMass,2)*RecoilMass - Power(EjectileMass,2)*RecoilMass - 3.*TargetMass*Power(RecoilMass,2) + Power(RecoilMass,3) + 
		      3.*Power(TargetMass,2)*LabEnergy - Power(EjectileMass,2)*LabEnergy - 4.*TargetMass*RecoilMass*LabEnergy + Power(RecoilMass,2)*LabEnergy - 2.*TargetMass*Power(LabEnergy,2) + 
		      Power(ProjectileMass,2)*(-TargetMass + RecoilMass + LabEnergy) - 2.*ProjectileMass*(Power(TargetMass,2) + Power(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)) - 
		      Sqrt(LabEnergy*(2.*RecoilMass + LabEnergy)*(Power(ProjectileMass,4) + Power((Power(TargetMass,2) - Power(EjectileMass,2) + Power(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)),2) - 
								  2.*Power(ProjectileMass,2)*(Power(TargetMass,2) + Power(EjectileMass,2) + Power(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)))));
  double denominator = 2.*(Power(TargetMass,2) + Power(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy));

  if(denominator == 0.) {
    return 0.;
  }

  return nominator/denominator;
}

//-------------------- geometric functions (solid angle and such) --------------------
double Physics::SinusTheta(double X, double Y, double Z, double Factor) {
  return Factor/Sqrt(1 + Power(Z,2)/(Power(X,2)+Power(Y,2)));
}

double Physics::SinusThetaZInt(double X, double Y, double Z1, double Z2, double Factor) {
  return - Factor * Sqrt(Power(X,2) + Power(Y,2)) * Log((Z1 + Sqrt(Power(X,2) + Power(Y,2) + Power(Z1,2)))/(Z2 + Sqrt(Power(X,2) + Power(Y,2) + Power(Z2,2))));
}

double Physics::SolidAngle(double x, double y, double z) {
  //solid angle at point x,y,z
  return y/Power(Power(x,2)+Power(y,2)+Power(z,2),3./2.);
}

double Physics::IntegratedSolidAngle(double x, double y, double z, double d, double w) {
  //integration of solid angle over a plane at position y from x-d/2 to x+d/2 and z-w/2 to z+w/2
  //division by d*w is not necessary (this would normalize the result to be independent from d and w)
  double x2 = Power(x,2);
  double y2 = Power(y,2);
  double z2 = Power(z,2);
  double d2 = Power(d,2);
  double w2 = Power(w,2);

  double xl = Power(x - d/2.,2);
  double xh = Power(x + d/2.,2);
  double zl = Power(z - w/2.,2);
  double zh = Power(z + w/2.,2);

  return (2*(ATan((4.*d*w*y)/(-(d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xl + zl) + (-w2 + Power(d - 2.*x,2) + 4.*(y2 + z2))*Sqrt(y2 + xh + zl)
			      + (-d2 + 4.*(x2 + y2) + Power(w - 2.*z,2))*Sqrt(y2 + xl + zh) + 4.*Sqrt((y2 + xl + zl)*(y2 + xh + zl)*(y2 + xl + zh))))
	     + ATan((4.*d*w*y)/((-d2 + 4.*(x2 + y2) + Power(w + 2.*z,2))*Sqrt(y2 + xh + zl) + (-w2 + Power(d + 2.*x,2) + 4.*(y2 + z2))*Sqrt(y2 + xl + zh)
				- (d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xh + zh) + 4.*Sqrt((y2 + xh + zl)*(y2 + xl + zh)*(y2 + xh + zh))))));
}

double Physics::SolidAngle(double* x, double* par) {
  //parameters: 0 = y, 1 = z-shift, 2 = multiplication factor
  //solid angle at point x[0],par[0],x[1]+par[1] multiplited with par[2]
  return par[2]*par[0]/Power(Power(x[0],2)+Power(par[0],2)+Power(x[1] - par[1],2),3./2.);
}

double Physics::IntegratedSolidAngle(double* x, double* par) {
  //x[0] = along strip, x[1] = along beam line
  //parameters: 0 = y, 1 = width in x, 2 = width in z, 3 = z-shift, 4 = multiplication factor
  //integration of solid angle over a plane at position par[0] from x-par[1]/2 to x+par[1]/2 and z-par[2]/2 to z+par[2]/2 
  //limit for par[1] and par[2] -> 0 is: par[0]/Power(Power(x[0],2)+Power(par[0],2)+Power(x[1],2),3./2.)
  //NOT NECESSARY:
  //1/(par[1]*par[2]) garantuees the normalization (integral over x and z from -inf to inf yields always 2 pi, independent of the pixel-size used)

  double z = x[1] - par[3];

  double x2 = Power(x[0],2);
  double y2 = Power(par[0],2);
  double z2 = Power(z,2);
  double d2 = Power(par[1],2);
  double w2 = Power(par[2],2);

  double xl = Power(x[0] - par[1]/2.,2);
  double xh = Power(x[0] + par[1]/2.,2);
  double zl = Power(z - par[2]/2.,2);
  double zh = Power(z + par[2]/2.,2);

  return par[4]*(2*(ATan((4.*par[1]*par[2]*par[0])/(-(d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xl + zl) + (-w2 + Power(par[1] - 2.*x[0],2) + 4.*(y2 + z2))*Sqrt(y2 + xh + zl)
								    + (-d2 + 4.*(x2 + y2) + Power(par[2] - 2.*z,2))*Sqrt(y2 + xl + zh) + 4.*Sqrt((y2 + xl + zl)*(y2 + xh + zl)*(y2 + xl + zh))))
				    + ATan((4.*par[1]*par[2]*par[0])/((-d2 + 4.*(x2 + y2) + Power(par[2] + 2.*z,2))*Sqrt(y2 + xh + zl) + (-w2 + Power(par[1] + 2.*x[0],2) + 4.*(y2 + z2))*Sqrt(y2 + xl + zh)
								      - (d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xh + zh) + 4.*Sqrt((y2 + xh + zl)*(y2 + xl + zh)*(y2 + xh + zh))))));
}

double Physics::BarrelSolidAngle(double* n, double* par) {
  //n[0] = pixel number along strip, n[1] = strip number
  //parameters: 0 = y, 1 = width in x, 2 = width in z, 3 = z-shift, 4 = barrel width
  //integration of solid angle over a plane at position par[0] from x-par[1]/2 to x+par[1]/2 and z-par[2]/2 to z+par[2]/2
  //in this case z is the middle of the strip and not the randomized z-position
  //NOT NECESSARY:
  //1/(par[1]*par[2]) garantuees the normalization (integral over x and z from -inf to inf yields always 2 pi, independent of the pixel-size used)

  double x = n[0]*par[1] - par[4]/2.;
  double z = n[1]*par[2] - par[3];

  double x2 = Power(x,2);
  double y2 = Power(par[0],2);
  double z2 = Power(z,2);
  double d2 = Power(par[1],2);
  double w2 = Power(par[2],2);

  double xl = Power(x - par[1]/2.,2);
  double xh = Power(x + par[1]/2.,2);
  double zl = Power(z - par[2]/2.,2);
  double zh = Power(z + par[2]/2.,2);

  return 1.*(2*(ATan((4.*par[1]*par[2]*par[0])/(-(d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xl + zl) + (-w2 + Power(par[1] - 2.*x,2) + 4.*(y2 + z2))*Sqrt(y2 + xh + zl)
								+ (-d2 + 4.*(x2 + y2) + Power(par[2] - 2.*z,2))*Sqrt(y2 + xl + zh) + 4.*Sqrt((y2 + xl + zl)*(y2 + xh + zl)*(y2 + xl + zh))))
				+ ATan((4.*par[1]*par[2]*par[0])/((-d2 + 4.*(x2 + y2) + Power(par[2] + 2.*z,2))*Sqrt(y2 + xh + zl) + (-w2 + Power(par[1] + 2.*x,2) + 4.*(y2 + z2))*Sqrt(y2 + xl + zh)
								  - (d2 + w2 - 4.*(x2 + y2 + z2))*Sqrt(y2 + xh + zh) + 4.*Sqrt((y2 + xh + zl)*(y2 + xl + zh)*(y2 + xh + zh))))));
}

double Physics::IdentifiedSolidAngle(double* x, double* par) {
  //x[0] = along strip, x[1] = along beam line
  //parameters: 0 = y, 1 = width in x, 2 = width in z, 3 = z-shift, 4 = multiplication factor, 5 = gap between delta-E and E-rest detectors, 6 = striplength, 7 = detectorlength + target gap
  //integration of solid angle over a plane at position par[0] from x-par[1]/2 to x+par[1]/2 and z-par[2]/2 to z+par[2]/2 
  //limit for par[1] and par[2] -> 0 is: par[0]/Power(Power(x[0],2)+Power(par[0],2)+Power(x[1],2),3./2.)

  //three edges of the detector are not covered by E-rest: the two sides and the edge away from the target
  //cover towards the sides is up to striplength/2 * (distance to beam) / (distance to beam + gap)
  //cover along beam line is up to (detectorlength+targetgap-Zshift) * (distance to beam) / (distance to beam + gap) - (targetgap-Zshift) (measure in detector coordinates)
  //or maximum-z * (distance to beam) / (distance to beam + gap) (measured in coordinates centered at target)
  //cout<<"--------------------"<<endl
  //    <<"x = "<<x[0]<<", "<<x[1]<<", par = "<<par[0]<<", "<<par[1]<<", "<<par[2]<<", "<<par[3]<<", "<<par[4]<<", "<<par[5]<<", "<<par[6]<<", "<<par[7]<<endl
  //    <<"--------------------"<<endl;

  //copy the parameters
  double p[8];
  for(int i = 0; i < 8; i++)
    {
      p[i] = par[i];
    }

  //correct for target shift
  //cout<<"x[1] = "<<x[1]<<", p[7] = "<<p[7]<<", p[3] = "<<p[3]<<" => ";
  x[1] -= p[3];
  p[7] -= p[3];
  //correction is done, now set the shift to zero, otherwise all other functions perform the correction again!
  p[3] = 0.;
  //cout<<"x[1] = "<<x[1]<<", p[7] = "<<p[7]<<", p[3] = "<<p[3]<<endl;

  //check whether the outer edge of the pixel is outside the covered area
  if(Abs(x[0])+p[1]/2. > p[6]/2.*p[0]/(p[0]+p[5]))
    {
      //cout<<"Abs(x[0])+p[1]/2. = Abs("<<x[0]<<")+"<<p[1]<<"/2. = "<<Abs(x[0])+p[1]/2.<<" > p[6]/2.*p[0]/(p[0]+p[5]) = "<<p[6]<<"/2.*"<<p[0]<<"/"<<(p[0]+p[5])<<" = "<<p[6]/2.*p[0]/(p[0]+p[5])<<endl;
      //if the inner edge is outside the covered area as well nothing is covered => solid angle is zero
      if(Abs(x[0])-p[1]/2. >= p[6]/2.*p[0]/(p[0]+p[5]))
	{
	  //cout<<"Abs(x[0])-p[1]/2. = Abs("<<x[0]<<")-"<<p[1]<<"/2. = "<<Abs(x[0])-p[1]/2.<<" >= p[6]/2.*p[0]/(p[0]+p[5]) = "<<p[6]<<"/2.*"<<p[0]<<"/"<<(p[0]+p[5])<<" = "<<p[6]/2.*p[0]/(p[0]+p[5])<<endl;
	  return 0.;
	}
      //calculate solid angle between inner edge and edge of covered area
      //to this end we simply set the width of the pixel to the distance between those edges and change the center of the pixel accordingly
      //the center should move to the middle between the edge of the covered area and the inner edge of the pixel
      //the edge of the covered area is of the same sign as x[0] and the inner edge is x[0]+w/2 if x[0]<0 and x[0]-w/2 otherwise
      //Sign(a,b) returns Abs(a) for positive b and -Abs(a) otherwise
      //to set the width we need the inner edge and thus the original center, to set the new center we need the original width (from par) => change the width first
      //cout<<"x[0] = "<<x[0]<<", p[1] = "<<p[1]<<" => ";
      p[1] = p[6]/2.*p[0]/(p[0]+p[5]) - (Abs(x[0])-p[1]/2.);
      x[0] = (Sign(p[6]/2.*p[0]/(p[0]+p[5])-par[1]/2., x[0]) + x[0])/2.;
      //cout<<"x[0] = "<<x[0]<<", p[1] = "<<p[1]<<endl;
    }

  //check whether the far edge of the pixel is outside the covered area
  if(Abs(x[1]-p[3])+p[2]/2. > (p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5]))
    {
      //cout<<"Abs(x[1]-p[3])+p[2]/2. = Abs("<<x[1]-p[3]<<")+"<<p[2]<<"/2. = "<<Abs(x[1]-p[3])+p[2]/2.<<" > (p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5]) = "<<(p[7]-Sign(p[3],x[1]))<<"*"<<p[0]<<"/"<<(p[0]+p[5])<<" = "<<(p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5])<<endl;
      //if the closer edge is outside the covered area as well nothing is covered => solid angle is zero
      if(Abs(x[1])-p[2]/2. >= (p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5]))
	{
	  //cout<<"Abs(x[1])-p[2]/2. = Abs("<<x[1]<<")-"<<p[2]<<"/2. = "<<Abs(x[1])-p[2]/2.<<" >= (p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5]) = "<<(p[7]-Sign(p[3],x[1]))<<"*"<<p[0]<<"/"<<(p[0]+p[5])<<" = "<<(p[7]-Sign(p[3],x[1]))*p[0]/(p[0]+p[5])<<endl;
	  return 0.;
	}
      //calculate solid angle between near edge and edge of covered area
      //same as before
      //cout<<"x[1] = "<<x[1]<<", p[2] = "<<p[2]<<" => ";
      p[2] = p[7]*p[0]/(p[0]+p[5]) - (Abs(x[1])-p[2]/2.);
      x[1] = (Sign(p[7]*p[0]/(p[0]+p[5])-par[2]/2., x[1]) + x[1])/2.;
      //cout<<"x[1] = "<<x[1]<<", p[2] = "<<p[2]<<endl;
    }

  //now the parameters are set right (or zero was returned) => just return the solid angle
  return IntegratedSolidAngle(x,p);
}

double Physics::CoveredPhiAngle(double Theta, double DistanceToBeam, double DetectorWidth, double GapWidth, double Shift) {
  //return value is the covered phi angle in percent (of 2 pi)
  //convert Theta from deg to rad
  Theta *= Pi()/180.;

  //factor 4 stems from 4 quadrants each with 2 halfs and divided by 2 (pi)
  if(Theta < Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/DistanceToBeam))
    {
      return 0.;
    }
  else if(Theta < Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
    {
      return 4*ATan2(DistanceToBeam, Sqrt(Power(GapWidth-Shift+DetectorWidth,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta))/Pi();
    }
  else if(Theta < Pi()/2. - ATan((GapWidth-Shift)/DistanceToBeam))
    {
      return 4*ATan2(DistanceToBeam, DetectorWidth/2.)/Pi();
    }
  else if(Theta < Pi()/2. - ATan((GapWidth-Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
    {
      return 4*(ATan2(DistanceToBeam, DetectorWidth/2.) - ATan2(DistanceToBeam, Sqrt(Power(GapWidth-Shift,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta)))/Pi();
    }
  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
    {
      return 0.;
    }
  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift)/DistanceToBeam))
    {
      return 4*(ATan2(DistanceToBeam, DetectorWidth/2.) - ATan2(DistanceToBeam, - Sqrt(Power(GapWidth+Shift,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta)))/Pi();
    }
  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
    {
      return 4*ATan2(DistanceToBeam, DetectorWidth/2.)/Pi();
    }
  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift+DetectorWidth)/DistanceToBeam))
    {
      return 4*ATan2(DistanceToBeam, - Sqrt(Power(GapWidth+Shift+DetectorWidth,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta))/Pi();
    }
  else
    {
      return 0.;
    }
}

double Physics::CoveredPhiAngle(double* x, double* par) {
  //return value is the covered phi angle in percent (of 2 pi)
  //parameter: 0 = DistanceToBeam, 1 = DetectorWidth, 2 = GapWidth, 3 = Shift
  //convert x[0] from deg to rad
  x[0] *= Pi()/180.;

  //factor 4 stems from 4 quadrants each with 2 halfs and divided by 2 (pi)
  if(x[0] < Pi()/2. - ATan((par[2]-par[3]+par[1])/par[0]))
    {
      return 0.;
    }
  else if(x[0] < Pi()/2. - ATan((par[2]-par[3]+par[1])/Sqrt(Power(par[0],2)+Power(par[1]/2.,2))))
    {
      return 4*ATan2(Sqrt(Power(par[2]-par[3]+par[1],2) - Power(par[0]/Tan(x[0]),2)) * Tan(x[0]), par[0])/Pi();
    }
  else if(x[0] < Pi()/2. - ATan((par[2]-par[3])/par[0]))
    {
      return 4*ATan2(par[1]/2., par[0])/Pi();
    }
  else if(x[0] < Pi()/2. - ATan((par[2]-par[3])/Sqrt(Power(par[0],2)+Power(par[1]/2.,2))))
    {
      return 4*(ATan2(par[1]/2., par[0]) - ATan2(Sqrt(Power(par[2]-par[3],2) - Power(par[0]/Tan(x[0]),2)) * Tan(x[0]), par[0]))/Pi();
    }
  else if(x[0] < Pi()/2. - ATan(-(par[2]+par[3])/Sqrt(Power(par[0],2)+Power(par[1]/2.,2))))
    {
      return 0.;
    }
  else if(x[0] < Pi()/2. - ATan(-(par[2]+par[3])/par[0]))
    {
      return 4*(ATan2(par[1]/2., par[0]) - ATan2(- Sqrt(Power(par[2]+par[3],2) - Power(par[0]/Tan(x[0]),2)) * Tan(x[0]), par[0]))/Pi();
    }
  else if(x[0] < Pi()/2. - ATan(-(par[2]+par[3]+par[1])/Sqrt(Power(par[0],2)+Power(par[1]/2.,2))))
    {
      return 4*ATan2(par[1]/2., par[0])/Pi();
    }
  else if(x[0] < Pi()/2. - ATan(-(par[2]+par[3]+par[1])/par[0]))
    {
      return 4*ATan2(- Sqrt(Power(par[2]+par[3]+par[1],2) - Power(par[0]/Tan(x[0]),2)) * Tan(x[0]), par[0])/Pi();
    }
  else
    {
      return 0.;
    }
}

//double Physics::CoveredPhiAngle(double ThetaLow, double ThetaHigh, double DistanceToBeam, double DetectorWidth, double GapWidth, double Shift)
//{
//  //CAUTION: the following code assumes that the difference between ThetaLow and ThetaHigh is smaller than the different ranges of evalutation
//  //this means that if ThetaLow is in the nth range ThetaHigh is either in the same range or in the n+1th
//  if(ThetaLow >= ThetaHigh)
//    {
//      cerr<<__PRETTY_FUNCTION__<<": Error, ThetaLow is greater or equal to ThetaHigh ("<<ThetaLow<<" >= "<<ThetaHigh<<")"<<endl;
//      return -1.;
//    }
//
//  //convert Theta from deg to rad
//  ThetaLow *= Pi()/180.;
//  ThetaHigh *= Pi()/180.;
//
//  //need to distinguish several possible cases
//  //range is below first detector, range is partly below first detector, range is on the low edge of the first detector,
//  //range is partly on the low edge and partly inside the first detector, range is insider the first detector
//  //range is partly on the high edge and partly inside the first detector, range is on the high edge of the first det.,
//  //range is partly on the high edge of the first det., range is between detectors and so on for the second detector
//
//  //easiest solutions are when the intervall is completly outside the detector or completly inside
//  //both below smallest theta or both above highest theta
//  if(ThetaHigh < Pi()/2. - ATan( (GapWidth-Shift+DetectorWidth)/DistanceToBeam) ||
//     ThetaLow  > Pi()/2. - ATan(-(GapWidth+Shift+DetectorWidth)/DistanceToBeam))
//    {
//      return 0.;
//    }
//  //both between the detectors
//  if(ThetaHigh < Pi()/2. - ATan(-(GapWidth-Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))) &&
//     ThetaLow  > Pi()/2. - ATan( (GapWidth-Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      return 0.;
//    }
//  //inside first or second detector
//  if((ThetaLow  > Pi()/2. - ATan( (GapWidth-Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))) &&
//      ThetaHigh < Pi()/2. - ATan( (GapWidth-Shift)/DistanceToBeam)) || 
//     (ThetaLow  > Pi()/2. - ATan(-(GapWidth-Shift)/DistanceToBeam) &&
//      ThetaHigh < Pi()/2. - ATan(-(GapWidth-Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2)))))
//    {
//      return (ThetaHigh-ThetaLow)*8*ATan(DistanceToBeam, DetectorWidth/2.);
//    }
//  //range is (at least partially) on the low edge
//  if(ThetaHigh < Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      //if the range is only partly on the edge of a detector and otherwise outside simply shrink it so that its fully on the edge
//      if(ThetaLow > Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/DistanceToBeam))
//	{
//	  ThetaLow = Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/DistanceToBeam);
//	}
//      return ;
//    }
//
//
//
//  //factor 8 stems from 4 quadrants each with 2 halfs
//  if(Theta < Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/DistanceToBeam))
//    {
//      return 0.;
//    }
//  else if(Theta < Pi()/2. - ATan((GapWidth-Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      return 8*ATan(DistanceToBeam, Sqrt(Power(GapWidth-Shift+DetectorWidth,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta));
//    }
//  else if(Theta < Pi()/2. - ATan((GapWidth-Shift)/DistanceToBeam))
//    {
//      return 8*ATan(DistanceToBeam, DetectorWidth/2.);
//    }
//  else if(Theta < Pi()/2. - ATan((GapWidth-Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      return 8*(ATan(DistanceToBeam, DetectorWidth/2.) - ATan(DistanceToBeam, Sqrt(Power(GapWidth-Shift,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta)));
//    }
//  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      return 0.;
//    }
//  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift)/DistanceToBeam))
//    {
//      return 8*(ATan(DistanceToBeam, DetectorWidth/2.) - ATan(DistanceToBeam, - Sqrt(Power(GapWidth+Shift,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta)));
//    }
//  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift+DetectorWidth)/Sqrt(Power(DistanceToBeam,2)+Power(DetectorWidth/2.,2))))
//    {
//      return 8*ATan(DistanceToBeam, DetectorWidth/2.);
//    }
//  else if(Theta < Pi()/2. - ATan(-(GapWidth+Shift+DetectorWidth)/DistanceToBeam))
//    {
//      return 8*ATan(DistanceToBeam, - Sqrt(Power(GapWidth+Shift+DetectorWidth,2) - Power(DistanceToBeam/Tan(Theta),2)) * Tan(Theta));
//    }
//  else
//    {
//      return 0.;
//    }
//}

double Physics::GammaEfficiency(double* x, double* par) {
  if(x[0] > 0)
    return par[0]*Log(x[0]) + par[1]*Log(x[0])/x[0] + par[2]*Power(Log(x[0]),2)/x[0] + par[3]*Power(Log(x[0]),4)/x[0] + par[4]*Power(Log(x[0]),5)/x[0];

  return 1.;
}

double Physics::GammaEfficiencyError(double* x, double* par) {
  //parameter 0 is the number of parameters that were fitted
  //parameters 1 to par[0]^2 are the entries of the variance-covariance matrix
  double result = 0.;

  if(par[0] > 5 || par[0] < 0)
    {
      cerr<<"Error, par[0] in GammaEfficiencyError should be between 0 and 5!"<<endl;

      return 0.;
    }

  if(x[0] <= 0.)
    {
      return 0.;
    }

  double partialDerivative[5] = {Log(x[0]), Log(x[0])/x[0], Power(Log(x[0]),2)/x[0], Power(Log(x[0]),4)/x[0], Power(Log(x[0]),5)/x[0]};

  for(int i = 0; i < par[0]; i++)
    {
      for(int j = 0; j < par[0]; j++)
	{
	  result += partialDerivative[i]*partialDerivative[j]*par[i*((int)par[0])+j+1];
	}
    }

  if(result >= 0.)
    {
      return Sqrt(result);
    }

  //cout<<endl<<"x[0] = "<<x[0]<<endl;
  //
  //for(int i = 0; i < par[0]; i++)
  //  {
  //    for(int j = 0; j < par[0]; j++)
  //	{
  //	  cout<<setw(16)<<partialDerivative[i]*partialDerivative[j]*par[i*((int)par[0])+j+1]<<" ";
  //	}
  //    cout<<endl;
  //  }
  //
  //for(int i = 0; i < par[0]; i++)
  //  {
  //    for(int j = 0; j < par[0]; j++)
  //	{
  //	  cout<<partialDerivative[i]<<"*"<<partialDerivative[j]<<"*"<<par[i*((int)par[0])+j+1]<<" ";
  //	}
  //    cout<<endl;
  //  }
  //
  //cout<<"===================="<<endl;

  return 0.;
}

double Physics::GammaEfficiencyLowerBound(double* x, double* par) {
  //parameters 0-4 are passed to GammaEfficiency and the remaining parameters to GammaEfficiencyError
  return GammaEfficiency(x,par)-GammaEfficiencyError(x,par+5);
}

double Physics::GammaEfficiencyUpperBound(double* x, double* par) {
  //parameters 0-4 are passed to GammaEfficiency and the remaining parameters to GammaEfficiencyError
  return GammaEfficiency(x,par)+GammaEfficiencyError(x,par+5);
}

//ratio between two gamma efficiencies (e.g. cluster to core :)
double Physics::GammaEfficiencyRatio(double* x, double* par) {
  //0-4:denominator parameters, 5-9: divisor parameters
  if(x[0] > 0)
    {
      if((par[5]*x[0] + par[6] + par[7]*Log(x[0]) + par[8]*Power(Log(x[0]),3) + par[9]*Power(Log(x[0]),4)) != 0)
	{
	  return (par[0]*x[0] + par[1] + par[2]*Log(x[0]) + par[3]*Power(Log(x[0]),3) + par[4]*Power(Log(x[0]),4))/
	    (par[5]*x[0] + par[6] + par[7]*Log(x[0]) + par[8]*Power(Log(x[0]),3) + par[9]*Power(Log(x[0]),4));
	}
    }

  return 1.;
}

double Physics::ReleaseCurve(double* x, double* par) {
  //parameter 0: number of components
  //parameter 1+n*5: rise time of nth component
  //parameter 2+n*5: fast release time of nth component
  //parameter 3+n*5: slow release time of nth component
  //parameter 4+n*5: ratio fast/slow of nth component
  //parameter 5+n*5: amplitude of nth component

  return ReleaseCurve(x[0], par);
}

double Physics::ReleaseCurve(double x, double* par) {
  //parameter 0: number of components
  //parameter 1+n*5: rise time of nth component
  //parameter 2+n*5: fast release time of nth component
  //parameter 3+n*5: slow release time of nth component
  //parameter 4+n*5: ratio fast/slow of nth component
  //parameter 5+n*5: amplitude of nth component

  double Result = 0;
  for(int i = 0; i < (int) par[0]; i++)
    {
      if(((par[4+i*5]*Power(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + (1.-par[4+i*5])*Power(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/Log(2.)) != 0)
	{
	  Result += par[5+i*5]*(1. - Exp(-Log(2.)*x/Abs(par[1+i*5]))) * (par[4+i*5]*Exp(-Log(2.)*x/Abs(par[2+i*5])) + 
								     (1.-par[4+i*5])*Exp(-Log(2.)*x/Abs(par[3+i*5])))
	    /((par[4+i*5]*Power(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + 
	       (1.-par[4+i*5])*Power(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/Log(2.));
	}
      else
	{
	  cerr<<"Error, normalization is zero: "<<((par[4+i*5]*Power(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + (1.-par[4+i*5])*Power(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/Log(2.))<<" = (("<<par[4+i*5]*Power(par[2+i*5],2)/(par[2+i*5]+par[1+i*5])<<" + "<<(1.-par[4+i*5])*Power(par[3+i*5],2)<<"/"<<(par[3+i*5]+par[1+i*5])<<")/Log(2.))"<<endl;
	  return 0.;
	}
    }

  return Result;
}

double Physics::IsoldeReleaseCurve(double* x, double* par) {
  //parameters 0-2: probability to get each 1-3th proton pulse
  //parameter 3: number of components
  //plus 5 parameters for each component => 5*(# of components)+4 parameters
  double Result = ReleaseCurve(x[0],par+3) + 
    par[0]*ReleaseCurve(x[0]+1200.,par+3) + 
    par[1]*ReleaseCurve(x[0]+2400.,par+3) + 
    par[2]*ReleaseCurve(x[0]+3600.,par+3) +
    ReleaseCurve(x[0]+4800.,par+3);

  if(x[0] < 1200.)
    {
      return Result;
    }

  for(int i = 0; i < 3; i++)
    {
      if(x[0] < 2400.+i*1200.)
	{
	  return (1.-par[i])*Result;
	}
    }

  return 0.;
}
