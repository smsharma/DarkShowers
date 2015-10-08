// originally from pythia8 SigmaHiggs.cc
// modified by Tim Lou (07/22/2014)

#include <iostream>

#include "higgsportal.hh"

using namespace std;

namespace Pythia8 {

  const double vh = 246.0;
  const double mt = 175.0;
  const double mh = 125.0;  
  const double lt2 = pow2(mt/vh) * 2.0 ;

// Initialize process.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceTheta::initConstants() {

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceTheta::calcWidth(bool) {

  if (id1Abs == 24) widNow = totalWidth;
  else widNow = 0;

}

  
void Sigma1gg2Theta::initProc() {

  particlePtr = particleDataPtr->particleDataEntryPtr(idTheta);

  mRes = particleDataPtr->m0(idTheta);

  GammaRes = particleDataPtr->mWidth(idTheta);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1gg2Theta::sigmaKin() {

  // Incoming width for gluons, gives colour factor of 1/8 * 1/8.

  //not on-shell
  //if ( mH <= 2* mRes ) 
  
  
  if ( mH <= 2*mRes ) 
    {
      cout<<sH<<" and "<<pow2(2*mRes)<<endl;
      sigma = 0;
      return;
    }
  

  double widthIn  = particleDataPtr->particleDataEntryPtr(25)
    ->resWidthChan( mH, 21, 21) / 64.;

  // widthIn = Higgs partial decay width to gluons
  // answer = width * 1/(s-mh^2)^2 * (yukawa)^2

  double extra_ME2 = pow2(lt2 * vh) * 
    3.0/(pow2(sH - pow2(mh)));

  // widthIn = |ME|^2 / (16*pi*m)
  double ME2 = widthIn* 16* M_PI * mH * extra_ME2;


  sigma = ME2 / (16 * M_PI * sH2) ;

  /*
  cout<< "width: "<<widthIn <<endl;
  cout<<"extra ME: "<<extra_ME<<endl;
  cout<<"ME2 : "<<ME2<<endl;
  cout<<"sigma  at "<<mH<<" is "<<sigma*1e10<<endl;
  */
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gg2Theta::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idTheta, idTheta);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1gg2Theta::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  return 1.;

}

} // end namespace Pythia8
