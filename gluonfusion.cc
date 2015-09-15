
// originally from pythia8 SigmaHiggs.cc
// modified by Tim Lou (07/22/2014)

#include <iostream>

#include "gluonfusion.hh"

using namespace std;

namespace Pythia8 {

  const double vh = 246.0;
  const double mt = 175.0;
  const double mh = 125.0;  
  const double lt2 = pow2(mt/vh) * 2.0 ;


//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceTheta::initConstants() {

  // ** TO DO ** 
  totalWidth = 10;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceTheta::calcWidth(bool) {

  if (id1Abs == 24) widNow = totalWidth;
  else widNow = 0;

}


//--------------------------------------------------------------------------

// Initialize process.
  
void Sigma1gg2Theta::initProc() {

  
  // Store H0, H1, H2 or A3 mass and width for propagator.

  // width = roughly 10 GeV
  // ** TODO **
  // ** need detailed calculations from pythong **

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

  //cout<<"sampling: "<<mH<<endl;
  double widthIn  = particleDataPtr->particleDataEntryPtr(25)
    ->resWidthChan( mH, 21, 21) / 64.;

  // widthIn = Higgs partial decay width to gluons
  // answer = width * 1/(s-mh^2)^2 * (yukawa)^2

  double extra_ME = lt2 * vh * 1.0/(sH - pow2(mh));
  
  //factor of two and * 3 due to convention, 3 = F-QCD  factor, 
  //2 = using m_stop instead of m_eta in python file

  //** TODO ** typically 0.1
  double wave_function = 0.1;
  double extra_factor = pow2(extra_ME) * sH * wave_function * 3.0/2.0;
  //cout<<"extra factor: "<<extra_factor<<endl;
  widthIn *= extra_factor;

  //cout<<"width In: "<<widthIn<<endl;


  // Set up Breit-Wigner.

  double width    = GammaRes;
  double sigBW    = 8. * M_PI/ ( pow2(sH - m2Res) + pow2(mH * width) );

  //cout<<"first piece "<<pow2(sH - m2Res)<<endl;
  //cout<<"second piece "<<pow2(mH * width)<<endl;

  double widthOut = width;

  //cout<<widthIn<<","<<sigBW<<","<<widthOut<<endl;

  // Done.
  sigma = widthIn * sigBW * widthOut;

  //cout<<sigma<<endl;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gg2Theta::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idTheta);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1gg2Theta::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  //WW channel only
  return 1.;

}

} // end namespace Pythia8
