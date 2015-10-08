// originally from pythia8 SigmaHiggs.cc
// modified by Tim Lou (07/22/2014)

#include <iostream>

#include "gluonportal.hh"

using namespace std;

namespace Pythia8 {

// Initialize process.
  
void Sigma2gg2ss::initProc() {

  particlePtr = particleDataPtr->particleDataEntryPtr(id_in);
  mRes = particleDataPtr->m0(id_in);
  GammaRes = particleDataPtr->mWidth(id_in);
  //if really small width
  if(GammaRes < 0.01)
    GammaRes = 0.01;
  
  m2Res    = mRes*mRes;
  cout<<"resonance mass: "<<mRes<<endl;
  GamMRat  = GammaRes / mRes;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2ss::sigmaKin() {

  // Incoming width for gluons, gives colour factor of 1/8 * 1/8.
  // spin and color factor
  double gluon_factor  = 1.0/(64.0*2.0*2.0);
  
  // matrix element = coupling / scale^2
  // scale = 100 GeV^2
  // coupling = 1.0
  // factor of 2 = polarization sum

  double ME2 = pow2(coupling/(4*M_PI) * sH/scale  )* 
    2* gluon_factor *
    mRes * GammaRes * 16 * M_PI /
    ( pow2(sH - m2Res) + pow2(sH * GamMRat) );


  
  /*
  cout<<"sampling: "<<sqrt(sH)<<endl;

  cout<<"first part: "<<pow2(coupling/(4*M_PI) * sH/scale  )<<endl;
  cout<<"second part: "<<    2* gluon_factor *
    mRes * GammaRes * 16 * M_PI<<endl;
  cout<<"third part: "<<     ( pow2(sH - m2Res) + pow2(sH * GamMRat) )<<endl;
  */

  cout<<"sampling: "<<sqrt(sH)<<" ME: "<<ME2<<endl;
  

  sigma = ME2;

  //cout<<sigma<<endl;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2ss::setIdColAcol() {

  // Flavours trivial.
  if(!anti)
    setId( 21, 21, id, (-1)*id);
  else
    setId( 21, 21, id, id);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2gg2ss::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  

  //WW channel only
  return 1.;

}

} // end namespace Pythia8
