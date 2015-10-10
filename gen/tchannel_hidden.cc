#include "Pythia8/Pythia.h"
#include "tchannel_hidden.hh"

void HiddenTChannel::initProc() {
  
  particlePtr = particleDataPtr->particleDataEntryPtr(idZprime);
  //get the mass of the Zprime
  mRes = particleDataPtr->m0(idZprime);
  m2Res  = mRes*mRes;

}

  
void HiddenTChannel::setIdColAcol() {

  setId( id1, id2, idHidden, -idHidden);
  
  // see if the decay product is a quark
  if (abs(idHidden) >= 1 && abs(idHidden) <= 6){
    //antiquark first
    if (id1 < 0)
      setColAcol( 0, 1, 1, 0, 2, 0, 0, 2);
    else
      setColAcol( 1, 0, 0, 1, 2, 0, 0, 2); 
  }
  
  else {
    //antiquark first
    if (id1 < 0)
      setColAcol( 0, 1, 1, 0, 0, 0, 0, 0);
    else
      setColAcol( 1, 0, 0, 1, 0, 0, 0, 0); 
  }
  
}

double HiddenTChannel::weightDecay( Event& process, int iResBeg, int iResEnd) {
  return 1.;
}


void HiddenTChannel::sigmaKin() {

  //spin trace for scalars
  //extra factor of two due to spin average
  double trace = pow2(tH);
  trace /= 4.0;
  
  //yukawa coupling strength
  double lambda=1.;
  //include the main resonance factor
  double ME2 = lambda*trace /
    ( pow2(tH - pow2(m_chi_tilde )) );

  ME2 /= 3.0;

  
  //the overall normalization is off at this moment
  sigma = ME2;
}
