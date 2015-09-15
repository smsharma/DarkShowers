// SigmaHiggs.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// Part of code written by Marc Montull, CERN summer student 2007.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Higgs process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef gluonfusion_H
#define gluonfusion_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {


// The ResonanceTheta class handles a F-squarkonium resonance
class ResonanceTheta : public ResonanceWidths {

public:

  // Constructor.
  ResonanceTheta(int idResIn=663) {
    initBasic(idResIn);
  }

private:

  double totalWidth;
 
  // Initialize constants.
  virtual void initConstants();
 
  // Calculate various common prefactors for the current mass.
  // Superfluous here, so skipped.
  //virtual void calcPreFac(bool = false);

  // Calculate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

 
// A derived class for g g -> H0 (SM) to resonance

class Sigma1gg2Theta : public Sigma1Process {

public:

  // Constructor.
  Sigma1gg2Theta(): idTheta(663) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const {return "g g -> h* -> Theta";}
  virtual int    code()       const {return 621;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return idTheta;}
  
private:

  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int idTheta;
  ParticleDataEntry* particlePtr;

};
 
} // end namespace Pythia8

#endif // Pythia8_SigmaHiggs_H
