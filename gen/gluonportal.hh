// SigmaHiggs.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// Part of code written by Marc Montull, CERN summer student 2007.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Higgs process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef gluonportal_H
#define gluonportal_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {


  //gluglu -> scalar scalar

class Sigma2gg2ss : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2ss(int id_out, int id_in=666, 
	      bool equal_anti = false,
	      double scale = 1000.0, 
	      double coupling= 1.0):
    id_in(id_in), anti(equal_anti),
    id(id_out), scale(scale), coupling(coupling) {}

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
  virtual string name()       const {return "g g -> H' -> s s";}
  virtual int    code()       const {return 10000;}
  virtual string inFlux()     const {return "gg";}
  virtual bool isSChannel() const {return true;}
  //virtual int    resonanceA() const {return id_in;}
  virtual int id3Mass()  const {return id;}
  virtual int id4Mass()  const {
    if(!anti) return -id;
    else return id;
      }
  virtual bool convertM2()  const {return true;}
private:

  bool anti;
  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int id_in;
  int id;
  double scale, coupling;
  ParticleDataEntry* particlePtr;

};
 
} // end namespace Pythia8

#endif // Pythia8_SigmaHiggs_H
