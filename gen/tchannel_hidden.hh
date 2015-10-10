#ifndef tchannel_hidden_H
#define tchannel_hidden_H

#include "Pythia8/SigmaProcess.h"

using namespace Pythia8;

// A derived class for q qbar -> Zprime -> hidden

class HiddenTChannel : public Sigma2Process {

public:

  // Constructor.
  // idZprime = custom pdgid of Zprime
  // idHidden = id of the particle Zprime decays into
  // Zprime decays to +idHidden and -idHidden
  HiddenTChannel(int idHidden, int idZprime=666, 
		 double width=0.01, double m_chi_tilde=1000.): 
    idHidden(idHidden), idZprime(idZprime), width(width), m_chi_tilde(m_chi_tilde) {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). Assumed flavour-independent so simple.
  virtual double sigmaHat() {
    return sigma;
  }

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Pythia8::Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual string name()       const 
  {return "q qbar -> t-channel -> hidden";}
  virtual int    code()       const {return 10000;}
  virtual string inFlux()     const {return "qqbarSame";}
  //virtual int    resonanceA() const {return idZprime;}
  virtual bool isSChannel() const {return false;}
  virtual int id3Mass()  const {return idHidden;}
  virtual int id4Mass()  const {return -idHidden;}
  
  virtual bool   convertM2() const  {return true;}
  
private:

  // Store flavour-specific process information and standard prefactor.
  int idHidden;
  int idZprime;
  double width;
  double mRes, GammaRes, m2Res, GamMRat, normTheta2qqbar, sigma;
  double alpha_dark, thetaWRat;
  double m_chi_tilde;
  // Pointer to properties of Theta, to access decay width.
  ParticleDataEntry* particlePtr;

};





#endif
  
