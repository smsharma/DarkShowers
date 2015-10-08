#ifndef __pythia_functions_h
#define __pythia_functions_h


// Pythia UserHook to apply parton level pt cut
class ParticlePTCut: public UserHooks {

public:

  double ptcut;
  int id;
  PseudoJet h;

  ParticlePTCut(double ptcut, int id=25): 
    ptcut(ptcut), id(id) {}

  virtual bool canVetoProcessLevel()
  {
    return true;
  }

  virtual bool doVetoProcessLevel(Pythia8::Event& evt)
  {
    for(int i = 0; i < evt.size(); ++i){
      Pythia8::Particle& p = evt[i];

      if(p.id()!= id)
        continue;
   
      h = PseudoJet(p.px(), p.py(), p.pz(), p.e());
 
      if (p.pT() < ptcut)
        return true;
    }
    return false;
  }

};

//intialize QCD events, for debugging only
void init_qcd(Pythia& pythia, string mode)
{

  pythia.readString("HardQCD:all = off");

  if (mode == "qcd" || mode == "gg"){
    pythia.readString("HardQCD:gg2gg = on");
    pythia.readString("HardQCD:qqbar2gg = on");
  }

  if (mode == "qcd" || mode == "qq"){
    pythia.readString("HardQCD:gg2qqbar = on");
    pythia.readString("HardQCD:qq2qq = on");
    pythia.readString("HardQCD:qqbar2qqbarNew = on");
  }

  if (mode == "qcd" || mode == "gq" || mode == "qg"){
    pythia.readString("HardQCD:qg2qg = on");    
  }
  
  cout<<"INFO: pp -> "<<mode<<" enabled "<<endl;
}

// Function to input all particles to Delphes
void Pythia_to_Delphes(DelphesFactory* factory, 
		       TObjArray* ary,
		       Pythia8::Event& evt,
		       vector<PseudoJet>* muons=NULL){

  // Loop over particles
  for(int i=0; i<evt.size(); ++i){
    Pythia8::Particle& p = evt[i];

    
    if(!p.isVisible() || (p.statusHepMC()) != 1)
	    continue;
    
    Candidate* can = factory->NewCandidate();
    can->PID = p.id();
    can->Status = p.statusHepMC();
    can->Charge = p.charge();
    can->Mass = p.m();

    can->Momentum.SetPxPyPzE
      (p.px(), p.py(), p.pz(), p.e());
    
    ary->Add(can);
  }
}



#endif

