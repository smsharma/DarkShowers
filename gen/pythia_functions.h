#ifndef __pythia_functions_h
#define __pythia_functions_h

#include "fastjet/ClusterSequence.hh"
#include "tchannel_hidden.hh"

// Delphes library
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

const double PI = 3.14159265358979323846;

using namespace fastjet;

string add_st(string input, double number){
  stringstream sout;
  sout<<number;
  return input + sout.str();
}

string to_st(double number){
  stringstream sout;
  sout<<number;
  return sout.str();
}

string to_st(int number){
  stringstream sout;
  sout<<number;
  return sout.str();
}

//initialize t_channel processes
void init_tchannel(Pythia& pythia,
		   double mphi=1000.,
		   double pt_cut=500)
{
  Sigma2Process* myprocess = new HiddenTChannel(4900101, 666,0.01,mphi);

  // apply a phase-space cut to speed up MC generation
  pythia.readString(add_st("PhaseSpace:pTHatMin = ",
				pt_cut));
  pythia.setSigmaPtr(myprocess);     
}

void init_hidden(Pythia& pythia, 
		 double mass=20.0,
		 double lambda=1.0,
		 double inv = 0.3,
		 bool run=true,
		 int Nc=2,
		 int NFf=2,
		 int NBf=0
		 )
{
  
  pythia.readString("HiddenValley:Ngauge  = " +
		    to_st(Nc) );

  if(run){
    pythia.readString("HiddenValley:Run = on");
    pythia.readString(add_st
		      ("HiddenValley:Lambda = ", lambda));
    
    cout<<"running coupling: conf. scale "<<lambda<<endl;

    double b = 11./3 *Nc - 1./6. * NBf - 1./3. * NFf;

    cout<<"alpha(TeV) = "
	<<2*PI / (b*log(1000/lambda))<<endl;
      
  }

  else{
    
    pythia.readString("HiddenValley:Run = off");
    //change FSR strength
    pythia.readString
      (add_st("HiddenValley:alphaFSR = ", lambda));
  }

  //decouple the heavy flavor state
  pythia.readString("4900001:m0 = 5000");
  pythia.readString("4900002:m0 = 5000");
  pythia.readString("4900003:m0 = 5000");
  pythia.readString("4900004:m0 = 5000");
  pythia.readString("4900005:m0 = 5000");
  pythia.readString("4900006:m0 = 5000");
  pythia.readString("4900011:m0 = 5000");  
  pythia.readString("4900012:m0 = 5000");
  pythia.readString("4900013:m0 = 5000");
  pythia.readString("4900014:m0 = 5000");
  pythia.readString("4900015:m0 = 5000");
  pythia.readString("4900016:m0 = 5000");
  
  //mass of scalar dark quark
  //fix width to be small
  pythia.readString(add_st("4900101:m0 = ", mass/2));
  pythia.readString
    (add_st("4900101:mWidth = ", mass/100));
  pythia.readString
    (add_st("4900101:mMin = ", mass/2 - mass/100));
  pythia.readString
    (add_st("4900101:mMax = ", mass/2 + mass/100));
  
  //this fixes the qv to be spin 1/2
  //change to 1 for spin zero
  //pythia8 has a weird option where spinFv
  //is correlated with spin of qv
  pythia.readString("HiddenValley:spinFv = 0");
  pythia.readString("HiddenValley:FSR = on");
  pythia.readString("HiddenValley:fragment = on");

  //fix mass of dark scalar mesons
  //spin 0 diagonal 
  pythia.readString(add_st("4900111:m0 = ", mass));
  //spin 1 diagonal
  pythia.readString(add_st("4900113:m0 = ", mass));

  //we mock DM production by having diagonal meson decay
  //into a pair of charged dark meson
  //hence each DM is slightly less than half of 
  //diagonal meson mass

  //spin 0 charged (DM)
  pythia.readString(add_st("4900211:m0 = ", mass/2.0-0.01));
  //spin 1 charged (DM)
  pythia.readString(add_st("4900213:m0 = ", mass/2.0-0.01));
 
  //stop showering when pt less than threshold
  pythia.readString
    (add_st("HiddenValley:pTminFSR = ", mass));
    
  
  //do one flavor showering
  //flavor running is achieved elsewhere
  pythia.readString
    (add_st("HiddenValley:nFlav = ", 1));
  
  //running for N bosonic flavor
  pythia.readString
    (add_st("HiddenValley:NBFlavRun = ", NBf));
  
  //running for N fermionic flavor
  pythia.readString
    (add_st("HiddenValley:NFFlavRun = ", NFf));

  cout << "Turning on direct production" << endl;

  pythia.readString("HiddenValley:gg2UvUvbar = on");


  //onMode bRatio meMode product1 product2
  /*
    pythia.readString("4900111:onechannel = 1 " 
    + to_st(1-inv)
    + " 91 21 21");
  */
  
  //coupling to up-type quark only
  //implies primary decay to charms from MFV
  //due to helicity suppression
  //ignore light flavor
  pythia.readString("4900111:onechannel = 1 " 
		    + to_st(1.0-inv)
		    + " 91 -3 3");

  //invisible ratio
  //proportion ~inv of the time 
  //the meson is invisible
  pythia.readString("4900111:addchannel = 1 " 
		    + to_st(inv)
		    + " 0 4900211 -4900211");

  //spin 1 meson decay
  //democratic to all flavors
  pythia.readString("4900113:onechannel = 1 " 
		    + to_st((1-inv)/5.)
		      + " 91 -1 1");

  pythia.readString("4900113:addchannel = 1 " 
		    + to_st((1-inv)/5.)
		      + " 91 -2 2");
  
  pythia.readString("4900113:addchannel = 1 " 
		    + to_st((1-inv)/5.)
		    + " 91 -3 3");

  pythia.readString("4900113:addchannel = 1 " 
		    + to_st((1-inv)/5.)
		    + " 91 -4 4");
  
  pythia.readString("4900113:addchannel = 1 " 
		    + to_st((1-inv)/5.)
		    + " 91 -5 5");
  
  //spin 1 invisible ratio
  pythia.readString("4900113:addchannel = 1 " 
		    + to_st(inv)
		    + " 0 4900213 -4900213");
  
  //can also decay into 3 gluons through loops
  //should be seriously suppressed, commented
  //out for future references
  
  /*
  pythia.readString("4900113:onechannel = 1 1 92 21 21 21");
  pythia.readString("4900113:onechannel = 1 1 92 21 21 21");
  */

  //1:3 ratio for producing spin1 vectors
  pythia.readString
    (add_st("HiddenValley:probVector = ", 0.75));
  
}

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
    can->PID = p.id() ;
    can->Status = p.statusHepMC();
    can->Charge = p.charge();
    can->Mass = p.m();

    can->Momentum.SetPxPyPzE
      (p.px(), p.py(), p.pz(), p.e());
    
    ary->Add(can);
  }
}


//class to estimate remaining time
class Timer{
  
 public:
  vector<time_t> time_list;
  vector<int> proc_list;
  //int day, hr, min, sec;
  int nevt, max_interval;

 Timer(int nevt, int max_interval = 20):
  nevt(nevt), max_interval(max_interval)
  {
    time_t now;
    time(&now);
    time_list.push_back(now);
    proc_list.push_back(0);
  }


  void update(int nproc){
    
    //update number of processed events
    time_t now;
    time(&now);

    time_list.push_back(now);
    proc_list.push_back(nproc);
    
  }

  friend ostream& operator<< (ostream &out, Timer& mytime);
};

ostream& operator<< (ostream &out, const PseudoJet& myjet){
  out << myjet.pt()<< ","
      << myjet.m()<< ","
      << myjet.eta()<< ","
      << myjet.phi()<< ",";
  return out;
}

ostream& operator<< (ostream &out, Timer& mytime){

  //if no info default to zero
  if(mytime.proc_list.size() <= 1){
    out<<"\r\033[K"<<"0.0% processed... ";
    return out;
  }
  
  //first get two timing info for comparison
  time_t now = mytime.time_list.back();
  time_t past = mytime.time_list[0];

  int proc = mytime.proc_list.back();
  int proc_past = mytime.proc_list[0];

  //see if there are more info
  if(mytime.time_list.size() > mytime.max_interval){

    past = mytime.time_list[ mytime.time_list.size() - 
			     mytime.max_interval ];

    proc_past = mytime.proc_list[ mytime.proc_list.size() - 
				  mytime.max_interval ];
  }

  //estimate percentage done
  double percent = round (proc*1000.0 /mytime.nevt)/10.0;
  
  //estimate time needed
  double rate = (proc-proc_past)*1.0 / (now-past);
  int sec_left = round( (mytime.nevt - proc) / (rate) );

  int min = sec_left/60 ;
  int hr = min/60;
  int day = hr/24;
  int sec = sec_left % 60;

  min = min % 60;
  hr = hr % 24;

  stringstream sout;
  sout<<"\r\033[K"<<std::setprecision(1)<<fixed<<percent<<"% processed... ";

  if(day > 0)
    sout<<day<<"d "
	<<setfill('0') << setw(2)<<hr<<"h "
	<<setfill('0') << setw(2)<<min<<"m "
	<<setfill('0') << setw(2)<<sec<<"s";

  else if(hr > 0)
    sout<<hr<<"h "
	<<setfill('0') << setw(2)<<min<<"m "
	<<setfill('0') << setw(2)<<sec<<"s";

  
  else if(min > 0)
    sout<<min<<"m "
	<<setfill('0') << setw(2)<<sec<<"s";
  
  else if(sec > 0)
    sout<<sec<<"s";
  
  out<<sout.str()<<flush;

  return out;
}

int get_nmeson(const Pythia8::Event& evt){
  
  int pdgid_rho = 4900213;
  int pdgid_meson = 4900211;

  // int pdgid_rho = 4900113;
  // int pdgid_meson = 4900111;

  int n_meson = 0;
  for(int i=0; i<evt.size(); ++i){
    const Pythia8::Particle& p = evt[i];


    if(((abs(p.id()) == pdgid_rho) ||
       (abs(p.id()) == pdgid_meson)) && (p.statusHepMC() == 1) ){
      // cout << "pdgid " << p.id() << " status " << p.statusHepMC() << endl;    
      n_meson ++;
    }

    // if(abs(p.id()) != pdgid_rho &&
    //    abs(p.id()) != pdgid_meson )
    //   continue;


    // if(abs(abs(p.statusHepMC())-83) > 1)
    //   continue;
    
    
  }
  cout<<"total particle number: "<<n_meson<<endl;
  return n_meson;

}


int get_glu(const Pythia8::Event& evt){
  
  int pdgid_glu = 4900991;

  int n_glu = 0;
  for(int i=0; i<evt.size(); ++i){
    const Pythia8::Particle& p = evt[i];

    if(abs(p.id()) != pdgid_glu)
      continue;


    n_glu ++;
  }
  //cout<<"total particle number: "<<n_meson<<endl;
  return n_glu;

}


double dot3(const PseudoJet& a, const PseudoJet& b){
return a.px() * b.px() + a.py() * b.py() + a.pz() * b.pz();
}


double get_dphijj(const PseudoJet& met, const vector<PseudoJet>& jets){
  
  double result = 999;
  for(int i=0; i<4 && i<jets.size(); ++i){
    double t_dphi = fabs(met.phi() - jets[i].phi());
    if(t_dphi > PI)
      t_dphi = 2*PI - t_dphi;
    result = result < t_dphi ? result : t_dphi;
  }
  return result;

}

double get_Mt(const PseudoJet& MEt,
	      const PseudoJet& jj){

  double ET= jj.mperp();
  double MT = sqrt(jj.m2()+ 2*
		   (MEt.pt()*jj.pt() -
		    dot3(jj, MEt)));

  return MT;
}


double get_Mt(const PseudoJet& MEt,
	      const PseudoJet& jet1_,
	      const PseudoJet& jet2_){

  return get_Mt(MEt, jet1_ + jet2_);
}




#endif

