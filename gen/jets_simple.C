// Tim Lou
// 10/04/2015

// Modified:
// 10/13/2015: Siddharth Mishra-Sharma

// C++ tools
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

// ROOT
#include "TROOT.h"
#include "TRandom.h"
#include "TApplication.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

// Delphes library
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

// FastJet?
// Duplicates in Delphes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

// Pythia
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/CombineMatchingInput.h" //  Matching not implemented yet

// Pythia8 library to perform t-channel production
#include "tchannel_hidden.hh"

//Pythia8 library to perform higher dimensional ops
//dim-7 GG-chichi production
//#include "gluonportal.hh"

#include "CmdLine/CmdLine.hh"

//extra tools, to be simplified
//#include "tools.h"

//to simplify code
//moving all unnecessary functions to this file
#include "pythia_functions.h"

//using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace std;  

int main(int argc, char** argv) {

  cout<<"Usage: -m (mode) -n (nevent = 100) -o (output) -pt_min (100) -mphi (1000) -metmin (0) -phimass (default=20) -alpha (dark confinement scale) -frag (fragmentation) -inv (invisible ratio) -v (verbose) -seed (0) -rehad (off)"<<endl;

  //parse input strings
  CmdLine cmdline(argc, argv);
  
  string mode = cmdline.value<string>("-m", "tchannel");
  double pt_min = cmdline.value<double>("-ptmin", 100);
  double met_min = cmdline.value<double>("-metmin", 800);
  double met_max = cmdline.value<double>("-metmax", 99999);
  double dphi_max = cmdline.value<double>("-dphimax", 999);

  string output = cmdline.value<string>("-o", "output");

  bool rehad = cmdline.present("-rehad");

  // Instantiate event-wide, object and info files
  // file_evt stores event wide variables
  // also last two line stores cxn, efficiency...etc
  ofstream file_evt;
  ofstream file_meta;

  //file_obj stores object info for each event
  //i.e. jets, leptons, dark mesons...etc
  ofstream file_obj;
  string input;

  if (mode == "lhe")
    input = cmdline.value<string>("-i");  
  try
  {
    file_evt.open((output + ".evt").c_str());
    file_obj.open((output + ".obj").c_str());
    file_meta.open((output + ".meta").c_str());
  }
  catch(...)
  {
    cerr<<"ERROR: cannot open "<<output<<", exiting..."<<endl;
    return 1;
  }

  int nEvent = cmdline.value<int>("-n", 100);
  int nAbort = 5;

  // Pythia generator
  Pythia pythia;

  // Initialization for LHC
  pythia.readString("Beams:eCM = " + 
	   to_st(cmdline.value<int>("-ECM", 13000)));

  // Custom Higgs pT cut
  /*
  ParticlePTCut* HiggsPTCut = new ParticlePTCut(pt_min,25);
  pythia.setUserHooksPtr(HiggsPTCut);  
  */

  bool m_lhe = false;

  // Check for verbose mode
  if(!cmdline.present("-v"))
    pythia.readString("Print:quiet = on");
  
  //hidden scalar production
  if (mode == "tchannel"){ 

    double mphi=
      cmdline.value<double>("-mphi", 1000.0);

    init_tchannel(pythia, mphi, 0.0);

    init_hidden(pythia,
		cmdline.value<double>("-phimass", 20.0),
		cmdline.value<double>("-alpha", 10), // Actually just the confinement scale, not alpha_dark
		cmdline.value<double>("-inv", 0.3)
		);  
  }  
  
  else if(mode == "lhe"){
    //read lhe file
    m_lhe=true;
    
    CombineMatchingInput combined;
    UserHooks* matching = combined.getHook(pythia);
    if (!matching) {
      cout<<"ERROR: cannot obtain matching pointer"<<endl;
      return 1;
    }
    
    pythia.setUserHooksPtr(matching);

    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Beams:LHEF = "+ input);
    pythia.readString("Beams:frameType = 4");
    pythia.readString("JetMatching:merge = on");
    pythia.readString("JetMatching:scheme = 1");
    pythia.readString("JetMatching:setMad = on");
    pythia.readString("JetMatching:jetAlgorithm = 2");
    pythia.readString("JetMatching:exclusive = 2");
    pythia.readString("JetMatching:nJetMax = " + to_st(cmdline.value<int>("-nmatch", 2)));
    cout<<"input: "<<"Beams:LHEF = "+input<<endl;
  } 
  
  else
  {
    cerr<<"ERROR: mode: " << mode << " not supported, exiting..." << endl;
    return 1;
  }
  
  // Set seed
  pythia.readString("Random:setSeed = on");
  // Fix random seed
  pythia.readString("Random:seed = " + 
		    to_st(cmdline.value<int>("-seed", 0)));

  // Rehadronization turned on
  if(rehad)
    pythia.readString("HadronLevel:all = off");
  
  // Initialize Pythia

  pythia.init();
  
  cout<<"INFO: pT cut: "<< pt_min <<endl;
  cout<<"INFO: MEt cut: "<< met_min <<endl;

  // Begin event loop.
  int iAbort = 0;
  
  //event level variable
  file_evt << "evt,Mt,Mjj,MEt,dphi,deta,"
	   << "nj,nl,ndark,ndiag,MR,R2,aT,Mjj_mc" 
	   << endl;


  //obj level variable
  file_obj << "evt,type,n,pt,m,eta,phi,ntrk"<<endl;

  // Access to pythia event
  Pythia8::Event& event = pythia.event;
  Pythia8::Event saved_event;

  int iEvent = 0;
  int iTotal = 0;

  // Start timer

  Timer mytime(nEvent);

  // Declare Delphes variables

  ExRootConfReader *config = new ExRootConfReader();
  config->ReadFile("delphes_card_CMS.tcl");
  Delphes *delphes = new Delphes("Delphes"); 
  delphes -> SetConfReader(config);
  DelphesFactory *factory = delphes->GetFactory();

  TObjArray* stable = 
    delphes->ExportArray("stableParticles");

  // Start root application mode

  gROOT->SetBatch();
  int appargc = 1;
  char appName[] = "Delphes";
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  // Initialize delphes

  delphes->InitTask();

  // Events to print
  int evt_print = 20;

  // Running simulation
  bool end=false;

  while ((!m_lhe && (iEvent < nEvent)) || 
	 (m_lhe && (iTotal < nEvent) && !end)) 
  {
    // Clear delphes
    delphes->Clear();

    // If rehadronization is turned on
    if(rehad) {
      
      // Renew an event
      if (iTotal % 5 == 0) {
      	while (!pythia.next()) {
	  
	  if(pythia.info.atEndOfFile()){
	    cout <<"Pythia reached end of file"<<endl;
	    end=true;
	    break;    
	  }

	  
	  if (++iAbort < nAbort) continue;
	  
	  cerr << "ERROR: Event generation aborted prematurely, owing to error!" << endl;
	  break;
	}
	
	saved_event = pythia.event;
      }

      else 
	pythia.event = saved_event;
	
      
      // Run hadronization
      pythia.forceHadronLevel();
    }
    else
    {
      // Tell pythia to run pythia.next()
      while (!pythia.next()) {

	  if(pythia.info.atEndOfFile()){
	    cout <<"Pythia reached end of file"<<endl;
	    end=true;
	    break;    
	  }
	  
	  if (++iAbort < nAbort) continue;
	  cerr<<"ERROR: Event generation aborted due to error!" 
	       <<endl;
	  break;
      }
    }


    // Increment tried events
    ++iTotal;

    // Now process through Delphes
    Pythia_to_Delphes(factory, stable, event);
    
    // Run delphes code
    delphes->ProcessTask();

    Candidate *can;

    const TObjArray* vMEt = delphes->ImportArray
      ("MissingET/momentum");
    
    can = (Candidate*) TIter(vMEt).Next() ;

    // Missing ET pointer must exist
    if(can == NULL){
      cout<<"ERROR: MET pointer not found!"<<endl;
      continue;
    }

    PseudoJet MEt
      (-can->Momentum.Px(),-can->Momentum.Py(), 
       0,can->Momentum.Pt());    
    
    // Now grab the jets
    const TObjArray* jets = delphes->ImportArray
      ("UniqueObjectFinder/jets");

    // Now grab muons
    const TObjArray* muons = delphes->ImportArray
      ("UniqueObjectFinder/muons");      

    // Now grab the electrons
    const TObjArray* electrons = delphes->ImportArray
      ("UniqueObjectFinder/electrons");      

    // Demand at least two jets above pt cut
    int njet = cmdline.value<int>("-njet", 2);

    if(njet < 0){
      cout<<"ERROR: cannot require negative jets"<<endl;
      njet=0;
    }

    //grab objects
    vector<PseudoJet> selected_jets;
    vector<int> jets_ntrk;
    vector<PseudoJet> selected_leptons;
    
    // Loop over muons and get information
    for(int i=0; i<muons->GetEntriesFast(); i++){
      
      Candidate* cmuon = (Candidate*) muons->At(i);
      if(fabs(cmuon->Momentum.Eta())>2.5)
        continue;
      
      PseudoJet cmuon_v(cmuon->Momentum.Px(), 
           cmuon->Momentum.Py(),
           cmuon->Momentum.Pz(),
           cmuon->Momentum.E());

      selected_leptons.push_back(cmuon_v);
    }    

    // Loop over electrons
    for(int i=0; i<electrons->GetEntriesFast(); i++){
      
      Candidate* celectron = (Candidate*) electrons->At(i);
      if(fabs(celectron->Momentum.Eta())>2.5)
        continue;
      
      PseudoJet celectron_v(celectron->Momentum.Px(), 
           celectron->Momentum.Py(),
           celectron->Momentum.Pz(),
           celectron->Momentum.E());

      selected_leptons.push_back(celectron_v);
    }    

    // Loop over jets and get information
    for(int i=0; i<jets->GetEntriesFast(); i++){
      
      Candidate* cjet = (Candidate*) jets->At(i);

      if(cjet->Momentum.Pt()<50.0)
	continue;
	
      if(fabs(cjet->Momentum.Eta())>3.0)
	continue;

      //get constituents
      TObjArray* jet_cands = cjet->GetCandidates(); 

      //count the number of charged
      //tracks within the jet
      int ntrk=0;

      for(int j=0; j<jet_cands->GetEntriesFast(); j++){
        Candidate* jet_cand = (Candidate*) jet_cands->At(j);
        if (jet_cand->Momentum.Pt() > 1.0 && 
	    jet_cand->Charge != 0)
        {
	  /*
          PseudoJet vcand(jet_cand->Momentum.Px(), 
          jet_cand->Momentum.Py(),
          jet_cand->Momentum.Pz(),
          jet_cand->Momentum.E());     
	  */
          ntrk++;          
        }
      }
      
      PseudoJet cjet_v(cjet->Momentum.Px(), 
		       cjet->Momentum.Py(),
		       cjet->Momentum.Pz(),
		       cjet->Momentum.E());

      //store the jets
      selected_jets.push_back(cjet_v);
      jets_ntrk.push_back(ntrk);
    }


    //demand njets > pt_min
    if(selected_jets.size() < njet || selected_jets[0].pt() < pt_min) 
      continue;

    if (get_dphijj(MEt, selected_jets) > dphi_max)
      continue;

    //print the jets
    for(int i=0; i<selected_jets.size(); i++){
      file_obj << iEvent << ",j,"
	       << i+1 << ","
	       << selected_jets[i] 
	       << jets_ntrk[i] <<endl;  
    }    

    //print the leptons
    for(int i=0; i<selected_leptons.size(); i++){
      file_obj << iEvent << ",l,"
	       << i+1 << ","
	       << selected_leptons[i] <<"1"
	       << endl;  
    }
    //print met
    file_obj << iEvent <<",met,"
	     << 1 << "," 
	     << MEt<<"0"<<endl;

    int n_diag=0;
    int n_dark=0;
    vector<PseudoJet> dark_mesons;
    PseudoJet inv_jet;

    // gather MC information
    for (int i=0; i<event.size(); ++i){
      Particle& p = event[i];

      if(p.id()==4900111)
        n_diag++;

      // Only get invisible pions
      if(p.id()!=4900211) continue;
      
      // We expect the next particle to be an anti-pion
      // if not, there is an error!
      if(i+1 >= event.size() || event[i+1].id() != -4900211)
	{
	  cout 
	    << "ERROR: invisible pions not coming in pairs!" 
	    << endl;
	  continue;
	}

      n_dark++;
      
      Particle& p_next = event[i+1];

      PseudoJet
	meson( p.px() + p_next.px(),
	       p.py() + p_next.py(),
	       p.pz() + p_next.pz(),
	       p.e() + p_next.e() );
      
      dark_mesons.push_back( meson );
      inv_jet += meson;
    }

    PseudoJet jj;
    //vector of first two jets
    //could have less than two jets
    for(int i=0; i<selected_jets.size() && i<2; ++i)
      jj += selected_jets[i];
    
    double Mt = get_Mt(MEt, jj);
    double Mjj = jj.m();
    double dphijj = get_dphijj(MEt, selected_jets);
    
    double deta = -1;
    if(selected_jets.size() >= 2)
      deta = fabs(selected_jets[0].eta()-
		  selected_jets[1].eta());
    
    file_evt<<iEvent<<","
	    << Mt <<","
	    << Mjj<<","
	    << MEt.pt()<<","    
	    << dphijj<<","
	    << deta<<","
	    << selected_jets.size()<<","
	    << selected_leptons.size()<<","
	    << n_dark<<","
	    << n_diag<<","
	    <<endl;

    ++iEvent;
    
    if(!m_lhe && (iEvent %10 ==0))
    {
      mytime.update(iEvent);
      cout<<mytime;
    }

    else if(m_lhe && (iTotal %10 ==0))
    {
      mytime.update(iTotal);
      cout<<mytime;
    }
  }

  if(!m_lhe)
    mytime.update(iEvent);

  else
    mytime.update(iTotal);

  cout<<mytime<<endl;
  cout<<iEvent<<" total events"<<endl;

  file_meta<<"nevt, npass, eff, total, pass, "
	   <<"ptcut, metcut, cxn, cxn_err"<<endl;
  file_meta<<iTotal<<","<<iEvent<<","
	   <<iEvent/double(iTotal)<<","
	   <<iTotal<<","
	   <<iEvent<<","
	   <<pt_min<<","
	   <<met_min<<","
	   <<pythia.info.sigmaGen()*1e9<<","
	   <<pythia.info.sigmaErr()*1e9<<endl;
  
  //clean up
  delphes->FinishTask();
  delete delphes;
  delete config;
  
  if(cmdline.present("-v"))
    pythia.stat();
  // Done.
  file_evt.close();
  file_obj.close();
  return 0;

}
