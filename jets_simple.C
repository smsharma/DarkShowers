// Tim Lou
// 10/04/2015

// c++ tools
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>

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

//fastjet??
//duplicates in delphes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "Pythia8/Pythia.h"
//matching not implemented yet
//#include "CombineMatchingInput.h"

//Pythia8 library to perform t-channel production
#include "tchannel_hidden.hh"

//Pythia8 library to perform higher dimensional ops
//dim-7 GG-chichi production
//#include "gluonportal.hh"

#include "CmdLine/CmdLine.hh"

//extra tools, to be simplified
#include "tools.h"

//to simplify code
//moving all unnecessary functions to this file
#include "pythia_functions.h"

//using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace std;  


//initialize t_channel processes
void init_tchannel(Pythia& pythia,
		   double m_chi_tilde=1000.)
{
  Sigma2Process* myprocess = new HiddenTChannel(4900101, 666,0.01,m_chi_tilde);

  // apply a phase-space cut to speed up MC generation
  pythia.readString("PhaseSpace:pTHatMin =  1000");
  pythia.setSigmaPtr(myprocess);     
}

void init_hidden(Pythia& pythia, 
		 double mass=20.0,
		 double frag=0.8,
		 double alpha=0.1,
		 double mZ = 1000,
		 double inv = 0.3,
		 bool run=false
		 )
{
  
  if(run){
    pythia.readString("HiddenValley:Run = on");
    pythia.readString(add_strings("HiddenValley:Lambda = ", alpha));
    cout<<"running enabled, alpha(1TeV) is: "<<alpha<<endl;
  }
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

    pythia.readString(add_strings("4900101:m0 = ", mass/2));
    pythia.readString
      (add_strings("4900101:mWidth = ", mass/100));
    pythia.readString
      (add_strings("4900101:mMin = ", mass/2 - mass/100));
    pythia.readString
      (add_strings("4900101:mMax = ", mass/2 + mass/100));

    pythia.readString("HiddenValley:spinqv = 0");
    pythia.readString("HiddenValley:FSR = on");
    pythia.readString("HiddenValley:fragment = on");


    pythia.readString(add_strings("4900111:m0 = ", mass));
    pythia.readString(add_strings("4900113:m0 = ", mass));

    pythia.readString(add_strings("4900211:m0 = ", mass/2-0.01));
    pythia.readString(add_strings("4900213:m0 = ", mass/2-0.01));

    pythia.readString(add_strings("HiddenValley:bmqv2 = ", frag));

    pythia.readString
      (add_strings("HiddenValley:pTminFSR = ", mass));

    
    pythia.readString
      (add_strings("HiddenValley:nFlav = ", 1));

    //onMode bRatio meMode product1 product2
    /*
    pythia.readString("4900111:onechannel = 1 " 
		      + to_st(1-inv)
		      + " 91 21 21");
    */

    
    pythia.readString("4900111:onechannel = 1 " 
		      + to_st((1-inv)/5.)
		      + " 91 -1 1");

    pythia.readString("4900111:addchannel = 1 " 
		      + to_st((1-inv)/5.)
		      + " 91 -2 2");

    pythia.readString("4900111:addchannel = 1 " 
		      + to_st((1-inv)/5.)
		      + " 91 -3 3");

    pythia.readString("4900111:addchannel = 1 " 
		      + to_st((1-inv)/5.)
		      + " 91 -4 4");

    pythia.readString("4900111:addchannel = 1 " 
		      + to_st((1-inv)/5.)
		      + " 91 -5 5");

    pythia.readString("4900111:addchannel = 1 " 
		      + to_st(inv)
		      + " 0 4900211 -4900211");

    pythia.readString("4900113:onechannel = 1 1 92 21 21 21");

    pythia.readString("4900113:onechannel = 1 1 92 21 21 21");

    //disable spin 1 mesons
    pythia.readString
      (add_strings("HiddenValley:probVector = ", 0.75));

    //change FSR strength
    pythia.readString
      (add_strings("HiddenValley:alphaFSR = ", alpha));
}


int main(int argc, char** argv) {

  cout<<"Usage: -m (mode -- source) -n (nevent = 100) -o (output.txt) -pt_min (200)  -met_min (0) -phimass (default=0.5) -alpha (dark FSR coupling) -frag (fragmentation) -inv (invisible ratio) -v (verbose) -seed (0) -rehad (off)"<<endl;

  
  //parse input strings
  CmdLine cmdline(argc, argv);
  
  string mode = cmdline.value<string>("-m", "t-channel");
  double pt_min = cmdline.value<double>("-ptmin", 100);
  double met_min = cmdline.value<double>("-metmin", 0);
  double met_max = cmdline.value<double>("-metmax", 99999);
  double dphi_max = cmdline.value<double>("-dphimax", 999);

  string output = cmdline.value<string>("-o", "output");

  bool rehad = cmdline.present("-rehad");

  // Initialize random numbers
  ran.SetSeed(cmdline.value<int>("-seed", 0));

  // Instantiate event-wide, object and info files
  ofstream file_evt;
  ofstream file_obj;
  ofstream file_meta;
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

  ParticlePTCut* HiggsPTCut = new ParticlePTCut(pt_min,25);
  pythia.setUserHooksPtr(HiggsPTCut);  

  bool m_lhe = false;

  // Check for verbose mode

  if(!cmdline.present("-v"))
    pythia.readString("Print:quiet = on");

  double mZ = cmdline.value<double>("-mZ", 1000);
  //hidden scalar production
  if (mode == "ss")
    init_hidden(pythia,
		cmdline.value<double>("-phimass", 15.0),
		cmdline.value<double>("-frag", 0.8),
		cmdline.value<double>("-alpha", 0.1),
		mZ,
		cmdline.value<double>("-inv", 0.3),
		cmdline.present("-run")
		);
  else if (mode == "tchannel"){ 
    init_hidden(pythia,
    cmdline.value<double>("-phimass", 15.0),
    cmdline.value<double>("-frag", 0.8),
    cmdline.value<double>("-alpha", 0.1),
    mZ,
    cmdline.value<double>("-inv", 0.3),
    cmdline.present("-run")
    );    
    init_tchannel(pythia,
    cmdline.value<double>("-mchitilde", 1000.0)
    );
  }  
  else if(mode == "qcd" ||
    mode == "qq" || 
    mode == "gq" ||
    mode == "qg" ||
    mode == "gg" )
  { 
    init_qcd(pythia, mode);
    pythia.readString("PhaseSpace:pTHatMin =  300");
  }
  else if(mode == "lhe"){
    //read lhe file
    m_lhe=true;
    
    // CombineMatchingInput combined;
    // UserHooks* matching = combined.getHook(pythia);
    // if (!matching) return 1;
    // pythia.setUserHooksPtr(matching);


    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = "+ input);

    // pythia.readString("JetMatching:merge = on");
    // pythia.readString("JetMatching:scheme = 1");
    // pythia.readString("JetMatching:doShowerKt = on");
    // pythia.readString("JetMatching:setMad = on");
    /*
    pythia.readString("JetMatching:qCut = 30.0000");
    pythia.readString("JetMatching:coneRadius = 1.0");
    pythia.readString("JetMatching:etaJetMax = 10.0");
    pythia.readString("JetMatching:nJetMax = 3");
    */

    /*    
    //pythia.readString("JetMatching:nJet = 2");
    //pythia.readString("JetMatching:doMerge = 1");
    pythia.readString("JetMatching:nQmatch = 5");

    pythia.readString("JetMatching:jetAlgorithm = 2");
    pythia.readString("JetMatching:slowJetPower = -1");
    pythia.readString("JetMatching:nQmatch = 5");
    */

    
    if (cmdline.present("-inclusive"))
    {
      cout<<"Inclusive mode, more jets allowed"<<endl;
      pythia.readString("JetMatching:exclusive = 0");
    }

    else 
    {
      cout<<"Exclusive mode, all jets must match"<<endl;
      pythia.readString("JetMatching:exclusive = 1");
    }
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
  pythia.readString("Random:seed = " + int_st(cmdline.value<int>("-seed", 0)));

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

  // file_evt << "evt,tag,met,dphi_higgs_met,dphi_jet1_met,dphi_jet2_met,dphi_jet1_jet2,dphi_min_jet_met,m_mu_mu,m_t,num_jets,num_reco_muons,"
  // <<"num_gen_muons,num_bs,num_charged,num_diag,num_dark,higgs_pt,nSub11,nSub12,nSub21,nSub22" << endl;


  // file_evt<<"evt, tag, met, dphi, vbf_m, "
  //   <<"hpt, hM, hMt, hMs, hMb, mc_hM, "
  //   <<"nj, nmu, mc_nmu, "
  //   <<"ndiag, ndark, higgs_pt_cut"
  //   <<endl;

  file_evt << "evt,mt,met,mj1,mj2,mjj,dphi_jj,dphi_min,eta_j1,eta_j2,deta_jj,num_tracks_jet,nSub11,nSub12,nSub21,nSub22" << endl;


  //object level variable
  file_obj<<"evt, type, n, pt, m, eta, phi"<<endl;


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

  while ((!m_lhe && (iEvent < nEvent)) || (m_lhe && (iTotal < nEvent))) 
  {
    // Clear delphes
    delphes->Clear();

    // If rehadronization is turned on
    if(rehad) {
      
      // Renew an event
      if (iTotal % 5 == 0) {
      	while (!pythia.next()) {
  	      if (++iAbort < nAbort) continue;
  	      cerr << "ERROR: Event generation aborted prematurely, owing to error!" << endl;
  	      break;
  	    }
  	
  	    saved_event = pythia.event;
      }

      else 
      {
	      pythia.event = saved_event;
      }
      
      // Run hadronization
      pythia.forceHadronLevel();
    }
    else
    {
      // Tell pythia to run pythia.next()
      while (!pythia.next()) {
	      if (++iAbort < nAbort) continue;
        cerr << "ERROR: Event generation aborted prematurely, owing to error!" << endl;
	      break;
      }
    }
    if(m_lhe && pythia.info.atEndOfFile())
      break;    

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

    // If met is too small or large, continue
    if ((can->Momentum.Pt() < met_min) || (can->Momentum.Pt() > met_max)){
      continue;
    }
    PseudoJet MEt(-can->Momentum.Px(),-can->Momentum.Py(), 0,can->Momentum.Pt());    
    
    // Now grab the jets
    const TObjArray* vjets = delphes->ImportArray
      ("FastJetFinder/jets");

    // Now grab muons
    const TObjArray* muons = delphes->ImportArray
      ("MuonIsolation/muons");      

    // Now grab the electrons
    const TObjArray* electrons = delphes->ImportArray
      ("UniqueObjectFinder/electrons");      

    // Demand at least two jets above pt cut
    int njet = cmdline.value<int>("-njet", 2);

    if(vjets->GetEntriesFast() < njet)
      continue;

    vector<PseudoJet> selected_jets;
    vector<int> njet_list;

    vector<PseudoJet> selected_muons;
    vector<int> nmuon_list;

    const float zMass = 91.1876;

    int num_gen_muons = 0;
    int num_b_mesons = 0;
    int num_charged = 0;
    int num_vis = 0;
    int num_invis = 0;

    // Loop over generator level particles and get muon information

    for(int i=0; i<event.size(); ++i)
    {
      Pythia8::Particle& p = event[i];

      // Save  objects for scatter plots

      if (p.pT() > 1.0)
      {
        if (p.isFinal() && p.isVisible())
        {
          PseudoJet vcand(p.px(), 
          p.py(),
          p.pz(),
          p.e());     
          file_obj<<print_obj(iEvent, num_vis, "vis", vcand);
          num_vis++;
        }

        if (p.isFinal() && !(p.isVisible()))
        {
          PseudoJet vcand(p.px(), 
          p.py(),
          p.pz(),
          p.e());     
          file_obj<<print_obj(iEvent, num_invis, "invis", vcand);
          num_invis++;
        }     

       if (p.isCharged() && p.isVisible())
        {
          PseudoJet vcand(p.px(), 
          p.py(),
          p.pz(),
          p.e());     
          file_obj<<print_obj(iEvent, num_charged, "charged", vcand);
          num_charged++;
        }

      }

      // Number of muons and B mesons

      if (fabs(p.id())==13)
      {
        if (p.isFinal())
          num_gen_muons++;
      }    
      if (((abs(p.id())>510) && (abs(p.id())<546)) || p.id() == 10511 || p.id() == 10521 || p.id() == 10513 || p.id() == 10523 || p.id() == 20513 || p.id() == 20523 || p.id() == 10531 || p.id() == 10533 || p.id() == 20533 || p.id() == 10541 || p.id() == 10543 || p.id() == 20543)
      {
        if (abs(p.status()) == 83 || abs(p.status()) == 84)
          num_b_mesons++;
      }  
    }

    int num_reco_muons = 0;

    // Loop over reco muons and get information

    for(int i=0; i<muons->GetEntriesFast(); i++){
      
      Candidate* cmuon = (Candidate*) muons->At(i);

      if(fabs(cmuon->Momentum.Eta())>2.5)
        continue;
      
      PseudoJet cmuon_v(cmuon->Momentum.Px(), 
           cmuon->Momentum.Py(),
           cmuon->Momentum.Pz(),
           cmuon->Momentum.E());
      num_reco_muons++;

      selected_muons.push_back(cmuon_v);
      nmuon_list.push_back(selected_muons.size());
    }    

    // Loop over jets and get information

    int num_tracks_jet = 0;
    // Get charged tracks in leading two tracks; otherwise if want all jets, 
    // i<vjets->GetEntriesFast()
    for(int i=0; i<2; i++){
      
      Candidate* cjet = (Candidate*) vjets->At(i);

      if(fabs(cjet->Momentum.Eta())>3.0)
	      continue;

      TObjArray* jet_cands = cjet->GetCandidates(); 
      for(int j=0; j<jet_cands->GetEntriesFast(); j++){
        Candidate* jet_cand = (Candidate*) jet_cands->At(j);
        // cout << jet_cand->Momentum.Pt() << endl;
        if (jet_cand->Momentum.Pt() > 2.0)
        {
          PseudoJet vcand(jet_cand->Momentum.Px(), 
          jet_cand->Momentum.Py(),
          jet_cand->Momentum.Pz(),
          jet_cand->Momentum.E());     
          file_obj<<print_obj(iEvent, num_tracks_jet, "jet_track", vcand);
          num_tracks_jet++;          
        }
      }


      PseudoJet cjet_v(cjet->Momentum.Px(), 
		       cjet->Momentum.Py(),
		       cjet->Momentum.Pz(),
		       cjet->Momentum.E());

      selected_jets.push_back(cjet_v);
      njet_list.push_back(selected_jets.size());
    }

    if(selected_jets.size() < njet)
      continue;

    if(selected_jets[0].pt() < pt_min)
      continue;    

    if (dphi(MEt, selected_jets) > dphi_max)
      continue;

    // Get the HT

    const TObjArray* vHT = delphes->ImportArray
      ("ScalarHT/energy");
    
    can = (Candidate*) TIter(vHT).Next() ;

    // Missing ET pointer must exist

    if(can == NULL){
      cout<<"ERROR: HT pointer not found!"<<endl;
      continue;
    }

    double HT = can->Momentum.M();  

    // Now grab the invisible guys

    vector<PseudoJet> unjet_input;
    PseudoJet all_inv;
    int num_diag=0;
    int num_dark=0;

    // Loop over particles

    for (int i=0; i<event.size(); ++i){
      Particle& p = event[i];

      if(p.id()==4900111)
        num_diag++;
      // Only get invisible pions
      if(p.id()!=4900211) continue;
      
      // We expect the next particle to be an anti-pion
      // if not, there is an error!

      if(i+1 >= event.size() || event[i+1].id() != -4900211)
      {
        cout << "ERROR: invisible pions not coming in pairs!" << endl;
        continue;
      }

      num_dark++;
      
      Particle& p_next = event[i+1];

      unjet_input.push_back( PseudoJet
           ( p.px() + p_next.px(),
             p.py() + p_next.py(),
             p.pz() + p_next.pz(),
             p.e() + p_next.e() ) );
      
      all_inv += unjet_input.back();
      charge(unjet_input.back(), 0);
    }
    // Jet reclustering

    // Reclustering radius
    double R=1.0;

    // Initialize jet algorithm
    // for invisible jets

    vector<PseudoJet>& reclustered_jets = selected_jets;
    
    JetDefinition jet_def(cambridge_algorithm, R);

    ClusterSequence cs(selected_jets, jet_def);
  
    reclustered_jets = sorted_by_pt(cs.inclusive_jets());


    // We have all the objects now; get and save variables
    
    double Mt = -1;
    double deta = 999;
    double drap = 999;
    double eta_max = 999;
    double dphi_jet1_met = 999;
    double dphi_jet2_met = 999;
    double dphi_jet1_jet2 = 999;
    double dphi_higgs_met = 999;
    double dphi_min_jet_met = 999;
    double m_mu_mu = 999;

    PseudoJet jet_higgs;

    double MjH = 999;
    double Mscaled = 999;
    double MMu2 = 999;
    double JetPt = 999;

    double nSub11 = 999;
    double nSub12 = 999;
    double nSub21 = 999;
    double nSub22 = 999;

    Candidate* cjetL = (Candidate*) vjets->At(0);
    Candidate* cjetSL = (Candidate*) vjets->At(1);

    nSub11 = cjetL->Tau[0];
    nSub12 = cjetL->Tau[1];

    nSub21 = cjetSL->Tau[0];
    nSub22 = cjetSL->Tau[1];    

    // Minimum deltaphi between gets and met

    dphi_min_jet_met = dphi(MEt, reclustered_jets);

    // Reconstruct higgs using two leading reclustered jets

    jet_higgs = reclustered_jets[0] + reclustered_jets[1];

  	Mt = MT(MEt, reclustered_jets[0], reclustered_jets[1]);

  	deta = fabs(reclustered_jets[0].eta()-reclustered_jets[1].eta());

  	drap = fabs(reclustered_jets[0].rap()-reclustered_jets[1].rap());

    dphi_jet1_met=dphi(selected_jets[0], MEt);
    dphi_jet2_met=dphi(selected_jets[1], MEt);
    dphi_jet1_jet2=dphi(selected_jets[0], selected_jets[1]);

    if (selected_muons.size() > 1){
      MMu2 = (selected_muons[0]+selected_muons[1]).m();
      //cout << selected_muons.size() << " " << MMu2 << endl;
    }

    string tag = "None";
    tag = "tchannel";

    dphi_higgs_met = dphi(jet_higgs, MEt);

    // cout << "jet1met: " << dphi_jet1_met << endl;
    // cout << "jet2met: " << dphi_jet2_met << endl << endl;

    // cout << std::min(dphi_jet1_met,dphi_jet2_met) << endl;

    if (num_reco_muons > 1)
      m_mu_mu = (selected_muons[0]+selected_muons[1]).m();
  
    // dphi_min_jet_met = std::min(dphi_jet1_met,dphi_jet2_met);


    // file_evt<<iEvent<<","<<tag<<","
    //   <<MEt.pt()<<","<<dphi_higgs_met<<","
    //   <<dphi_jet1_met<<","
    //   <<dphi_jet2_met<<","<<dphi_jet1_jet2<<","
    //   <<dphi_min_jet_met<<","      
    //   <<m_mu_mu<<","
    //   <<Mt<<","
    //   <<selected_jets.size()<<","
    //   <<selected_muons.size()<<","
    //   <<num_gen_muons<<","
    //   <<num_b_mesons<<","
    //   <<num_tracks_jet<<","
    //   <<num_diag<<","
    //   <<num_dark<<","
    //   <<HiggsPTCut->h.pt()<<","
    //   <<nSub11<<","
    //   <<nSub12<<","
    //   <<nSub21<<","
    //   <<nSub22<<endl;


    file_evt<<iEvent<<","
    << Mt <<","
    <<MEt.pt()<<","    
    << reclustered_jets[0].m()<<","
    << reclustered_jets[1].m()<<","
    << (reclustered_jets[0]+reclustered_jets[1]).m()<<","
    <<dphi_jet1_jet2<<","
    <<dphi_min_jet_met<<","
    <<reclustered_jets[0].eta()<<","
    <<reclustered_jets[1].eta()<<","
    <<fabs(reclustered_jets[0].eta()-reclustered_jets[1].eta())<<","
    <<num_tracks_jet<<","
    <<nSub11<<","
    <<nSub12<<","
    <<nSub21<<","
    <<nSub22<<endl;    

    //print met
    file_obj<<print_obj(iEvent, 1, "met", MEt);

    //print all the jets
    for(int i=0; i<selected_jets.size(); ++i)
      file_obj<<print_obj(iEvent, i, "jet", selected_jets[i]);

    //print all the muons
    for(int i=0; i<selected_muons.size(); ++i)
      file_obj<<print_obj(iEvent, i, "muon", selected_muons[i]);
        
    
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

  file_meta<<"nevt, npass, eff, cxn, cxn_err"<<endl;
  file_meta<<iTotal<<","<<iEvent<<","
     <<iEvent/double(iTotal)<<","
     <<pythia.info.sigmaGen()*1e9<<","
     <<pythia.info.sigmaErr()*1e9<<endl;
  
  file_meta<<"# cxn in pb"<<endl;

  file_evt<<"# total event = "<<iTotal<<endl;
  file_evt<<"# pass event = "<<iEvent<<endl;
  file_evt<<"# pt_min = "<<pt_min<<endl;
  file_evt<<"# trigger efficiency = "<<iEvent/double(iTotal)<<endl;
  file_evt<<"# cxn = "<<pythia.info.sigmaGen()*1e9<<endl;
  file_evt<<"# cxn_err = "<<pythia.info.sigmaErr()*1e9<<endl;
  

  //clean up
  delphes->FinishTask();
  delete delphes;
  delete config;
  
  if(cmdline.present("-v"))
    pythia.stat();
  // Done.
  file_evt.close();
  file_obj.close();
  file_meta.close();
  return 0;

}
