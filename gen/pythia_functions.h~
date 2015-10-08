#ifndef __tools_h
#define __tools_h

#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>
#include "matrix_inverse.h"
#include <complex>
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

double PI = 3.14159265358979323846;
TRandom3 ran;

typedef std::complex<double> cplx;


//invert matrix and return determinant
double Invert_Matrix(double* matrix, int rank)
{
  MatrixClass<double> M(rank, rank);
  
  //dump matrix entries to external class
  for(int i=0; i<rank; i++)
  for(int j=0; j<rank; j++)
    M.setvalue(i,j,matrix[i*rank+j]);
  
  //invert matrix
  double matrix_determinant=M.invert();
  assert(matrix_determinant!=0);
  
  bool state;
  //dump values back to matrix
  for(int i=0; i<rank; i++)
  for(int j=0; j<rank; j++)
    {
      M.getvalue(i,j,matrix[i*rank+j], state);
      //make sure there is no error
    }

  return matrix_determinant;
}

//function to apply matrix
double apply3(double m[3][3], double (&b)[3]){
  
  double a[3]={b[0], b[1], b[2]};

  for(int i=0; i<3; ++i){
    
    b[i] = 0;
    
    for(int j=0; j<3; ++j)
      b[i] += m[i][j]*a[j];
    
  }
  return 0;
  
}


double normalize(double (&v)[3]){
  double norm = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  norm = 1.0/sqrt(norm);
  v[0]*= norm; v[1]*= norm; v[2]*= norm;
  return 0;
}


bool zero(double v[3]){
  return v[0]==0 && v[1]==0 && v[2]==0;
}

double EigenValues3(double m[3][3], double (&eigen)[3],
        double (&vec)[3][3])
{
  double p1 = m[0][1]*m[0][1] +
    m[0][2]*m[0][2] + m[1][2]*m[1][2];
  
  //if diagonal
  if(p1 == 0){
    eigen[0] = m[0][0];
    eigen[1] = m[1][1];
    eigen[2] = m[2][2]; 
  }
  
  else {
    double q = (m[0][0]+m[1][1]+m[2][2])/3.0;
    double p2 = (m[0][0]-q)*(m[0][0]-q) + 
      (m[1][1]-q)*(m[1][1]-q) + 
      (m[2][2]-q)*(m[2][2]-q) + 2*p1;
    double p = sqrt(p2/6.0);

    double B[3][3];

    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++){
  B[i][j] = m[i][j];
  
  if(i==j)
    B[i][j] -= q;
  B[i][j]/= p;
      }
    double r = Invert_Matrix(&B[0][0], 3)/2.0;
        
    double phi=0;
    if(r<=-1)
      phi = PI/3.0;
    else if(r<1)
      phi = acos(r)/3.0;
    
    eigen[0] = q + 2*p*cos(phi);
    eigen[2] = q + 2*p*cos(phi + 2*PI/3.0);
    eigen[1] = 3*q - eigen[0] - eigen[2];
  }  

  // various projections
  double A[3][3][3];
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      
      for(int k=0; k<3; k++){
  A[k][i][j] = m[i][j];

  if(i==j)
    A[k][i][j] -= eigen[k];
      }
    }

  //first try a vector
  double v[3]={ran.Gaus(), ran.Gaus(), ran.Gaus()};

  int index=0;
  for(; index<2; index++){
    apply3(A[index], v);
    if(zero(v)){
      cout<<"ERROR: degenerate problems"<<endl;
    }
  }
  
  //found first eigenvector
  normalize(v);
  vec[2][0]=v[0];
  vec[2][1]=v[1];
  vec[2][2]=v[2];
  
  double v1[3]={ran.Gaus(), ran.Gaus(), ran.Gaus()};
  double dot = v1[0]*v[0] + v1[1]*v[1] + v1[2]*v[2];
  v1[0] -= dot*v[0];
  v1[1] -= dot*v[1];
  v1[2] -= dot*v[2];

  apply3(A[0], v1);
  if(zero(v1)){
    cout<<"ERROR: degenerate problems"<<endl;
  }

  normalize(v1);
  vec[1][0]=v1[0];
  vec[1][1]=v1[1];
  vec[1][2]=v1[2];

  //get the remaining one by cross product
  vec[0][0]=v[1]*v1[2] - v[2]*v1[1];
  vec[0][1]=v[2]*v1[0] - v[0]*v1[2];
  vec[0][2]=v[0]*v1[1] - v[1]*v1[0];
  return 0;
  
}

//return root for quadratic equation ax^2 + bx + c = 0
double quad_root(double a, double b, double c, bool plus=true){

  double dis= b*b - 4*a*c;

  if(dis<0){
    cout<<"WARNING: imaginary root, taking real part"<<endl;
    dis=0;
  }
    
  double rN1 = (-b + sqrt(dis))/(2*a);
  double rN2 = (-b - sqrt(dis))/(2*a);

  if (rN1 < rN2) swap(rN1, rN2);

  if (plus) return rN1;
  else return rN2;
}

// covariance matrix for a jet
class jet_cov{
  
 public:
  double mat[4][4];
  double matL[2][2];
  // Changed indives from 3 to 4
  double mat3[4][4];
  double matL3[2][2];
  double cov[4][4];
  double covL[2][2];
  double cov3[3][3];
  double eigen[3];
  double eigenv[3][3];
  double covL3[2][2];
  double sum2[4][4];
  double sum[4];
  double N;
  double store_N;
  bool b_compute_cov;

  vector<PseudoJet> constit;

  jet_cov(){
    
    for(int i=0; i<4; i++){
      sum[i]=0;
    }

    for(int i=0; i<3; i++){
      eigen[i]=0;
    }

    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++){
      sum2[i][j]=0;
      mat[i][j]=0;
      cov[i][j]=0;      
    }

    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      mat3[i][j]=0;
      cov3[i][j]=0;      
      eigenv[i][j]=0;
    }

    for(int i=0; i<2; i++)
    for(int j=0; j<2; j++){
      matL[i][j]=0;
      covL[i][j]=0;            
      matL3[i][j]=0;
      covL3[i][j]=0;            
    }

    N=0;
    b_compute_cov=false;
  }

  void scale(double factor=1.0){

    double inv_factor = 1.0/factor;
    
    for(int i=0; i<4; i++){
      sum[i]*= (factor);
    }

    for(int i=0; i<3; i++){
      eigen[i]*= (factor);
    }

    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++){
      sum2[i][j]*= (factor*factor);
      mat[i][j]*= (inv_factor*inv_factor);
      cov[i][j]*= (factor*factor);      
    }

    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      mat3[i][j]*= (inv_factor*inv_factor);
      cov3[i][j]*= (factor*factor);      
    }

    for(int i=0; i<2; i++)
    for(int j=0; j<2; j++){
      matL[i][j]*= (inv_factor*inv_factor);
      covL[i][j]*= (factor*factor);            
      matL3[i][j]*= (inv_factor*inv_factor);
      covL3[i][j]*= (factor*factor);            
    }
    
  }

  //jet_cov(Candidate* can);  

  void fill(double px, double py, double pz, double e){
    
    double input[4]={px, py, pz, e};

    for(int i=0; i<4; i++){
      sum[i] += input[i];
    }

    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++){
  sum2[i][j] += input[i]*input[j];
      }
    
    constit.push_back(PseudoJet(px, py, pz, e));

    N+=1.0;
  }
  
  void fill(Candidate* c1){
    const TLorentzVector& p1 = c1->Momentum;
    fill(p1.Px(), p1.Py(), p1.Pz(), p1.E());
  }
  
  void fill(const PseudoJet& c1){
    fill(c1.px(), c1.py(), c1.pz(), c1.e());
  }
  
  bool compute_cov() {

    if(N<2){
      cout<<"ERROR: N too small for covariance matrix"<<endl;
      return false;
    }
    
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++){
  mat[i][j] = sum2[i][j]/N - sum[i]*sum[j]/(N*N);
  mat[i][j] *= (N / (N-1));
  cov[i][j] = mat[i][j];

  if(i<3 && j<3){
    mat3[i][j]=mat[i][j];
    cov3[i][j]=cov[i][j];   
  }

      }

    double det = Invert_Matrix(&mat[0][0], 4);

    matL[0][0]=mat[3][3]; matL[0][1]=mat[3][2];
    matL[1][0]=mat[2][3]; matL[1][1]=mat[2][2];

    covL[0][0]=mat[3][3]; covL[0][1]=mat[3][2];
    covL[1][0]=mat[2][3]; covL[1][1]=mat[2][2];
    
    double detL = Invert_Matrix(&covL[0][0], 2);

    //now do the 3D ones
    double det3 = Invert_Matrix(&mat3[0][0], 3);

    matL3[0][0]=mat3[3][3]; matL3[0][1]=mat3[3][2];
    matL3[1][0]=mat3[2][3]; matL3[1][1]=mat3[2][2];

    covL3[0][0]=mat3[3][3]; covL3[0][1]=mat3[3][2];
    covL3[1][0]=mat3[2][3]; covL3[1][1]=mat3[2][2];

    //find the eigenvalues for the covariance matrix;
    EigenValues3(cov3, eigen, eigenv);

    
    /*
    cout<<"cov3: "<<endl;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
  cout<<cov3[i][j]<<", ";
      }
    cout<<endl;
    }

    cout<<"eigenvalues: "<<eigen[0]
  <<", "<<eigen[1]
  <<", "<<eigen[2]<<endl;
    
    cout<<"eigenvectors: "<<endl;

    for(int k=0; k<3; k++){
      cout<<"vec"<<k<<": ";
      for(int i=0; i<3; i++){
  cout<<eigenv[k][i]<<", ";
      }
      cout<<endl;
    }
    cout<<endl;
    */

    //take square root 
    for(int i=0; i<3; i++){
      eigen[i] = sqrt(eigen[i]);
    }
    
    double detL3 = Invert_Matrix(&covL3[0][0], 2);
    
    if(det==0 || detL==0 || det3==0 || detL3==0){
      cout<<"ERROR: covariance matrix not invertible!"<<endl;
      return false;
    }

    b_compute_cov=true;
    return true;
  }

  double mat_dot3(double px1, double py1, double pz1,
      double px2, double py2, double pz2){
    
    double input1[3]={px1, py1, pz1};
    double input2[3]={px2, py2, pz2};

    if(!b_compute_cov)
      compute_cov();

    //if there are errors
    if(!b_compute_cov)
      return 0;

    double result =0;
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++){
  result += input1[i]*mat3[i][j]*input2[j]; 
      }
    return result; 
    
  }

  double mat_dot3
    (const PseudoJet& jet1, const PseudoJet& jet2){
    return mat_dot3(jet1.px(), jet1.py(), jet1.pz(), 
        jet2.px(), jet2.py(), jet2.pz()); 
  }

  //compute the dot product
  double operator()
    (double px1, double py1, double pz1, double e1,
     double px2, double py2, double pz2, double e2, bool bcov=false){
    double input1[4]={px1, py1, pz1, e2};
    double input2[4]={px2, py2, pz2, e2};
    
    if(!b_compute_cov)
      compute_cov();
  
    //if there are errors
    if(!b_compute_cov)
      return 0;
    
    double result =0;
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++){
  if(bcov)
    result += input1[i]*cov[i][j]*input2[j];  
  else
    result += input1[i]*mat[i][j]*input2[j];  
      }
    return result; 
  }

  double operator()
    (const PseudoJet& jet1, const PseudoJet& jet2, 
     bool bcov=false){
    return (*this)(jet1.px(), jet1.py(), jet1.pz(), jet1.e(), 
       jet2.px(), jet2.py(), jet2.pz(), jet2.e(),
       bcov);
  }

  void fill_jet(Candidate* jet, bool flip=false){
    
    TObjArray* constit = jet->GetCandidates();
    //recluster first
    vector<PseudoJet> cons;
    for (int i=0; i<constit->GetEntriesFast(); ++i){
      
      Candidate* c1 = (Candidate*) constit->At(i);
      
      TLorentzVector p1 = c1->Momentum;

      if(flip){
  p1[0] *= -1.0;
  p1[1] *= -1.0;
      }

      cons.push_back(PseudoJet(p1.Px(), p1.Py(),
             p1.Pz(), p1.E()));
      fill(c1);
    }
    
    /*
    //cluster into little jets
    JetDefinition jet_def(kt_algorithm, 0.2);
    ClusterSequence cs(cons, jet_def);

    vector<PseudoJet> subjets = 
      sorted_by_pt(cs.inclusive_jets());
    
    //now fill
    for(int i=0; i<subjets.size(); ++i){
      fill(subjets[i]);
    }
    */

  }

  void app_matL(double a, double b,
    double& ra, double& rb) const{
    ra = matL[0][0]*a + matL[0][1]*b;
    rb = matL[1][0]*a + matL[1][1]*b;
  }

  void app_covL(double a, double b,
    double& ra, double& rb) const{
    ra = covL[0][0]*a + covL[0][1]*b;
    rb = covL[1][0]*a + covL[1][1]*b;
  }

  double scale_covL(double a, double b,
    double ra, double rb) const{
    return covL[0][0]*a*ra + covL[0][1]*a*rb +
      covL[1][0]*b*ra + covL[1][1]*b*rb;
  }

  double scale_matL(double a, double b,
    double ra, double rb) const{
    return matL[0][0]*a*ra + matL[0][1]*a*rb +
      matL[1][0]*b*ra + matL[1][1]*b*rb;
  }


  //return the energy component of MEt
  double get_E(const PseudoJet& jet) {
    
    //we need to evaluate the integral
    //~dy exp{-iyJ} Prod dx Prob(x) exp{iyx} sum(|x|)
    // = dy f^{N-1}(y) M(y) exp{-yJ} /
    // dy f^{N}(y) M(y) exp{-yJ}
    // 
    // where f(y) = int dx Prob(x)exp{iyx}
    // so first we compute f(y)
    //
    // f(y) = sum_constit exp{-y x}
    
    return sqrt(jet.modp2());
    
  }


  PseudoJet get_invjet(const PseudoJet& jet, 
           const PseudoJet& met){
    
    if(!b_compute_cov)
      compute_cov();

    PseudoJet p(jet), m(met.px(), met.py(), 0, 0);
    
    double eff_n=1.0;
    p*= (eff_n/N);
    scale(eff_n);

    double pz_mod = p.pz() + 
      (mat3[2][0]*p.px() + mat3[2][1]*p.py())/mat3[2][2];
    double STm = (mat3[2][0]*m.px() + mat3[2][1]*m.py()); 
    double a = mat3[2][2]*pz_mod*pz_mod - mat_dot3(p, p);
    //double b = -STm*pz_mod;
    double c = mat_dot3(m, m)- STm*STm/mat3[2][2];
    
    cout<<"a,c: "<<c<<", "<<a<<endl;
    double rN = sqrt(-c/a);
    //double rN = quad_root(a, b, c, true);
    //double rNn = quad_root(a, b, c, false);
    cout<<"roots found: "<<rN<<endl;//", "<<rNn<<endl;

    //make the mass zero for now
    double kz = rN*pz_mod - STm/mat3[2][2];


    //now compute the mass
    //compute integral of sqrt( N*(N-1)sigma^-1 ) x + k )^2 average over gaussian(x)

    double k0 = sqrt(m.modp2() + kz*kz);
    double avg=0;
    if(rN <1){
      cout<<"small rN, mass set to zero"<<endl;
    }
      
    else{
      
      double jet_e = 
  get_E(PseudoJet(m.px(), m.py(), kz, 0));
      
      
      double x_fact = sqrt((rN-1)*rN);
      int n_int = 500;
      
      PseudoJet j(m.px(), m.py(), kz, 0);
      for(int i=0; i<n_int; ++i){

  double j0 = ran.Gaus()*eigen[0]*x_fact;
  double j1 = ran.Gaus()*eigen[1]*x_fact;
  double j2 = ran.Gaus()*eigen[2]*x_fact;
  
  PseudoJet avg_j
    (j0*eigenv[0][0]+j1*eigenv[1][0]+j2*eigenv[2][0],
     j0*eigenv[0][1]+j1*eigenv[1][1]+j2*eigenv[2][1],
     j0*eigenv[0][2]+j1*eigenv[1][2]+j2*eigenv[2][2],
     0 );

  double p2=(avg_j+j).modp2();
  if(p2 < 0){
    cout<<"negative norm error!!"<<endl;
    p2=fabs(p2);
  }

  avg+= sqrt(p2);
      }
      avg /= n_int;
    }

    double mass= sqrt(avg*avg-k0*k0);
    cout<<"mass: "<<mass<<endl;
    
    //k0 = sqrt(k0*k0+mass*mass);
    cout<<"result: "
  <<m.px()<<", "<<m.py()<<", "<<kz<<", "<<k0<<endl;

    cout<<"jet: "
  <<jet.px()<<", "<<jet.py()<<", "
  <<jet.pz()<<", "<<jet.e()<<endl;

    store_N = rN;
    
    return PseudoJet(m.px(), m.py(), kz, k0);
  }

  PseudoJet get_invjet(const PseudoJet& met){
    
    return get_invjet
      (PseudoJet(sum[0], sum[1], sum[2], sum[3]), met);
    
  }


  
  PseudoJet get_invjet4(const PseudoJet& jet, 
      const PseudoJet& met){
    
    if(!b_compute_cov)
      compute_cov();
    
    PseudoJet p(jet), m(met.px(), met.py(), 0, 0);
    p*= (1.0/N);

    double S0p = mat[3][0]*p.px() + mat[3][1]*p.py();
    double Szp = mat[2][0]*p.px() + mat[2][1]*p.py();

    double S0m = mat[3][0]*m.px() + mat[3][1]*m.py();
    double Szm = mat[2][0]*m.px() + mat[2][1]*m.py();
    
    double p0_mod, pz_mod;
    app_covL(S0p, Szp, p0_mod, pz_mod);
    p0_mod = p.E() + p0_mod;
    pz_mod = p.pz() + pz_mod;

    double a = scale_matL(p0_mod, pz_mod, p0_mod, pz_mod)-
      (*this)(p, p);
    double c = (*this)(m, m) - 
      scale_covL(S0m, Szm, S0m, Szm);

    cout<<"a, c: "<<a<<", "<<c<<endl;
    double rN = sqrt(-c/a);
    cout<<"roots found: "<<rN<<endl;
    store_N = rN;

    double k0, kz;
    app_covL(S0m, Szm, k0, kz);
    k0 = rN*p0_mod - k0;
    kz = rN*pz_mod - kz;
    

    cout<<"jet: "
  <<p.px()<<", "<<p.py()<<", "<<p.pz()<<", "<<p.e()<<endl;
    cout<<"result: "
  <<m.px()<<", "<<m.py()<<", "<<kz<<", "<<k0<<endl;

    if(k0 < sqrt(m.modp2() + kz*kz)){
      k0 = sqrt(m.modp2() + kz*kz);
      cout<<"mod k0: "<<k0<<endl;
    }
    return PseudoJet(m.px(), m.py(), kz, k0);
  }
  
  double get_mjj(const PseudoJet& jet, 
     const PseudoJet& met){
    
    return (get_invjet4(jet, met)+jet).m();
  }


  double get_mjj(const PseudoJet& met){

    return get_mjj
      (PseudoJet(sum[0], sum[1], sum[2], sum[3]), met);

  }



};




double dot3(const PseudoJet& a, const PseudoJet& b){
return a.px() * b.px() + a.py() * b.py() + a.pz() * b.pz();
}

class wdouble{

public:
double val, w;
  
wdouble(double val, double w=1.0):
val(val), w(w){}
  
  bool operator<(const wdouble& other) const{
return val<other.val;
}
    
};


double quad_sum(double a, double b){
  return sqrt(a*a + b*b);
}

void swap(double& a, double& b){
  double tempa = a;
  double tempb = b;
  
  a=tempb;
  b=tempa;
}


double max(double n1, double n2){
  if(n1>n2) return n1;
  else return n2;
}

double min(double n1, double n2){
  if(n1>n2) return n2;
  else return n1;
}

double fix_phi(double phi){
  while(phi<0)
    phi+=2*PI;
  while(phi > 2*PI)
    phi-=2*PI;
  return phi;
}

double rand_angle(){
  return ran.Rndm()*2*PI;    
}


void reset_PtEtaPhiM(PseudoJet& input, 
         double pt,
         double eta, 
         double phi, 
         double mass=0){
  double px = pt*cos(phi);
  double py = pt*sin(phi);

  double pz = pt*sinh(eta);

  double e = sqrt(pt*pt + pz*pz + mass*mass);

  input.reset_momentum(px, py, pz, e);  
}

PseudoJet reverse(const PseudoJet& j){
  return PseudoJet(-j.px(), -j.py(), -j.pz(), j.E());
}

double invertm2(const PseudoJet& j, const PseudoJet& inv){

  //require the same boost factor
  double gamma = j.E()/sqrt(j.modp2());
  double invE = sqrt(inv.modp2())*gamma;
  return sqrt(invE*invE - inv.modp2());

  if(dot3(j, inv) < 0)
    return 0;
  
  //first find out m/pt ratio
  double mpt = j.m()/j.pt();
  
  //then get total mass
  double m = (j+inv).pt() *mpt;

  /*
  cout<<"original mass: "<<j.m()<<endl;
  cout<<"new mass: "<<m<<endl;
  cout<<"jet p: "<<j.px()<<","<<j.py()<<","<<j.pz()<<endl;
  cout<<"inv p: "<<inv.px()<<","<<inv.py()<<","<<inv.pz()<<endl;
  */

  //now invert
  double A = m*m - j.m2() + 2*dot3(j, inv);
  double E2 = j.E()*j.E();
  double B = inv.modp2();
  
  double result2 = -2*sqrt(A*E2 + B*E2 + E2*E2) + A + 2*E2;
  if (result2 < 0) 
    result2 = 0;

  PseudoJet test;
  reset_PtEtaPhiM(test, inv.pt(), 
      inv.eta(), inv.phi(), sqrt(result2));

  cout<<"inv: "<<test.px()
      <<","<<test.py()
      <<","<<test.pz()
      <<","<<test.e()<<endl;

  cout<<"jet: "<<j.px()
      <<","<<j.py()
      <<","<<j.pz()
      <<","<<j.e()<<endl;
  cout<<endl;

  return sqrt(result2);
}



double average(const vector<wdouble>& input){
  double result =0;
  double wsum = 0;
  for(int i=0; i<input.size(); ++i)
    {
      result += input[i].val * input[i].w;
      wsum += input[i].w;
    }

  result/= wsum;
  return result;
}

double percentile(const vector<wdouble>& input, double p=0.5, bool v= false){
  
  if(input.size()==0)
    return -1;

  if(input.size()==1)
    return input[0].val;


  double wsum = 0;

  for(int i=0; i<input.size(); ++i)
      wsum+=input[i].w;
    
  double normal=1.0/wsum;

  wsum = 0;
  double wsum_before = 0;
  for(int i=0; i<input.size(); ++i)
      {
  wsum_before = wsum;
  wsum+=input[i].w * normal;
  
  if(v)
    cout<<"at i:"<<i<<","<<std::scientific
        <<input[i].val<<","
        <<input[i].w<<","<<normal<<endl;
  
  if(wsum > p){

    if(i==0)
      return input[i].val;
    
    if(v){
      cout<<"returning: "
    <<input[i-1].val<<", "
    <<(input[i].w * normal)<<","
    <<(input[i].val-input[i-1].val)<<","
    <<(p-wsum_before)<<endl;
    }
      
    
    return input[i-1].val + (p-wsum_before)
      *(input[i].val-input[i-1].val)/(input[i].w * normal);
  }
  
      }
  return input[input.size()-1].val;
  
}


double stdev(const vector<wdouble>& input){
  
  double val2 = 0;
  double val = 0;
  double wsum = 0;
  double w2sum = 0;

  for(int i=0; i<input.size(); ++i)
    {
      val += input[i].val * input[i].w;
      val2 += input[i].val * input[i].val * input[i].w;
      wsum += input[i].w;
      w2sum += input[i].w * input[i].w;
    }
  
  return sqrt((val2*wsum - val*val)/(wsum*wsum - w2sum))
    * wsum/(val);
  
}

double average(const vector<double>& input){
  double result =0;
  for(int i=0; i<input.size(); ++i)
    result += input[i];

  result/= input.size();
  return result;
}

double median(const vector<double>& input, bool bsort=true){
  
  vector<double> temp(input);

  if(bsort)
    sort(temp.begin(), temp.end());

  return temp[input.size()/2];
}

double width(const vector<double>& input, bool bsort=true){
  
  vector<double> temp(input);

  if(bsort)
    sort(temp.begin(), temp.end());

  return (temp[3*input.size()/4] - temp[input.size()/4])
    /temp[input.size()/2];
}


double dphi(const PseudoJet& met, const vector<PseudoJet>& jets){
  
  double result = 999;
  for(int i=0; i<2 && i<jets.size(); ++i){
    double t_dphi = fabs(met.phi() - jets[i].phi());
    if(t_dphi > PI)
      t_dphi = 2*PI - t_dphi;
    result = result < t_dphi ? result : t_dphi;
  }
  return result;
}

double dphi(const PseudoJet& met, const PseudoJet& jet){
  
  double t_dphi = fabs(met.phi() - jet.phi());
    if(t_dphi > PI)
      t_dphi = 2*PI - t_dphi;
    
  return t_dphi;
}

double dphi(double phi1, double phi2){
  
  double t_dphi = fabs(phi1 - phi2);
    if(t_dphi > PI)
      t_dphi = 2*PI - t_dphi;
    
  return t_dphi;
}

//convert x,y to phi
double convert_phi(double init_phi=0, double width=1.0){
  
  /*
  double x = fabs(gauss(gen));
  double y = gauss(gen);
  
  double angle = init_phi+atan2(x,y);
  */

  double angle = init_phi + ran.Gaus()*width;

  while(angle > 2*PI)
    angle -= 2*PI;

  while(angle < 0)
    angle += 2*PI;

  return angle;
   
}


double HT(const vector<PseudoJet>& jets){

  double result = 0;
  for(int i=0; i<jets.size(); ++i){
    result += jets[i].mt();
  }

  return result;
}

double MJ(const vector<PseudoJet>& jets){

  double result = 0;
  for(int i=0; i<jets.size(); ++i){
    result += jets[i].m();
  }

  return result;
}

double invert(double (& input)[2][2]){

  double a = input[0][0];
  double b = input[0][1];
  double c = input[1][0];
  double d = input[1][1];


  double det = 
    input[0][0]*input[1][1] -
    input[0][1]*input[1][0];

  if(det==0) return 0;

  double inv = 1.0/det;
  swap(input[0][0], input[1][1]);
  
  input[0][0]*=inv;
  input[0][1]*=-inv;
  input[1][0]*=-inv;
  input[1][1]*=inv;

  /*
  cout<<input[0][0]*a + input[0][1]*c <<", "
      <<input[0][0]*b + input[0][1]*d <<", "<<endl
      <<input[1][0]*a + input[1][1]*c <<", "
      <<input[1][0]*b + input[1][1]*d <<", "<<endl;
  */

  return det;
}


// timer

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



//functions for jets
void charge(PseudoJet& jet, double c){
  jet.set_user_index((int) (c*3.0));
}

double charge(const PseudoJet& jet){
  return jet.user_index()/3.0;
}



double charge(const vector<PseudoJet>& constit){
  
  double jet_charge = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit.size(); ++i){
    double cpt = constit[i].pt();
    jet_charge += charge(constit[i])* cpt;
    jet_sumpt += cpt;
  }

  jet_charge /= jet_sumpt;

  return jet_charge;
}

double charge(const TObjArray* constit){
  
  double jet_charge = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){

    Candidate* c1 = (Candidate*) constit->At(i);
    const TLorentzVector& p1 = c1->Momentum;

    double cpt = p1.Pt();
    jet_charge += c1->Charge * cpt;
    jet_sumpt += cpt;
  }

  jet_charge /= jet_sumpt;

  return jet_charge;
}

//compute met vector
PseudoJet met(const vector<PseudoJet>& input){

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  
  for (int i=0; i<input.size(); ++i){
    Px -= input[i].px();
    Py -= input[i].py();
    Pz -= input[i].pz();
  }

  return PseudoJet(Px, Py, Pz, sqrt(Px*Px+Py*Py));

}

double charge_ratio(const vector<PseudoJet>& constit){
  
  double N_charge=0;

  for (int i=0; i<constit.size(); ++i){
    if(charge(constit[i]) != 0)
      N_charge += 1.0;
  }  

  if (constit.size() > 0)
    return N_charge / constit.size();

  return 0;
}


double charge_ratio(TObjArray* constit){
  
  double N_charge=0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){
    if(((Candidate*) constit->At(i))->Charge != 0)
      N_charge += 1.0;
  }  

  if (constit->GetEntriesFast() > 0)
    return N_charge / constit->GetEntriesFast();

  return 0;
}

double dipole(const vector<PseudoJet>& constit){
  
  double jet_dipole = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit.size(); ++i){
    for (int j=i; j<constit.size(); ++j){
      
      double cpt1 = constit[i].pt();
      double cpt2 = constit[j].pt();

      jet_dipole += charge(constit[i])*charge(constit[j])*
  cpt1*cpt2;

      jet_sumpt += cpt1*cpt2;

    }}

  if(jet_sumpt > 0)
    jet_dipole /= jet_sumpt;

  return jet_dipole;
}

double dipole(const TObjArray* constit){
  
  double jet_dipole = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){
    for (int j=i; j<constit->GetEntriesFast(); ++j){

      Candidate* c1 = (Candidate*) constit->At(i);
      const TLorentzVector& p1 = c1->Momentum;

      Candidate* c2 = (Candidate*) constit->At(j);
      const TLorentzVector& p2 = c2->Momentum;      

      double cpt1 = p1.Pt();
      double cpt2 = p2.Pt();

      jet_dipole += c1->Charge * c2->Charge*
  cpt1*cpt2;

      jet_sumpt += cpt1*cpt2;

    }}

  if(jet_sumpt > 0)
    jet_dipole /= jet_sumpt;

  return jet_dipole;
}




double tripole(const vector<PseudoJet>& constit){
  
  double jet_tripole = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit.size(); ++i){
    for (int j=i; j<constit.size(); ++j){
      for (int k=j; k<constit.size(); ++k){
      
      double cpt1 = constit[i].pt();
      double cpt2 = constit[j].pt();
      double cpt3 = constit[k].pt();

      jet_tripole += 
  charge(constit[i])*
  charge(constit[j])*
  charge(constit[k])*
  cpt1*cpt2*cpt3;

      jet_sumpt += cpt1*cpt2*cpt3;

      }}}

  if(jet_sumpt > 0)
    jet_tripole /= jet_sumpt;

  return jet_tripole;
}


double tripole(const TObjArray* constit){
  
  double jet_tripole = 0;
  double jet_sumpt = 0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){
    for (int j=i; j<constit->GetEntriesFast(); ++j){
      for (int k=j; k<constit->GetEntriesFast(); ++k){
      
  Candidate* c1 = (Candidate*) constit->At(i);
  const TLorentzVector& p1 = c1->Momentum;
  
  Candidate* c2 = (Candidate*) constit->At(j);
  const TLorentzVector& p2 = c2->Momentum;      

  Candidate* c3 = (Candidate*) constit->At(k);
  const TLorentzVector& p3 = c3->Momentum;      
  
  double cpt1 = p1.Pt();
  double cpt2 = p2.Pt();
  double cpt3 = p3.Pt();
  
  jet_tripole += 
    c1->Charge*
    c2->Charge*
    c3->Charge*
    cpt1*cpt2*cpt3;
  
  jet_sumpt += cpt1*cpt2*cpt3;
  
      }}}

  if(jet_sumpt > 0)
    jet_tripole /= jet_sumpt;

  return jet_tripole;
}

double corr1(const TObjArray* constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;
  double Rf = pow(R, beta);
  
  for (int i=0; i<constit->GetEntriesFast(); ++i){
      
    Candidate* c1 = (Candidate*) constit->At(i);
    const TLorentzVector& p1 = c1->Momentum;
    double cpt1 = p1.Pt();
    
    jet_sumpt += cpt1;

  }

  return jet_sumpt*Rf;
}

double corr2(const TObjArray* constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){
    for (int j=i; j<constit->GetEntriesFast(); ++j){
      
      Candidate* c1 = (Candidate*) constit->At(i);
      const TLorentzVector& p1 = c1->Momentum;

      Candidate* c2 = (Candidate*) constit->At(j);
      const TLorentzVector& p2 = c2->Momentum;      

      double cpt1 = p1.Pt();
      double cpt2 = p2.Pt();
      double dR12 = p1.DeltaR(p2);

      jet_sumpt += cpt1*cpt2*pow(dR12, beta);

    }}

  return jet_sumpt;
}

double corr3(const TObjArray* constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;

  for (int i=0; i<constit->GetEntriesFast(); ++i){
    for (int j=i; j<constit->GetEntriesFast(); ++j){
      for (int k=j; k<constit->GetEntriesFast(); ++k){
      
  Candidate* c1 = (Candidate*) constit->At(i);
  const TLorentzVector& p1 = c1->Momentum;

  Candidate* c2 = (Candidate*) constit->At(j);
  const TLorentzVector& p2 = c2->Momentum;      

  Candidate* c3 = (Candidate*) constit->At(k);
  const TLorentzVector& p3 = c3->Momentum;      
  
  double cpt1 = p1.Pt();
  double cpt2 = p2.Pt();
  double cpt3 = p3.Pt();
  
  double dR12 = p1.DeltaR(p2);
  double dR23 = p2.DeltaR(p3);
  double dR31 = p3.DeltaR(p1);

  double dR = dR12*dR23*dR31;
  
  jet_sumpt += cpt1*cpt2*cpt3* pow(dR, beta);

      }}}

  return jet_sumpt;
}


double corr1(const vector<PseudoJet>& constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;
  double Rf = pow(R, beta);

  for (int i=0; i<constit.size(); ++i){
      
      double cpt1 = constit[i].pt();

      jet_sumpt += cpt1;

    }

  return jet_sumpt*Rf;
}

double corr2(const vector<PseudoJet>& constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;

  for (int i=0; i<constit.size(); ++i){
    for (int j=i; j<constit.size(); ++j){
      
      double cpt1 = constit[i].pt();
      double cpt2 = constit[j].pt();
      double dR12 = constit[i].delta_R(constit[j]);

      jet_sumpt += cpt1*cpt2*pow(dR12, beta);

    }}

  return jet_sumpt;
}

double corr3(const vector<PseudoJet>& constit, 
       double beta=1.0, double R=0.5){
  
  double jet_sumpt = 0;

  for (int i=0; i<constit.size(); ++i){
    for (int j=i; j<constit.size(); ++j){
      for (int k=j; k<constit.size(); ++k){
      
  double cpt1 = constit[i].pt();
  double cpt2 = constit[j].pt();
  double cpt3 = constit[k].pt();
  
  double dR12 = constit[i].delta_R(constit[j]);
  double dR23 = constit[j].delta_R(constit[k]);
  double dR31 = constit[k].delta_R(constit[i]);

  double dR = dR12*dR23*dR31;
  
  jet_sumpt += cpt1*cpt2*cpt3* pow(dR, beta);

      }}}

  return jet_sumpt;
}

string to_st(double num){
  stringstream sout;
  sout<<num;
  return sout.str();
}

string int_st(int num){
  stringstream sout;
  sout<<num;
  return sout.str();
}

// function to add string to numbers
const char* add_st(string st, double num){
  stringstream sout;
  sout<<st<<num;
  return sout.str().c_str();
}

// function to add string to numbers
const char* add_fb(double id, string st, double num){
  stringstream sout;
  sout<<id<<st<<num;
  return sout.str().c_str();
}


// Function to print object information

string print_obj(int evt, int n, string type, const PseudoJet& v)
{
  
  stringstream sout;

  if(type=="met")
    sout << evt << "," << type << "," << 1 << "," << v.pt() << "," << 0 << "," << 0 << "," << v.phi() << endl;  
  else
    sout << evt << "," << type << "," << n << "," << v.pt() << "," << v.m() << "," << v.eta() << "," << v.phi() << endl;  
  return sout.str();
}


// MT with two jets

double MT(const PseudoJet& MEt, const PseudoJet& jet1_, const PseudoJet& jet2_)
{
  PseudoJet jj = jet1_ + jet2_;
  double ET = jj.mperp();
  double MT = sqrt(jj.m2() + 2*(MEt.pt()*jj.pt() - dot3(jj, MEt)));
  return MT;
}

// MT with one jet

double MT(const PseudoJet& MEt, const PseudoJet& jet)
{
  double ET = jet.mperp();
  double MT = sqrt(jet.m2() + 2*(MEt.pt()*jet.pt() - dot3(jet, MEt)));

  return MT;
}

double Mcol(const PseudoJet& MEt, const PseudoJet& jet1_, const PseudoJet& jet2_)
{
  
  PseudoJet jet1(jet1_), jet2(jet2_);

  reset_PtEtaPhiM(jet1, 1./cosh(jet1.eta()), jet1.eta(), jet1.phi());
  reset_PtEtaPhiM(jet2, 1./cosh(jet2.eta()), jet2.eta(), jet2.phi());
  
  double mat[2][2] =
    { {jet1.px(), jet2.px()} ,
      {jet1.py(), jet2.py()} };
  
  // Matrix must be invertible
  if(invert(mat) == 0)
    {
      cout<<"ERROR: not invertible!"<<endl;
      return -1;
    }
  

  double c1 = 
    mat[0][0]*MEt.px() + mat[0][1]*MEt.py();
  
  double c2 = 
    mat[1][0]*MEt.px() + mat[1][1]*MEt.py();

  jet1*= c1;
  jet2*= c2;
  
  reset_PtEtaPhiM(jet1, jet1.pt(), jet1.eta(), jet1.phi());
  reset_PtEtaPhiM(jet2, jet2.pt(), jet2.eta(), jet2.phi());
  
  return (jet1+jet2+ jet1_ + jet2_).m();  
}

string add_strings(string st1, double num1)
{
  string added = st1 + to_st(num1);
  return added;
}

#endif

