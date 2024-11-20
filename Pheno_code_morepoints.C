#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"
#include <string>
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include <stdlib.h> 
//#include "T1F.h"                                                                                                                                                                     
#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standart Library of Complex Numbers */
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "RConfigure.h"

#ifdef R__HAS_MATHMORE
#include "Math/MultiRootFinder.h"
#endif
#include "Math/WrappedMultiTF1.h"
#include "TF2.h"
#include "TError.h"
 //#ifndef ROOT_Math_GSLMultiRootSolver
 #define ROOT_Math_GSLMultiRootSolver
 
// #include "gsl/gsl_vector.h"
// #include "gsl/gsl_matrix.h"
// #include "gsl/gsl_multiroots.h"
// #include "gsl/gsl_blas.h"
 //#include "GSLMultiRootFunctionWrapper.h"
 #include <Math/RootFinderAlgorithms.h>
 #include "Math/IFunction.h"
 #include "Math/Error.h"
 
 #include <vector>
 #include <string>
 #include <cassert>
    #include "TRandom.h"
    
// #include <GSLMultiRootSolver.h>
 using namespace std;
void SetPlotStyle();

//Hgg and Hgz inputs from Gilbert
float GF = 0.0000116639; //Fermi constant
float MW = 80. ; 
float alpha = 1/137;
float Qplus = 1;
float Qplusplus = 2;
float sw = TMath::Sqrt(0.23);
float cw = TMath::Sqrt(1 - sw*sw);

double coefgammagamma(double MH);
TComplex Fgg0(double rta);
TComplex GgZ0(double rta);
	  
TComplex FF1(double rta);
TComplex F0(double rta);
TComplex F12(double rta) ;
TComplex A0(double tau, double lam, bool condition);
TComplex A12(double tau, double lam, bool condition) ;
TComplex A1(double tau, double lam, bool condition) ;
TComplex I1(double tau, double lam) ;
TComplex I2(double tau, double lam) ;
TComplex A1prime(double tau, double lam) ;

void DrawAtlas(float xshift=0, float yshift=0.);

bool unitarity2corr(double kappa, double l0, double l1, double l2, double l3, double l4) ;
bool BFBold(double l0, double l1, double l2, double l3, double l4) ;
bool RKcombinedBFBunitold(double kappa, double l0, double l1, double l2, double l3, double l4);
bool BFBnew(double l0, double l1, double l2, double l3, double l4);
bool RKcombinedBFBunit(double kappa, std::vector<double> lambdas);

void  RKComputeMasses2_Hpp(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;
void  muRKComputeMasses2_Hpp(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;
void  sinamuRKComputeMasses2_Hpp(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;
void  sinamuRKComputeMasses2_Hpp2(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;
void  sinamuRKComputeMasses2_Hpp3(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;

void  test_RKComputeMasses2_Hpp(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) ;

void  sinamumh0(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp) ;


void Solve(vector<double>& roots, int nest, TF1* f, double s, double e);
void exampleMultiRoot(const char * algo , int printlevel ) ;
using namespace ROOT::Math;

float FixedDoublyChargedHiggsMass = 220.;

int option_fixed =2;
int option_vary_lambda =1;
//vary lambda [-1,1] option_vary_lambda = 1;
//vary lambda [-2,2] option_vary_lambda = 2;
//vary lambda [-3,3] option_vary_lambda = 3;
//vary lambda [-4,4] option_vary_lambda = 4;
//vary lambda [-5,5] option_vary_lambda = 5;

//fixing only FixedDoublyChargedHiggsMass and h0=125 option_fixed = 1;
//fixing FixedDoublyChargedHiggsMass and FixedSingeChargedHiggsMass and h0=125 option_fixed = 2;
//fixing FixedDoublyChargedHiggsMass and A0 and h0=125 option_fixed = 3;
//fixing h0=125 and A0 option_fixed = 4;




void Pheno_code(){ 
SetPlotStyle();




std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!---------------------------------SearchHpp ----------------------------------------!!!!!!!!!!!!!!!"<<std::endl;

option_fixed = 2;

//fixing only FixedDoublyChargedHiggsMass and h0=125 option_fixed = 1;
//fixing FixedDoublyChargedHiggsMass and FixedSingeChargedHiggsMass and h0=125 option_fixed = 2;
//fixing FixedDoublyChargedHiggsMass and A0 and h0=125 option_fixed = 3;
//fixing h0=125 and A0 option_fixed = 4;


/*

///////////////////--------------sina0.99------------//////////////
//RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.99,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);



///////////////////--------------sina0.1------------//////////////
//RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);

RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);



///////////////////--------------sina0.01------------//////////////
//RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);


///////////////////--------------sina0.001------------//////////////
//RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);


///////////////////--------------sina0.0001------------//////////////
//RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);
*/

///////////////////--------------sina0.00001------------//////////////

//RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 50. ,1., 50., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 100. ,1., 100., 1., option_fixed, option_vary_lambda,0);

///what Adam needs sinalpha 8*1e-4////////////////////
RKComputeMasses2_Hpp(0.0015,0.1,246,0.1, 0.2, 0.3, 0.4, 220. ,1., 220., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);



//RKComputeMasses2_Hpp(0.0001,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.001,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.01,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
//RKComputeMasses2_Hpp(0.1,0.1,246,0.1, 0.2, 0.3, 0.4, 200. ,1., 200., 1., option_fixed, option_vary_lambda,0);
/*
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 300. ,1., 300., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 400. ,1., 400., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 500. ,1., 500., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 600. ,1., 600., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 700. ,1., 700., 1., option_fixed, option_vary_lambda,0);
RKComputeMasses2_Hpp(0.00001,0.1,246,0.1, 0.2, 0.3, 0.4, 800. ,1., 800., 1., option_fixed, option_vary_lambda,0);
*/

std::cout<<"!!!!!!!!!!!!!-------------------------------- ----------------------------------------!!!!!!!!!!!!!!!!"<<std::endl;

}

double coefgammagamma(double MH) {
double formula = GF*(1/137)* (1/137)*pow(MH,3)/128/TMath::Sqrt(2)/pow(TMath::Pi(),3);
return formula;

}



//unitarity2corr is the new unitarity constraints which have to be summed with the BFB constraints (BFBold should be for only lambda3 < 0 while BFBnew is for every lambda3)
//RKcombinedBFBunit combines unitarity and BFBnew. RKcombinedBFBunitold should be identical to RKcombinedBFBunitold (only for cross-checks and archiving purposes)


bool unitarity2corr(double kappa, double l0, double l1, double l2, double l3, double l4) {
 bool  tmpbool =   l1 <= kappa*TMath::Pi() && l1 >= -(kappa*TMath::Pi()) && l1 + l4 <= kappa*TMath::Pi() &&  l1 + l4 >= -(kappa*TMath::Pi()) && l1 + (3*l4)/2 <= kappa*TMath::Pi() &&     l1 + (3*l4)/2 >= -(kappa*TMath::Pi()) && 10/2 <= kappa*TMath::Pi() && l0/2 >= -(kappa*TMath::Pi()) &&     2*l2 <= kappa*TMath::Pi() &&  2*l2 >= -(kappa*TMath::Pi()) && 2*(l2 + l3) <= kappa*TMath::Pi() &&  2*(l2 + l3) >= -(kappa*TMath::Pi()) &&    (l0 + 4*l2 + 8*l3 +   sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/4 <=  kappa*TMath::Pi() &&    (l0 + 4*l2 + 8*l3 +  sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/  4 >=     -(kappa*TMath::Pi()) && (l0 + 4*l2 + 8*l3 -         sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/4 <=  kappa*TMath::Pi() &&    (l0 + 4*l2 + 8*l3 -  sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/  4 >=     -(kappa*TMath::Pi()) && (3*l0 + 16*l2 + 12*l3 + sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/4 <=  kappa*TMath::Pi() &&    (3*l0 + 16*l2 + 12*l3 +  sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/ 4 >= -(kappa*TMath::Pi()) &&    (3*l0 + 16*l2 + 12*l3 - sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) +  24*(2*l1 + l4)*(2*l1 + l4)))/ 4 <= kappa*TMath::Pi() &&  (3*l0 + 16*l2 + 12*l3 -  sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/ 4 >= -(kappa*TMath::Pi()) && l1 - l4/2 <= kappa*TMath::Pi() &&    l1 - l4/2 >= -(kappa*TMath::Pi()) &&  l2 - l3 <= kappa*TMath::Pi() && 2 * l2 - l3 >= -kappa*TMath::Pi();
 
return tmpbool;}


bool BFBold(double l0, double l1, double l2, double l3, double l4) {
bool tmpbool = l0 >= 0&& l2 + l3 >= 0&& l2 + l3/2 >= 0 &&  l1 + sqrt(l0*(l2 + l3)) >= 0&& l1 + sqrt(l0*(l2 + l3/2)) >= 0&&  l1 + l4 + sqrt(l0*(l2 + l3)) >= 0&&  l1 + l4 + sqrt(l0*(l2 + l3/2)) >= 0;
return tmpbool;}

bool RKcombinedBFBunitold(double kappa, double l0, double l1, double l2, double l3, double l4) {
bool tmpbool =  unitarity2corr(kappa, l0, l1, l2, l3, l4) &&  BFBold(l0, l1, l2, l3, l4); 
return tmpbool;}

bool BFBnew(double l0, double l1, double l2, double l3, double l4) {
bool tmpbool0= l0 >= 0&& l2 + l3 >= 0&& l2 + l3/2 >= 0;
bool tmpbool1= l1 + sqrt(l0*(l2 + l3)) >= 0&& sqrt(l0)*l3 <= sqrt((l2 + l3)*l4*l4)&& l1 + l4 + sqrt(l0*(l2 + l3)) >= 0;
bool tmpbool2= sqrt(l0)*l3 >=  sqrt((l2 + l3)*l4*l4)&& -sqrt(l0*(l2 + l3/2) * (1 - l4*l4/2/l0/l3)) <=  l1 + l4/2;
bool tmpbool3 = tmpbool1 || tmpbool2;
return tmpbool0 + tmpbool3;
}

bool RKcombinedBFBunit(double kappa, std::vector<double> lambdas) {  
  double l0 = lambdas[0]; double l1 = lambdas[1]; double l2 = lambdas[2]; 
  double l3 = lambdas[3]; double l4 = lambdas[4]; 

bool tmpbool  =  unitarity2corr(kappa, l0, l1, l2, l3, l4) &&  BFBnew(l0, l1, l2, l3, l4);
return tmpbool;
}
  

double mHpp = 220.;


void  RKComputeMasses2_Hpp(double sina, double vt, double vd, double l1, double l2,double  l3, double l4, double mHpp, double dmHpp, double mHp, double dmHp, int option_fixed, int option_vary_lambda, int ratio_value) {
 ostringstream mysstream; 
      ofstream eff_file_e; 
      ofstream eff_file_condition_lam; 
   string mHpp_char = to_string(mHpp);
   string sina_char = to_string(sina);
   string vt_char = to_string(vt);
   string vd_char = to_string(vd);

std::cout<<"mHpp_char = "<<mHpp_char<<std::endl;
std::string message = std::string("pp_HppHmm_H++mass_200_800") +std::string("_sina_") + sina_char+std::string("_vt_")+vt_char+std::string("_vd_")+vd_char+std::string("_fixedH++_fixedh0_fixedH+_generateHpHm.dat");
std::string message2 = std::string("lam1_H++mass_200_800") +std::string("_sina_") + sina_char+std::string("_vt_")+vt_char+std::string("_vd_")+vd_char+std::string("_fixedH++_fixedh0_fixedH+_generateHpHm.dat");

 string line = mysstream.str();
  eff_file_e.open (message,ios::app);
  eff_file_condition_lam.open (message2,ios::app);

double mu = 0.00000001;
double l0 = 0.55;
 l1=0;
double FixedDoublyChargedHiggsMass  = mHpp;


cout << setprecision(10);
eff_file_e<<"import model Type-II_seesaw_mix_UFO"<<std::endl;
eff_file_e<<"define EXCL = H ha h"<<std::endl;
eff_file_e<<"define Hpp = H++ H-- "<<std::endl;
eff_file_e<<"define Hp = H+ H- "<<std::endl;
eff_file_e<<"generate p p > Hp Hp /EXCL"<<std::endl;
eff_file_e<<"output pptoHpHm_AnaPoints"<<std::endl;
eff_file_e<<"set LHC 13"<<std::endl;
eff_file_e<<"set bwcutoff 1000"<<std::endl;
eff_file_e<<"set cut_decays F"<<std::endl;
eff_file_e<<"set pdfwgt F"<<std::endl;
eff_file_e<<"set use_syst F"<<std::endl;
eff_file_e<<"set fixed_ren_scale F"<<std::endl;
eff_file_e<<"set fixed_fac_scale F"<<std::endl;


const char * algo = 0;
int printlevel = 1;

bool tmptmHpp = 0; bool tmpbool = false; 

double mHdpm2=mHpp*mHpp;
double mA2 = 0;
double mHpm2 = mHp*mHp;
double mh02  = 0;
double mH02  = 0;

int contor=0;
int contor_found=0;
int contor_found1=0;
int contor_found2=0;
int contor_found3=0;
int contor_found4=0;
int mp_contor_found=0;
int mu_contor_found=0;

double A = 1;
double B = 1;
double CC = 1;
double sA = 1; 
double cA = 1; 

  
    vector<double> sina_vec_mh02;
  vector<double> sina_vec1_mh02;
  vector<double> sina_vec2_mh02;
  vector<double> sina_vec3_mh02;
  vector<double> sina_vec4_mh02;

   vector<double> sina_vec_sina;
  vector<double> sina_vec_vt;
  vector<double> sina_vec_vd;
  vector<double> sina_vec_l0;
  vector<double> sina_vec_l1;
  vector<double> sina_vec_l2;
  vector<double> sina_vec_l3;
  vector<double> sina_vec_l4;
  vector<double> sina_vec_mu;
  vector<double> sina_vec_mHdpm;
  vector<double> sina_vec_mA2;
  vector<double> sina_vec_mHpm2;
  vector<double> sina_vec_mH02;
  
     vector<double> sina_vec1_sina;
  vector<double> sina_vec1_vt;
  vector<double> sina_vec1_vd;
  vector<double> sina_vec1_l0;
  vector<double> sina_vec1_l1;
  vector<double> sina_vec1_l2;
  vector<double> sina_vec1_l3;
  vector<double> sina_vec1_l4;
  vector<double> sina_vec1_mu;
  vector<double> sina_vec1_mHdpm;
  vector<double> sina_vec1_mA2;
  vector<double> sina_vec1_mHpm2;
  vector<double> sina_vec1_mH02;
  
     vector<double> sina_vec2_sina;
  vector<double> sina_vec2_vt;
  vector<double> sina_vec2_vd;
  vector<double> sina_vec2_l0;
  vector<double> sina_vec2_l1;
  vector<double> sina_vec2_l2;
  vector<double> sina_vec2_l3;
  vector<double> sina_vec2_l4;
  vector<double> sina_vec2_mu;
  vector<double> sina_vec2_mHdpm;
  vector<double> sina_vec2_mA2;
  vector<double> sina_vec2_mHpm2;
  vector<double> sina_vec2_mH02;
     vector<double> sina_vec3_sina;
  vector<double> sina_vec3_vt;
  vector<double> sina_vec3_vd;
  vector<double> sina_vec3_l0;
  vector<double> sina_vec3_l1;
  vector<double> sina_vec3_l2;
  vector<double> sina_vec3_l3;
  vector<double> sina_vec3_l4;
  vector<double> sina_vec3_mu;
  vector<double> sina_vec3_mHdpm;
  vector<double> sina_vec3_mA2;
  vector<double> sina_vec3_mHpm2;
  vector<double> sina_vec3_mH02;
  
     vector<double> sina_vec4_sina;
  vector<double> sina_vec4_vt;
  vector<double> sina_vec4_vd;
  vector<double> sina_vec4_l0;
  vector<double> sina_vec4_l1;
  vector<double> sina_vec4_l2;
  vector<double> sina_vec4_l3;
  vector<double> sina_vec4_l4;
  vector<double> sina_vec4_mu;
  vector<double> sina_vec4_mHdpm;
  vector<double> sina_vec4_mA2;
  vector<double> sina_vec4_mHpm2;
  vector<double> sina_vec4_mH02;
  
  vector<double> mp_vec_sina;
  vector<double> mp_vec_vt;
  vector<double> mp_vec_vd;
  vector<double> mp_vec_l0;
  vector<double> mp_vec_l1;
  vector<double> mp_vec_l2;
  vector<double> mp_vec_l3;
  vector<double> mp_vec_l4;
  vector<double> mp_vec_mu;
  vector<double> mp_vec_mHdpm;
  vector<double> mp_vec_mA2;
  vector<double> mp_vec_mHpm2;
  vector<double> mp_vec_mH02;
  vector<double> mp_vec_mh02;
  
    vector<double> mu_vec_sina;
  vector<double> mu_vec_vt;
  vector<double> mu_vec_vd;
  vector<double> mu_vec_l0;
  vector<double> mu_vec_l1;
  vector<double> mu_vec_l2;
  vector<double> mu_vec_l3;
  vector<double> mu_vec_l4;
  vector<double> mu_vec_mu;
  vector<double> mu_vec_mHdpm;
  vector<double> mu_vec_mA2;
  vector<double> mu_vec_mHpm2;
  vector<double> mu_vec_mH02;
  vector<double> mu_vec_mh02;
  


sA=sina;
cA = sqrt(1-sina*sina);
mh02=125.*125.;





  TH2F*  mp_l1l4 = new TH2F("l14","l14",100,-2,2,100,-2,2);
  TH2F*  mp_l1l2 = new TH2F("l12","l12",100,-2,2,100,-2,2);
  TH2F*  mp_l1l3 = new TH2F("l13","l13",100,-2,2,100,-2,2);
  TH2F*  mp_l2l3 = new TH2F("l23","l23",100,-2,2,100,-2,2);
  TH2F*  mp_l2l4 = new TH2F("l24","l24",100,-2,2,100,-2,2);
  TH2F*  mp_l3l4 = new TH2F("l34","l34",100,-2,2,100,-2,2);
  
  TH2F*  mp_l1mu = new TH2F("mp_l1mu","mp_l1mu",100,0,1000,100,0,100);
  TH2F*  mp_l2mu = new TH2F("mp_l2mu","mp_l2mu",100,0,1000,100,0,100);
  TH2F*  mp_l3mu = new TH2F("mp_l3mu","mp_l3mu",100,0,1000,100,0,100);
  TH2F*  mp_l4mu = new TH2F("mp_l4mu","mp_l4mu",100,0,1000,100,0,100);
  TH2F*  mp_l0mu = new TH2F("mp_l0mu","mp_l0mu",100,0,100,100,0,100);
  
  
  TH2F*  mp_l1sina = new TH2F("mp_l1sina","mp_l1sina",100,-2,2,100,0,50);
  TH2F*  mp_l2sina = new TH2F("mp_l2sina","mp_l2sina",100,-2,2,100,0,50);
  TH2F*  mp_l3sina = new TH2F("mp_l3sina","mp_l3sina",100,-2,2,100,0,50);
  TH2F*  mp_l4sina = new TH2F("mp_l4sina","mp_l4sina",100,-2,2,100,0,50);
  TH2F*  mp_l0sina = new TH2F("mp_l0sina","mp_l0sina",100,-2,2,100,0,50);


  TH2F*  mp_l1mHpp = new TH2F("mp_l1mHpp","mp_l1mHpp",100,-2,2,100,0,1000);
  TH2F*  mp_l2mHpp = new TH2F("mp_l2mHpp","mp_l2mHpp",100,-2,2,100,0,1000);
  TH2F*  mp_l3mHpp = new TH2F("mp_l3mHpp","mp_l3mHpp",100,-2,2,100,0,1000);
  TH2F*  mp_l4mHpp = new TH2F("mp_l4mHpp","mp_l4mHpp",100,-2,2,100,0,1000);
  TH2F*  mp_l0mHpp = new TH2F("mp_l0mHpp","mp_l0mHpp",100,-2,2,100,0,1000);

  TH2F*  mp_l1mHp = new TH2F("mp_l1mHp","mp_l1mHp",100,-2,2,100,0,1000);
  TH2F*  mp_l2mHp = new TH2F("mp_l2mHp","mp_l2mHp",100,-2,2,100,0,1000);
  TH2F*  mp_l3mHp = new TH2F("mp_l3mHp","mp_l3mHp",100,-2,2,100,0,1000);
  TH2F*  mp_l4mHp = new TH2F("mp_l4mHp","mp_l4mHp",100,-2,2,100,0,1000);
  TH2F*  mp_l0mHp = new TH2F("mp_l0mHp","mp_l0mHp",100,-2,2,100,0,1000);
  
  
  TH2F*  mp_l1mA0 = new TH2F("mp_l1mA0","mp_l1mA0",100,-2,2,100,0,1000);
  TH2F*  mp_l2mA0 = new TH2F("mp_l2mA0","mp_l2mA0",100,-2,2,100,0,1000);
  TH2F*  mp_l3mA0 = new TH2F("mp_l3mA0","mp_l3mA0",100,-2,2,100,0,1000);
  TH2F*  mp_l4mA0 = new TH2F("mp_l4mA0","mp_l4mA0",100,-2,2,100,0,1000);
  TH2F*  mp_l0mA0 = new TH2F("mp_l0mA0","mp_l0mA0",100,-2,2,100,0,1000);
  
  
  TH2F*  mp_l1mH0 = new TH2F("mp_l1mH0","mp_l1mH0",100,-2,2,100,0,1000);
  TH2F*  mp_l2mH0 = new TH2F("mp_l2mH0","mp_l2mH0",100,-2,2,100,0,1000);
  TH2F*  mp_l3mH0 = new TH2F("mp_l3mH0","mp_l3mH0",100,-2,2,100,0,1000);
  TH2F*  mp_l4mH0 = new TH2F("mp_l4mH0","mp_l4mH0",100,-2,2,100,0,1000);
  TH2F*  mp_l0mH0 = new TH2F("mp_l0mH0","mp_l0mH0",100,-2,2,100,0,1000);
  
  
  TH2F*  mp_l1mh0 = new TH2F("mp_l1mh0","mp_l1mh0",100,-2,2,100,0,1000);
  TH2F*  mp_l2mh0 = new TH2F("mp_l2mh0","mp_l2mh0",100,-2,2,100,0,1000);
  TH2F*  mp_l3mh0 = new TH2F("mp_l3mh0","mp_l3mh0",100,-2,2,100,0,1000);
  TH2F*  mp_l4mh0 = new TH2F("mp_l4mh0","mp_l4mh0",100,-2,2,100,0,1000);
  TH2F*  mp_l0mh0 = new TH2F("mp_l0mh0","mp_l0mh0",100,-2,2,100,0,1000);
 
 
  TH1F*  mp_Hppmu = new TH1F("mp_Hppmu","mp_Hppmu",100,0,1000);
  TH1F*  mp_Hpmu = new TH1F("mp_Hpmu","mp_Hpmu",100,0,1000);
  TH1F*  mp_A0mu = new TH1F("mp_A0mu","mp_A0mu",100,0,1000);
  TH1F*  mp_H0mu = new TH1F("mp_H0mu","mp_H0mu",100,0,1000);
  TH1F*  mp_h0mu = new TH1F("mp_h0mu","mp_h0mu",100,0,1000);

float mu_new;
int mass=0;

////coment this for and also a brakeat at the end if you want just a fixed H++///////////
//for(int massHpp=0; massHpp<=600; massHpp = massHpp + 20.){

//FixedDoublyChargedHiggsMass  = mHpp + massHpp;
FixedDoublyChargedHiggsMass=220.;
if(option_fixed ==2 ){
////fixing FixedDoublyChargedHiggsMass and FixedSingeChargedHiggsMass
  mp_contor_found=0;

while ( mp_contor_found<1  ) {
sA=sina;
cA = sqrt(1-sA*sA);
mh02=125.011244065*125.011244065;//125.*125.;

double  FixedA0Mass =1.;
double FixedSingeChargedHiggsMass=1. ;//



 TF1 *frand6 = new TF1("x","x",-100,100);

  TF1 *frand5 = new TF1("x","x",-5,5);
  TF1 *frand4 = new TF1("x","x",-4,4);
  TF1 *frand3 = new TF1("x","x",-3,3);
  TF1 *frand2 = new TF1("x","x",-2,2);
  TF1 *frand1 = new TF1("x","x",-1,1);


vd = 246.;




sA=(double) sina;

cA = (double) sqrt(1-sina*sina);


mHpm2=FixedSingeChargedHiggsMass*FixedSingeChargedHiggsMass;


double kappa = 8.;
 bool  tmpbool; 
double condition =0;
double condition_mu_mu_new=0;
double condition_lam =0;


double Deltam= FixedDoublyChargedHiggsMass - FixedSingeChargedHiggsMass;
double Deltam_lower=1;
double Deltam_upper=1;

if(FixedDoublyChargedHiggsMass > (1/2.)*vd) {

Deltam_lower = FixedDoublyChargedHiggsMass - 1/2.*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd) ;
Deltam_upper =  FixedDoublyChargedHiggsMass - 1/2.*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass - vd*vd) ;
std::cout<<"FixedDoublyChargedHiggsMass > (1/2)*vd for FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<std::endl;
std::cout<<"Deltam_lower = "<<Deltam_lower<<std::endl;
std::cout<<"Deltam_upper = "<<Deltam_upper<<std::endl;
std::cout<<"1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd) = "<<1/2.*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd)<<std::endl;
std::cout<<"FixedDoublyChargedHiggsMass - 1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass - vd*vd) = "<<FixedDoublyChargedHiggsMass - 1/2.*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass - vd*vd)<<std::endl;
std::cout<<"FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<std::endl;
std::cout<<"vd = "<<vd<<std::endl;

}

if(FixedDoublyChargedHiggsMass < (1/2.)*vd) {


Deltam_lower = FixedDoublyChargedHiggsMass - 1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd) ;
Deltam_upper =  FixedDoublyChargedHiggsMass + 1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd);
std::cout<<"FixedDoublyChargedHiggsMass < (1/2)*vd for FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<std::endl;
std::cout<<"Deltam_lower = "<<Deltam_lower<<std::endl;
std::cout<<"Deltam_upper = "<<Deltam_upper<<std::endl;
std::cout<<"1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd) = "<<1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass + vd*vd)<<std::endl;
std::cout<<"FixedDoublyChargedHiggsMass - 1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass - vd*vd) = "<<FixedDoublyChargedHiggsMass - 1/2*sqrt(4*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass - vd*vd)<<std::endl;
std::cout<<"FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<std::endl;
std::cout<<"vd = "<<vd<<std::endl;
}
std::cout<<"std::max(-20, Deltam_lower)"<<std::max<double>(-20, Deltam_lower)<<std::endl;
std::cout<<"std::min(20,Deltam_upper)"<<std::min<double>(20,Deltam_upper)<<std::endl;

//for(int mass=0;mass<=40;mass++){
for(float mass=std::max<float>(-20, Deltam_lower);mass<=std::min<float>(20,Deltam_upper);mass=mass+0.2){
if(FixedDoublyChargedHiggsMass<200.) FixedSingeChargedHiggsMass= sqrt(mh02 + FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass)/sqrt(2)  + (double) mass/5.;
if(FixedDoublyChargedHiggsMass>=200.) {
FixedSingeChargedHiggsMass =  FixedDoublyChargedHiggsMass + mass;
std::cout<<"FixedSingeChargedHiggsMass after for(int mass=std::max(-20, Deltam_lower);mass<=std::min(-20, Deltam_lower);mass++){ is "<<FixedSingeChargedHiggsMass<<std::endl;
}
//FixedSingeChargedHiggsMass=400.;
FixedA0Mass=  FixedDoublyChargedHiggsMass - 100. + (double) mass;
mHdpm2=FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass;
mHpm2 = FixedSingeChargedHiggsMass*FixedSingeChargedHiggsMass;

   string mHpp_char = to_string(FixedDoublyChargedHiggsMass);
sA = 0.0015;


mu = (double) ((2*vt*FixedDoublyChargedHiggsMass*FixedDoublyChargedHiggsMass  + 2*l3*vt*vt*vt - vd*vd*((4*vt*FixedSingeChargedHiggsMass*FixedSingeChargedHiggsMass)/(vd*vd+2*vt*vt)) )/(-sqrt(2)*vd*vd));
//std::cout<<" FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<std::endl;
//std::cout<<" l3 = "<<l3<<std::endl;
//std::cout<<" vt = "<<vt<<std::endl;
//std::cout<<" vd = "<<vd<<std::endl;
//std::cout<<" FixedSingeChargedHiggsMass = "<<FixedSingeChargedHiggsMass<<std::endl;

mh02 = 125.*125.;

l4 = (2*sqrt(2)*mu - (4*vt*FixedSingeChargedHiggsMass*FixedSingeChargedHiggsMass)/(vd*vd+2*vt*vt) )/vt;

l1 = (double) (-l4 + (sqrt(2)*mu)/vt + (sA*(-((mu*vd)/(sqrt(2)*vt*vt)) + (mh02 - 2*(l2 + l3)*vt*vt)/(vd*vt)))/sqrt(1 - sA*sA));


l0 = (2.*(pow((vd * ( -sqrt(2.)* mu + (l1 + l4) *vt)),2) + (sqrt(2.)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt)*mh02 - pow(mh02,2)))/(((sqrt(2.)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt) - mh02)*pow(vd,2));

condition = (-(sA*vd) + 2*cA*vt)*(-2*mh02 + (l1 + l4)*pow(vd,2) +4*(l2 + l3)*pow(vt,2));

kappa = 8.;
tmpbool =   l1 <= kappa*TMath::Pi() && l1 >= -(kappa*TMath::Pi()) && l1 + l4 <= kappa*TMath::Pi() &&  l1 + l4 >= -(kappa*TMath::Pi()) && l1 + (3*l4)/2 <= kappa*TMath::Pi() &&     l1 + (3*l4)/2 >= -(kappa*TMath::Pi()) && 10/2 <= kappa*TMath::Pi() && l0/2 >= -(kappa*TMath::Pi()) &&     2*l2 <= kappa*TMath::Pi() &&  2*l2 >= -(kappa*TMath::Pi()) && 2*(l2 + l3) <= kappa*TMath::Pi() &&  2*(l2 + l3) >= -(kappa*TMath::Pi()) &&    (l0 + 4*l2 + 8*l3 +   sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/4 <=  kappa*TMath::Pi() &&    (l0 + 4*l2 + 8*l3 +  sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/  4 >=     -(kappa*TMath::Pi()) && (l0 + 4*l2 + 8*l3 -         sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/4 <=  kappa*TMath::Pi() &&    (l0 + 4*l2 + 8*l3 -  sqrt((l0 - 4*l2 - 8*l3)*(l0 - 4*l2 - 8*l3) + 16*l4*l4))/  4 >=     -(kappa*TMath::Pi()) && (3*l0 + 16*l2 + 12*l3 + sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/4 <=  kappa*TMath::Pi() &&    (3*l0 + 16*l2 + 12*l3 +  sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/ 4 >= -(kappa*TMath::Pi()) &&    (3*l0 + 16*l2 + 12*l3 - sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) +  24*(2*l1 + l4)*(2*l1 + l4)))/ 4 <= kappa*TMath::Pi() &&  (3*l0 + 16*l2 + 12*l3 -  sqrt((-3*l0 + 16*l2 + 12*l3)*(-3*l0 + 16*l2 + 12*l3) + 24*(2*l1 + l4)*(2*l1 + l4)))/ 4 >= -(kappa*TMath::Pi()) && l1 - l4/2 <= kappa*TMath::Pi() &&    l1 - l4/2 >= -(kappa*TMath::Pi()) &&  l2 - l3 <= kappa*TMath::Pi() && 2 * l2 - l3 >= -kappa*TMath::Pi();

 mu_new = (sqrt(2.)*vt*(-(cA*sA*vt*(-2.*mh02 + (l1 + l4)*pow(vd,2) + 4.*(l2 + l3)*pow(vt,2))) + vd*(-2.*(l1 + l4)*pow(vt,2) + pow(sA,2)*(mh02 + 2.*(l1 - l2 - l3 + l4)*pow(vt,2)))))/ (vd*(pow(sA,2)*pow(vd,2) - 4.*pow(cA,2)*pow(vt,2)));
std::cout<<" mu/mu_new = "<<mu/mu_new<<" !!!!!!!!!!!!!sina = "<<sA<<std::endl;

double sAvalue =  sA - 0.0001;
double isA = 0.000001;

////comment this for fixed sinalpha
/*
while( (abs(l0)>1.05 || abs(l1)>1.05 || abs(l2)>1.05|| abs(l3)>1.05|| abs(l4)>1.05)) {
condition_lam=0;
std::cout<<"@@@@@@@abs(l0)>1.05 || abs(l1)>1.05 || abs(l2)>1.05|| abs(l3)>1.08|| abs(l4)>1.05@@@@@@@@@  "<<std::endl;
std::cout<<"Deltam_lower = "<<Deltam_lower<<" Deltam_upper = "<<Deltam_upper<<std::endl;
std::cout<<"FixedDoublyChargedHiggsMass = "<<FixedDoublyChargedHiggsMass<<" FixedSingeChargedHiggsMass = "<<FixedSingeChargedHiggsMass<<std::endl;
std::cout<<"sina = "<<sA<<" vt = "<<vt<<" vd = "<<vd<<std::endl;
std::cout<<"l0 = "<<l0<<" l1 = "<<l1<<" l2 = "<<l2<<" l3 = "<<l3<<" l4 = "<<l4<<" mu = "<<mu<<" mu_new = "<<mu_new<<std::endl;


std::cout<<"End of When Caluculating!!!!!!!"<<std::endl;
l2 = frand1->GetRandom(); 
l3 = frand1->GetRandom();


l4 = (2*sqrt(2)*mu - (4*vt*FixedSingeChargedHiggsMass*FixedSingeChargedHiggsMass)/(vd*vd+2*vt*vt) )/vt;

l1 = (double) (-l4 + (sqrt(2)*mu)/vt + (sA*(-((mu*vd)/(sqrt(2)*vt*vt)) + (mh02 - 2*(l2 + l3)*vt*vt)/(vd*vt)))/sqrt(1 - sA*sA));


l0 = (2.*(pow((vd * ( -sqrt(2.)* mu + (l1 + l4) *vt)),2) + (sqrt(2.)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt)*mh02 - pow(mh02,2)))/(((sqrt(2.)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt) - mh02)*pow(vd,2));



sA =  sAvalue + isA;
isA = isA + 0.000005;

}
*/
condition_lam=1.;
std::cout<<"When Caluculating!!!!!!!"<<std::endl;

std::cout<<"sina = "<<sA<<" vt = "<<vt<<" vd = "<<vd<<std::endl;
std::cout<<"l0 = "<<l0<<" l1 = "<<l1<<" l2 = "<<l2<<" l3 = "<<l3<<" l4 = "<<l4<<" mu = "<<mu<<" mu_new = "<<mu_new<<std::endl;
std::cout<<"End of When Caluculating!!!!!!!"<<std::endl;
if(abs(mu/mu_new-1)<0.0001) condition_mu_mu_new=1;
if(condition>0) std::cout<<"condition>0"<<std::endl;
else std::cout<<"condition<0"<<std::endl;
if(condition_mu_mu_new>0) std::cout<<"condition_mu_mu_new>0"<<std::endl;
else std::cout<<"condition_mu_mu_new<0 ratio"<<mu/mu_new<<std::endl;
if(condition_lam>0) std::cout<<"condition_lam>0"<<std::endl;
else std::cout<<"condition_lam<0"<<std::endl;

if(condition>0 && condition_mu_mu_new>0 && condition_lam>0) {
if(tmpbool){

A = l0/2. * vd*vd;
B = (vd * ( -sqrt(2.)* mu + (l1 + l4) *vt));
CC = (sqrt(2.)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt);
mA2 = mu* (vd*vd +  4. *vt*vt)/sqrt(2)/vt;
mH02  = 0.5 *(  l0/2. * vd*vd + (sqrt(2)* mu *vd*vd + 4*(l2 + l3)* pow(vt,3))/(2.*vt) + sqrt( ( l0/2. * vd*vd - (sqrt(2.)* mu *vd*vd + 4.*(l2 + l3)* pow(vt,3))/(2.*vt))*( l0/2. * vd*vd - (sqrt(2.)* mu *vd*vd + 4.*(l2 + l3)* pow(vt,3))/(2.*vt)) + 4. *(vd * ( -sqrt(2.)* mu + (l1 + l4) *vt))*(vd * ( -sqrt(2.)* mu + (l1 + l4) *vt))));

cA = sqrt(1 - sA*sA); 

std::cout<<"mHdpm2 = "<<mHdpm2<<" mA2 = "<<mA2<<" mHpm2 = "<<mHpm2<<" mh02 = "<<mh02<<" mH02 = "<<" sA = "<<sA<<std::endl;

if(option_fixed==2 && mHdpm2>0. && mA2>0. && mHpm2>0. && mh02>0. && mH02>0. ) {
std::cout<<"FOUND!!!!!!!"<<std::endl;
std::cout<<"A = "<<A<<std::endl;
std::cout<<"B = "<<B<<std::endl;
std::cout<<"CC = "<<CC<<std::endl;
std::cout<<" vt = "<<vt<<" vd = "<<vd<<std::endl;
cout << setprecision(10);


std::cout<<" mHdpm = "<<sqrt(mHdpm2)<<std::endl;
std::cout<<" mHpm2 = "<<sqrt(mHpm2)<<std::endl;
std::cout<<" mA2 = "<<sqrt(mA2)<<std::endl;
std::cout<<" mH02 = "<<sqrt(mH02)<<std::endl;
std::cout<<" mh02 = "<<sqrt(mh02)<<std::endl;
std::cout<<" sina = "<<sA<<std::endl;
std::cout<<" mu = "<<mu<<std::endl;
std::cout<<" l0 = "<<l0<<std::endl;
std::cout<<" l1 = "<<l1<<std::endl;
std::cout<<" l2 = "<<l2<<std::endl;
std::cout<<" l3 = "<<l3<<std::endl;
std::cout<<" l4 = "<<l4<<std::endl;

   string mp_contor_found_char = to_string(mp_contor_found);

eff_file_e<<"launch--name=pp_HpHm_H++"+mHpp_char+"_"+mp_contor_found_char<<std::endl;
eff_file_e<<"set wh__2 auto"<<std::endl;
eff_file_e<<"set wha auto"<<std::endl;
eff_file_e<<"set whp auto"<<std::endl;
eff_file_e<<"set whpp auto"<<std::endl;

eff_file_e<<"set mHpp "<<sqrt(mHdpm2)<<std::endl;
eff_file_e<<"set mHp "<<sqrt(mHpm2)<<std::endl;
eff_file_e<<"set mha "<<sqrt(mA2)<<std::endl;
eff_file_e<<"set mh__2 "<<sqrt(mH02)<<std::endl;
eff_file_e<<"set mh "<<sqrt(mh02)<<std::endl;
eff_file_e<<"set vt "<<vt<<std::endl;
eff_file_e<<"set sinalpha "<<sA<<std::endl;
eff_file_e<<"set mu "<<mu<<std::endl;
eff_file_e<<"set lam "<<l0<<std::endl;
eff_file_e<<"set lam1 "<<l1<<std::endl;
eff_file_e<<"set lam2 "<<l2<<std::endl;
eff_file_e<<"set lam3 "<<l3<<std::endl;
eff_file_e<<"set lam4 "<<l4<<std::endl;
eff_file_e<<std::endl;


 mp_l1l4->Fill(l1,l4);
 mp_l1l2->Fill(l1,l2);
 mp_l1l3->Fill(l1,l3);
 mp_l2l3->Fill(l2,l3);
 mp_l2l4->Fill(l2,l4);
 mp_l3l4->Fill(l3,l4);

 mp_l1mu->Fill(sqrt(mHdpm2),mu*100);
 mp_l2mu->Fill(sqrt(mHpm2),mu*100);
 mp_l3mu->Fill(sqrt(mH02),mu*100);
 mp_l4mu->Fill(sqrt(mh02),mu*100);
 mp_l0mu->Fill(sA*100,mu*100);


 mp_l1sina->Fill(l1,sA*1000);
 mp_l2sina->Fill(l2,sA*1000);
 mp_l3sina->Fill(l3,sA*1000);
 mp_l4sina->Fill(l4,sA*1000);
 mp_l0sina->Fill(l0,sA*1000);

 mp_l1mHpp->Fill(l1,sqrt(mHdpm2));
 mp_l2mHpp->Fill(l2,sqrt(mHdpm2));
 mp_l3mHpp->Fill(l3,sqrt(mHdpm2));
 mp_l4mHpp->Fill(l4,sqrt(mHdpm2));
 mp_l0mHpp->Fill(l0,sqrt(mHdpm2));

 mp_l1mHp->Fill(l1,sqrt(mHpm2));
 mp_l2mHp->Fill(l2,sqrt(mHpm2));
 mp_l3mHp->Fill(l3,sqrt(mHpm2));
 mp_l4mHp->Fill(l4,sqrt(mHpm2));
 mp_l0mHp->Fill(l0,sqrt(mHpm2));


 mp_l1mA0->Fill(l1,sqrt(mA2));
 mp_l2mA0->Fill(l2,sqrt(mA2));
 mp_l3mA0->Fill(l3,sqrt(mA2));
 mp_l4mA0->Fill(l4,sqrt(mA2));
 mp_l0mA0->Fill(l0,sqrt(mA2));

 mp_l1mH0->Fill(l1,sqrt(mH02));
 mp_l2mH0->Fill(l2,sqrt(mH02));
 mp_l3mH0->Fill(l3,sqrt(mH02));
 mp_l4mH0->Fill(l4,sqrt(mH02));
 mp_l0mH0->Fill(l0,sqrt(mH02));


 mp_l1mh0->Fill(l1,sqrt(mh02));
 mp_l2mh0->Fill(l2,sqrt(mh02));
 mp_l3mh0->Fill(l3,sqrt(mh02));
 mp_l4mh0->Fill(l4,sqrt(mh02));
 mp_l0mh0->Fill(l0,sqrt(mh02));







   mp_vec_sina.push_back(sA*100); mp_vec_vt.push_back(vt); mp_vec_vd.push_back(vd); mp_vec_l0.push_back(l0); mp_vec_l1.push_back(l1); mp_vec_l2.push_back(l2); mp_vec_l3.push_back(l3); mp_vec_l4.push_back(l4); mp_vec_mu.push_back(mu*100);  mp_vec_mHdpm.push_back(sqrt(mHdpm2)); mp_vec_mA2.push_back(sqrt(mA2));  mp_vec_mHpm2.push_back(sqrt(mHpm2));   mp_vec_mH02.push_back(sqrt(mH02)); mp_vec_mh02.push_back(sqrt(mh02));
}

 }
//////////comment here the breaket for a fixed one single mass of H++///////////////
//}

mp_contor_found++;
}



} //end of while
} //end of for
}
TCanvas *mp_mass = new TCanvas ("mp_mass for vt=0.1 and mH0=125 and |lambda|<=5 and sina in [0,0.001]","mp_mass for vt=0.1 and mH0=125 and |lambda|<=5 and sina in [0,0.001]",950, 950);

 //mp_mass -> SetLogy() ;// SetLogy in current pad
   TLegend *legend_mp_mass = new TLegend(0.85,0.65,0.95,0.95);
    legend_mp_mass->SetBorderSize(0);
    legend_mp_mass->SetTextFont(42);
    legend_mp_mass->SetTextSize(0.03);
    legend_mp_mass->SetLineColor(0);
    legend_mp_mass->SetLineColor(0);


    
TH1F*  mp_mHdpm = new TH1F("mHpp","mHpp", mp_vec_mHdpm.size(),0,mp_vec_mHdpm.size());
TH1F*  mp_mA2 = new TH1F("mA","mA",mp_vec_mA2.size(),0,mp_vec_mA2.size());
TH1F*  mp_mHpm2 = new TH1F("mHp","mHp",mp_vec_mHpm2.size(),0,mp_vec_mHpm2.size());
TH1F*  mp_mH02 = new TH1F("mH0","mH0",mp_vec_mH02.size(),0,mp_vec_mH02.size());
TH1F*  mp_mh02 = new TH1F("mh0","mh0",mp_vec_mh02.size(),0,mp_vec_mh02.size());
TH1F*  mp_sina = new TH1F("sina","sina",mp_vec_sina.size(),0,mp_vec_sina.size());
TH1F*  mp_mu = new TH1F("mu","mu",mp_vec_mu.size(),0,mp_vec_mu.size());

if(ratio_value == 0) {

       mp_sina->SetBinContent(0,0);
       mp_sina->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_sina.size();i++) {
       mp_sina->SetBinContent(i+1,mp_vec_sina[i]);
       mp_sina->SetBinError(i+1,0);

       } 
       
       mp_mu->SetBinContent(0,0);
       mp_mu->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mu.size();i++) {
       mp_mu->SetBinContent(i+1,mp_vec_mu[i]);
       mp_mu->SetBinError(i+1,0);

       } 
       mp_mHdpm->SetBinContent(0,0);
       mp_mHdpm->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mHdpm.size();i++) {
       mp_mHdpm->SetBinContent(i+1,mp_vec_mHdpm[i]);
       mp_mHdpm->SetBinError(i+1,0);

       } 
     

       mp_mA2->SetBinContent(0,0);
       mp_mA2->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mA2.size();i++) {
       mp_mA2->SetBinContent(i+1,mp_vec_mA2[i]);
       mp_mA2->SetBinError(i+1,0);

       }      

       mp_mHpm2->SetBinContent(0,0);
       mp_mHpm2->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mHpm2.size();i++) {
       mp_mHpm2->SetBinContent(i+1,mp_vec_mHpm2[i]);
       mp_mHpm2->SetBinError(i+1,0);

       } 
       

       mp_mH02->SetBinContent(0,0);
       mp_mH02->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mH02.size();i++) {
       mp_mH02->SetBinContent(i+1,mp_vec_mH02[i]);
       mp_mH02->SetBinError(i+1,0);

       } 

       mp_mh02->SetBinContent(0,0);
       mp_mh02->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mh02.size();i++) {
       mp_mh02->SetBinContent(i+1,mp_vec_mh02[i]);
       mp_mh02->SetBinError(i+1,0);

       } 
       }
       
       
if(ratio_value == 1) {

       mp_mHdpm->SetBinContent(0,0);
       mp_mHdpm->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mHdpm.size();i++) {
       mp_mHdpm->SetBinContent(i+1,mp_vec_mHdpm[i]/mp_vec_mHdpm[i]);
       mp_mHdpm->SetBinError(i+1,0);

       } 
     

       mp_mA2->SetBinContent(0,0);
       mp_mA2->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mA2.size();i++) {
       mp_mA2->SetBinContent(i+1,mp_vec_mA2[i]/mp_vec_mHdpm[i]);
       mp_mA2->SetBinError(i+1,0);

       }      

       mp_mHpm2->SetBinContent(0,0);
       mp_mHpm2->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mHpm2.size();i++) {
       mp_mHpm2->SetBinContent(i+1,mp_vec_mHpm2[i]/mp_vec_mHdpm[i]);
       mp_mHpm2->SetBinError(i+1,0);

       } 
       

       mp_mH02->SetBinContent(0,0);
       mp_mH02->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mH02.size();i++) {
       mp_mH02->SetBinContent(i+1,mp_vec_mH02[i]/mp_vec_mHdpm[i]);
       mp_mH02->SetBinError(i+1,0);

       } 

       mp_mh02->SetBinContent(0,0);
       mp_mh02->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_mh02.size();i++) {
       mp_mh02->SetBinContent(i+1,mp_vec_mh02[i]/mp_vec_mHdpm[i]);
       mp_mh02->SetBinError(i+1,0);

       } 
       }
  if(ratio_value == 0){     
legend_mp_mass->AddEntry(mp_mHdpm,"mHpp");
legend_mp_mass->AddEntry(mp_mA2,"mA0");
legend_mp_mass->AddEntry(mp_mHpm2,"mHp");
legend_mp_mass->AddEntry(mp_mH02,"mH0");
legend_mp_mass->AddEntry(mp_mh02,"mh0");
legend_mp_mass->AddEntry(mp_sina,"sina*100");
legend_mp_mass->AddEntry(mp_mu,"mu*100");

}
  if(ratio_value == 1){     
legend_mp_mass->AddEntry(mp_mHdpm,"mHpp/mHpp");
legend_mp_mass->AddEntry(mp_mA2,"mA0/mHpp");
legend_mp_mass->AddEntry(mp_mHpm2,"mHp/mHpp");
legend_mp_mass->AddEntry(mp_mH02,"mH0/mHpp");
legend_mp_mass->AddEntry(mp_mh02,"mh0/mHpp");
}

mp_mHdpm->GetXaxis()->SetTitle("Variations;"); 
if(ratio_value == 0) mp_mHdpm->GetYaxis()->SetTitle("M (GeV)"); 
if(ratio_value == 1) mp_mHdpm->GetYaxis()->SetTitle("M/mHpp"); 

if(ratio_value == 0) mp_mHdpm->GetYaxis()->SetRangeUser(-10,FixedDoublyChargedHiggsMass+500); 
if(ratio_value == 1) mp_mHdpm->GetYaxis()->SetRangeUser(0,1.5);//FixedDoublyChargedHiggsMass+500); 

mp_mHdpm->SetMarkerColorAlpha(kBlue, 0.35);
mp_mA2->SetMarkerColorAlpha(kOrange, 0.35);
mp_mHpm2->SetMarkerColorAlpha(kRed, 0.35);
mp_mH02->SetMarkerColorAlpha(kGreen, 0.35);
mp_mh02->SetMarkerColorAlpha(kMagenta, 0.35);
mp_sina->SetMarkerColorAlpha(kBlack, 0.35);
mp_mu->SetMarkerColorAlpha(kYellow, 0.35);

mp_mHdpm->SetMarkerSize(1.5);
mp_mA2->SetMarkerSize(1.5);
mp_mHpm2->SetMarkerSize(1.5);
mp_mH02->SetMarkerSize(1.5);
mp_mh02->SetMarkerSize(1.5);
mp_sina->SetMarkerSize(1.5);
mp_mu->SetMarkerSize(1.5);


mp_mHdpm->Draw("P0");
mp_mA2->Draw("P0same");
mp_mHpm2->Draw("P0same");
mp_mH02->Draw("P0same");
mp_mh02->Draw("P0same");
//mp_sina->Draw("P0same");
//mp_mu->Draw("P0same");


legend_mp_mass->Draw("same");

//vary lambda [-1,1] option_vary_lambda = 1;
//vary lambda [-2,2] option_vary_lambda = 2;
//vary lambda [-3,3] option_vary_lambda = 3;
//vary lambda [-4,4] option_vary_lambda = 4;
//vary lambda [-5,5] option_vary_lambda = 5;

//fixing only FixedDoublyChargedHiggsMass and h0=125 option_fixed = 1;
//fixing FixedDoublyChargedHiggsMass and FixedSingeChargedHiggsMass and h0=125 option_fixed = 2;
//fixing FixedDoublyChargedHiggsMass and A0 and h0=125 option_fixed = 3;
//fixing h0=125 and A0 option_fixed = 4;
string tip = "";


 if(option_vary_lambda == 1 && option_fixed == 1)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_and_abslambda1";
  if(option_vary_lambda == 2 && option_fixed == 1)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_and_abslambda2";
   if(option_vary_lambda == 3 && option_fixed == 1)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_and_abslambda3";
    if(option_vary_lambda == 4 && option_fixed == 1)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_and_abslambda4";
     if(option_vary_lambda == 5 && option_fixed == 1)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_and_abslambda5";
     
  if(option_vary_lambda == 1 && option_fixed == 2)  tip = "MpformH0125andFixedDoublyChargedHiggsMassandFixedSingeChargedHiggsMassandabslambda1";
  if(option_vary_lambda == 2 && option_fixed == 2)  tip = "MpformH0125andFixedDoublyChargedHiggsMassandFixedSingeChargedHiggsMassandabslambda2";
   if(option_vary_lambda == 3 && option_fixed == 2)  tip = "MpformH0125andFixedDoublyChargedHiggsMassandFixedSingeChargedHiggsMassandabslambda3";
    if(option_vary_lambda == 4 && option_fixed == 2)  tip = "MpformH0125andFixedDoublyChargedHiggsMassandFixedSingeChargedHiggsMassandabslambda4";
     if(option_vary_lambda == 5 && option_fixed == 2)  tip = "MpformH0125andFixedDoublyChargedHiggsMassandFixedSingeChargedHiggsMassandabslambda5";
     
   if(option_vary_lambda == 1 && option_fixed == 3)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_andFixedA0_and_abslambda1";
  if(option_vary_lambda == 2 && option_fixed == 3)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_andFixedA0_and_abslambda2";
   if(option_vary_lambda == 3 && option_fixed == 3)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_andFixedA0_and_abslambda3";
    if(option_vary_lambda == 4 && option_fixed == 3)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_andFixedA0_and_abslambda4";
     if(option_vary_lambda == 5 && option_fixed == 3)  tip = "Mp_for_mH0125_and_FixedDoublyChargedHiggsMass_andFixedA0_and_abslambda5";
     
    if(option_vary_lambda == 1 && option_fixed == 4)  tip = "Mp_for_mH0125_and_varyDoublyChargedHiggsMass_andFixedA0_and_abslambda1";
  if(option_vary_lambda == 2 && option_fixed == 4)  tip = "Mp_for_mH0125_and_varyDoublyChargedHiggsMass_andFixedA0_and_abslambda2";
   if(option_vary_lambda == 3 && option_fixed == 4)  tip = "Mp_for_mH0125_and_varyDoublyChargedHiggsMass_andFixedA0_and_abslambda3";
    if(option_vary_lambda == 4 && option_fixed == 4)  tip = "Mp_for_mH0125_and_varyDoublyChargedHiggsMass_andFixedA0_and_abslambda4";
     if(option_vary_lambda == 5 && option_fixed == 4)  tip = "Mp_for_mH0125_and_varyDoublyChargedHiggsMass_andFixedA0_and_abslambda5";
     
                     
 if(ratio_value == 0) mysstream << tip <<"TargerMass"<<FixedDoublyChargedHiggsMass<<"kappa8.eps";
 if(ratio_value == 1) mysstream << tip <<"TargerMass"<<FixedDoublyChargedHiggsMass<<"andvt"<<vt<<"andvd" << vd<<"kappa8_ratioMHpp.eps";

 


mp_mass->Print(line.c_str());


TCanvas *sina_mass = new TCanvas ("sina for vary mu","sina for vary mu",950, 950);

 //sina_mass -> SetLogy() ;// SetLogy in current pad
   TLegend *legend_sina_mass = new TLegend(0.75,0.75,0.95,0.95);
    legend_sina_mass->SetBorderSize(0);
    legend_sina_mass->SetTextFont(42);
    legend_sina_mass->SetTextSize(0.023);
    legend_sina_mass->SetLineColor(0);
    legend_sina_mass->SetLineColor(0);


TH1F*  sina_sina_mass = new TH1F("sina for vt=0.1 and mH0=125 and |lambda|<=5 and sina in [0,0.001]","sina for vt=0.1 and mH0=125 and |lambda|<=5 and sina in [0,0.001]",mp_vec_sina.size(),0,mp_vec_sina.size());

       sina_sina_mass->SetBinContent(0,0);
       sina_sina_mass->SetBinError(0,0);
       for(unsigned i=0;i<mp_vec_sina.size();i++) {
       sina_sina_mass->SetBinContent(i+1,mp_vec_sina[i]);
       sina_sina_mass->SetBinError(i+1,0);

       } 
 

legend_sina_mass->AddEntry(sina_sina_mass,"sina");

sina_sina_mass->GetXaxis()->SetTitle("mu"); 
sina_sina_mass->GetYaxis()->SetTitle("sina"); 
sina_sina_mass->GetYaxis()->SetRange(0,1);
sina_sina_mass->SetMarkerColorAlpha(kOrange, 0.35);

sina_sina_mass->SetMarkerSize(1.5);


sina_sina_mass->Draw("P0");


legend_sina_mass->Draw("same");

//vary lambda [-1,1] option_vary_lambda = 1;
//vary lambda [-2,2] option_vary_lambda = 2;
//vary lambda [-3,3] option_vary_lambda = 3;
//vary lambda [-4,4] option_vary_lambda = 4;
//vary lambda [-5,5] option_vary_lambda = 5;

//fixing only FixedDoublyChargedHiggsMass and h0=125 option_fixed = 1;
//fixing FixedDoublyChargedHiggsMass and FixedSingeChargedHiggsMass and h0=125 option_fixed = 2;
//fixing FixedDoublyChargedHiggsMass and A0 and h0=125 option_fixed = 3;
//fixing h0=125 and A0 option_fixed = 4;
string tip2 = "";

 ostringstream mysstream2; 
 if(option_fixed == 1)  tip2 = "sina";
 if(option_fixed == 2)  tip2 = "sina";
 if(option_fixed == 3)  tip2 = "sina";
 if(option_fixed == 4)  tip2 = "sina";

                     

 mysstream2 << "vt="<<vt<<"vd=" << vd;
 string line2 = mysstream2.str();
 
  TPaveText *t=new TPaveText(0.1,0.8,0.6,0.85,"brNDC");
    t->AddText(line2.c_str());
   //t->Draw();

 string line3 = mysstream2.str()+".eps";


sina_mass->Print(line3.c_str());



TCanvas *graphmuHpp = new TCanvas ("graphmuHpp","graphmuHpp",950, 950);
   TLegend *legend_mu_mp = new TLegend(0.85,0.65,0.95,0.95);
    legend_mu_mp->SetBorderSize(0);
    legend_mu_mp->SetTextFont(42);
    legend_mu_mp->SetTextSize(0.03);
    legend_mu_mp->SetLineColor(0);
    legend_mu_mp->SetLineColor(0);



   Int_t n = mp_vec_mu.size();
   Double_t x[n], y[n];
   for (Int_t i=0; i<n; i++) {
      x[i] = mp_vec_mu[i];
      y[i] = mp_vec_l0[i];
   }
   TGraph *gr0 = new TGraph (n, x, y);
   gr0->SetMarkerColor(4);    gr0->SetLineColor(5);

   gr0->SetMarkerStyle(20);
   
   Double_t x1[n], y1[n];
   for (Int_t i=0; i<n; i++) {
      x1[i] = mp_vec_mu[i];
      y1[i] = mp_vec_l1[i];
   }
   TGraph *gr1 = new TGraph (n, x1, y1);
   gr1->SetMarkerColor(1);    gr1->SetLineColor(1);

   gr1->SetMarkerStyle(20);   


   Double_t x2[n], y2[n];
   for (Int_t i=0; i<n; i++) {
      x2[i] = mp_vec_mu[i];
      y2[i] = mp_vec_l2[i];
   }
   TGraph *gr2 = new TGraph (n, x2, y2);
   gr2->SetMarkerColor(2);    gr2->SetLineColor(2);

   gr2->SetMarkerStyle(20); 

   Double_t x3[n], y3[n];
   for (Int_t i=0; i<n; i++) {
      x3[i] = mp_vec_mu[i];
      y3[i] = mp_vec_l3[i];
   }
   TGraph *gr3 = new TGraph (n, x3, y3);
   gr3->SetMarkerColor(3);   gr3->SetLineColor(3);

   gr3->SetMarkerStyle(20); 
   
   Double_t x4[n], y4[n];
   for (Int_t i=0; i<n; i++) {
      x4[i] = mp_vec_mu[i];
      y4[i] = mp_vec_l4[i];
   }
   TGraph *gr4 = new TGraph (n, x4, y4);
   gr4->SetMarkerColor(5);
   gr4->SetLineColor(5);
   gr4->SetMarkerStyle(20); 
         
      
    TMultiGraph  *mg  = new TMultiGraph();
   mg->Add(gr0);
   mg->Add(gr1);
   mg->Add(gr2);
   mg->Add(gr3);
   mg->Add(gr4);
mg->GetYaxis()->SetTitle("Lambdas"); 
mg->GetXaxis()->SetTitle("mu*100"); 

mg->GetYaxis()->SetRangeUser(-10,FixedDoublyChargedHiggsMass+500); 

  
legend_mu_mp->AddEntry(gr0,"l0");
legend_mu_mp->AddEntry(gr1,"l1");
legend_mu_mp->AddEntry(gr2,"l2");
legend_mu_mp->AddEntry(gr3,"l3");
legend_mu_mp->AddEntry(gr4,"l4");

   mg->Draw("ALP");
 legend_mu_mp->Draw("same");  
  
  

graphmuHpp->Print(line2.c_str());




}
  
 #include <iostream>
#include <vector>
void Solve(vector<double>& roots, int nest, TF1* f, double s, double e)
{
  double r = f->GetX(0, s, e);
  cout << "(" << nest << ") start=" << s << " end=" << e << " root=" << r << endl;
  roots.push_back(r);
  if (r != s && r !=e) {
    Solve(roots, nest+1, f, s, r);
    Solve(roots, nest+1, f, r, e);
  }
}
    
void SetPlotStyle()
{

  // use plain black on white colors
  Int_t icol=0; // WHITE
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetFrameFillColor(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetStatColor(icol);
  //gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);

  // set margin sizes
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  //double tsize=0.05;
  double tsize=0.04;
  gStyle->SetTextFont(font);

  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");
  
  gStyle->SetLabelSize(tsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(tsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(tsize,"z");
  gStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //gStyle->SetErrorX(0.001);
  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);




 gROOT->ForceStyle();

}
