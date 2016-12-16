/*
  vector<double> Fn1/2; // lookup table for Fn based on different gamma, for gamma from -1 to -99 should go to absolute+100
  vector<double> Fn1/2_d; // lookup table for Fn_derivative of gamma, based on different gamma, for gamma from -1 to -99 should go to absolute+100
 * Add subroutine GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn) under PRFCluster
 */
#if !defined(PRFCluster_H)
#define  PRFCluster_H

#define CRI(a) ((a)==0?"BIC":(a)==1?"AIC":"AICc")

#define NAME "MASS-PRF"	//Program name
#define FULLNAME "Model Averaged Site Selection - Poission Random Field"
#define VERSION    "1.3"
#define LASTUPDATE "December 16th, 2016"
#define FUNCTION "Estimate selection intensity for each single site in coding sequences by using polymorphism and divergence data."
#define REFERENCE "Zi-Ming Zhao, Ning Li, Zhang Zhang and Jeffrey P. Townsend. (2016) Regions within coding gene sequences experience diverse intensities of natural selection inferred from polymorphism and divergence."



#define IR_H 50      //Initial value of r
#define IR_L -50
#define ER 0.001   //f(r)<=ER, this r is the one we are trying to find
#define MinDx 0.4   //f(r)<=MinDx, use the optimal gamma meeting this criteria even though ER is not reached
#define JMAX 200  //Maximum number of iteration
#define INF 50 // Maximum gamma
#define NINF -50 // Minimum gamma

//#define Model_Number 10000 //Choose models randomly, how many models are chosen.


#include "base.h"
#include <map>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif
#include <time.h>
#include <thread>

using namespace std;

class CandidateModels {
 public:
  double CW; //criterion and weight
  double LogLikelihood;
  
  long pos_start, pos_end;
  double p0, pc;
  long cs, ce;
  
  double InL0, InL;
  double AIC0, AIC;
  double AICc0, AICc;
  double BIC0, BIC;
  
  CandidateModels() {}
  
  //For all possible candidate models
  CandidateModels(double CW, long pos_start, long pos_end, long cs, long ce, double p0, double pc,double LogLikelihood){
    this->CW = CW;
    this->pos_start = pos_start;		this->pos_end = pos_end;
    this->cs = cs;		this->ce = ce;
    this->p0 = p0;		this->pc = pc;
    this->LogLikelihood = LogLikelihood;
  }
  
  //For best-fit selected models
  CandidateModels(long pos_start, long pos_end, long cs, long ce, double p0, double pc, double InL0, double InL,double AIC0, double AIC, double AICc0, double AICc, double BIC0, double BIC) {
    this->pos_start = pos_start;		this->pos_end = pos_end;
    this->cs = cs;		this->ce = ce;
    this->p0 = p0;		this->pc = pc;
    this->InL0 = InL0;		this->InL = InL;
    this->AIC0 = AIC0;		this->AIC = AIC;
    this->AICc0 = AICc0;		this->AICc = AICc;
    this->BIC0 = BIC0;		this->BIC = BIC;
  }  
};

class CI {
 public:
  float weight;
  float p;
  float LogLikelihood;
  
  CI() {}

  CI(float weight, float p, float LogLikelihood){
    this->weight = weight;		
    this->p = p;
    this->LogLikelihood=LogLikelihood;
  }
};

struct SiteModels{
  long pos;
  vector<CI> sms;
};

struct more_than_CI
{
    bool operator() (const CI& struct1, const CI& struct2)
    {
       if ((struct1.weight > struct2.weight) || ((struct1.weight == struct2.weight) && (struct1.p > struct2.p)))
    	   return true;
       else return false;
    }
};

struct pthread_stochastic_args
{
  struct SiteModels *pr;
  struct SiteModels *dr;
  long N;
  long species_n;
  long i;
};

class rModels {
 public:
  double weight;
  double r;
  
  rModels(){}
  
  rModels(double weight, double r){
    this->weight=weight;
    this->r=r;
  }
};


class PRFCluster: public Base {	
 public:	
  PRFCluster();
  ~PRFCluster();	
  
  //Main function
  int Run(int argc, const char*argv[]);
  //Lookup table from Mathematica to get the integration for gamma calculation
  vector<double> GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn);
  //Clustering using maximum likelihood
  int RunML(vector<string> pol_seq, vector<string> div_seq);
  

  //AIC/BIC based on sub sequence
  int ClusterSubSeq(int pos_start, int pos_end, char symbol='S', int flag_seq=0, struct SiteModels *p=NULL);
  void Time(int time_start);
  
  int init(long N);
  int ModelAveraging(long pos_start, long pos_end, long cs, long ce, double p0, double pc, double min_cri,struct SiteModels *p);
  int CI_MA(struct SiteModels *p,long N);
  void CI_MA_threaded(struct SiteModels *pointer,long i);

  int SitePRF(int species_n, long N);
  int SiteNI(long N);

  void SiteSpecificDivergentTime(long species_n, struct SiteModels* psp, struct SiteModels* psd, int model_ps, int model_ds, int i, vector<double> &div);
  void DivergentTime(long species_n, long N, vector<double> &div);
  //double ClusterPRF(long species_n, double pr, double dr);
	

  int r_stochastic(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n);
  void r_stochastic_threaded(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n, long i, std::ostringstream* myout);


  int BubbleSort(struct SiteModels *p,long site);
  long RandomModel_Num(struct SiteModels *p,long site);
  int CI_UpLow_r_stochastic(long site);
  int CI_UpLow_r_thread(long site,double min_weight, vector<rModels> vec_rModels);
  
  int r_exact(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n);
  void r_exact_threaded(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n, long i, std::ostringstream* myout);
  
  double ModelSpecific_r_PRF(double p_pr, double p_dr,long species_n, double site_divergent_time);
  int CI_UpLow_r_exact(long site,double min_weight);
  int ModelAveraged_r(long site,double min_weight);
  int ModelAveraged_r_thread(long site,double min_weight,vector<rModels> vec_rModels);

  int CI_UpLow_r(long site,double min_weight);


  int output(long N);

 protected:
  
  double BinomialProb(long n, long i);

  
  //int nw (from pos_start to pos_end, all of the nonsynonymous)
  //int nc (from cs to ce, all of the nonsynonymous)
  double getp0pc_MK(int pos_start, int pos_end, int cs, int ce, double &p0, double &pc, int nw, int nc);


  //flag_N_pol=0 means there is all A/T/G/C in the polymophism data, otherwise there is Ns or other ambiguous nucleotide
  //flag_N_div=0 means there is all A/T/G/C in the divergence data, otherwise there is Ns or other ambiguous nucleotide
  int flag_N_pol;
  int flag_N_div;
  string getPolSysRep(vector<string> seq);
  string getDivSysRep(vector<string> pol, vector<string> div);
  
  //Consider the sites other than A, T, G or C 
  int ReplaceCodon4PolSys(vector <string>& codon, vector <string>& codon_other);
  int ReplaceCodon4DivSys(vector <string>& codon, vector <string>& codon_other, string& div_codon);
  char ReplaceSite(char codon_other_symbol, int pointer_site[], string symbol);

  
  long getDifference(string seq, int pos_start, int pos_end, char symbol='1');
 
  
  //Parse parameters
  bool parseParameter(int argc, const char* argv[]);
  //Show help information
  void showHelpInfo();
  
  
 public:
  vector<double> Fn1; // lookup table for Fn based on different gamma
  vector<double> Fn1_d; // lookup table for Fn_derivative based on different gamma

  vector<double> Fn2; // lookup table for Fn based on different gamma
  vector<double> Fn2_d; // lookup table for Fn derivative based on different gamma

  //A vector of selected models
  vector<CandidateModels> vec_SelectedModels;
  //A vector of all candidate models
  vector<CandidateModels> vec_AllModels;
  
  //Store all the models from pr*dr 
  //Instead of p, it is estimated r.
  //Just for one site each time
  vector<rModels> vec_rModels;
  
  //Lower CI
  vector<double> vec_lower_r;
  //Upper CI
  vector<double> vec_upper_r;

   
  
  //rate by model selection
  vector<double> vec_MS_rate;
  //H=Hot Spot, C=Cold Spot, default=-
  vector<char> vec_spot;
  //rate by model averaging
  vector<double> vec_MA_rate;
  //Lower CI
  vector<double> vec_lower_rate;
  //Upper CI
  vector<double> vec_upper_rate;
  
  //Store the divergence time and r for each site.
  vector<double> vec_time;
  vector<double> vec_r;
  
  
  vector<double> divergent_time;
  vector<double> divergent_time_sums;
  vector<double> divergent_time_weights;
  
  vector<double> rWeightSums;
  vector<double> pWeightSums;
  vector<double> RandomModel_NumFastInit(struct SiteModels *p,long site);
  long RandomModel_NumFast(vector <double> &pWeightSums);

  //ps/pr/ds/dr>1 && flag_found==0, it means it could not reject the null model
  //It means when cs=0, ce=N-1, then cri==cri0, and cri could never be less than cri0 (0-N-1,para=0)
  int flag_found_pr;
  int flag_found_dr;
  int flag_found_ps;
  int flag_found_ds;
  
  int modelAveraged_p_gamma;

  vector<CandidateModels> vec_SelectedModels_ps;
  vector<CandidateModels> vec_SelectedModels_pr;
  vector<CandidateModels> vec_SelectedModels_ds;
  vector<CandidateModels> vec_SelectedModels_dr;
  
  vector<CandidateModels> vec_AllModels_ps;
  vector<CandidateModels> vec_AllModels_pr;
  vector<CandidateModels> vec_AllModels_ds;
  vector<CandidateModels> vec_AllModels_dr;
  
  
  vector<double> vec_MS_rate_ps;
  vector<double> vec_MS_rate_pr;
  vector<double> vec_MS_rate_ds;
  vector<double> vec_MS_rate_dr;

  vector<double> vec_MA_rate_ps;
  vector<double> vec_MA_rate_pr;
  vector<double> vec_MA_rate_ds;
  vector<double> vec_MA_rate_dr;
  
  vector<double> vec_lower_rate_ps;
  vector<double> vec_lower_rate_pr;
  vector<double> vec_lower_rate_ds;
  vector<double> vec_lower_rate_dr;
  
  vector<double> vec_upper_rate_ps;
  vector<double> vec_upper_rate_pr;
  vector<double> vec_upper_rate_ds;
  vector<double> vec_upper_rate_dr;

  //For estimating Neutrality Index
  vector<double> vec_NI;

 protected:
  //Polymorphic sequences: within species
  string pol_seqfile;
  vector<string> pol_seqname;
  vector<string> pol_seq;
  string pol_codon_consensus;
  
  //Divergent sequences: between species
  string div_seqfile;
  vector<string> div_seqname;
  vector<string> div_seq;	
  string div_codon_consensus;

  int output_format_num;
  int input_format_num;
  int species_num;
  int Model_Number;
  int model_num_flag;
  int verbose;
  //Choose the genetic code for this species
  int genetic_code;

  //Criterion, 0=BIC, 1=AIC, 2=AICc, 3=LRT
  int criterion_type;

  //Confidence interval for model averaging, default=0.95
  double confidence_interval;
  //Only model selection, default=0
  int MS_only;


  //Calculate 95% confidence intervals for Model Averaging, default 0
  int ci_ma;
  //Estimate selection coefficient, default 1
  int r_estimate;
  //Calculate 95% confidence intervals for r, default 1
  int ci_r;
  //Algorithm for calculating 95% confidence interval for r, default 0
  int ci_r_exact;
  

  //Default is 0 without showing the clustering results of synonymous from polymorphism and divergency sequences
  //User could choose 1 to show the clustering results of synonymous.
  int Sys_cluster;
  
  //User could input their own species divergence time
  double Div_time;
  int divtime_flag;
  int site_specific_flag;
  //Default is 1 which means to replace ambiguous nucleotide (N, R, Y, etc.) with the most frequently used nucleotide in other sequences
  //Otherwise is 0, see this codon as a gap
  int Nuc_replace;


  //Defaul is 0, do not estimate Neutrality Index=(PR/DR)/(PS/DS)
  int NI_estimate;

 private:
};

#endif


