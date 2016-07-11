/*
To change: for direct input of RS,
for krunching, need to increase/decrease the prob according due to krunching. 
Output format option: default amino acid, alternative nucleotide.
 *  To Remove for the final version: in Run, add srand(1234);
* ! Added in Run, srand(1234); // to fix random number generator, to make sure the program is repeatable; need to remove for the final version of the program
12/13/2013 To modify: time function; CIs_PRF, Site_PRF, to split to three SitePRF, gamma, CIs_gamma; output amino acid, nucleotide to two functions.
01/06/2014 Gauss–Seidel method, instead of Newton's method

! 01/09/2014 Removed gamma min_dx output in SitePRF
! 01/09/2014 Fixed bug: parameter k in ClusterSubSeq; null model, para=1, all others, para=2.Since parameters represent p0 and pc, estimated based on the model of cs and ce.
  AIC0=AICc0=BIC0=-2*InL0; will be changed to   AIC0=AICc0=BIC0=-2*InL0+2;
      if (cs==pos_start && ce==pos_end) para = 0.0; to   if (cs==pos_start && ce==pos_end) para = 1;
      else if(pos_end==ce || cs==pos_start) para=1.0; to else if(pos_end==ce || cs==pos_start) para=2;

! 01/06/2014 Modification: use gamma model average, still keep SitePRF, but model averaged gamma based on all gamma models is only calculation when CI_r=1, and model averaged gamma will replace the gamma based on model averaged pr and dr.
! 01/06/2014 Modification: gamma calcuation in SitePRF and CIr_PRF, use r_optimal when flag_root==false and dx_min<MinDx; decide the value for the criteria MinDx for dx_min, set it as 0.4 based on the output

*! 12/13/2013 Modification:  IR_H 50 IR_L -50 instead 10, -10, in .h
*! 12/13/2013 Modification: use SitePRF for gamma calculation, using the model averaged pr and dr. Remove the model medium gamma or model averaged gamma. [further revised 01/06/2014]
*! 12/13/2013 Modification: divergent_time is used directly in the output, instead of using vec_time[i]=divergent_time; in SitePRF.
*! 12/13/2013 Modification: add the option: pr=0 & dr=0, NULL, vec_r[i]=-66; in SitePRF and CIs_PRF.

! 12/10/2013 Fixed Bug - CIr_exact inside CIr_stochastic, wrong, it should be as original, since it is for the specific site only, not all sites.
! 12/09/2013 Model medium gamma, instead of model average gamma

! 12/04/2013 Modification - .h add INF and NINF for gamma values. #define INF 50 // Maximum gamma  #define NINF -50 // Minimum gamma
! 12/04/2013 Fixed Bug - polymorphism species number for the input consensus. add the parameter -SN, and species_num for the case of -ic, input consensus sequences.
! 12/04/2013 Revision - 95% CIr, for the borders, use discrete way to get lower, upper and model averaged gamma. vec_upper_r[site] = vec_rModels[i].r; averaged_gamma+=vec_rModels[i].r *(uci-rWeightSums[i-1]);
! 12/04/2013 Fixed Bug - 95% CIr, add twice for the lower_model, add twice for the uci (flag_lower==1 && i>lower_model && rWeightSums[i-1]<=uci && rWeightSums[i]>uci), not i>=lower_model, or rWeightSums[i]>=uci

! 12/03/2013 Modification - added the option of input format, using consensus sequences directly. input_format_num in the parseParameters, add readFastaConsensus in base.cpp and base.h; add options of input_format_num in Run and RunML
! 12/03/2013 Fixed Bug - Stop running SitePRF in RunML; Add vec_r[i]=-0 for the case pr=0 dr=0 together with 95%CIr in CIr_stochastic;
! 12/03/2013 Fixed Bug - CIs_PRF, model averaged gamma for infinite positive/negative gamma. pr=0, r=50, dr=0, r=-50; Use the maximum and minimum gamma in the lookup table, 50 instead of 99. The gamma value is used in the model averaged calculation. Keep the original code CIr_stochastic, if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0); not if(pr[i].sms[model_pr].p==0.0 or dr[i].sms[model_dr].p==0.0)
! 12/03/2013 Fixed Bug - CI_UpLow_rc, for averaged gamma, for the border of lower and upper 95% CI gamma, added the next or previous model twice. Change code from {//averaged_gamma+=vec_rModels[i-1].weight*vec_rModels[i-1].r+vec_upper_r[site] *(uci-rWeightSums[i-1]);} to { averaged_gamma+=vec_upper_r[site] *(uci-rWeightSums[i-1]);}
! 12/03/2013 Fixed Bug - ModelAveraging, use i instead of j in keeping the loglikelihood value for each model j for the site i. Change Code from {//double site_likelihood=vec_AllModels[i].LogLikelihood;} to { double site_likelihood=vec_AllModels[j].LogLikelihood;}

12/03/2013 Not Bug - RandomModel_NumFastInit LogLikelihood sum, not weight sum. - change back, right originally, should use weight, but not loglikelihood, 0 to 1

! 12/02/0213 Fixed Bug - log likelihood in CIr_exact and CIr_stochastic. BinomialProb already performed log for the likelihood, so use the log likelihood directly in the gamma weight calculation, no log any more. Change From {double prob_tmp=log(pr[i].sms[model_pr].LogLikelihood)+log(dr[i].sms[model_dr].LogLikelihood);} To {double log_prob_tmp=pr[i].sms[model_pr].LogLikelihood+dr[i].sms[model_dr].LogLikelihood;}


! 11/26/2013 Modification - always run ModelAveraging in ClusterSubSeq, since site specific models are summarized in ModelAveraging.
! 11/26/2013 Fixed Bug - gamma weight calculation. AIC weight recalculation for gamma using both polymorphism and divergence likelihood, but not the mutation rate;
     !added likelihood in the struct CI in PRFCluster.h and its assignment in ModelAveraging, in class CandidateModels in PRFCluster.h and its assignment in ClusterSubSeq
	 !Changed CIr_stochastic and CIr_exact by using likelihood, instead of mutation rate. double prob_tmp=log(pr[i].sms[model_pr].p)+log(dr[i].sms[model_dr].p);
! 11/26/2013 Fixed Bug - model averaged gamma, it is calculated based on summarized gamma weight from 0.025 to 0.975, so the final model_averaged gamma should divide 0.95 in CI_UpLow_rc.
! 11/26/2013 Fixed Bug - lower CI gamma and model averaged gamma. In CI_UpLow_rc, for the first model with the weight over 0.025 (lower 95%), I am calculating y intercept by using both model one and model two, and further get exactly 95% CI 0.025 border (interpolation) for lower CI gamma value, and model averaged gamma;
! 11/26/2013 Modification - simply using CIr_exact(dr, N) in CIr_stochastic; Inside CIr_stochastic, use CIr_exact while the model number less than N_Random, instead of 100


! 11/19/2013 Fixed bug - gamma AIC weight the number of parameters k, should be 4 instead of 2 in CIr_exact, CIr_stochastic, -2log(Likelihood)+2k, as 2 parameters for ClusterSubSeq.
*! 11/16/2013 Compiling issues with OS X 10.9 system
* !!! Fixed bug due to the new OS X 10.9 system [error: variable length array of non-POD element type 'struct SiteModels'], use a very large number instead of the parameter N for the gene length, for keeping all models for each gene site. struct SiteModels sm_div[10000]; struct SiteModels sm_pol[10000]; if (N>10000) { cout<<"The length of gene exceeds 10000, Revise the SiteModels upboundary array size!"<<endl; throw 1;}
* !!! Fixed compiling issue due to new OS X 10.9 system: use c++ compiler instead of gcc compiler for the new updated MAC OS X 10.9 system. Compile Command: c++ -O3 -o cMAC-PRF_test cMACprf.cpp cPRFCluster.cpp base.cpp -w
*
* 11/16/2013: model average gamma and 95% CI gamma. BUG: The model averaged gamma should not be the gamma corresponding to the model averaged probability; it should be gamma averaged for all gamma weights.  Model averaged gamma is the sum of all model specific gamma multiplying gamma AIC weight calculated from the two likelihood, polymorphism and divergence probablity.
* Fixed Bug - gamma weight calculation in CIr_stochastic and CIr_exact, AIC calculation should use pr and dr probability instead of weight. double prob_tmp=log(pr[i].sms[j].weight)+log(dr[i].sms[k].weight);  double weight_tmp=-2*prob_tmp+2*parameters;
* Fixed Bug - model averaged gamma: deleted the function of SitePRF to calculate the model averaged gamma; both model averaged gamma and 95% CI gamma in CI_UpLow_r are calculated in CI_UpLow_r.
*
*
*
* 10/20/2013 Add an option of having amino acid as the default output, add mutation status in the output, revise the output position starting from 0 to starting from 1.
* ! Revision - added an option for output format using amino acid as default: changed to output amino acid (Average three codon nucleotide sites when output) as default -o 0, nucleotide output as an alternative option -o 1. Codes changed in parseParameter and showHelpInfo
	* ! Added output_format_num in output(), for amino acid default output, and alternative nucleotide output; Bug fixed: div_codon_consensus contains "_" in the sequence in addition to "*RS"; so program quit using else throw if not seeing *RS in output().
	* ! Added output_flag in parseParameter() and showHelpInfo(); Added output_format_num in PRFCluster.h
	* ! Converted the start and end positions for the clusters of ps, pr, ds, dr, the pos_start and pos_end; and cs, ce in each cluster according to the the option output_format_num of amino acid or nucleotide output in output(), {cout<<vec_SelectedModels_dr[i].pos_start<<" ~ "<<vec_SelectedModels_dr[i].pos_end;} {cout<<"\tcs= "<<vec_SelectedModels_ps[i].cs<<"\tce= "<<vec_SelectedModels_ps[i].ce;}
* ! Revision - added the output file format by adding two columns of the PolymorphismMutationStatus and DivergenceMutationStatus as Replacement (R) or Synonymous (S)
* ! Fixed Bug - one position deviation for the output. Converted output positions to the right position in the sequence for clusters' pos_start, pos_end, cs, ce, and site position, plus one to make it the real position in the sequence.


v1.2: Major revisions
1) For time and memory saving, Sort instead of bubble sort in CIr_stochastic, sort(dr[i].sms.begin(), dr[i].sms.end(), more_than_CI()); Added Struct more_than_CI in cPRFCluster.h
2) Calculate for each site of nucleotide, instead of three sites in one codon, for purpose of krunching.
* ! SitePRF, CI_MA,CIr_stochastic, CIr_exact, SiteModels, ClusterSubSeq, ModelAveraging, Changed from N/CODONSIZE to N, to calculate each nucleotide site, instead of only the first nucleotide for each codon
3) For time and memory saving, Stochastic probabilities for gamma. Added CI_UpLow_r in .cpp and .h, use it instead of CI_UpLow_r_stochastic and CI_UpLow_r_exact separately.
* ! CI_UpLow_r: rank by gamma, partial sum of weight list, get the summed weight 2.5% and 97.5% - position of the model for the gamma CI
* ! Added int cPRFCluster::CI_UpLow_r(long site,double min_weight), to exchange codes, CI_UpLow_r_exact and CI_UpLow_r_stochastic; to make them consistent
4) Remind the requirement of lookup table in the running folder. Quit program immediately with errors if not seeing the LookupTable. Added in GammaLookupTable, if (!myfileFn2) throw "Error in opening LambdaCILookupTable...\n";
5) Added time monitor for each major step
6) Changed, CIr_stochastic, replace RandomModel_Num with RandomModel_NumFast
7) Save memory. Changed from double to float for struct CI
*

v1.1: Major revision of C++ codes from MACPRF version 1.0:
1) In PRFCluster::Run in 'PRFCluster.cpp', add codes to open and read the lookup tables for the specific number of sequences in the polymorphism data (n), with the subroutine PRFCluster::GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn)'
2) Add a subroutine 'PRFCluster::GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn)' in 'MACprf.cpp', to extract vectors Fn1,Fn2,Fn1_d, Fn2_d for calculating values of gx1, gx2, gx1_d, gx2_d for gamma calculation.
3) Revise codes accordingly for gamma calculation in the subroutines 'PRFCluster::SitePRF' and 'PRFCluster::CIs_PRF' in 'MACprf.cpp'; for gamma values within the gamma range but not recorded exactly as the interval in the Lookup Table, corresponding values are calculated linearly based on the fraction of the values between two consecutive gamma available in the Lookup Table.  
4) In PRFCluster.h, add vectors Fn1,Fn2,Fn1_d, Fn2_d and GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn)

Major bugs fixed from version 1 in subroutines of 'SitePRF' and 'CIs_PRF' in 'PRFCluster.cpp'  for the gamma calculation: 
1) Add initial_r=rtn; and rtn=initial_r inside the 'for' loop, since once the code jumps out from the inside JMAX loop, it can get to the original gamma with the assigned gamma interval increment in the rtn<=IR_H loop. Without the assignment of initial_r, the bug is caused by the new calculated gamma passing from JMAX loop outside the gamma scope between IR_L and IR_H to quit without visiting all potential gamma from IR_L to IR_H with the assigned gamma interval, or to get stuck in the dead loop with gamma always within the IR_L and IR_H scope without the assigned gamma interval increment. [Note: IR_L means the lower boundary of gamma; IR_H means the upper boundary of gamma]
2) Add r_optimal and min_dx, and a series of other related parameters for gamma calculation, since an optimal gamma with the smallest dx can be passed and saved if no gamma value qualified the criteria of ER (Expected Error 0.001). 

Updated by Zi-Ming Zhao, 10/20/2013
*
 */
#include "PRFCluster.h"


PRFCluster::PRFCluster() {
  flag_N_pol=0;
  flag_N_div=0;

  confidence_interval=0.95;
  confidence_interval=(1.0-confidence_interval)/2.0;

  divergent_time=0.0;
  
  output_format_num=0;
  input_format_num=0;
  species_num=0;
  flag_found_pr=0;
  flag_found_dr=0;
  flag_found_ps=0;
  flag_found_ds=0;

  //Initiate parameters
  pol_seqfile = "";
  div_seqfile = "";
  genetic_code = 1;
  criterion_type = 0;
  Sys_cluster=0;
  MS_only=0;
  ci_ma=0;
  Div_time=0.0;
  r_estimate=1;
  ci_r=1;
  ci_r_exact=0;
  Nuc_replace=1;
  NI_estimate=0;
}

PRFCluster::~PRFCluster() {
  pol_seq.clear();
  pol_seqname.clear();
  div_seq.clear();
  div_seqname.clear();


  vec_rModels.clear();
  vec_lower_r.clear();
  vec_upper_r.clear();


  vec_time.clear();
  vec_r.clear();
  vec_NI.clear();

  vec_SelectedModels.clear();
  vec_AllModels.clear();
  vec_MA_rate.clear();
  vec_lower_rate.clear();
  vec_upper_rate.clear();


  vec_lower_rate_ps.clear();
  vec_lower_rate_pr.clear();
  vec_lower_rate_ds.clear();
  vec_lower_rate_dr.clear();
  
  vec_upper_rate_ps.clear();
  vec_upper_rate_pr.clear();
  vec_upper_rate_ds.clear();
  vec_upper_rate_dr.clear();

  vec_MA_rate_ps.clear();
  vec_MA_rate_pr.clear();
  vec_MA_rate_ds.clear();
  vec_MA_rate_dr.clear();
  
  vec_MS_rate_ps.clear();
  vec_MS_rate_pr.clear();
  vec_MS_rate_ds.clear();
  vec_MS_rate_dr.clear();
  
  vec_SelectedModels_ps.clear();
  vec_SelectedModels_pr.clear();
  vec_SelectedModels_ds.clear();
  vec_SelectedModels_dr.clear();
  
  vec_AllModels_ps.clear();
  vec_AllModels_pr.clear();
  vec_AllModels_ds.clear();
  vec_AllModels_dr.clear();
}

// Log likelihood of Bernoulli distribution. BernoulliProb= (i/n)^i*[(n-i)/n]^(n-i); prob=Log(BernoulliProb)=i*log(i/n)+(n-i)*log[(n-i)/n]
double PRFCluster::BinomialProb(long n, long i) {
  double prob = 0.0;
  prob += (i==0)?0:i*log(double(i)/n);
  prob += (n-i==0)?0:(n-i)*log(double(n-i)/n);
  return prob;
}

double PRFCluster::getp0pc_MK(int pos_start, int pos_end, int cs, int ce, double &p0, double &pc, int nw, int nc) {  
  //p0 means the whole sequence except the center part.
  int non_cent_len=pos_end-pos_start-ce+cs;
  if(cs==pos_start && ce==pos_end){
    p0=0.0;
  }else{
    p0=(double)(nw-nc)/non_cent_len;
  }

  //pc means the center part.
  int cent_len=ce-cs+1;
  pc=(double)nc/cent_len;

  return 1;
}


int PRFCluster::Run(int argc, const char*argv[]) {	
  int i, flag=1;
  ////srand(1234); // to fix random number generator, to make sure the program is repeatable; need to remove for the final version of the program

  try {		
    //Parse input parameters
    if (parseParameter(argc, argv)!=1) throw 1;
    
    //Check input file names
    if(pol_seqfile=="" || div_seqfile=="") throw 1;
    
    cout<<endl<<NAME<<", Version: "<<VERSION<<" ["<<LASTUPDATE<<"]"<<endl;	
    cout<<"Reference: "<<REFERENCE<<endl<<endl;
    
    static time_t time_start = time(NULL);
    long p_species_n=0;
    //Read input seqeunces: polymorphic and divergent sequences
    if (input_format_num==0){

    if (readFasta(pol_seqfile, pol_seqname, pol_seq)!=1) throw "Error in reading polymorphic sequences.";
    if (readFasta(div_seqfile, div_seqname, div_seq)!=1) throw "Error in reading divergent sequences.";
    p_species_n=pol_seq.size();
    }

    else if (input_format_num==1){
        if (readFastaConsensus(pol_seqfile, pol_seqname, pol_seq)!=1) throw "Error in reading polymorphic sequences.";
        if (readFastaConsensus(div_seqfile, div_seqname, div_seq)!=1) throw "Error in reading divergent sequences.";
        if (species_num!=0) {
        	p_species_n=species_num;
        }
        else {

        	cout<< "Error! The number of polymorphism sequence number has to be input when using consensus sequence!\n";
        	throw 1;
        }

    }

    //Print the sequences and their names
    cout<<endl<<"Criterion used: "<<CRI(criterion_type)<<endl<<endl;

    cout<<endl<<"Polymorphism data:"<<endl;
    cout<<endl<<">";
    for(int i=0; i<pol_seqname.size();i++){
      cout<<pol_seqname[i].c_str()<<"/";
    }
    cout<<endl;
    for(int i=0; i<pol_seq.size();i++){
      cout<<pol_seq[i].c_str()<<endl;
    }
    cout<<endl;
    cout<<endl<<"Divergence data:"<<endl;
    cout<<endl<<">";
    for(int i=0; i<div_seqname.size();i++){
      cout<<div_seqname[i].c_str()<<"/";
    }
    cout<<endl;
    for(int i=0; i<div_seq.size();i++){
      cout<<div_seq[i].c_str()<<endl;
    }
    cout<<endl<<endl;		    
    
    
    //Check sequences' length
    for (i=0; i<pol_seq.size(); i++) {
      if (div_seq[0].length()!=pol_seq[i].length()) throw "Input sequences are not equal in length.";
    }
    if(pol_seq[0].size()%3!=0) cout<<"Warning: the length of sequences can be divided by 3 (codon size).\n";
    	//throw "Error: the length of sequences can be divided by 3 (codon size).";


    //Open and Read the lookup table for numeric integration for the specific number of species
     //special attention to gamma=0,purposely put the corresponding numeric integration to 0; treat separately for this case.
       cout<<endl<<"Open and Read the lookup table for numeric integration:"<<endl;

       cout<<endl<<"The number of species in the polymorphism data:"<<p_species_n<<endl;
       int tmp_n=p_species_n; //the number of species in the polymorphism data, according to the line in the lookup table
       if (p_species_n>100 and p_species_n<=1500)
       {
    	   tmp_n=98+int(p_species_n/50); //for n>100, the interval is 50 in the lookup table
       }
       else if (p_species_n>1500)
       {
    	   cout<<"Error! The number of sequences in the polymorphism data is over the range in the Lookup Table 1500."<<endl;
       }
       cout<<"The row to be extracted in the lookup table is "<<tmp_n<<endl;
       

  
pol_codon_consensus = getPolSysRep(pol_seq);
 div_codon_consensus = getDivSysRep(pol_seq, div_seq);
 long species_n=pol_seq.size();

 if (input_format_num==1){
   pol_codon_consensus = pol_seq[0];
   div_codon_consensus = div_seq[0];
   species_n=species_num;
 }

 cout<<"Polymorphism: "<<pol_codon_consensus<<endl;
 cout<<"Divergence:   "<<div_codon_consensus<<endl;
 cout<<endl;




  


       //test the sucess of assignment of values in the vectors of Fn and Fn_d
       //cout<<"Test Fn1[0]: "<<Fn1[0]<<endl;
       //cout<<"Test Fn1[86]: "<<Fn1[86]<<endl;
       //cout<<"Test Fn1_d[0]: "<<Fn1_d[0]<<endl;

    //Run the main function
   // RunML(pol_seq, div_seq);
    
    //Display on screen
    cout<<endl<<"Mission accomplished. (Time elapsed: ";
    time_t t = time(NULL)-time_start;
    int h=t/3600, m=(t%3600)/60, s=t-(t/60)*60;
    if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
    else   cout<<m<<":"<<s<<")"<<endl;	
    
  }	
  catch (const char* e) {
    cout<<e<<endl;		
    flag = 0;
  }
  catch (...) {		
    flag = 0;
  }
  
  return flag;
}
/*************************************************************/





/*************************************************************/
/*
Subfunction - open and read the lookup table, to extract the values for the certain species number.
Input: the lookup table file name, p_specie_n converted tmp_n;
Output: vector<double> Fn
*/
vector<double> PRFCluster::GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn){

   double t_value; // Temporary value
   string t_Fn2; // Temporary Fn
   int line_n=0; // the number of the item in the file
   double t_r=-199; //the Temporary gamma value
   int ii=0;

   ifstream myfileFn2(input_f_name.c_str());
   if (!myfileFn2) throw "Error in opening GammaLookupTable...\n";


   if (myfileFn2.is_open())
   {
	   //t_r=-199;
	   line_n=0;
	   ii=0;

	   while ( myfileFn2.good() )
		  {
		   if (line_n>tmp_n*202 and line_n<(tmp_n+1)*202) // the number 202 is based on the number of gamma values in the lookup table, from -50 to 50 in the interval of 0.5, plus the n description column.
			  {
				myfileFn2 >>t_value;
				t_r=(ii-100.0)/2.0;
				//cout<<"***r"<<t_r<<":"<<t_value<<"***\t"; // the values in the table for specific gamma
				Fn[ii]=t_value;
				ii+=1;

			  }
				else
				{
					myfileFn2 >>t_Fn2;
				}
				line_n=line_n+1;
			   }

	 myfileFn2.close();
   }
   else
	   {
	   cout << "Unable to open file"<<endl;
	   }
   cout <<endl<<"Vector size after extraction: "<<Fn.size()<<endl;
   //cout<<"Test Fn[0]: "<<Fn[0]<<endl;
   return Fn;
}



int PRFCluster::init(long N){
  vec_time.resize(N,0.0);
  vec_r.resize(N,0.0);
  vec_NI.resize(N,0.0);

  vec_MS_rate_ps.resize(N,0.0);
  vec_MS_rate_pr.resize(N,0.0);
  vec_MS_rate_ds.resize(N,0.0);
  vec_MS_rate_dr.resize(N,0.0);

  vec_MA_rate_ps.resize(N,0.0);
  vec_MA_rate_pr.resize(N,0.0);
  vec_MA_rate_ds.resize(N,0.0);
  vec_MA_rate_dr.resize(N,0.0);
  
  vec_lower_rate_ps.resize(N,0.0);
  vec_lower_rate_pr.resize(N,0.0);
  vec_lower_rate_ds.resize(N,0.0);
  vec_lower_rate_dr.resize(N,0.0);

  vec_upper_rate_ps.resize(N,0.0);
  vec_upper_rate_pr.resize(N,0.0);
  vec_upper_rate_ds.resize(N,0.0);
  vec_upper_rate_dr.resize(N,0.0);

  vec_SelectedModels_ps.clear();
  vec_SelectedModels_pr.clear();
  vec_SelectedModels_ds.clear();
  vec_SelectedModels_dr.clear();

  vec_AllModels_ps.clear();
  vec_AllModels_pr.clear();
  vec_AllModels_ds.clear();
  vec_AllModels_dr.clear();


  vec_SelectedModels.clear();
  vec_AllModels.clear();
  vec_MA_rate.resize(N,0.0);
  vec_lower_rate.resize(N,0.0);
  vec_upper_rate.resize(N,0.0);


  vec_rModels.clear();
  vec_lower_r.resize(N,0.0);
  vec_upper_r.resize(N,0.0);


  return 1;
}
/*************************************************************/



/*************************************************************/
int PRFCluster::output(long N){
  
  cout<<endl<<"//Results based on model selection: "<<endl;
  if(Sys_cluster==1){
    cout<<endl<<"Clusters from Polymorphism Synonymous:"<<endl;
    if(vec_SelectedModels_ps.size()==0){
      cout<<"Note: PS=1 or 0. There is not enough information for clustering!"<<endl<<endl;
    }else if(vec_SelectedModels_ps.size()==1 && vec_SelectedModels_ps[0].pos_start==vec_SelectedModels_ps[0].cs && vec_SelectedModels_ps[0].pos_end==vec_SelectedModels_ps[0].ce){
      cout<<"Note: There is no cluster in this sequences for Polymorphism Synonymous"<<endl<<endl;
    }else{
      for(long i=0; i<vec_SelectedModels_ps.size(); i++){

      	//Converted the start and end positions for the cluster, the pos_start, pos_end, cs, ce according to the the option output_format_num of amino acid or nucleotide output in output()
    	  if (output_format_num==1)
    	{
    		cout<<(vec_SelectedModels_ps[i].pos_start+1)<<" nucleotide ~ "<<(vec_SelectedModels_ps[i].pos_end+1)<<" nucleotide";
    		cout<<"\tcs= "<<(vec_SelectedModels_ps[i].cs+1)<<"\tce= "<<(vec_SelectedModels_ps[i].ce+1);
    	}
    	else if (output_format_num==0)
    	{
    		cout<<(vec_SelectedModels_ps[i].pos_start/3+1)<<" amino acid ~ "<<(vec_SelectedModels_ps[i].pos_end/3+1)<<" amino acid";
    		cout<<"\tcs= "<<(vec_SelectedModels_ps[i].cs/3+1)<<"\tce= "<<(vec_SelectedModels_ps[i].ce/3+1);
    	}
		//cout<<vec_SelectedModels_ps[i].pos_start<<" ~ "<<vec_SelectedModels_ps[i].pos_end;
		//cout<<"\tcs= "<<vec_SelectedModels_ps[i].cs<<"\tce= "<<vec_SelectedModels_ps[i].ce;
		cout<<endl;

		cout<<"InL0= "<<vec_SelectedModels_ps[i].InL0<<"\tInL= "<<vec_SelectedModels_ps[i].InL;
		cout<<"\tAIC0= "<<vec_SelectedModels_ps[i].AIC0<<"\tAIC= "<<vec_SelectedModels_ps[i].AIC;
		cout<<"\tAICc0= "<<vec_SelectedModels_ps[i].AICc0<<"\tAICc= "<<vec_SelectedModels_ps[i].AICc;
		cout<<"\tBIC0= "<<vec_SelectedModels_ps[i].BIC0<<"\tBIC= "<<vec_SelectedModels_ps[i].BIC;
		cout<<endl;

		cout<<"P0ps= "<<vec_SelectedModels_ps[i].p0<<"\tPcps= "<<vec_SelectedModels_ps[i].pc;
		cout<<endl<<endl;
      }
    }
  }

  cout<<endl<<"Clusters from Polymorphism Replacement:"<<endl;
  if(vec_SelectedModels_pr.size()==0){
    cout<<"Note: PR=1 or 0. There is not enough information for clustering!"<<endl<<endl;
  }else if(vec_SelectedModels_pr.size()==1 && vec_SelectedModels_pr[0].pos_start==vec_SelectedModels_pr[0].cs && vec_SelectedModels_pr[0].pos_end==vec_SelectedModels_pr[0].ce){
    cout<<"Note: There is no cluster in this Sequences for Polymorphism non-Synonymous"<<endl<<endl;
  }else{
    for(long i=0; i<vec_SelectedModels_pr.size(); i++){

    	//Converted the start and end positions for the cluster, the pos_start, pos_end, cs, ce according to the the option output_format_num of amino acid or nucleotide output in output()
    	if (output_format_num==1)
    	{
    		cout<<(vec_SelectedModels_pr[i].pos_start+1)<<" nucleotide ~ "<<(vec_SelectedModels_pr[i].pos_end+1)<<" nucleotide";
    		cout<<"\tcs= "<<(vec_SelectedModels_pr[i].cs+1)<<"\tce= "<<(vec_SelectedModels_pr[i].ce+1);
    	}
    	else if (output_format_num==0)
    	{
    		cout<<(vec_SelectedModels_pr[i].pos_start/3+1)<<" amino acid ~ "<<(vec_SelectedModels_pr[i].pos_end/3+1)<<" amino acid";
    		cout<<"\tcs= "<<(vec_SelectedModels_pr[i].cs/3+1)<<"\tce= "<<(vec_SelectedModels_pr[i].ce/3+1);
    	}
      //cout<<vec_SelectedModels_pr[i].pos_start<<" ~ "<<vec_SelectedModels_pr[i].pos_end;
      //cout<<"\tcs= "<<vec_SelectedModels_pr[i].cs<<"\tce= "<<vec_SelectedModels_pr[i].ce;
      cout<<endl;
      
      cout<<"InL0= "<<vec_SelectedModels_pr[i].InL0<<"\tInL= "<<vec_SelectedModels_pr[i].InL;
      cout<<"\tAIC0= "<<vec_SelectedModels_pr[i].AIC0<<"\tAIC= "<<vec_SelectedModels_pr[i].AIC;
      cout<<"\tAICc0= "<<vec_SelectedModels_pr[i].AICc0<<"\tAICc= "<<vec_SelectedModels_pr[i].AICc;
      cout<<"\tBIC0= "<<vec_SelectedModels_pr[i].BIC0<<"\tBIC= "<<vec_SelectedModels_pr[i].BIC;
      cout<<endl;

      cout<<"\tP0pr= "<<vec_SelectedModels_pr[i].p0<<"\tPcpr= "<<vec_SelectedModels_pr[i].pc;
      cout<<endl<<endl;
    }
  }


  if(Sys_cluster==1){
    cout<<endl<<"Clusters from Divergence Synonymous:"<<endl;
    if(vec_SelectedModels_ds.size()==0){
      cout<<"Note: DS=1 or 0. There is not enough information for clustering!"<<endl<<endl;
    }else if(vec_SelectedModels_ds.size()==1 && vec_SelectedModels_ds[0].pos_start==vec_SelectedModels_ds[0].cs && vec_SelectedModels_ds[0].pos_end==vec_SelectedModels_ds[0].ce){
      cout<<"Note: There is no cluster in this Sequences for Divergence Synonymous"<<endl<<endl;
    }else{
      for(long i=0; i<vec_SelectedModels_ds.size(); i++){

     	//Converted the start and end positions for the cluster, the pos_start, pos_end, cs, ce according to the the option output_format_num of amino acid or nucleotide output in output()
    	  if (output_format_num==1)
      	{
      		cout<<(vec_SelectedModels_ds[i].pos_start+1)<<" nucleotide ~ "<<(vec_SelectedModels_ds[i].pos_end+1)<<" nucleotide";
      		cout<<"\tcs= "<<(vec_SelectedModels_ds[i].cs+1)<<"\tce= "<<(vec_SelectedModels_ds[i].ce+1);
      	}
      	else if (output_format_num==0)
      	{
      		cout<<(vec_SelectedModels_ds[i].pos_start/3+1)<<" amino acid ~ "<<(vec_SelectedModels_ds[i].pos_end/3+1)<<" amino acid";
      		cout<<"\tcs= "<<(vec_SelectedModels_ds[i].cs/3+1)<<"\tce= "<<(vec_SelectedModels_ds[i].ce/3+1);
      	}
		//cout<<vec_SelectedModels_ds[i].pos_start<<" ~ "<<vec_SelectedModels_ds[i].pos_end;
		//cout<<"\tcs= "<<vec_SelectedModels_ds[i].cs<<"\tce= "<<vec_SelectedModels_ds[i].ce;
		cout<<endl;

		cout<<"InL0= "<<vec_SelectedModels_ds[i].InL0<<"\tInL= "<<vec_SelectedModels_ds[i].InL;
		cout<<"\tAIC0= "<<vec_SelectedModels_ds[i].AIC0<<"\tAIC= "<<vec_SelectedModels_ds[i].AIC;
		cout<<"\tAICc0= "<<vec_SelectedModels_ds[i].AICc0<<"\tAICc= "<<vec_SelectedModels_ds[i].AICc;
		cout<<"\tBIC0= "<<vec_SelectedModels_ds[i].BIC0<<"\tBIC= "<<vec_SelectedModels_ds[i].BIC;
		cout<<endl;

		cout<<"P0ds= "<<vec_SelectedModels_ds[i].p0<<"\tPcds= "<<vec_SelectedModels_ds[i].pc;
		cout<<endl<<endl;
      }
    }
  }

  cout<<endl<<"Clusters from Divergence Replacement:"<<endl;
  if(vec_SelectedModels_dr.size()==0){
    cout<<"Note: DR=1 or 0. There is not enough information for clustering!"<<endl<<endl;
  }else if(vec_SelectedModels_dr.size()==1 && vec_SelectedModels_dr[0].pos_start==vec_SelectedModels_dr[0].cs && vec_SelectedModels_dr[0].pos_end==vec_SelectedModels_dr[0].ce){
    cout<<"Note: There is no cluster in this Sequences for Divergence Replacement"<<endl<<endl;
  }else{
    for(long i=0; i<vec_SelectedModels_dr.size(); i++){
    	//Converted the start and end positions for the cluster, the pos_start, pos_end, cs, ce according to the the option output_format_num of amino acid or nucleotide output in output()
    	if (output_format_num==1)
    	{
    		cout<<(vec_SelectedModels_dr[i].pos_start+1)<<" nucleotide ~ "<<(vec_SelectedModels_dr[i].pos_end+1)<<" nucleotide";
    		cout<<"\tcs= "<<(vec_SelectedModels_dr[i].cs+1)<<"\tce= "<<(vec_SelectedModels_dr[i].ce+1);
    	}
    	else if (output_format_num==0)
    	{
    		cout<<(vec_SelectedModels_dr[i].pos_start/3+1)<<" amino acid ~ "<<(vec_SelectedModels_dr[i].pos_end/3+1)<<" amino acid";
    		cout<<"\tcs= "<<(vec_SelectedModels_dr[i].cs/3+1)<<"\tce= "<<(vec_SelectedModels_dr[i].ce/3+1);
    	}
      //cout<<vec_SelectedModels_dr[i].pos_start<<" ~ "<<vec_SelectedModels_dr[i].pos_end;
      //cout<<"\tcs= "<<vec_SelectedModels_dr[i].cs<<"\tce= "<<vec_SelectedModels_dr[i].ce;
      cout<<endl;

      cout<<"InL0= "<<vec_SelectedModels_dr[i].InL0<<"\tInL= "<<vec_SelectedModels_dr[i].InL;
      cout<<"\tAIC0= "<<vec_SelectedModels_dr[i].AIC0<<"\tAIC= "<<vec_SelectedModels_dr[i].AIC;
      cout<<"\tAICc0= "<<vec_SelectedModels_dr[i].AICc0<<"\tAICc= "<<vec_SelectedModels_dr[i].AICc;
      cout<<"\tBIC0= "<<vec_SelectedModels_dr[i].BIC0<<"\tBIC= "<<vec_SelectedModels_dr[i].BIC;
      cout<<endl;
      
      cout<<"P0dr= "<<vec_SelectedModels_dr[i].p0<<"\tPcdr= "<<vec_SelectedModels_dr[i].pc;
      cout<<endl<<endl;
    }
  }

  if (N%3!=0) { cout<< "Warning: the length of sequence can not be divided by 3 (codon size)."<<endl;}

  //Print the results of Model Averaging
  if (MS_only==0) {
    cout<<endl<<"//Results based on model averaging: "<<endl;
    cout.setf(ios::left);
    int width=15;

    //Output the headings for the data table containing the values of gamma and other relatives.
    cout.width(width); cout<<"Position\t";
    if(Sys_cluster==1){
      cout.width(width); cout<<"MS_PolSys\t";
      cout.width(width); cout<<"MA_PolSys\t";
      if(ci_ma==1){
	cout.width(width); cout<<"Lower_CI_PolSys\t";
	cout.width(width); cout<<"Upper_CI_PolSys\t";
      }
    }
    cout.width(width); cout<<"MS_PolRep\t";
    cout.width(width); cout<<"MA_PolRep\t";

    if(ci_ma==1){
      cout.width(width); cout<<"Lower_CI_PolRep\t";
      cout.width(width); cout<<"Upper_CI_PolRep\t";
    }

    if(Sys_cluster==1){
      cout.width(width); cout<<"MS_DivSys\t";
      cout.width(width); cout<<"MA_DivSys\t";
      if(ci_ma==1){
	cout.width(width); cout<<"Lower_CI_DivSys\t";
	cout.width(width); cout<<"Upper_CI_DivSys\t";
      }
    }

    cout.width(width); cout<<"MS_DivRep\t";
    cout.width(width); cout<<"MA_DivRep";
    if(ci_ma==1){
      cout.width(width); cout<<"\tLower_CI_DivRep\t";
      cout.width(width); cout<<"Upper_CI_DivRep";
    }

    if(r_estimate==1){
      cout.width(width); cout<<"\tDivergentTime\t";
      cout.width(width); cout<<"Gamma";
      if(ci_r==1){
	cout.width(width); cout<<"\tLower_CI_Gamma\t";
	cout.width(width); cout<<"Upper_CI_Gamma";
      }
    }
    if(NI_estimate==1){
      cout.width(width); cout<<"\tNI";
    }
    cout.width(width); cout<<"\tPolymorphismMutationStatus\t";
    cout.width(width); cout<<"DivergenceMutationStatus";
    cout<<endl;


    //Output the values for the data table containing the values of gamma and other relatives.
	//nucleotide/codon format output
	if (output_format_num==1)
	{
	    for(long i=0; i<N; i++) {
			  cout.width(width);cout<<i+1<<"\t";
    	      //Output synonymous polymorphism model selection and model average probabilities, and confidence interval if ci_ma==1
			  if(Sys_cluster==1){
				cout.width(width);cout<<vec_MS_rate_ps[i]<<"\t";
				cout.width(width);cout<<vec_MA_rate_ps[i]<<"\t";
				if(ci_ma==1){
				  cout.width(width);cout<<vec_lower_rate_ps[i]<<"\t";
				  cout.width(width);cout<<vec_upper_rate_ps[i]<<"\t";
				}
			  }
    	      //Output replacment polymorphism model selection and model average probabilities, and confidence interval if ci_ma==1
			  cout.width(width);cout<<vec_MS_rate_pr[i]<<"\t";
			  cout.width(width);cout<<vec_MA_rate_pr[i]<<"\t";
			  if(ci_ma==1){
				cout.width(width);cout<<vec_lower_rate_pr[i]<<"\t";
				cout.width(width);cout<<vec_upper_rate_pr[i]<<"\t";
			  }
    	      //Output synonymous replacment model selection and model average probabilities, and confidence interval if ci_ma==1
			  if(Sys_cluster==1){
				cout.width(width);cout<<vec_MS_rate_ds[i]<<"\t";
				cout.width(width);cout<<vec_MA_rate_ds[i]<<"\t";
				if(ci_ma==1){
				  cout.width(width);cout<<vec_lower_rate_ds[i]<<"\t";
				  cout.width(width);cout<<vec_upper_rate_ds[i]<<"\t";
				}
			  }
    	      //Output divergence replacement  model selection and model average probabilities, and confidence interval if ci_ma==1
			  cout.width(width);cout<<vec_MS_rate_dr[i]<<"\t";
			  cout.width(width);cout<<vec_MA_rate_dr[i];
			  if(ci_ma==1){
				cout<<"\t";cout.width(width);cout<<vec_lower_rate_dr[i]<<"\t";
				cout.width(width);cout<<vec_upper_rate_dr[i];
			  }

    	      //Output gamma value and confidence interval of gamma if ci_r==1
			  if(r_estimate==1){
				cout<<"\t";cout.width(width);cout<<divergent_time<<"\t";
				if(vec_r[i]==50){
				  cout.width(width);cout<<INF;
				}else if (vec_r[i]==-50){
				  cout.width(width);cout<<NINF;
				}

				else if (vec_r[i]==0 || vec_r[i]==-66){
				  cout.width(width);cout<<"NULL";
				}else{
				  cout.width(width);cout<<vec_r[i];
				}

				if(ci_r==1){
				  if(vec_lower_r[i]==50){
					cout<<"\t";cout.width(width);cout<<INF<<"\t";
				  }else if(vec_lower_r[i]==-50){
					cout<<"\t";cout.width(width);cout<<NINF<<"\t";
				  }

				  else if(vec_lower_r[i]==0 || vec_lower_r[i]==-66){
					cout<<"\t";cout.width(width);cout<<"NULL"<<"\t";
				  }else{
					cout<<"\t";cout.width(width);cout<<vec_lower_r[i]<<"\t";
				  }

				  if(vec_upper_r[i]==50){
				   cout.width(width);cout<<INF;
				  }else if(vec_upper_r[i]==-50){
					cout.width(width);cout<<NINF;
				  }

				  else if(vec_upper_r[i]==0 || vec_upper_r[i]==-66){
					cout.width(width);cout<<"NULL";
				  }else{
					cout.width(width);cout<<vec_upper_r[i];
				  }
				}//ci_r
			  }//r_estimate


			  //Estimate Neutrality Index
			  if(NI_estimate==1){
				if(vec_NI[i]==99){
				  cout<<"\t";cout.width(width);cout<<"INF";
				}else if(vec_NI[i]==-99){
				  cout<<"\t";cout.width(width);cout<<"N-INF";
				}

				else if(vec_NI[i]==-66){
				  cout<<"\t";cout.width(width);cout<<"NULL";
				}else{
				  cout<<"\t";cout.width(width);cout<<vec_NI[i];
				}
			  }//NI

			  //Output the Replacement or Synonymous status of the mutation of polymorphism and divergence.
			  cout<<"\t";cout.width(width);cout<<pol_codon_consensus[i]<<"\t";
			  cout.width(width);cout<<div_codon_consensus[i];
			  cout<<endl;

		}
	}



	//amino acid format output
	else if (output_format_num==0)
	{
	    for(long i=0; i<N/3; i++) {
    	      cout.width(width);cout<<i+1<<"\t";
    	      //Output synonymous polymorphism model selection and model average probabilities, and confidence interval if ci_ma==1
    	      if(Sys_cluster==1){
    			cout.width(width);cout<<(vec_MS_rate_ps[i*3]+vec_MS_rate_ps[i*3+1]+vec_MS_rate_ps[i*3+2])/3<<"\t";
    			cout.width(width);cout<<(vec_MA_rate_ps[i*3]+vec_MA_rate_ps[i*3+1]+vec_MA_rate_ps[i*3+2])/3<<"\t";
    			if(ci_ma==1){
    			  cout.width(width);cout<<(vec_lower_rate_ps[i*3]+vec_lower_rate_ps[i*3+1]+vec_lower_rate_ps[i*3+2])/3<<"\t";
    			  cout.width(width);cout<<(vec_upper_rate_ps[i*3]+vec_upper_rate_ps[i*3+1]+vec_upper_rate_ps[i*3+2])/3<<"\t";
    			}
    	      }
    	      //Output replacment polymorphism model selection and model average probabilities, and confidence interval if ci_ma==1
    	      cout.width(width);cout<<(vec_MS_rate_pr[i*3]+vec_MS_rate_pr[i*3+1]+vec_MS_rate_pr[i*3+2])/3<<"\t";
    	      cout.width(width);cout<<(vec_MA_rate_pr[i*3]+vec_MA_rate_pr[i*3+1]+vec_MA_rate_pr[i*3+2])/3<<"\t";
    	      if(ci_ma==1){
    			cout.width(width);cout<<(vec_lower_rate_pr[i*3]+vec_lower_rate_pr[i*3+1]+vec_lower_rate_pr[i*3+2])/3<<"\t";
    			cout.width(width);cout<<(vec_upper_rate_pr[i*3]+vec_upper_rate_pr[i*3+1]+vec_upper_rate_pr[i*3+2])/3<<"\t";
    	      }

    	      //Output synonymous replacment model selection and model average probabilities, and confidence interval if ci_ma==1
    	      if(Sys_cluster==1){
    			cout.width(width);cout<<(vec_MS_rate_ds[i*3]+vec_MS_rate_ds[i*3+1]+vec_MS_rate_ds[i*3+2])/3<<"\t";
    			cout.width(width);cout<<(vec_MA_rate_ds[i*3]+vec_MA_rate_ds[i*3+1]+vec_MA_rate_ds[i*3+2])/3<<"\t";
    			if(ci_ma==1){
    			  cout.width(width);cout<<(vec_lower_rate_ds[i*3]+vec_lower_rate_ds[i*3+1]+vec_lower_rate_ds[i*3+2])/3<<"\t";
    			  cout.width(width);cout<<(vec_upper_rate_ds[i*3]+vec_upper_rate_ds[i*3+1]+vec_upper_rate_ds[i*3+2])/3<<"\t";
    			}
    	      }
    	      //Output divergence replacement  model selection and model average probabilities, and confidence interval if ci_ma==1
    	      cout.width(width);cout<<(vec_MS_rate_dr[i*3]+vec_MS_rate_dr[i*3+1]+vec_MS_rate_dr[i*3+2])/3<<"\t";
    	      cout.width(width);cout<<(vec_MA_rate_dr[i*3]+vec_MA_rate_dr[i*3+1]+vec_MA_rate_dr[i*3+2])/3;
    	      if(ci_ma==1){
    			cout<<"\t";cout.width(width);cout<<(vec_lower_rate_dr[i*3]+vec_lower_rate_dr[i*3+1]+vec_lower_rate_dr[i*3]+2)/3<<"\t";
    			cout.width(width);cout<<(vec_upper_rate_dr[i*3]+vec_upper_rate_dr[i*3+1]+vec_upper_rate_dr[i*3+2])/3;
    	      }

    	      //Output gamma value and confidence interval of gamma if ci_r==1
    	      if(r_estimate==1){
    			cout<<"\t";cout.width(width);cout<<divergent_time<<"\t";
    			if(vec_r[i*3]==50 or vec_r[i*3+1]==50 or vec_r[i*3+2]==50){
    			  cout.width(width);cout<<"INF";
    			}else if (vec_r[i*3]==-50 or vec_r[i*3+1]==-50 or vec_r[i*3+2]==-50){
    			  cout.width(width);cout<<"N-INF";
    			}else if (vec_r[i*3]==0 || vec_r[i*3]==-66 or vec_r[i*3+1]==0 || vec_r[i*3+1]==-66  or vec_r[i*3+2]==0 || vec_r[i*3+2]==-66 ){
    			  cout.width(width);cout<<"NULL";
    			}else{
    			  cout.width(width);cout<<(vec_r[i*3]+vec_r[i*3+1]+vec_r[i*3+2])/3;
    			}

    			if(ci_r==1){
    			  if(vec_lower_r[i*3]==50 or vec_lower_r[i*3+1]==50 or vec_lower_r[i*3+2]==50){
    				cout<<"\t";cout.width(width);cout<<"INF"<<"\t";
    			  }else if(vec_lower_r[i*3]==-50 or vec_lower_r[i*3+1]==-50 or vec_lower_r[i*3+2]==-50){
    				cout<<"\t";cout.width(width);cout<<"N-INF"<<"\t";
    			  }else if(vec_lower_r[i*3]==0 || vec_lower_r[i*3]==-66 or vec_lower_r[i*3+1]==0 || vec_lower_r[i*3+1]==-66 or vec_lower_r[i*3+2]==0 || vec_lower_r[i*3+2]==-66){
    				cout<<"\t";cout.width(width);cout<<"NULL"<<"\t";
    			  }else{
    				cout<<"\t";cout.width(width);cout<<(vec_lower_r[i*3]+vec_lower_r[i*3+1]+vec_lower_r[i*3+2])/3<<"\t";
    			  }

    			  if(vec_upper_r[i*3]==50 or vec_upper_r[i*3+1]==50 or vec_upper_r[i*3+2]==50){
    			   cout.width(width);cout<<"INF";
    			  }else if(vec_upper_r[i*3]==-50 or vec_upper_r[i*3+1]==-50 or vec_upper_r[i*3+2]==-50){
    				cout.width(width);cout<<"N-INF";
    			  }else if(vec_upper_r[i*3]==0 || vec_upper_r[i*3]==-66 or vec_upper_r[i*3+1]==0 || vec_upper_r[i*3+1]==-66 or vec_upper_r[i*3+2]==0 || vec_upper_r[i*3+2]==-66){
    				cout.width(width);cout<<"NULL";
    			  }else{
    				cout.width(width);cout<<(vec_upper_r[i*3]+vec_upper_r[i*3+1]+vec_upper_r[i*3+2])/3;
    			  }
    			}//ci_r
    	      }//r_estimate

    	      //Estimate Neutrality Index
    	      if(NI_estimate==1){
    			if(vec_NI[i*3]==99 or vec_NI[i*3+1]==99 or vec_NI[i*3+2]==99){
    			  cout<<"\t";cout.width(width);cout<<"INF";
    			}else if(vec_NI[i*3]==-99 or vec_NI[i*3+1]==-99 or vec_NI[i*3+2]==-99){
    			  cout<<"\t";cout.width(width);cout<<"N-INF";
    			}else if(vec_NI[i*3]==-66 or vec_NI[i*3+1]==-66 or vec_NI[i*3+2]==-66){
    			  cout<<"\t";cout.width(width);cout<<"NULL";
    			}else{
    			  cout<<"\t";cout.width(width);cout<<(vec_NI[i*3]+vec_NI[i*3+1]+vec_NI[i*3+2])/3;
    			}
    	      }//NI

    	      //Output the Replacement or Synonymous status of the mutation of polymorphism and divergence.

    	      if (pol_codon_consensus[i*3]=='R' or pol_codon_consensus[i*3+1]=='R' or pol_codon_consensus[i*3+2]=='R') { cout<<"\tR\t";  }
    	      else if (pol_codon_consensus[i*3]=='S' or pol_codon_consensus[i*3+1]=='S' or pol_codon_consensus[i*3+2]=='S') {	cout<<"\tS\t";  }
    	      else { cout<<"\t*\t"; }
    	      //(pol_codon_consensus[i*3]=='*' and pol_codon_consensus[i*3+1]=='*' and pol_codon_consensus[i*3+2]=='*')
    	      /*else { cout <<endl<<"Error in polymorphism codon consensus output, the three codon positions are not *, R or S! "<<pol_codon_consensus[i*3]<<pol_codon_consensus[i*3+1]<<pol_codon_consensus[i*3+2]<<endl;
    	      throw 1;}
    	      */

    	      if (div_codon_consensus[i*3]=='R' or div_codon_consensus[i*3+1]=='R' or div_codon_consensus[i*3+2]=='R') { cout.width(width);cout<<'R'; }
    	      else if (div_codon_consensus[i*3]=='S' or div_codon_consensus[i*3+1]=='S' or div_codon_consensus[i*3+2]=='S') { cout.width(width);cout<<'S'; }
    	      else { cout.width(width);cout<<'*'; }
    	      //div_codon_consensus contains "_" in the sequence in addition to "*RS"; so program quit using else throw if not seeing *RS in output().
    	      /*
    	      if (div_codon_consensus[i*3]=='*' and div_codon_consensus[i*3+1]=='*' and div_codon_consensus[i*3+2]=='*') { cout.width(width);cout<<'*'; }
    	      else { cout <<endl<<"Error in divergence codon consensus output, the three codon positions are not *, R or S! "<<div_codon_consensus[i*3]<<div_codon_consensus[i*3+1]<<div_codon_consensus[i*3+2]<<endl;
    	      throw 1;}
    	      */

    	      //cout<<"\t";cout.width(width);cout<<pol_codon_consensus[i*3]<<"\t";
    	      //cout.width(width);cout<<div_codon_consensus[i*3];
    	      cout<<endl;
    	}
    }

	else { throw 1;}



    cout<<endl<<"Abbreviation: MS=Model Selection; MA=Model Averaging; CI=Confidence Interval; ps=Polymorphism Synonymous; pr=Polymorphism Replacement; ds=Divergence Synonymous; dr=Divergence Replacement; Gamma=N*s (Gamma: scaled selection coefficient (selection intensity); N: haploid effective population size; s: selection coefficient); NI=Neutrality Index (NI=(pr/dr)/(ps/ds), >1 Negative selection, <1 Positive selection); INF=Infinite; N-INF=Negative Infinite; NULL=Not enough information for this site"<<endl;
  }

  cout<<endl<<"#End of clustering"<<endl<<endl;

  return 1;
}
/*************************************************************/



/*************************************************************/
int PRFCluster::RunML(vector<string> pol_seq, vector<string> div_seq) {

  cout<<"****** Start RunML"<<endl;
  time_t time_startRunML = time(NULL); // Record the start time
  long N=pol_seq[0].length();
  init(N);

  pol_codon_consensus = getPolSysRep(pol_seq);
  div_codon_consensus = getDivSysRep(pol_seq, div_seq);
  long species_n=pol_seq.size();

  if (input_format_num==1){
	  pol_codon_consensus = pol_seq[0];
	  div_codon_consensus = div_seq[0];
	  species_n=species_num;
  }


  cout<<"Polymorphism: "<<pol_codon_consensus<<endl;
  cout<<"Divergence:   "<<div_codon_consensus<<endl;
  cout<<endl;
  if(flag_N_pol==1){
    cout<<endl<<"*****************"<<endl;
    cout<<"Note: Nucleotide other than A, T, G or C is in the polymorphism data sequences! Please see details in manual document."<<endl;
    if(Nuc_replace==1){
      cout<<"      It is substituted by the most frequent used nucleotide in polymorphism sequence at this site!"<<endl;
    }else{
      cout<<"      This codon is seen as gap!"<<endl;
    }
    cout<<"*****************"<<endl<<endl;
  }
  if(flag_N_div==1){
    cout<<endl<<"*****************"<<endl;
    cout<<"Note: Nucleotide other than A, T, G or C is in the divergence data sequence! Please see details in manual document."<<endl;
    if(Nuc_replace==1){
      cout<<"      It is substituted by the most frequent used nucleotide in polymorphism sequence at this site!"<<endl;
    }else{
      cout<<"      This codon is seen as gap!"<<endl;
    }
    cout<<"*****************"<<endl<<endl;
  }


  //Print the #Div & #Pol
  double ps=0.0;
  double ds=0.0;
  double pr=0.0;
  double dr=0.0;

  ps=getDifference(pol_codon_consensus,0,N-1,'S');
  ds=getDifference(div_codon_consensus,0,N-1,'S');
  pr=getDifference(pol_codon_consensus,0,N-1,'R');
  dr=getDifference(div_codon_consensus,0,N-1,'R');

  cout<<"The gene length: "<<N<<" bp"<<endl;
  cout<<"PS: "<<ps<<endl;
  cout<<"DS: "<<ds<<endl;
  cout<<"PR: "<<pr<<endl;
  cout<<"DR: "<<dr<<endl;


  //Estimate the divergence time between species using the information from whole sequence.

  if(Div_time==0.0){
    divergent_time=DivergentTime(species_n, N);
    if(ps==0.0 || ds==0.0){
      cout<<endl<<"*************"<<endl;
      cout<<"Warning: There is not enough information from sysnonymous sites to estimate divergence time! The estimation of r is not right!"<<endl;
      cout<<"Note:    Please give the empirical divergence time, run the program to estimate r."<<endl;
      cout<<"*************"<<endl;
    }else{
      cout<<endl<<"*************"<<endl;
      cout<<"Note: Divergence time is estimated basing on this sequence which might be not or even close to species divergence time!"<<endl;
      cout<<"      Prefer to use species divergence time to estimated r!"<<endl;
      cout<<"*************"<<endl;
    }
  }else{
    divergent_time=Div_time;
  }

  cout<<"Species_n="<<species_n<<endl;
  cout<<"Divergent Time: "<<divergent_time<<endl;
  cout<<endl;
  cout<<"Before Clustering (Total Time elapsed: ";
  time_t t2total = time(NULL)-time_startRunML;
  int h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  //Use ps/pr/ds/dr information to find the cluster from the sequences
  //flag_seq==0 means using pol_seq
  //flag_seq==1 means using div_seq
  int flag_seq=0;

  //// Fixed bug due to the new OS X 10.9 system [error: variable length array of non-POD element type 'struct SiteModels']
  ////Solution: use a very large number instead of the parameter N for the gene length, for keeping all models for each gene site, to make sure the number is larger than the gene length.
 // struct SiteModels sm_pol[N];
  struct SiteModels sm_pol[10000];
  if (N>10000) { cout<<"The length of the gene exceeds 10000, Revise the SiteModels upper-boundary array size!"<<endl; throw 1;}

  if(Sys_cluster==1){
    //Initialize for PS
    vec_SelectedModels.clear();
    vec_MS_rate.clear();
    vec_MA_rate.clear();
    vec_lower_rate.clear();
    vec_upper_rate.clear();

    vec_MS_rate.resize(N,0.0);
    vec_MA_rate.resize(N,0.0);
    vec_lower_rate.resize(N,0.0);
    vec_upper_rate.resize(N,0.0);
    
    flag_seq=0;

    cout<<"****** Start Clustering Polymorphism Sysnonymous"<<endl;
    time_t time_start1 = time(NULL); // Record the start time

    ClusterSubSeq(0, N-1,'S',flag_seq,sm_pol);

    cout<<"End Clustering Polymorphism Sysnonymous (Time elapsed: ";
    time_t t2 = time(NULL)-time_start1;
    h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
    if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
    else   cout<<m<<":"<<s<<")"<<endl;

    cout<<"End Clustering Polymorphism Sysnonymous (Total Time elapsed: ";
    t2total = time(NULL)-time_startRunML;
    h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
    if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
    else   cout<<m<<":"<<s<<")"<<endl;


    vec_SelectedModels_ps=vec_SelectedModels;
    vec_MS_rate_ps=vec_MS_rate;
    vec_MA_rate_ps=vec_MA_rate;
    if(MS_only==0 && ci_ma==1 && ps>1){
      CI_MA(sm_pol,N);
    }
    vec_lower_rate_ps=vec_lower_rate;
    vec_upper_rate_ps=vec_upper_rate;
  }



  //Initialize for PR
  vec_SelectedModels.clear();
  vec_MS_rate.clear();
  vec_MA_rate.clear();
  vec_lower_rate.clear();
  vec_upper_rate.clear();

  vec_MS_rate.resize(N,0.0);
  vec_MA_rate.resize(N,0.0);
  vec_lower_rate.resize(N,0.0);
  vec_upper_rate.resize(N,0.0);

  flag_seq=0;
  cout<<"****** Start Clustering Polymorphism Replacement"<<endl;
  time_t time_start1 = time(NULL); // Record the start time
  ClusterSubSeq(0, N-1,'R',flag_seq,sm_pol);
  cout<<"End Clustering Polymorphism Replacement (Time elapsed: ";
  time_t t2 = time(NULL)-time_start1;
  h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  cout<<"End Clustering Polymorphism Replacement (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;


  vec_SelectedModels_pr=vec_SelectedModels;
  vec_MS_rate_pr=vec_MS_rate;
  vec_MA_rate_pr=vec_MA_rate;
  if(MS_only==0 && ci_ma==1 && pr>1){
    CI_MA(sm_pol,N);
  }
  vec_lower_rate_pr=vec_lower_rate;
  vec_upper_rate_pr=vec_upper_rate;
  
  //// Fixed bug due to the new OS X 10.9 system [error: variable length array of non-POD element type 'struct SiteModels']
  ////Solution: use a very large number instead of the parameter N for the gene length, for keeping all models for each gene site, to make sure the number is larger than the gene length.
 // struct SiteModels sm_div[N];
  struct SiteModels sm_div[10000];
  if (N>10000) { cout<<"The length of the gene exceeds 10000, Revise the SiteModels upper-boundary array size!"<<endl; throw 1;}

  if(Sys_cluster==1){
    //Initialize for DS
    vec_SelectedModels.clear();
    vec_MS_rate.clear();
    vec_MA_rate.clear();
    vec_lower_rate.clear();
    vec_upper_rate.clear();

    vec_MS_rate.resize(N,0.0);
    vec_MA_rate.resize(N,0.0);
    vec_lower_rate.resize(N,0.0);
    vec_upper_rate.resize(N,0.0);

    flag_seq=1;


    ClusterSubSeq(0, N-1,'S',flag_seq,sm_div);
    vec_SelectedModels_ds=vec_SelectedModels;
    vec_MS_rate_ds=vec_MS_rate;
    vec_MA_rate_ds=vec_MA_rate;
    if(MS_only==0 && ci_ma==1 && ds >1){
      CI_MA(sm_div,N);
    }
    vec_lower_rate_ds=vec_lower_rate;
    vec_upper_rate_ds=vec_upper_rate;
  }

  //Initialize for DR
  vec_SelectedModels.clear();
  vec_MS_rate.clear();
  vec_MA_rate.clear();
  vec_lower_rate.clear();
  vec_upper_rate.clear();

  vec_MS_rate.resize(N,0.0);
  vec_MA_rate.resize(N,0.0);
  vec_lower_rate.resize(N,0.0);
  vec_upper_rate.resize(N,0.0);

  flag_seq=1;

  cout<<"****** Start Clustering Divergence Replacement"<<endl;
  time_start1 = time(NULL); // Record the start time
  ClusterSubSeq(0, N-1,'R',flag_seq,sm_div);
  cout<<"End Clustering Divergence Replacement (Time elapsed: ";
  t2 = time(NULL)-time_start1;
  h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  cout<<"End Clustering Divergence Replacement (Total Time elapsed: ";
   t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  vec_SelectedModels_dr=vec_SelectedModels;
  vec_MS_rate_dr=vec_MS_rate;
  vec_MA_rate_dr=vec_MA_rate;
  if(MS_only==0 && ci_ma==1 && dr>1){
    CI_MA(sm_div,N);
  }
  vec_lower_rate_dr=vec_lower_rate;
  vec_upper_rate_dr=vec_upper_rate;


  //Do model averaging and estimate gamma
  if (MS_only==0 && r_estimate==1) {

	  //If pr & dr ==1 or 0, r couldn't be estimated.

    if(pr>1.0 || dr>1.0){

      cout<<"****** Start SitePRF using model averaged pr and dr.To get model averaged gamma, ci_r has to be 1."<<endl;
      time_start1 = time(NULL); // Record the start time
      SitePRF(species_n,N);
      cout<<"End SitePRF (Time elapsed: ";
      t2 = time(NULL)-time_start1;
      h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
      else   cout<<m<<":"<<s<<")"<<endl;

      cout<<"End SitePRF (Total Time elapsed: ";
       t2total = time(NULL)-time_startRunML;
      h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
      else   cout<<m<<":"<<s<<")"<<endl;

    }else{
      cout<<endl<<"*************"<<endl;
      cout<<"Note: There are not enough nonsynonymous sites for estimating gamma!"<<endl;
      cout<<"*************"<<endl;
    }


    if(pr>1.0  && dr >1.0 && ci_r==1){
      if(ci_r_exact==1){
		  cout<<"****** Start CIr_exact"<<endl;
		  time_start1 = time(NULL); // Record the start time

		CIr_exact(sm_pol,sm_div,N,species_n);
		cout<<"End CIr_exact (Time elapsed: ";
		t2 = time(NULL)-time_start1;
		h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
		if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
		else   cout<<m<<":"<<s<<")"<<endl;
	      cout<<"End CIr_exact (Total Time elapsed: ";
	       t2total = time(NULL)-time_startRunML;
	      h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
	      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
	      else   cout<<m<<":"<<s<<")"<<endl;

      }else{
          cout<<"****** Start CIr_stochastic"<<endl;
          time_start1 = time(NULL); // Record the start time
			CIr_stochastic(sm_pol,sm_div,N,species_n);
			cout<<"End CIr_stochastic (Time elapsed: ";
			t2 = time(NULL)-time_start1;
			h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
			if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
			else   cout<<m<<":"<<s<<")"<<endl;
		      cout<<"End CIr_stochastic (Total Time elapsed: ";
		       t2total = time(NULL)-time_startRunML;
		      h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
		      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
		      else   cout<<m<<":"<<s<<")"<<endl;
      }
    }else if(ci_r==1){
      cout<<endl<<"*************"<<endl;
      cout<<"Note: There are not enough nonsynonymous sites for estimating confidence interval of gamma!"<<endl;
      cout<<"*************"<<endl;
    }


  }
  time_t t_endAfterCIr = time(NULL)-time_startRunML;
  h=t_endAfterCIr/3600, m=(t_endAfterCIr%3600)/60, s=t_endAfterCIr-(t_endAfterCIr/60)*60;
  cout<<"End CIr (Total Time elapsed: ";
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

    
  if(MS_only==1){
    cout<<endl<<"*************"<<endl;
    cout<<"Warning:"<<endl<<"In terms of NO model averaging, it won't estimate selection coefficient (gamma) and its confidence intervals. Please check tutorial for more details!"<<endl;
    cout<<"*************"<<endl;
  }

  if(MS_only==0 && Sys_cluster==1 && NI_estimate==1){
    SiteNI(N);
  }else{
    NI_estimate=0;
    cout<<endl<<"**************"<<endl;
    cout<<"Warning: "<<endl<<"For estimating Neutrality Index (NI), model-averaging and synonymous cluster will be needed. Please set '-s 1 -m 0 -NI 1'. "<<endl;
    cout<<"*************"<<endl;
  }
  time_t t_endBeforeOutput = time(NULL)-time_startRunML;
  h=t_endBeforeOutput/3600, m=(t_endBeforeOutput%3600)/60, s=t_endBeforeOutput-(t_endBeforeOutput/60)*60;
  cout<<"Before Output (Total Time elapsed: ";
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  output(N);

  time_t t_endRunML = time(NULL)-time_startRunML;
  h=t_endRunML/3600, m=(t_endRunML%3600)/60, s=t_endRunML-(t_endRunML/60)*60;
  cout<<"End RunML (Total Time elapsed: ";
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  return 1;
}


string PRFCluster::getPolSysRep(vector<string> seq) {
  long i, j, k;
  long seq_length = seq[0].length();
  string cons_seq = "";
  
  for (i=0; i<seq_length; i+=CODONSIZE) {	
    //Get all of the codons from all of the sequences
    vector<string> codon;
    //Add all of the codons with other letters than A/T/G/C
    vector<string> codon_other;
    for(j=0; j<seq.size(); j++){
      string word=seq[j].substr(i,CODONSIZE);
      if((word[0]=='A' ||word[0]=='T'||word[0]=='G'||word[0]=='C') && (word[1]=='A' ||word[1]=='T'||word[1]=='G'||word[1]=='C') && (word[2]=='A' ||word[2]=='T'||word[2]=='G'||word[2]=='C')){
	codon.push_back(seq[j].substr(i,CODONSIZE));
      }else{
	codon_other.push_back(word);
      }
    }


    //If some of the codons are with other than A, T, G or C, they have to be replaced.
    if(codon_other.size()>0){
      flag_N_pol=1;

      //Replace the ambiguous nucleotide
      if(Nuc_replace==1){
	//Find the most frequently used nucleotide symbol (A, T, G or C) for site0, site1, site2
	//And replace the sites for each codon that are not A, T, G or C
	ReplaceCodon4PolSys(codon,codon_other);
	
	//Merge codon array and codon_other array
	for(j=0;j<codon_other.size();j++){
	  codon.push_back(codon_other[j]);
	}
      }else{
	//Ambiguous nucleotide is seen as gap
	cons_seq += '*';
	cons_seq += '*';
	cons_seq += '*';
	continue;
      }
	
    }


    //Find the R & S for the sequences
    //R: replacement (nonsynonymous), S: sysnonymous
    bool isSys = true;		
    bool isSameNuc[CODONSIZE];
    isSameNuc[0] = isSameNuc[1] = isSameNuc[2] = true;
    for(j=0; j<codon.size()-1; j++) {
      if (isSys==true) isSys = (isSys && getAminoAcid(codon[j])==getAminoAcid(codon[j+1]));
      if (isSameNuc[0]==true) isSameNuc[0] = (isSameNuc[0] && codon[j][0]==codon[j+1][0]);
      if (isSameNuc[1]==true) isSameNuc[1] = (isSameNuc[1] && codon[j][1]==codon[j+1][1]);
      if (isSameNuc[2]==true) isSameNuc[2] = (isSameNuc[2] && codon[j][2]==codon[j+1][2]);					     
    }
    
    
    //S: synonymous; R: replacement
    string SysOrRep = "S";
    if (isSys==false) SysOrRep = "R";
    
    string tmp = "";
    for(k=0; k<CODONSIZE; k++) {			
      if (isSameNuc[k]==true) tmp += "*";
      else tmp += SysOrRep;
    }    
    cons_seq += tmp;
  }
  
  return cons_seq;
}

int PRFCluster::ReplaceCodon4PolSys(vector <string>& codon, vector <string>& codon_other){
  //To store the number of A, T, G, C in other sequences
  string symbol="ATGC";
  int site0[4]={0,0,0,0};
  int site1[4]={0,0,0,0};
  int site2[4]={0,0,0,0};
  for(int j=0; j<codon.size(); j++) {
    if(codon[j][0]=='A') site0[0]++;
    if(codon[j][0]=='T') site0[1]++;
    if(codon[j][0]=='G') site0[2]++;
    if(codon[j][0]=='C') site0[3]++;
    
    if(codon[j][1]=='A') site1[0]++;
    if(codon[j][1]=='T') site1[1]++;
    if(codon[j][1]=='G') site1[2]++;
    if(codon[j][1]=='C') site1[3]++;
    
    if(codon[j][2]=='A') site2[0]++;
    if(codon[j][2]=='T') site2[1]++;
    if(codon[j][2]=='G') site2[2]++;
    if(codon[j][2]=='C') site2[3]++;
  }

  //Replace each codon for three sites that are not A, T, G or C
  for(int j=0;j<codon_other.size();j++){
    char word=codon_other[j][0];
    if(codon_other[j][0]!='A' && codon_other[j][0]!='T' && codon_other[j][0]!='G' && codon_other[j][0]!='C'){
      word=ReplaceSite(codon_other[j][0],site0,symbol);
      codon_other[j][0]=word;
    }

    word=codon_other[j][1];
    if(codon_other[j][1]!='A' && codon_other[j][1]!='T' && codon_other[j][1]!='G' && codon_other[j][1]!='C'){
      word=ReplaceSite(codon_other[j][1],site1,symbol);
      codon_other[j][1]=word;
    }

    word=codon_other[j][2];
    if(codon_other[j][2]!='A' && codon_other[j][2]!='T' && codon_other[j][2]!='G' && codon_other[j][2]!='C'){
      word=ReplaceSite(codon_other[j][2],site2,symbol);
      codon_other[j][2]=word;
    }
  }

  return 1;
}


char PRFCluster::ReplaceSite(char codon_other_symbol, int pointer_site[], string symbol){
  char word=codon_other_symbol;

  //A, T, G or C
  if(codon_other_symbol=='N'){
    long site_tmp=0;
    for(int j=1;j<4;j++){
      if(pointer_site[j]>=pointer_site[site_tmp]) site_tmp=j;
    }
    word=symbol[site_tmp];
    return(word);
  }

  //A or G
  if(codon_other_symbol=='R'){
    if(pointer_site[0]>pointer_site[2]){
      word=symbol[0];
    }else{
      word=symbol[2];
    }
    return(word);
  }

  //C or T
  if(codon_other_symbol=='Y'){
    if(pointer_site[1]>pointer_site[3]){
      word=symbol[1];
    }else{
      word=symbol[3];
    } 
    return(word);
  }

  //A or C
  if(codon_other_symbol=='M'){
    if(pointer_site[0]>pointer_site[3]){
      word=symbol[0];
    }else{
      word=symbol[3];
    } 
    return(word);
  }
  
  //G or T
  if(codon_other_symbol=='K'){
    if(pointer_site[1]>pointer_site[2]){
      word=symbol[1];
    }else{
      word=symbol[2];
    } 
    return(word);
  }

  //C or G
  if(codon_other_symbol=='S'){
    if(pointer_site[3]>pointer_site[2]){
      word=symbol[3];
    }else{
      word=symbol[2];
    } 
    return(word);
  }

  //A or T
  if(codon_other_symbol=='W'){
    if(pointer_site[0]>pointer_site[1]){
      word=symbol[0];
    }else{
      word=symbol[1];
    } 
    return(word);
  }

  //A, C or T
  if(codon_other_symbol=='H'){
    if(pointer_site[0]>pointer_site[1]){
      if(pointer_site[0]>pointer_site[3]){
	word=symbol[0];
      }else{
	word=symbol[3];
      }
    }else{
      if(pointer_site[1]>pointer_site[3]){
	word=symbol[1];
      }else{
	word=symbol[3];
      }
    }
    return(word);
  }

  //C, G or T
  if(codon_other_symbol=='B'){
    if(pointer_site[1]>pointer_site[2]){
      if(pointer_site[1]>pointer_site[3]){
	word=symbol[1];
      }else{
	word=symbol[3];
      } 
    }else{
      if(pointer_site[2]>pointer_site[3]){
	word=symbol[2];
      }else{
	word=symbol[3];
      } 
    }  
    return(word);
  }

  //A, C or G
  if(codon_other_symbol=='V'){
    if(pointer_site[0]>pointer_site[2]){
      if(pointer_site[0]>pointer_site[3]){
	word=symbol[0];
      }else{
	word=symbol[3];
      } 
    }else{
      if(pointer_site[2]>pointer_site[3]){
	word=symbol[2];
      }else{
	word=symbol[3];
      } 
    }  
    return(word);
  }

  //A, G or T
  if(codon_other_symbol=='D'){
    if(pointer_site[0]>pointer_site[1]){
      if(pointer_site[0]>pointer_site[2]){
	word=symbol[0];
      }else{
	word=symbol[2];
      } 
    }else{
      if(pointer_site[1]>pointer_site[2]){
	word=symbol[1];
      }else{
	word=symbol[2];
      } 
    }
    return(word);
  }

  return (word);
}

string PRFCluster::getDivSysRep(vector<string> pol, vector<string> div) {
  long i, j, k;
  long seq_length = div[0].length(); 
  string cons_seq = "";

  for (i=0; i<seq_length; i+=CODONSIZE) {

    string div_codon = div[0].substr(i,CODONSIZE);

    //Get all of the polymorphism codon
    vector<string> codon;
    //Add all of the codons without A, T, G or C                                                                                                                    
    vector<string> codon_other;
    for(j=0; j<pol.size(); j++){
      string word=pol[j].substr(i,CODONSIZE);
      if((word[0]=='A' ||word[0]=='T'||word[0]=='G'||word[0]=='C') && (word[1]=='A' ||word[1]=='T'||word[1]=='G'||word[1]=='C') && (word[2]=='A' ||word[2]=='T'||word[2]=='G'||word[2]=='C')){
	codon.push_back(pol[j].substr(i,CODONSIZE));
      }else{
        codon_other.push_back(word);
      }
    }

    
    //Consider about the nucleotide other than A, T, G or C
    if(codon_other.size()>0 || ((div_codon[0]!='A' && div_codon[0]!='T' && div_codon[0]!='G' && div_codon[0]!='C') || (div_codon[1]!='A' && div_codon[1]!='T' && div_codon[1]!='G' && div_codon[1]!='C') || (div_codon[2]!='A' && div_codon[2]!='T' && div_codon[2]!='G' && div_codon[2]!='C'))){
      
      //Replace the ambiguous nucleotide
      if(Nuc_replace==1){
	//Find the most frequently used nucleotide symbol (A, T, G or C) for site0, site1, site2
	//And replace the nucleotide which is not A, T, G or C with proper nucleotide (A, T, G or C)
	ReplaceCodon4DivSys(codon,codon_other,div_codon);
	
	if(codon_other.size()>0){
	  //Merge codon array and codon_other array
	  for(j=0;j<codon_other.size();j++){
	    codon.push_back(codon_other[j]);
	  }
	}
      }else{
	//Ambiguous nucleotide is seen as gap
        cons_seq += '*';
        cons_seq += '*';
        cons_seq += '*';
        continue;
      }

    }
      

    //Find the R & S for divergence sequence
    bool isDivCodon = true;
    for(j=0; j<codon.size() && isDivCodon==true; j++) {
      if (div_codon==codon[j]) {
	isDivCodon = false;
      }
    }

    bool isSys = true;
    bool isPolNuc[CODONSIZE];
    isPolNuc[0] = isPolNuc[1] = isPolNuc[2] = false;
    bool isSameNuc[CODONSIZE];
    isSameNuc[0] = isSameNuc[1] = isSameNuc[2] = true;
    if(isDivCodon==true) {
      for(j=0; j<codon.size(); j++) {
	if (isSys==true) isSys = (isSys && getAminoAcid(div_codon)==getAminoAcid(codon[j]));
	if (isPolNuc[0]==false) isPolNuc[0] = (isPolNuc[0] || div_codon[0]==codon[j][0]);
	if (isPolNuc[1]==false) isPolNuc[1] = (isPolNuc[1] || div_codon[1]==codon[j][1]);
	if (isPolNuc[2]==false) isPolNuc[2] = (isPolNuc[2] || div_codon[2]==codon[j][2]);
	if(j!=codon.size()-1){
	  if (isSameNuc[0]==true) isSameNuc[0] = (isSameNuc[0] && codon[j][0]==codon[j+1][0]);
	  if (isSameNuc[1]==true) isSameNuc[1] = (isSameNuc[1] && codon[j][1]==codon[j+1][1]);
	  if (isSameNuc[2]==true) isSameNuc[2] = (isSameNuc[2] && codon[j][2]==codon[j+1][2]);
	}
      }
    }
    
    string tmp = "";
    string SysOrRep = "S";
    if (isSys==false) SysOrRep = "R";
    
    if (isDivCodon==false) {
      tmp += "***";
    }
    else {
      for(k=0; k<CODONSIZE; k++) {
	if (isPolNuc[k]==true) tmp += "*";
	else if(isSameNuc[k]==true) tmp += SysOrRep;
	else tmp += '-';
      }
    }
    cons_seq += tmp;
  }

  return cons_seq;
}


int PRFCluster::ReplaceCodon4DivSys(vector <string>& codon, vector <string>& codon_other, string& div_codon){
  //The number of A, T, G, C in other sequences                                                                                                        
  string symbol="ATGC";
  int site0[4]={0,0,0,0};
  int site1[4]={0,0,0,0};
  int site2[4]={0,0,0,0};
  for(int j=0; j<codon.size(); j++) {
    if(codon[j][0]=='A') site0[0]++;
    if(codon[j][0]=='T') site0[1]++;
    if(codon[j][0]=='G') site0[2]++;
    if(codon[j][0]=='C') site0[3]++;

    if(codon[j][1]=='A') site1[0]++;
    if(codon[j][1]=='T') site1[1]++;
    if(codon[j][1]=='G') site1[2]++;
    if(codon[j][1]=='C') site1[3]++;

    if(codon[j][2]=='A') site2[0]++;
    if(codon[j][2]=='T') site2[1]++;
    if(codon[j][2]=='G') site2[2]++;
    if(codon[j][2]=='C') site2[3]++;
  }

  //Replace each codon for three sites that are not A, T, G or C
  if(codon_other.size()>0){
    for(int j=0;j<codon_other.size();j++){
      char word=codon_other[j][0];
      if(codon_other[j][0]!='A' && codon_other[j][0]!='T' && codon_other[j][0]!='G' && codon_other[j][0]!='C'){
	word=ReplaceSite(codon_other[j][0],site0,symbol);
	codon_other[j][0]=word;
      }
      
      word=codon_other[j][1];
      if(codon_other[j][1]!='A' && codon_other[j][1]!='T' && codon_other[j][1]!='G' && codon_other[j][1]!='C'){
	word=ReplaceSite(codon_other[j][1],site1,symbol);
	codon_other[j][1]=word;
      }
      
      word=codon_other[j][2];
      if(codon_other[j][2]!='A' && codon_other[j][2]!='T' && codon_other[j][2]!='G' && codon_other[j][2]!='C'){
	word=ReplaceSite(codon_other[j][2],site2,symbol);
	codon_other[j][2]=word;
      }
    }
  }

  //Replace the codon from divergence sequence if A, T, G or C exists in the codon
  if(div_codon[0]!='A' && div_codon[0]!='T' && div_codon[0]!='G' && div_codon[0]!='C'){
    char word=ReplaceSite(div_codon[0],site0,symbol);
    div_codon[0]=word;
    flag_N_div=1;
  }
  if(div_codon[1]!='A' && div_codon[1]!='T' && div_codon[1]!='G' && div_codon[1]!='C'){
    char word=ReplaceSite(div_codon[1],site1,symbol);
    div_codon[1]=word;
    flag_N_div=1;
  }
  if(div_codon[2]!='A' && div_codon[2]!='T' && div_codon[2]!='G' && div_codon[2]!='C'){
    char word=ReplaceSite(div_codon[2],site2,symbol);
    div_codon[2]=word;
    flag_N_div=1;
  }


  return 1;
}

long PRFCluster::getDifference(string seq, int pos_start, int pos_end, char symbol) {
  long i, n = 0;
  for (i=pos_start; i<=pos_end; i++) {
    if (seq[i]==symbol) n++;
  }
  return n;
}


int PRFCluster::ClusterSubSeq(int pos_start, int pos_end,char symbol,int flag_seq, struct SiteModels *pointer) {
  long N = pos_end - pos_start + 1;
  
  long symbol_n = 0;
  if(flag_seq==0){
    symbol_n=getDifference(pol_codon_consensus, pos_start, pos_end, symbol);
  }else if(flag_seq==1){
    symbol_n=getDifference(div_codon_consensus, pos_start, pos_end, symbol);
  }
  if (N==0 || N==1 || symbol_n==0 || symbol_n==1) return 1;
  double para=1.0;
  double InL0 = BinomialProb(N, symbol_n);
  double AIC0,AICc0,BIC0;
  AIC0=AICc0=BIC0=-2*InL0+2*para; //the parameter is p0 and pc, estimated from the models based on cs and ce, no pc for the null model, so para=1
  
  double InL = InL0;
  double AIC = AIC0;
  double AICc = AICc0;
  double BIC = BIC0;
	
  double min_cri = 10000000;
  

  long cs, ce, cs_max, ce_max;
  double p0_max,pc_max;
  //Found the lowest AIC/BIC
  int found=0;
  double para_min;
  

  vec_AllModels.clear();
  for (cs=pos_start; cs<=pos_end; cs+=1) {
    for(ce=cs; ce<=pos_end; ce+=1) {

      if (cs==pos_start && ce==pos_end) para = 1.0; //null model, without pc, only p0, thus para=1
      else if(pos_end==ce || cs==pos_start) para=2.0; //have both parameters p0 and pc, thus para=2
      else para=2.0;//have both parameters p0 and pc, thus para=2
      
      
      double symbol_cn = 0.0;
      if(flag_seq==0){
	symbol_cn = getDifference(pol_codon_consensus, cs, ce, symbol);
      }else if(flag_seq==1){
	symbol_cn = getDifference(div_codon_consensus, cs, ce, symbol);
      }
      
      //Center & Non-center region
      double InL_tmp = BinomialProb(ce-cs+1, symbol_cn) + BinomialProb(pos_end-ce+cs-pos_start, symbol_n-symbol_cn); //Log likelihood
      double AIC_tmp  = -2*InL_tmp + 2*para;
      double AICc_tmp = AIC_tmp;
      if (N-para-1>0.0) AICc_tmp += 2*para*(para+1)/(N-para-1);
      else AICc_tmp = 2*AIC_tmp;
      double BIC_tmp = -2*InL_tmp + para*log(double(N));
      
      double cri, cri0;
      //BIC
      if (criterion_type==0){
	cri = BIC_tmp;
	cri0=BIC;
      }
      //AIC
      else if (criterion_type==1){
	cri = AIC_tmp;
	cri0=AIC;
      }
      //AICc
      else if (criterion_type==2){
	cri = AICc_tmp;
	cri0=AICc;
      }
      
      
      if(min_cri > cri) min_cri=cri;
      
      double p0, pc;
      getp0pc_MK(pos_start, pos_end, cs, ce, p0, pc, symbol_n, symbol_cn);
      CandidateModels tmp_CM(cri, pos_start, pos_end, cs, ce, p0, pc,InL_tmp);
      vec_AllModels.push_back(tmp_CM);

      //cout<<"Model AIC weight:\t"<<AIC_tmp<<"\tLogLikelihood:\t"<<InL_tmp<<endl;

      
      //Evaluate the cluster by the criterion
      if (cri <= cri0) {
	if ( (cs-pos_start>3 && pos_end-ce>3 && ce-cs>3) ||
	     (cs==pos_start && pos_end-ce>3 && ce-cs>3) ||
	     (cs-pos_start>3 && pos_end==ce && ce-cs>3)
	     ) {
	  found = 1;						
	  cs_max = cs;
	  ce_max = ce;
	  p0_max=p0;
	  pc_max=pc;
	  
	  InL = InL_tmp;						
	  AIC = AIC_tmp;
	  AICc = AICc_tmp;
	  BIC = BIC_tmp;
	  para_min = para;
	  
	}
      }		
    }		
  }

  cout<<"The total models in ClusterSubSeq: "<<vec_AllModels.size()<<endl;
  if (found==0){
    //If it could not reject the null model, then keep the null model.
    if((symbol=='S' && flag_seq==0 && flag_found_ps==0) || (symbol=='R' && flag_seq==0 && flag_found_pr==0) || (symbol=='S' && flag_seq==1 && flag_found_ds==0) || (symbol=='R' && flag_seq==1 && flag_found_dr==0)){
      double p_tmp=(double)symbol_n/(double)N;
      CandidateModels nullmodel(0, N-1, 0, N-1, p_tmp,p_tmp, InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
      vec_SelectedModels.push_back(nullmodel);
      if(MS_only==0) ModelAveraging(0, N-1, 0, N-1, p_tmp,p_tmp, min_cri,pointer);
      ////ModelAveraging(0, N-1, 0, N-1, p_tmp,p_tmp, min_cri,pointer);

    }
    return 1;
  }else{
    if(symbol=='S' && flag_seq==0) flag_found_ps++;
    if(symbol=='R' && flag_seq==0) flag_found_pr++;
    if(symbol=='S' && flag_seq==1) flag_found_ds++;
    if(symbol=='R' && flag_seq==1) flag_found_dr++;
  }
  
  

  CandidateModels selectedmodel(pos_start, pos_end, cs_max, ce_max, p0_max, pc_max, InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
  vec_SelectedModels.push_back(selectedmodel);
  
  if(MS_only==0) ModelAveraging(pos_start, pos_end, cs_max, ce_max, p0_max, pc_max, min_cri,pointer);
  ////ModelAveraging(pos_start, pos_end, cs_max, ce_max, p0_max, pc_max, min_cri,pointer);
  /* Divide and Conquer */
  if (ce_max!=pos_end || cs_max!=pos_start) {				
    if (cs_max>pos_start+3) ClusterSubSeq(pos_start, cs_max-1,symbol,flag_seq,pointer);
    if (ce_max<pos_end-3) ClusterSubSeq(ce_max+1, pos_end,symbol,flag_seq,pointer);
    if (cs_max<ce_max-3) ClusterSubSeq(cs_max, ce_max, symbol,flag_seq,pointer);
  }  
  return 1;
}


int PRFCluster::ModelAveraging(long pos_start, long pos_end, long cs_max, long ce_max, double p0_max, double pc_max, double min_cri, struct SiteModels *pointer){
	long i, j;

	for (i=pos_start; i<=pos_end; i++) {
		vec_MS_rate[i] = p0_max;
	}
	for (i=cs_max; i<=ce_max; i++) {
		vec_MS_rate[i] = pc_max;
	}

	/* Model Averaging */
	double all_weight = 0.0;


	for (i=0; i<vec_AllModels.size(); i++) {
		vec_AllModels[i].CW = exp(-0.5*(vec_AllModels[i].CW - min_cri));
		all_weight += vec_AllModels[i].CW;
	}

	for (i=0; i<vec_AllModels.size(); i++) {
		vec_AllModels[i].CW = vec_AllModels[i].CW/all_weight;
	}


	for (i=pos_start; i<=pos_end; i++) {
		//Model Averaging and CI for this site
		double rate = 0.0;
		//probability and weight for all possible models
		vector<CI> CIs;
		CIs.clear();

		for (j=0; j<vec_AllModels.size(); j++) {
			float site_weight=0.0;
			site_weight = vec_AllModels[j].CW;
			float site_rate=0.0;
			site_rate = vec_AllModels[j].pc;
			float site_LogLikelihood=0.0;
			site_LogLikelihood=vec_AllModels[j].LogLikelihood;
			//cout<<"Model "<<j<<"\tLogLikelihood:\t"<<site_LogLikelihood<<endl;


			if (i<vec_AllModels[j].cs || i>vec_AllModels[j].ce) {
				site_rate = vec_AllModels[j].p0;
			}

			CI tmp_cis(site_weight, site_rate,site_LogLikelihood);
			CIs.push_back(tmp_cis);


			rate += site_rate*site_weight;
		}
		vec_MA_rate[i] = rate;


		pointer[i].sms.clear();
		pointer[i].sms=CIs;
		pointer[i].pos=i;


		CIs.clear();
	}
	return 1;
}

int PRFCluster::CI_MA(struct SiteModels *pointer,long N){
  for(long i=0;i<N;i++){
    //probability and weight for all possible models
    vector<CI> CIs;
    CIs.clear();
    CIs=pointer[i].sms;
    double lower=0.0;
    double upper=0.0;


    while (lower<confidence_interval || upper<confidence_interval) {
      long min_pos=0, max_pos=0;
      long j = 0;
      //flag_lower==0, max_pos is in front of min_pos
      //flag_lower==1, max_pos is behinde of min_pos and min_pos will be deleted.
      int flag_lower=0;
      while (j<CIs.size()) {
        if(CIs[j].p<CIs[min_pos].p) min_pos = j;
        if(CIs[j].p>CIs[max_pos].p) max_pos = j;
        j++;
      }

      if(lower<confidence_interval){
        if(max_pos>min_pos){
          flag_lower=1;
        }

        lower += CIs[min_pos].weight;
        vec_lower_rate[i] = CIs[min_pos].p;

        CIs.erase(CIs.begin() + min_pos);
      }

      if(upper<confidence_interval){
        if(flag_lower==1 && max_pos!=0){
          max_pos--;
        }
        upper += CIs[max_pos].weight;
        vec_upper_rate[i] = CIs[max_pos].p;

        CIs.erase(CIs.begin() + max_pos);
      }
      flag_lower=0;

      //For some sites,one weight is significant hight than others, so lower=upper
      if(CIs.size()==0){
        if(upper > confidence_interval){
          vec_lower_rate[i]=vec_upper_rate[i];

        }else if(lower > confidence_interval){
          vec_upper_rate[i]=vec_lower_rate[i];

        }
        break;
      }
      if(CIs.size()==1){
        if(lower<confidence_interval && upper<confidence_interval){
          vec_lower_rate[i]=CIs[0].p;
          vec_upper_rate[i]=CIs[0].p;

          break;
        }
      }    
    }             
    CIs.clear();
  }
  return 1;
}

double PRFCluster::DivergentTime(long species_n, long N){
  double ps=0.0;
  double ds=0.0;

  ps=getDifference(pol_codon_consensus,0,N-1,'S');
  ds=getDifference(div_codon_consensus,0,N-1,'S');

  if(ps==0.0){
    ps=0.0001;
  }
  if(ds==0.0){
    ds=0.0001;
  }


  //Estimate divergent time                                                                                                         
  double tmp=0.0;
  double j;
  for(j=1.0;j<=species_n-1;j++){
    tmp=tmp+(1/j);
  }
  double t=(ds/ps)*(tmp)-1-(1/species_n);
  return(t);
}


int PRFCluster::SiteNI(long N){
  for(int i=0;i<N;i++){
    if(vec_MA_rate_pr[i]==0){
      vec_NI[i]=99;
    }else if(vec_MA_rate_dr[i]==0){
      vec_NI[i]=-99;
    }else if(vec_MA_rate_ps[i]==0 || vec_MA_rate_ds[i]==0){
      vec_NI[i]=-66;
    }else{
      vec_NI[i]=(vec_MA_rate_pr[i]/vec_MA_rate_dr[i])/(vec_MA_rate_ps[i]/vec_MA_rate_ds[i]);
    }
  }
  return 1;
}





//Calculate the 95% confidence interval for gamma stochastically, randomly pick dr and pr, calculate gamma using CIs_PRF, and get gamma weight; rank gamma, and sum weight, get the 95% weight boundary gamma.
int PRFCluster::CIr_stochastic(struct SiteModels *pr, struct SiteModels *dr, long N,long species_n){
	int parameters=4; // p0,pc for dr and pr, rate-within-cluster, rate-without-cluster for polymorphism and divergence
	for(long i=0;i<N;i++){
		vec_rModels.clear();

		//Each time get one pr & dr, estimate r
		long end_num=N_Random;
		double min_weight=1000000;

		if(pr[i].sms.size()<=N_Random && dr[i].sms.size()<=N_Random){
			//Use the exact way to calculate the CIs for gamma
			//Get all the models
			cout<<"Polymorphism Model size: "<<pr[i].sms.size()<<";\tDivergence Model Size: "<<dr[i].sms.size()<<endl;
			cout<<"******Use exact algorithm to estimate gamma CI since model size smaller than defined.******"<<endl;
			CIr_exact(pr,dr,N,species_n);
		}

		else{
			//Sort the models by weight (pr/dr)
			//BubbleSort(pr,i);
			//BubbleSort(dr,i);

			sort(dr[i].sms.begin(), dr[i].sms.end(), more_than_CI());
			sort(pr[i].sms.begin(), pr[i].sms.end(), more_than_CI());
			end_num=N_Random;

			//Get 10000 models from pr & dr,respectively
			//If all the selected models from pr & dr are zero, flag_mutation=0, won't calculate CIs for r.
			int flag_mutation=1;
			int repeat_num=0;
			vector<double> PpWeightSums; //polymorphism weight sum
			vector<double> DpWeightSums; //divergence weight sum

			//vector with SummedWeight
			PpWeightSums=RandomModel_NumFastInit(pr,i);
			DpWeightSums=RandomModel_NumFastInit(dr,i);


			for(long j=1;j<=end_num;j++){
				//long model_pr=RandomModel_Num(pr,i);
				//long model_dr=RandomModel_Num(dr,i);
				long model_dr=RandomModel_NumFast(DpWeightSums); // create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
				long model_pr=RandomModel_NumFast(PpWeightSums); // create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2

				//cout<<"Pr num=="<<model_pr<<"===Dr num=="<<model_dr<<"=="<<pr[i].sms[model_pr].weight<<"\t"<<dr[i].sms[model_dr].weight<<endl;

				//Re-generate random sample if both pr and dr=0
				if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0){
					j--;
					if(repeat_num==N_Random){
						cout<<"For this site "<<i<<", there is no chance to mutate!"<<endl;
						flag_mutation=0;
						break;
					}else{
						repeat_num++;
					}
					continue;
				}else{

					//cout<<model_pr<<"\t"<<pr[i].sms[model_pr].weight<<"\t"<<pr[i].sms[model_pr].p<<endl;
					double new_r=CIs_PRF(pr[i].sms[model_pr].p,dr[i].sms[model_dr].p,species_n);
					//double prob_tmp=log(pr[i].sms[model_pr].LogLikelihood)+log(dr[i].sms[model_dr].LogLikelihood);
					double log_prob_tmp=pr[i].sms[model_pr].LogLikelihood+dr[i].sms[model_dr].LogLikelihood; //already log likelihood, no log any more
					//cout<<"Test pr dr LogLikelihood: "<<pr[i].sms[model_pr].LogLikelihood<<"\t"<<dr[i].sms[model_dr].LogLikelihood<<endl;
					////cout<<"Model\t"<<j<<"\tGamma:\t"<<new_r<<"\tPR rate\t"<<pr[i].sms[model_pr].p<<"\tDR rate\t"<<dr[i].sms[model_dr].p<<"\tpr/dr\t"<<pr[i].sms[model_pr].p/dr[i].sms[model_dr].p<<endl;
					////cout<<"\tdr_model:\t"<<model_dr<<"\tpr_model:\t"<<model_pr<<"\tpr LogLikelihood:\t"<<pr[i].sms[model_pr].LogLikelihood<<"\tdr LogLikelihood:\t"<<dr[i].sms[model_dr].LogLikelihood<<endl;

					double weight_tmp=-2*log_prob_tmp+2*parameters;
					if(weight_tmp<min_weight) min_weight=weight_tmp;
					rModels tmp_rm(weight_tmp,new_r);
					//rModels tmp_rm(0,new_r);
					vec_rModels.push_back(tmp_rm);
				}
			}

			//Find CI for r
			if(flag_mutation==1){
				//CI_UpLow_r_stochastic(i);
				CI_UpLow_r(i,min_weight);

			}else{
				////cout<<"Both pr and dr are 0. gamma was not calculated for model averaged gamma, 95% CI gamma!"<<endl;
				vec_r[i]=-0;
				vec_lower_r[i]=0;
				vec_upper_r[i]=0;

			}
			vec_rModels.clear();
		}
	}
	return 1;
}

//Calculate the 95% confidence interval for gamma exactly, randomly pick dr and pr, calculate gamma using CIs_PRF, and get gamma weight; rank gamma, and sum weight, get the 95% weight boundary gamma.
int PRFCluster::CIr_exact(struct SiteModels *pr, struct SiteModels *dr, long N,long species_n){
	int parameters=4; // p0,pc for dr and pr, rate-within-cluster, rate-without-cluster for polymorphism and divergence
	double min_weight=1000000;
	for(long i=0;i<N;i++){
		//Get all the models
		vec_rModels.clear();
		for(long j=0;j<pr[i].sms.size();j++){
			for(long k=0;k<dr[i].sms.size();k++){
				if(pr[i].sms[j].p==0.0 && dr[i].sms[k].p==0.0){
					continue;
				}
				double new_r=CIs_PRF(pr[i].sms[j].p,dr[i].sms[k].p,species_n);
				//cout<<"Test pr dr LogLikelihood: "<<pr[i].sms[j].LogLikelihood<<"\t"<<dr[i].sms[j].LogLikelihood<<endl;
				double log_prob_tmp=pr[i].sms[j].LogLikelihood+dr[i].sms[k].LogLikelihood;
				double weight_tmp=-2*log_prob_tmp+2*parameters;
				if(weight_tmp<min_weight) min_weight=weight_tmp;
				rModels tmp_rm(weight_tmp,new_r);
				vec_rModels.push_back(tmp_rm);
			}
		}

		//Caculate the CIs for r
		if(vec_rModels.size()==0){
			vec_lower_r[i]=0;
			vec_upper_r[i]=0;

		}else{
			//CI_UpLow_r_exact(i,min_weight);
			CI_UpLow_r(i,min_weight);
		}
		vec_rModels.clear();
	}
	return 1;
}

/***************************************************
* Function: Get a random model by choosing a value from 0 to 1, and select the clustering model for p with the weight above the random value;
* Function: create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
* Input Parameter: SiteModels; site
* Output:model_num
* Return Value: model_num
***************************************************/

vector<double> PRFCluster::RandomModel_NumFastInit(struct SiteModels *p,long site) {
  int pss = p[site].sms.size();
  CI *cip = &(p[site].sms[0]);

  pWeightSums.clear();
  pWeightSums.push_back(cip->weight);
  ++cip;
  for(long i = 1; i < pss; ++i, ++cip) {
    //if (cip->weight == 0.0) cerr << "null weight " << site << " " << i << endl;
    pWeightSums.push_back(pWeightSums[i-1] + cip->weight);
  }
  return (pWeightSums);
}

long PRFCluster::RandomModel_NumFast(vector <double> &pWeightSums){
  double rand_tmp=rand();
  //RAND_MAX is from system (2147483647)
  double rand_num=rand_tmp/RAND_MAX;
  if(rand_num<=pWeightSums[0]){
      return 0;
  }
  long lo = 0, hi = pWeightSums.size();
  long mid = (lo + hi)/2, nmid;
  while (1) {
    if (pWeightSums[mid] > rand_num)
      hi = mid;
    else
      lo = mid;
    nmid = (lo + hi)/2;
    if (mid == nmid) break;
    mid = nmid;
  }
  //cerr << "blew fast model num " << lo << "  " << hi <<  " " << rand_num << ": " << pWeightSums[lo] << " " << pWeightSums[hi] << endl;

  // not exactly the same condition as the original (which i believe is a little broken).
  if(!(hi == 0 || pWeightSums[hi-1] <= rand_num && rand_num <= pWeightSums[hi])) {
    cerr << "blew fast model num " << lo << "  " << hi <<  " " << rand_num << ": " << pWeightSums[lo] << " " << pWeightSums[hi-1] << " " << pWeightSums[hi] << endl;
    exit(1);
  }
  return(hi);

}

long PRFCluster::RandomModel_Num(struct SiteModels *p,long site){
  long model_num=-1;
  //randomize();
  double rand_tmp=rand();
  //RAND_MAX is from system (2147483647)
  double rand_num=rand_tmp/RAND_MAX;
  //double rand_tmp=rand()%(N_Random*N_Random);
  //double rand_num=rand_tmp/(double)(N_Random*N_Random);
  double sum=0.0;
  for(long i=0;i<p[site].sms.size();i++){
    if(rand_num==0.0){
      model_num=i;
      break;
    }
    if(rand_num>sum && rand_num <= p[site].sms[i].weight+sum){
      model_num=i;
      break;
    }else{
      sum=p[site].sms[i].weight+sum;
    }
  }

  return(model_num);
}


// calculating gamma using model averaged polymorphism and divergence probablity directly, alternatively model averaged gamma can be the sum of all model specific gamma multiplying gamma AIC weight calculated from the two likelihood.
// Calculate gamma for site site
int PRFCluster::SitePRF(int species_n,long N){

  //Estimate r for each site one by one
  for(long i=0; i<N; i++){

    //vec_time[i]=divergent_time;
    if(vec_MA_rate_pr[i]==0 && vec_MA_rate_dr[i]==0){
    	vec_r[i]=-66;
    	continue;
    }

    else if(vec_MA_rate_pr[i]==0 && vec_MA_rate_dr[i]!=0){
      //Under infinite positive selection
      vec_r[i]=50;
      continue;
    }else if(vec_MA_rate_dr[i]==0 && vec_MA_rate_pr[i]!=0){
      //Under infinite negative selection
      vec_r[i]=-50;
      continue;
    }

    double tmp=vec_MA_rate_dr[i]/vec_MA_rate_pr[i];
    double f, df; //f and df are corresponding to Equation (12) and (13) in the paper
    double rtn,dx;
    int gamma=0;
    double gx1,gx1_d,gx2,gx2_d;//gx is F(n); gx_d is the first derivation of r for F(n)
    bool flag_root=false;

   double initial_r=0.0; //bug1: to get the original gamma in the rtn<=IR_H loop if break from the inside JMAX loop, without this, the bug is caused by the new calculated gamma passing from JMAX loop not following the rules in IR_H, and can be outside the gamma boudary of IR_L and IR_H.

   //bug2: keep the optimal gamma with the smallest dx if no gamma value qualified the criteria ER
   double min_dx=49;
   double r_optimal=-66;
   double f_optimal= 49;
   double df_optimal=49;
   double tmp_optimal=49;

    //For each initial gamma, iterating till hit JMAX=200 or the criteria.
    //cout<<endl<<endl<<"******Start SitePRF for Site "<<i<< "******"<<endl;
    for(rtn=IR_L;rtn<=IR_H;rtn+=0.1){
      if(flag_root==true){
		break;
      }
      initial_r=rtn;
      //cout<<endl<<"****** Restart cycle for Gamma "<<rtn<<endl;

      for(int j=1; j<=JMAX; j++){
    	  gx1=0.0; //initialization. x=0; g(x=0)=0
          gx1_d=0.0;
          gx2=0.0;
          gx2_d=0.0;
          gamma=int (rtn);
          double fraction=rtn-gamma;
          //cout<<endl<<"******Iteration: "<<j<<"  Start Gamma: "<<rtn<<" r_LookupTable: "<<gamma<<" r_fraction: "<<fraction;
          //For gamma outside the Lookup table range, break the loop
          if (gamma>50 or gamma<-50) {
        	  break;
          }
         gamma= (gamma+50)*2; //gamma in the table are from -50 to 50

         //Get Fn, Fn' values from Lookup table, the vector of Fn and Fn_r_derivative
         //consider separately for fraction, solving the bugs of failing to calculated gx
          if (1>fraction >= 0.5)
          {
              gamma=gamma+1;
              fraction=fraction-0.5;
        	  gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
              gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
              gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
              gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
          }
          else if (0<fraction<0.5)
          {
        	  gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
              gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
              gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
              gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
          }
          else if (fraction==0.0) //for gamma=50, there is no gamma+1
		  {

			  gx1=Fn1[gamma];
			  gx1_d=Fn1_d[gamma];
			  gx2=Fn2[gamma];
			  gx2_d=Fn2_d[gamma];
		  }

          else if (-1<fraction <=-0.5)
          {
              gamma=gamma-2;
              fraction=fabs(fraction+1);
        	  gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
              gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
              gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
              gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
          }

          else if (-0.5<fraction < 0)
          {
              gamma=gamma-1;
              fraction=fabs(fraction+0.5);
        	  gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
              gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
              gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
              gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
          }
          else
          {
        	  break;
          }

          //Test: potential bugs, values in the vectors Fn or values in gx are not assigned or calculated correctly, and the values are zero.
          //cout<<endl<<" Fn1[gamma]:"<<Fn1[gamma]<<" Fn2[gamma]:"<<Fn2[gamma]<<" Fn1_d[gamma]:"<<Fn1_d[gamma]<<" Fn2_d[gamma]:"<<Fn2_d[gamma];
          //cout<<endl<<" gx1:"<<gx1<<" gx2:"<<gx2<<" gx1_d:"<<gx1_d<<" gx2_d:"<<gx2_d;

          //the integration of x is divided into two parts for f and df, since the parameter tmp (Prd/Prp) is changing.
          f=divergent_time+gx1-tmp*gx2; //f(gamma)
          df=gx1_d-tmp*gx2_d; // f(gamma) derivative
          dx=f/df; //f(r0)/f(r0)'
          rtn -= dx; // iteration for gamma using Newton's method, r1=r0-f(r0)/f(r0)'

          //Keep the minimal dx, and optimal gamma
          if (fabs(dx)<fabs(min_dx)) {

        	  min_dx=fabs(dx);
        	  r_optimal=rtn;
        	  f_optimal=f;
        	  df_optimal=df;
        	  tmp_optimal=tmp;

          }
          //cout<<" Prd/Prp: "<<tmp<<" fr: "<<f<<" fr_d: "<<df<<" dx: "<<dx<<" New Gamma: "<<rtn<<endl;
		if(fabs(dx)<ER) //criteria for quitting the iterations for gamma. ER is the error for the exact gamma, ER=0.001.
		{
		  vec_r[i]=rtn;

		  //cout<<endl<<"******Successfully Estimated gamma value!!!!******"<<endl;
		  //cout<<" Prd/Prp: "<<tmp<<" fr: "<<f<<" fr_d: "<<df<<" dx: "<<dx<<" New Gamma: "<<rtn<<endl<<endl;
		  flag_root=true;
		  break;
		}
      }

      rtn=initial_r;
    }
    //cout<<"Site: "<<i*CODONSIZE<<" Gamma value from SitePRF: "<<vec_r[i*CODONSIZE]<<endl;
    //If after all iterations, no gamma meets the criteria ER, still keeps the optimal gamma if min_dx<MinDx.
    if (flag_root==false and min_dx<=MinDx){
    	vec_r[i]=r_optimal;
    }


    if(flag_root==false and min_dx>MinDx){
      vec_r[i]=-66; // solved bug2:for No gamma calculated meeting the criteria of ER (Expected Error), keep the gamma with the smallest dx as the optimal to be output.

      //cout<<endl<<"****** Site: "<<i<<" Optimal Gamma was kept in SitePRF!!!******"<<endl;
      //cout<<" Prd/Prp optimal: "<<tmp_optimal<<" fr_optimal: "<<f_optimal<<" fr_d_optimal: "<<df_optimal<<" min_dx:"<<min_dx<<" New Gamma: "<<r_optimal<<endl<<endl;
    }

    //cout<<"Gamma: "<<vec_r[i]<<"\toptimal_r: "<<r_optimal<<"\tmin_dx: "<<min_dx<<endl;
  }

  return 1;
}


/***************************************************
* Function: Calculate gamma for the specific model from polymorphism and divergence probability/likelihood
* Input Parameter: polymorphism probability, divergence probability, the number of samples
* Output: gamma
***************************************************/

double PRFCluster::CIs_PRF (double p_pr, double p_dr,long species_n){
  double new_r;

  if(p_pr==0.0 && p_dr==0.0){
    return(-66);
  }

  //for pr=0, very strong positive selection, use the maximum gamma in the lookup table.
  if(p_pr==0.0){
    return(50);
  }
  //for dr=0, very strong negative selection, use the minimum gamma in the lookup table.
  if(p_dr==0.0){
    return(-50);
  }
//cout<<endl<<"****** Start Confidence Interval PRF - gamma calculation******"<<endl;

  double tmp=p_dr/p_pr;
  double df,f;
   double rtn,dx;
   int gamma=0;
   double gx1,gx1_d,gx2,gx2_d;//gx is F(n); gx_d is the first derivation of r for F(n)
   bool flag_root=false;

  double initial_r=0.0; //bug1: to get the original gamma in the rtn<=IR_H loop if break from the inside JMAX loop, without this, the bug is caused by the new calculated gamma passing from JMAX loop not following the rules in IR_H, and can be outside the gamma boudary of IR_L and IR_H.

  //bug2: keep the optimal gamma with the smallest dx if no gamma value qualified the criteria ER
  double min_dx=49;
  double r_optimal=-66;
  double f_optimal= 49;
  double df_optimal=49;
  double tmp_optimal=49;

   //For each initial gamma, iterating till hit JMAX=200 or the criteria.
   for(rtn=IR_L;rtn<=IR_H;rtn+=0.1){
     if(flag_root==true){
		break;
     }
     initial_r=rtn;
     ////cout<<endl<<"****** Restart cycle for Gamma "<<rtn<<endl;

     for(int j=1; j<=JMAX; j++){

         //cout<<endl<<"******Iteration: "<<j<<"  Start Gamma: "<<rtn;

   	     gx1=0.0; //initialization. x=0; g(x=0)=0
         gx1_d=0.0;
         gx2=0.0;
         gx2_d=0.0;
         gamma=int (rtn);
         double fraction=rtn-gamma;

         //cout<<" r_LookupTable: "<<gamma<<" r_fraction: "<<fraction;
         //For gamma outside the Lookup table range, break the loop
         if (gamma>50 or gamma<-50) {
       	  break;
         }
        gamma= (gamma+50)*2; //gamma in the table are from -50 to 50

        //Get Fn, Fn' values from Lookup table, the vector of Fn and Fn_r_derivative
        //consider separately for fraction, solving the bugs of failing to calculated gx
         if (1>fraction >= 0.5)
         {
             gamma=gamma+1;
             fraction=fraction-0.5;
       	     gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
             gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
             gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
             gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
         }
         else if (0<fraction<0.5)
         {
       	     gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
             gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
             gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
             gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
         }
         else if (fraction==0.0) //for gamma=50, there is no gamma+1

		 {

			 gx1=Fn1[gamma];
			 gx1_d=Fn1_d[gamma];
			 gx2=Fn2[gamma];
			 gx2_d=Fn2_d[gamma];
		 }

         else if (-1<fraction <=-0.5)
         {
             gamma=gamma-2;
             fraction=fabs(fraction+1);
       	     gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
             gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
             gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
             gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
         }

         else if (-0.5<fraction < 0)
         {
             gamma=gamma-1;
             fraction=fabs(fraction+0.5);
       	     gx1=Fn1[gamma]+fraction*2*(Fn1[gamma+1]-Fn1[gamma]);
             gx1_d=Fn1_d[gamma]+fraction*2*(Fn1_d[gamma+1]-Fn1_d[gamma]);
             gx2=Fn2[gamma]+fraction*2*(Fn2[gamma+1]-Fn2[gamma]);
             gx2_d=Fn2_d[gamma]+fraction*2*(Fn2_d[gamma+1]-Fn2_d[gamma]);
         }
         else
         {
       	  break;
         }

         //Test: potential bugs, values in the vectors Fn or values in gx are not assigned or calculated correctly, and the values are zero.
         //cout<<endl<<" Fn1[gamma]:"<<Fn1[gamma]<<" Fn2[gamma]:"<<Fn2[gamma]<<" Fn1_d[gamma]:"<<Fn1_d[gamma]<<" Fn2_d[gamma]:"<<Fn2_d[gamma];
         //cout<<endl<<" gx1:"<<gx1<<" gx2:"<<gx2<<" gx1_d:"<<gx1_d<<" gx2_d:"<<gx2_d;

         //the integration of x is divided into two parts for f and df, since the parameter tmp (Prd/Prp) is changing.
         f=divergent_time+gx1-tmp*gx2; //f(gamma)
         df=gx1_d-tmp*gx2_d; // f(gamma) derivative
         dx=f/df; //f(r0)/f(r0)'
         rtn -= dx; // iteration for gamma using Newton's method, r1=r0-f(r0)/f(r0)'

         //Keep the minimal dx, and optimal gamma
         if (fabs(dx)<fabs(min_dx)) {

       	  min_dx=fabs(dx);
       	  r_optimal=rtn;
       	  f_optimal=f;
       	  df_optimal=df;
       	  tmp_optimal=tmp;

         }

        if(fabs(dx)<ER) //criteria for quitting the iterations for gamma. ER is the error for the exact gamma, ER=0.001. fabs is the absolute value for dx.
		{
		  new_r=rtn;
		  flag_root=true;
		  break;
		}
    }
    rtn=initial_r;
  }

   if (flag_root==false and min_dx<=MinDx){
	   new_r=r_optimal;
   }

   if(flag_root==false and min_dx>MinDx){
	  new_r=-66;
  }

return (new_r);
}




// Model medium gamma instead of model averaged gamma
/***************************************************
* Function: Calculate gamma and 95% CI gamma. rank by gamma, create a partial summed weight models first, then find the model with the weight according to the CI
* Function: Use the same one for both Stochastic and Exact to find the gamma CI
* Input Parameter: SiteModels; site; min_weight_c
* Output: lower and upper gamma
* Return Value: int 1
***************************************************/
// this is an expedient lie...
inline bool operator<(const rModels& a, const rModels& b)
{
	return a.r < b.r;
}
int PRFCluster::CI_UpLow_r(long site,double min_weight) {
	double all_weight=0.0;
	long lower_model=-1;
	long upper_model=-1;
	rWeightSums.clear();
    int flag_lower=0;
    double lower=0.00, upper=0.00;

    sort(vec_rModels.begin(), vec_rModels.end()); // sort by gamma values incrementally
  //Calculate weight for each model and get a summed weight for gamma
  for(long i=0;i<vec_rModels.size();i++){
	vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
	all_weight+=vec_rModels[i].weight;
  }
  vec_rModels[0].weight=vec_rModels[0].weight/all_weight;
  rWeightSums.push_back(vec_rModels[0].weight);

  	long last_item=vec_rModels.size()-1;
  	double lci=confidence_interval;
  	double uci=1-confidence_interval;
  	double averaged_gamma=0; // average only for weights between 95% CI
  	double medium_gamma=0; // model medium gamma
  	double medium=0.5;

  	//lower_model, for the case the first item is already above the CI; record the lower model as 0; use interpolation, and the model one and model two to get the y axis intercept, to get lci gamma.
  	if (rWeightSums[0]>=lci) {
		lower_model=0;
		flag_lower=1;
		double yinterp=0;

		double r1=vec_rModels[0].r;
		double w1=vec_rModels[0].weight;
		double r2=vec_rModels[1].r;
		double w2=vec_rModels[1].weight;
		//Get the y axis intercept when weight=0, the gamma
		double g0=((w2+w1)*r1-w1*r2)/w2; //get the y intercept for gamma while weight equals 0; the assumption is a linear relationship between weight and gamma value
		////vec_lower_r[site]= (g0+r1)*lci/w1;//the Lower 95% CI gamma, use interpolation
		////averaged_gamma=vec_lower_r[site]*(w1-lci); //gamma, use interpolation, lci gamma*weight plus the first model gamma*weight
		vec_lower_r[site]= r1;
		averaged_gamma=r1*(w1-lci);

		////cout<<"***Lower CI model: "<< lower_model<<"\t SummedWeight: "<<rWeightSums[lower_model]<<"y intercept:\t"<<g0<<"lowerCI gamma:\t"<<vec_lower_r[site]<<"ModelAveragedGamma:\t"<<averaged_gamma<<endl<<endl;

		//To Be Deleted
		//vec_lower_r[site]= vec_rModels[lower_model].r; //the Lower 95% CI gamma, not use interpopation, since another point required for that, simply use the first model as the lower 95% CI model
  		//averaged_gamma=vec_rModels[lower_model].weight*vec_rModels[lower_model].r; //Get the gamma if more than the lower 95% CI

  	}

  	////Added: If the first model is above the medium, record the medium_gamma
  	if (rWeightSums[0]>=medium) {
		medium_gamma=vec_rModels[0].r;
  	}


  	//***Debug the model averaged gamma and lower CI gamma values
  	////cout<<"ModelID:\t"<<"ModelAveragedGamma:\t"<<"lowerCI:\t"<< "UpperCI:\t"<<"gamma[i]:\t"<< "weight[i]:\t"<<"gamma[i-1]:\t"<<"SummedWeight"<<endl;

  	for(long i=1;i<vec_rModels.size();i++){
  		vec_rModels[i].weight=vec_rModels[i].weight/all_weight;
  		double tmp_weight=vec_rModels[i].weight;
  		double SummedWeight=rWeightSums[i-1] + tmp_weight;
  		//if (tmp_weight == 0.0) cerr << "null weight " << site << " " << i << endl;
  		rWeightSums.push_back(SummedWeight);

  		//Get the lower CI gamma; record the model number for the rWeightSums at the two borders of the lower 95% CI lci 0.025, lower_model
  		if (rWeightSums[i-1]<lci && rWeightSums[i]>=lci) {
  			double lci_interpolation_rl=vec_rModels[i-1].r+(vec_rModels[i].r-vec_rModels[i-1].r)/( rWeightSums[i]- rWeightSums[i-1])*(lci-rWeightSums[i-1]); //linear interpolation of lower CI gamma lci 0.025 being between (i-1) and i
  			////vec_lower_r[site] = lci_interpolation_rl;
  			lower_model=i;
  			flag_lower=1;
  			//the first averaged gamma is the part between lci --- right trapezoid; the r[i-1]*weight plus the lci.r*weight
  			////averaged_gamma=(rWeightSums[i]-lci)*vec_lower_r[site]; // model averaged gamma start from the lci
  			vec_lower_r[site] = vec_rModels[i].r;
  			averaged_gamma=(rWeightSums[i]-lci)*vec_rModels[i].r;

  			//cout<<"***Lower CI model id:\t"<< lower_model<<"\tLowerCI gamma:\t"<<lci_interpolation_rl<<"\tSummedWeight:\t"<<rWeightSums[lower_model]<<endl<<endl;
  		}

  		//getting model averaged gamma between 95% CI
  		if (flag_lower==1 && i>lower_model && rWeightSums[i]<=uci){
  			averaged_gamma+=vec_rModels[i].weight*vec_rModels[i].r;
  		}
  		//Get the upper CI gamma; record the model number for the rWeightSums at the two sides of the upper 95% CI uci 0.975, upper_model
  		if (flag_lower==1 && i>lower_model && rWeightSums[i-1]<=uci && rWeightSums[i]>uci) {
  			upper_model=i;
  			////double uci_interpolation_ru=vec_rModels[i-1].r+(vec_rModels[i].r-vec_rModels[i-1].r)/( rWeightSums[i]- rWeightSums[i-1])*( uci-rWeightSums[i-1]); //linear interpolation of upper CI gamma
  			////vec_upper_r[site] = uci_interpolation_ru;
  			//bug: already add the previous one
  			//the final part of averaged gamma, the border of upper 95% CI. the r[i-1]*weight plus the uci.r*weight;
  			//averaged_gamma+=vec_rModels[i-1].weight*vec_rModels[i-1].r+vec_upper_r[site] *(uci-rWeightSums[i-1]);
  			////averaged_gamma+=vec_upper_r[site] *(uci-rWeightSums[i-1]);
  			vec_upper_r[site] = vec_rModels[upper_model].r;
  			averaged_gamma+=vec_rModels[i].r *(uci-rWeightSums[i-1]);
  		}

  		////Added: medium_gamma=vec_rModels[i].r;
  		if (rWeightSums[i-1]<medium && rWeightSums[i]>=medium){
  			medium_gamma=vec_rModels[i].r;
  		}


  		//***Debug the model averaged gamma and lower CI gamma values
  		////cout<<i<<"\t"<< averaged_gamma<<"\t"<< vec_lower_r[site]<<"\t"<< vec_upper_r[site]<<"\t"<< vec_rModels[i].r<<"\t"<<tmp_weight<<"\t"<<vec_rModels[i-1].r<<"\t"<<rWeightSums[i]<<endl;
  	}


 /* 	if (medium_gamma!=0) {vec_r[site]=medium_gamma;}
  	else {
  		cout<<"Warning: model medium gamma equals 0, and converted to -66"<<endl;
  		vec_r[site]=-66;}
*/
	if (averaged_gamma!=0) {
		vec_r[site]=averaged_gamma/(uci-lci); // model averaged gamma only takes from lci to uci, so it should divide (uci-lci).
	}
	else {
		cout<<"Model averaged gamma: NULL"<<endl;
				vec_r[site]=-66;

	}

  	if (vec_r[site]!=-66 )
  	{
  		//Make sure the LowerCI and UpperCI is two sides of value of model average; quit the program if model averaged gamma is not between 95% CI gamma

  		if ( CONVERT<int>(vec_lower_r[site])>CONVERT<int>(vec_r[site]) || CONVERT<int>(vec_upper_r[site]) < CONVERT<int>(vec_r[site]) || flag_lower!=1 || upper_model<lower_model)
  		{
  			cout<<endl<<endl<<"Warning: Problem in getting proper values of model medium gamma, lower and upper 95% CI gamma from all models!!!!"<<endl;
  			cout<<"Site: "<<site<<"\tEstimatedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<"\tTotalModelNO: "<<last_item<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
  			////cout<<"Total number of models: "<<last_item<<" Model 1: "<<vec_rModels[0].r<<"\tSumWeight: "<<rWeightSums[0]<<"\tModel 2: "<<vec_rModels[1].r<<"\tSumWeight: "<<rWeightSums[1]<<"\tModel "<<last_item-1<<": "<<vec_rModels[last_item-1].r<<"\tSumWeight: "<<rWeightSums[last_item-1]<<"; Model"<<last_item<<": "<<vec_rModels[last_item].r<<"\tSumWeight: "<<rWeightSums[last_item]<<endl<<"Lower95% weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
  			throw "Error in getting the right 95% confidence interval gamma in CI_UpLow_rc!";
  		}

  	}
  	else
  	{
  		vec_lower_r[site]=-66;
  		vec_upper_r[site]=-66;
  		cout<<"Warning in getting proper gamma and 95% CI gamma: ";
  		cout<<"Site: "<<site<<"\tEstimatedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<endl;
  		//throw "Error in getting the right gamma in CI_UpLow_rc!";

  	}

  	/* to be deleted; model averaged gamma
  	if (averaged_gamma!=0) {vec_r[site]=averaged_gamma/(uci-lci);}
  	else {
  		cout<<"Warning: averaged gamma equals 0, and converted to -66"<<endl;
  		vec_r[site]=-66;}

  	if (vec_r[site]!=-66 )
  	{
  		//Make sure the LowerCI and UpperCI is two sides of value of model average; quit the program if model averaged gamma is not between 95% CI gamma

  		if ( CONVERT<int>(vec_lower_r[site])>CONVERT<int>(vec_r[site]) || CONVERT<int>(vec_upper_r[site]) < CONVERT<int>(vec_r[site]) || flag_lower!=1 || upper_model<lower_model)
  		{
  			cout<<endl<<endl<<"Warning: Problem in getting proper values of model averaged gamma, lower and upper 95% CI gamma from all models!!!!"<<endl;
  			cout<<"Site: "<<site<<"\tAveragedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<"\tTotalModelNO: "<<last_item<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
  			////cout<<"Total number of models: "<<last_item<<" Model 1: "<<vec_rModels[0].r<<"\tSumWeight: "<<rWeightSums[0]<<"\tModel 2: "<<vec_rModels[1].r<<"\tSumWeight: "<<rWeightSums[1]<<"\tModel "<<last_item-1<<": "<<vec_rModels[last_item-1].r<<"\tSumWeight: "<<rWeightSums[last_item-1]<<"; Model"<<last_item<<": "<<vec_rModels[last_item].r<<"\tSumWeight: "<<rWeightSums[last_item]<<endl<<"Lower95% weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
  			throw "Error in getting the right 95% confidence interval gamma in CI_UpLow_rc!";
  		}

  	}
  	else
  	{
  		vec_lower_r[site]=-66;
  		vec_upper_r[site]=-66;
  		cout<<"Warning in getting proper gamma and 95% CI gamma: ";
  		cout<<"Site: "<<site<<"\tAveragedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<endl;
  		//throw "Error in getting the right gamma in CI_UpLow_rc!";

  	}
*/


  return 1;

}


/***************************************************
* Function: parse Parameters in the running command line
* Input Parameter: int, const char* []
* Output: Parse the input parameters
* Return Value: bool 
***************************************************/
bool PRFCluster::parseParameter(int argc, const char* argv[]) {
  bool flag = true;
  int i;
  string temp;
	
  try {	
    if (argc==2) {
      temp = stringtoUpper(argv[1]);
      if (temp=="-H") showHelpInfo();
      else throw 1;
    }
    else if (argc!=5 && argc!=7 && argc!=9 && argc!=11 && argc!=13 && argc!=15 && argc!=17 && argc!=19 && argc!=21 && argc!=23 && argc!=25 && argc!=27 && argc!=29 && argc!=31 && argc!=33) {
      throw 1;			
    }		
    else {			
      //parse parameters
      int pol_num_flag=0, input_consensus_flag=0, output_flag=0, pol_flag=0, div_flag=0, code_flag=0, criterion_flag=0,ms_flag=0, synonymous_flag=0,divtime_flag=0, ci_ma_flag=0,r_flag=0,ci_r_flag=0,ci_method_flag=0, nuc_replace_flag=0, ni_flag=0;
      for (i=1; i<argc; i++) {				
	temp = stringtoUpper(argv[i]);
	//Polymorphic input fasta file
	if (temp=="-P" && (i+1)<argc && pol_flag==0) {
	  pol_seqfile = argv[++i];						
	  pol_flag++;
	}
	//Divergent input fasta file
	else if (temp=="-D" && (i+1)<argc && div_flag==0) {
	  div_seqfile = argv[++i];
	  div_flag++;
	}
	//Option of using consensus polymorphism and divergence sequences directly.
	else if (temp=="-IC" && (i+1)<argc && input_consensus_flag==0) {
		  int format = CONVERT<int>(argv[++i]);
		  if (format==0) {

			  input_consensus_flag=0;
			  input_format_num=0;
		  }
		  else if (format==1)
		  {
			  input_consensus_flag=1;
			  input_format_num=1;
		  }
		  else {
		    throw 1;
		  }

	}

	//Polymorphism sequence number - the number of samples, required for IC input consensus.
	else if (temp=="-SN" && (i+1)<argc && pol_num_flag==0) {
		species_num = CONVERT<int>(argv[++i]);
		pol_num_flag++;
		if (species_num<1){
			cout<< "Error! The polymorphism sequence number should be no less than one!\n";
			throw 1;
		}
	}


    //Choice of the output format, amino acid or nucleotide level output
	else if (temp=="-O" && (i+1)<argc && output_flag==0) {

		  int format = CONVERT<int>(argv[++i]);
		  if (format==0) {
			  output_flag=0;
			  output_format_num=0;
		  }
		  else if (format==1)
		  {
			  output_flag=1;
			  output_format_num=1;
		  }
		  else {
		    throw 1;
		  }
	}

	//BIC/AIC/AICc/LRT
	else if (temp=="-C" && (i+1)<argc && criterion_flag==0) {
	  int num = CONVERT<int>(argv[++i]);
	  if (num>=0 && num<=3) {
	    criterion_type = num;
	    criterion_flag=1;
	  }
	  else {
	    throw 1;
	  }
	}
	//Genetic code
	else if (temp=="-G" && (i+1)<argc && code_flag==0) {
	  int num = CONVERT<int>(argv[++i]);
	  if (num>0 &&  num<24) {
	    genetic_code = num;
	    code_flag=1;
	  }
	  else {
	    throw 1;
	  }
	}
	//Model selection and model averaging
	else if (temp=="-M" && (i+1)<argc &&  ms_flag==0) {
	  int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
	    MS_only = num;
	    ms_flag=1;
	  }
	  else {
	    throw 1;
	  }
	}
	//Show the clustering results from synonymous sites
	else if (temp=="-S" && (i+1)<argc && synonymous_flag==0){
	  int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
	    Sys_cluster = num;
	    synonymous_flag=1;
	  }else{
	    throw 1;
	  }
	}
	//Species divergence time from users
	else if (temp=="-T" && (i+1)<argc && divtime_flag==0){
          double num=CONVERT<double>(argv[++i]);
	  Div_time = num;
	  divtime_flag=1;
        }
	//Confidence Intervals for Model averaging
	else if(temp=="-CI_M" && (i+1)<argc && ci_ma_flag==0){
	  int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
            ci_ma=num;
	    ci_ma_flag=1;
          }else{
            throw 1;
          }
	}
	//Estimate gamma
	else if(temp=="-R" && (i+1)<argc && r_flag==0){
          int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
            r_estimate=num;
	    r_flag=1;
          }else{
            throw 1;
          }
        }
	//Confidence Intervals of gamma
	else if(temp=="-CI_R" && (i+1)<argc && ci_r_flag==0){
          int num=CONVERT<int>(argv[++i]);	  
	  if(num==0 || num==1){
            ci_r=num;
	    ci_r_flag=1;
          }else{
            throw 1;
          }
        }
	//Exact algorithm for esimating CIs for gamma
	else if(temp=="-EXACT" && (i+1)<argc && ci_method_flag==0){
          int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
            ci_r_exact=num;
	    ci_method_flag=1;
          }else{
            throw 1;
          }
        }

	//Replace or see as gap in the ambiguous nucleotide site 
	else if(temp=="-N" && (i+1)<argc && nuc_replace_flag==0){
	  int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
	    Nuc_replace=num;
	    nuc_replace_flag=1;
	  }else{
	    throw 1;
	  }
	}

	//Estimate Neutrality Index
	else if(temp=="-NI" && (i+1)<argc && ni_flag==0){
	  int num=CONVERT<int>(argv[++i]);
	  if(num==0 || num==1){
	    NI_estimate=num;
	    ni_flag=1;
	  }else{
	    throw 1;
	  }
	}
        else{
          throw 1;
        }


      }			
    }
  }
  catch (...) {
    cout<<"Error in input parameter(s)."<<endl;
    cout<<"Type -h for help information."<<endl;
    cout<<NAME<<", Version: "<<VERSION<<" [Last Update: "<<LASTUPDATE<<"]"<<endl;
    flag = false;		
  }
  
  return flag;
}

//Information and help when type ./cMAC-PRF -h
void PRFCluster::showHelpInfo() {
  cout<<"***********************************************************************"<<endl;
  cout<<NAME<<", Version: "<<VERSION<<" ["<<LASTUPDATE<<"]"<<endl;	
  cout<<"Function: "<<FUNCTION<<endl;
  cout<<"Usage: "<<NAME<<" [OPTIONS]"<<endl;	
  cout<<"***********************************************************************"<<endl<<endl;

  cout<<"OPTIONS:"<<endl;	

  cout<<"  -p\tInput file name for polymorphic sequences [string, required]"<<endl;
  cout<<"  -o\tChoice of output format [integer, optional], {0: amino acid level output || 1: nucleotide level output, default=0}"<<endl;
  cout<<"  -ic\tChoice of input format for divergence and polymorphism sequences [integer, optional], {0: DNA sequences || 1: consensus sequences, default=0}"<<endl;
  cout<<"  -d\tInput file name for divergent sequences [string, required]"<<endl;	
  cout<<"  -sn\tInput the number of polymorphism sequences [integer, required when -ic=1]"<<endl;
  cout<<"  -c\tCriterion used for clustering [integer, optional], {0:BIC || 1:AIC || 2:AICc || 3:LRT}, default = 0"<<endl;
  cout<<"  -g\tGenetic code used for sequences [integer, optional], {1:standard}, default = 1"<<endl;
  cout<<"    \tMore information about the genetic codes can be found at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"<<endl;
  cout<<"  -m\tModel selection and model averaging  [integer, optional], {0: use both model selection and model averaging || 1: use only model selection}, default = 0"<<endl;
  cout<<"  -ci_m\tCalculate 95% confidence intervals for results of model averaging [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 0"<<endl;
  cout<<"  -s\tShow clustering results of synonymous sites from Polymorphism and Divergent sequences [integer, optional], {0: without clustering results of synonymous sites || 1: with clustering results of synonymous and replacement sites}, default = 0"<<endl;
  cout<<"  -r\tEstimate selection coefficient for each site [integer, optional], {0: NOT estimate selection coefficient || 1: estimate selection coefficient}, default=1"<<endl;
  cout<<"  -ci_r\tCalculate 95% confidence intervals for selection coefficient [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 1"<<endl;
  cout<<"  -exact\tAlgorithm for calculating 95% confidence intervals for selection coefficient [integer, optional], {0: use stochastic algorithm || 1: use exact algorithm}, default = 0"<<endl;
  cout<<"  -t\tSpecies divergence time [Default: estimate species divergence time from the sequences,optional]."<<endl;
  cout<<"  -n\tNucleotide is replaced or seen as gap when it is not A, T, G or C in the sequences [integer, optional], {0: see it as gap || 1: replace this nucleotide with the most frequently used nucleotide in other sequences}, default = 1"<<endl;
  cout<<"  -NI\tEstimate the Neutrality Index for each site [integer, optional], {0: NOT estimate Neutrality Index || 1: estimate Neutrality Index}, default=0"<<endl;
  cout<<"  -h\tShow help information"<<endl;	
  cout<<endl;
	
  cout<<"COPYRIGHT & LICENSE:"<<endl;
  cout<<NAME<<" is distributed as open-source software and licensed under the GNU General Public License "<<endl;
  cout<<"(Version 3; http://www.gnu.org/licenses/gpl.txt), in the hope that it will be useful, but "<<endl;
  cout<<"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A "<<endl;
  cout<<"PARTICULAR PURPOSE. See the GNU General Public License for more details."<<endl;
  cout<<endl;
  
  cout<<"REFERENCE:"<<endl;	
  cout<<REFERENCE<<endl;
  cout<<endl;

  cout<<"SEE ALSO:"<<endl;
  cout<<"For more information, please see <http://www.yale.edu/townsend/software.html>."<<endl;
  cout<<endl;
  
  cout<<"CONTACT:"<<endl;
  cout<<"Please send suggestions or bugs to Ning Li (lin.cmb@gmail.com), Zi-Ming Zhao (ziming.gt@gmail.com) and Jeffrey Townsend (Jeffrey.Townsend@Yale.edu)."<<endl;
  cout<<endl;

}


/* To be delete
 *
 * // To be deleted

int PRFCluster::CI_UpLow_r_stochastic(long site){
  //Sort all the models by r, descending
  for(long i=0;i<vec_rModels.size();i++){
    for(long j=i+1;j<vec_rModels.size();j++){
      if(vec_rModels[j].r > vec_rModels[i].r){
        double tmp_r=0.0;
	double tmp_w=0.0;

        tmp_r=vec_rModels[j].r;
        vec_rModels[j].r=vec_rModels[i].r;
        vec_rModels[i].r=tmp_r;

	tmp_w=vec_rModels[j].weight;
	vec_rModels[j].weight=vec_rModels[i].weight;
	vec_rModels[i].weight=tmp_w;
      }
    }
  }

  long models_num=vec_rModels.size();
  long upper_num=(int)(confidence_interval*models_num)-1;
  long lower_num=models_num-(int)(confidence_interval*models_num);
  vec_lower_r[site*CODONSIZE]=vec_rModels[lower_num].r;
  vec_lower_r[site*CODONSIZE+1]=vec_rModels[lower_num].r;
  vec_lower_r[site*CODONSIZE+2]=vec_rModels[lower_num].r;

  vec_upper_r[site*CODONSIZE]=vec_rModels[upper_num].r;
  vec_upper_r[site*CODONSIZE+1]=vec_rModels[upper_num].r;
  vec_upper_r[site*CODONSIZE+2]=vec_rModels[upper_num].r;

  return 1;
}
// To be deleted
int PRFCluster::CI_UpLow_r_exact(long site,double min_weight){
  double all_weight=0.0;
  for(long i=0;i<vec_rModels.size();i++){
    vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
    all_weight+=vec_rModels[i].weight;
  }
  for(long i=0;i<vec_rModels.size();i++){
    vec_rModels[i].weight=vec_rModels[i].weight/all_weight;
  }

  //sort by weight
  double lower=0.0, upper=0.0;
  long j;


  while (lower<confidence_interval || upper<confidence_interval) {
    long min_pos=0, max_pos=0;
    j = 0;
    int flag_lower=0;
    while (j<vec_rModels.size()) {
      if(vec_rModels[j].r<vec_rModels[min_pos].r) min_pos = j;
      if(vec_rModels[j].r>vec_rModels[max_pos].r) max_pos = j;
      j++;
    }

    if(lower<confidence_interval){
      if(max_pos>min_pos){
        flag_lower=1;
      }
      lower += vec_rModels[min_pos].weight;
      vec_lower_r[site*CODONSIZE] = vec_rModels[min_pos].r;
      vec_lower_r[site*CODONSIZE+1] = vec_rModels[min_pos].r;
      vec_lower_r[site*CODONSIZE+2] = vec_rModels[min_pos].r;
      vec_rModels.erase(vec_rModels.begin() + min_pos);
    }

    if(upper<confidence_interval){
      if(flag_lower==1 && max_pos!=0){
        max_pos--;
      }
      upper += vec_rModels[max_pos].weight;
      vec_upper_r[site*CODONSIZE] = vec_rModels[max_pos].r;
      vec_upper_r[site*CODONSIZE+1] = vec_rModels[max_pos].r;
      vec_upper_r[site*CODONSIZE+2] = vec_rModels[max_pos].r;
      vec_rModels.erase(vec_rModels.begin() + max_pos);
    }
    flag_lower=0;

    if(vec_rModels.size()==0){
      if(upper > confidence_interval){
        vec_lower_r[site*CODONSIZE]=vec_upper_r[site*CODONSIZE];
	vec_lower_r[site*CODONSIZE+1]=vec_upper_r[site*CODONSIZE];
	vec_lower_r[site*CODONSIZE+2]=vec_upper_r[site*CODONSIZE];
      }else if(lower > confidence_interval){
        vec_upper_r[site*CODONSIZE]=vec_lower_r[site*CODONSIZE];
	vec_upper_r[site*CODONSIZE+1]=vec_lower_r[site*CODONSIZE];
	vec_upper_r[site*CODONSIZE+2]=vec_lower_r[site*CODONSIZE];
      }
      break;
    }
    if(vec_rModels.size()==1){
      if(lower<confidence_interval && upper<confidence_interval){
        vec_lower_r[site*CODONSIZE]=vec_rModels[0].r;
	vec_lower_r[site*CODONSIZE+1]=vec_rModels[0].r;
	vec_lower_r[site*CODONSIZE+2]=vec_rModels[0].r;

        vec_upper_r[site*CODONSIZE]=vec_rModels[0].r;
	vec_upper_r[site*CODONSIZE+1]=vec_rModels[0].r;
	vec_upper_r[site*CODONSIZE+2]=vec_rModels[0].r;
        break;
      }
    }

  }

  return 1;
}



//Bubble sort all the models based on weight, very slow and take a lot of memory
  int PRFCluster::BubbleSort(struct SiteModels *p, long site){

  //Descending
  for(long i=0;i<p[site].sms.size();i++){
    for(long j=i+1;j<p[site].sms.size();j++){
      if((p[site].sms[j].weight > p[site].sms[i].weight) || ((p[site].sms[j].weight==p[site].sms[i].weight) && (p[site].sms[j].p>p[site].sms[i].p))){
	double tmp_w=0.0;
	double tmp_p=0.0;

	tmp_w=p[site].sms[j].weight;
	p[site].sms[j].weight=p[site].sms[i].weight;
	p[site].sms[i].weight=tmp_w;

	tmp_p=p[site].sms[j].p;
        p[site].sms[j].p=p[site].sms[i].p;
        p[site].sms[i].p=tmp_p;
      }
    }

  }
  return 1;
}
 */
