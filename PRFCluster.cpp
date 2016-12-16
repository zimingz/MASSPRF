/* Updated by Zi-Ming Zhao, 07/11/2016 */
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
#include "PRFCluster.h"
#ifdef __linux__
unsigned int NUM_CPU = get_nprocs_conf();
#else
unsigned int NUM_CPU = std::thread::hardware_concurrency();
#endif

PRFCluster::PRFCluster() {
  flag_N_pol=0;
  flag_N_div=0;

  if (NUM_CPU == 0) NUM_CPU = 8;
  confidence_interval=0.95;
  confidence_interval=(1.0-confidence_interval)/2.0;

  
  site_specific_flag = 0;
  divtime_flag = 0;
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
  modelAveraged_p_gamma=0;
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
/***************************************************
* Function: Log likelihood of Bernoulli distribution. BernoulliProb= (i/n)^i*[(n-i)/n]^(n-i); prob=Log(BernoulliProb)=i*log(i/n)+(n-i)*log[(n-i)/n]
* Input Parameter:
* Output:
***************************************************/
double PRFCluster::BinomialProb(long n, long i) {
  double prob = 0.0;
  prob += (i==0)?0:i*log(double(i)/n);
  prob += (n-i==0)?0:(n-i)*log(double(n-i)/n);
  return prob;
}
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::Run(int argc, const char*argv[]) {	
  int i, flag=1;
  ////srand(1234); // to fix random number generator, to make sure the program is repeatable; need to remove for the final version of the program

  try {		
    //Parse input parameters
    if (parseParameter(argc, argv)!=1) throw 1;
    if (site_specific_flag && !(Sys_cluster))
    {
      cout << "Can't do site_specific divergence time calculation without silent clustering -s option" <<endl;
      throw 1;
    }
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

   if (verbose==1){
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
  }

    //Check sequences' length
  for (i=0; i<pol_seq.size(); i++) {
    if (div_seq[0].length()!=pol_seq[i].length()) throw "Input sequences are not equal in length.";
  }
  if(pol_seq[0].size()%3!=0) cout<<"Warning: the length of sequences can not be divided by 3 (codon size).\n";
    	//throw "Error: the length of sequences can not be divided by 3 (codon size).";
    //Open and Read the lookup table for numeric integration for the specific number of species
     //special attention to gamma=0,purposely put the corresponding numeric integration to 0; treat separately for this case.
  if (verbose==1){
   cout<<endl<<"Open and Read the lookup table for numeric integration:"<<endl;
   cout<<endl<<"The number of species in the polymorphism data:"<<p_species_n<<endl;
 }
       int tmp_n=p_species_n; //the number of species in the polymorphism data, according to the line in the lookup table
       if (p_species_n>100 and p_species_n<=1500)
       {
    	   tmp_n=98+int(p_species_n/50); //for n>100, the interval is 50 in the lookup table
      }
      else if (p_species_n>1500)
      {
        cout<<"Error! The number of sequences in the polymorphism data is over the range in the Lookup Table 1500."<<endl;
      }
       //cout<<"The row to be extracted in the lookup table is "<<tmp_n<<endl;

       //Initialize the vectors,the size of the vector is based on the number of gamma values in the lookup table, from -50 to 50 in the interval of 0.5.
      Fn1.resize(201,0.0);
      Fn1_d.resize(201,0.0);
      Fn2.resize(201,0.0);
      Fn2_d.resize(201,0.0);

      cout<<endl<<"Extract LookupTable_gx1_n_gamma_integration_50_v9.dat"<<endl;
      Fn1=GammaLookupTable(tmp_n,"LookupTable_gx1_n_gamma_integration_50_v9.dat",Fn1);

      cout<<endl<<"Extract LookupTable_gx1_derivative_n_gamma_50_v9.dat"<<endl;
      Fn1_d=GammaLookupTable(tmp_n,"LookupTable_gx1_derivative_n_gamma_50_v9.dat",Fn1_d);

      cout<<endl<<"Extract LookupTable_gx2_n_gamma_integration_50_v9.dat"<<endl;
      Fn2=GammaLookupTable(tmp_n,"LookupTable_gx2_n_gamma_integration_50_v9.dat", Fn2);

      cout<<endl<<"Extract LookupTable_gx2_derivative_n_gamma_50_v9.dat"<<endl;
      Fn2_d=GammaLookupTable(tmp_n,"LookupTable_gx2_derivative_n_gamma_50_v9.dat",Fn2_d);

       //test the sucess of assignment of values in the vectors of Fn and Fn_d
       //cout<<"Test Fn1[0]: "<<Fn1[0]<<endl;
       //cout<<"Test Fn1[86]: "<<Fn1[86]<<endl;
       //cout<<"Test Fn1_d[0]: "<<Fn1_d[0]<<endl;

    //Run the main function
      RunML(pol_seq, div_seq);

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

/***************************************************
* Function: open and read the lookup table, to extract the values for the certain species number.
* Input Parameter:the lookup table file name, p_specie_n converted tmp_n;
* Output: vector<double> Fn
***************************************************/
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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


if (site_specific_flag && vec_SelectedModels_ds.size()!=0 && vec_SelectedModels_ps.size()!=0)
{
  for (long i=0;i<N;i++)
  {
    divergent_time[i] = divergent_time_sums[i]/divergent_time_weights[i];//model averaged divergence time for site i
  }
}
else{
cout<<"Warning in site specific divergence time calculation; -ssd. There is no silent clustering for PS or DS."<<endl;
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
	if (modelAveraged_p_gamma==0) {cout<<endl<<"//Results based on model averaging of gamma using different models: "<<endl;}
	if (modelAveraged_p_gamma==1) {cout<<endl<<"//Results based on gamma calculated from model averaged probablity of polymorphism and divergence replacement: "<<endl;}
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
    cout<<"\t";cout.width(width);cout<<divergent_time[i]<<"\t";
    if(vec_r[i]==50){
      cout.width(width);cout<<INF;
    }else if (vec_r[i]==-50){
      cout.width(width);cout<<NINF;
    }

    else if (vec_r[i]==0 || vec_r[i]==-66){
      cout.width(width);cout<<"NULL";
    }
    else if (vec_r[i]==-99){
      cout.width(width);cout<<"NA";
    }
    else{
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
     }
     else if(vec_lower_r[i]==-99){
       cout<<"\t";cout.width(width);cout<<"NA"<<"\t";
     }
     else{
       cout<<"\t";cout.width(width);cout<<vec_lower_r[i]<<"\t";
     }

     if(vec_upper_r[i]==50){
       cout.width(width);cout<<INF;
     }else if(vec_upper_r[i]==-50){
       cout.width(width);cout<<NINF;
     }

     else if(vec_upper_r[i]==0 || vec_upper_r[i]==-66){
       cout.width(width);cout<<"NULL";
     }
     else if(vec_upper_r[i]==-99){
       cout.width(width);cout<<"NA";
     }

     else{
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
         cout<<"\t";cout.width(width);cout<<divergent_time[i*3]<<"\t";
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
       cout<<endl<<"Abbreviation: MS=Model Selection; MA=Model Averaging; CI=Confidence Interval; ps=Polymorphism Synonymous; pr=Polymorphism Replacement; ds=Divergence Synonymous; dr=Divergence Replacement; Gamma=N*s (Gamma: scaled selection coefficient (selection intensity); N: haploid effective population size; s: selection coefficient); NI=Neutrality Index (NI=(pr/dr)/(ps/ds), >1 Negative selection, <1 Positive selection); INF=Infinite; N-INF=Negative Infinite;"<<endl;
       cout<<"NULL (Not enough information for this site or both pr=0 and dr=0); NA - gamma value is not available from our calculation."<<endl;
     }
     cout<<endl<<"#End of clustering"<<endl<<endl;
     return 1;
   }

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::RunML(vector<string> pol_seq, vector<string> div_seq) {
  cout<<"****** Start RunML"<<endl;
  time_t time_startRunML = time(NULL); // Record the start time
  long N=pol_seq[0].length();
  init(N);
  int q;
  //if model number is not input from the user, default 10000 is used, else, use the user input Model_Number.
  if (model_num_flag==0) {
   Model_Number=10000;
 }

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

cout<<"Number of Individuals in Polymorphism="<<species_n<<endl;
/*  cout<<"Divergent Time: "<<divergent_time<<endl;
  cout<<endl;
*/
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
struct SiteModels rep_pol[10000];
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

    q = ClusterSubSeq(0, N-1,'S',flag_seq,sm_pol);
    if (q == -1)
    {
      cout << "Error clustering (PS is too low!)" << endl;
      throw 1;
    }
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
  cout<<"****** Start Model Averaging Polymorphism Sysnonymous"<<endl;
      time_t time_start1 = time(NULL); // Record the start time
      CI_MA(sm_pol,N);
      cout<<"End Model Averaging Polymorphism Sysnonymous (Time elapsed: ";
        time_t t2 = time(NULL)-time_start1;
        h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
        if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End Model Averaging Polymorphism Sysnonymous (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
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
  q=ClusterSubSeq(0, N-1,'R',flag_seq,rep_pol);
  if (q == -1)
  {
    cout << "Error clustering (PR is too low!)" << endl;
    throw 1;
  }
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
 cout<<"****** Start Model Averaging Polymorphism Replacement"<<endl;
	  time_t time_start111 = time(NULL); // Record the start time
	  CI_MA(rep_pol,N);
	  cout<<"End Model Averaging Polymorphism Replacement (Time elapsed: ";
     time_t t22 = time(NULL)-time_start111;
     h=t22/3600, m=(t22%3600)/60, s=t22-(t22/60)*60;
     if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End Model Averaging Polymorphism Replacement (Total Time elapsed: ";
 t2total = time(NULL)-time_startRunML;
 h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
 if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
}
vec_lower_rate_pr=vec_lower_rate;
vec_upper_rate_pr=vec_upper_rate;

  //// Fixed bug due to the new OS X 10.9 system [error: variable length array of non-POD element type 'struct SiteModels']
  ////Solution: use a very large number instead of the parameter N for the gene length, for keeping all models for each gene site, to make sure the number is larger than the gene length.
 // struct SiteModels sm_div[N];
struct SiteModels sm_div[10000];
struct SiteModels rep_div[10000];
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
  q=ClusterSubSeq(0, N-1,'S',flag_seq,sm_div);
  if (q == -1)
  {
    cout << "Error clustering (DS is too low!)" << endl;
    throw 1;
  }
  vec_SelectedModels_ds=vec_SelectedModels;
  vec_MS_rate_ds=vec_MS_rate;
  vec_MA_rate_ds=vec_MA_rate;
  if(MS_only==0 && ci_ma==1 && ds >1){
    cout<<"****** Start Model Averaging Divergence Silent"<<endl;
    	  time_t time_start111 = time(NULL); // Record the start time
        CI_MA(sm_div,N);
        cout<<"EndModel Averaging Divergence Silent (Time elapsed: ";
         time_t t22 = time(NULL)-time_start111;
         h=t22/3600, m=(t22%3600)/60, s=t22-(t22/60)*60;
         if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End Model Averaging Divergence Silent (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
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
  q= ClusterSubSeq(0, N-1,'R',flag_seq,rep_div);
  if (q == -1)
  {
    cout << "Error clustering (DR is too low!)" << endl;
    throw 1;
  }
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

  cout<<"****** Start Model Averaging Polymorphism Replacement"<<endl;
    time_t time_start111 = time(NULL); // Record the start time

    CI_MA(rep_div,N);

    cout<<"End Model Averaging Polymorphism Replacement (Time elapsed: ";
      time_t t22 = time(NULL)-time_start111;
      h=t22/3600, m=(t22%3600)/60, s=t22-(t22/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End Model Averaging Divergence Replacement (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;

}

vec_lower_rate_dr=vec_lower_rate;
vec_upper_rate_dr=vec_upper_rate;

if(!divtime_flag){
  DivergentTime(species_n, N, divergent_time);
  for (long i = 0; i < N; i++)
  {
    divergent_time_sums.push_back(0);
    divergent_time_weights.push_back(0);
  }
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
}

else{
  divergent_time=vector<double>(N,Div_time);
}




  //Do model averaging and estimate gamma
if (MS_only==0 && r_estimate==1) {
    //No matter ci_r=0 or 1, we always run it, but if no ci for gamma, do ModelAverage_r only, but not CI_UpLow_r.
    //if(pr>1.0  && dr >1.0 && ci_r==1){
  if(pr>1.0  && dr >1.0){
    if(ci_r_exact==1){
      cout<<"****** Start r_exact"<<endl;
		time_start1 = time(NULL); // Record the start time
		r_exact(rep_pol,rep_div,sm_pol,sm_div,N,species_n);
		cout<<"End r_exact (Time elapsed: ";
      t2 = time(NULL)-time_start1;
      h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End r_exact (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
}
else{
  cout<<"****** Start r_stochastic"<<endl;
		time_start1 = time(NULL); // Record the start time
		r_stochastic(rep_pol,rep_div,sm_pol,sm_div,N,species_n);
		cout<<"End r_stochastic (Time elapsed: ";
      t2 = time(NULL)-time_start1;
      h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
cout<<"End r_stochastic (Total Time elapsed: ";
  t2total = time(NULL)-time_startRunML;
  h=t2total/3600, m=(t2total%3600)/60, s=t2total-(t2total/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
else   cout<<m<<":"<<s<<")"<<endl;
}
}else{
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
  //add an option of outputing gamma value estimated by model averaged pr and dr.
if (modelAveraged_p_gamma==1) {
  SitePRF(species_n,N);
}

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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
string PRFCluster::getPolSysRep(vector<string> seq) {
cout<<"Step-getPolSysRep:"<<endl;
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
cout<<"Finish getPolSysRep: "<<endl;
cout<<"consensus seq: "<<cons_seq<<endl;
return cons_seq;
}
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::ReplaceCodon4PolSys(vector <string>& codon, vector <string>& codon_other){
cout<<"Step_ReplaceCodon4PolSys: "<<endl;
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
long PRFCluster::getDifference(string seq, int pos_start, int pos_end, char symbol) {
  long i, n = 0;
  for (i=pos_start; i<=pos_end; i++) {
    if (seq[i]==symbol) n++;
  }
  return n;
}

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::ClusterSubSeq(int pos_start, int pos_end,char symbol,int flag_seq, struct SiteModels *pointer) {
  long N = pos_end - pos_start + 1;
  long symbol_n = 0;
  if(flag_seq==0){
    symbol_n=getDifference(pol_codon_consensus, pos_start, pos_end, symbol);
  }else if(flag_seq==1){
    symbol_n=getDifference(div_codon_consensus, pos_start, pos_end, symbol);
  }
  if (N==0 || N==1 || symbol_n==0 || symbol_n==1) return -1;
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
        if ( (cs-pos_start>1 && pos_end-ce>1 && ce-cs>1) ||
          (cs==pos_start && pos_end-ce>1 && ce-cs>1) ||
          (cs-pos_start>1 && pos_end==ce && ce-cs>1)
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
}
else{
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
  if (cs_max>pos_start+1) ClusterSubSeq(pos_start, cs_max-1,symbol,flag_seq,pointer);
  if (ce_max<pos_end-1) ClusterSubSeq(ce_max+1, pos_end,symbol,flag_seq,pointer);
  if (cs_max<ce_max-1) ClusterSubSeq(cs_max, ce_max, symbol,flag_seq,pointer);
}  
return 1;
}

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::CI_MA(struct SiteModels *pointer,long N){
  vector<thread> MA_threads;
  int index=0;

  while(index<N)
  {
    int num_threads=0;
    for (int i=0;(i<NUM_CPU && index<N);i++)
    {
     MA_threads.push_back(std::thread(&PRFCluster::CI_MA_threaded, 
      this,pointer,index));
     index++;
     num_threads++;
   }
   for (long i=0;i<num_threads;i++){
    MA_threads[i].join();
  }
 	MA_threads.clear();
}
return 1;
}
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
void PRFCluster::CI_MA_threaded(struct SiteModels *pointer,long i){
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
}
/***************************************************
* Function: Get Site Specific divergence time for each site i under each model of ps and ds.
* Input Parameter:
* Output:
***************************************************/
void PRFCluster::SiteSpecificDivergentTime(long species_n, struct SiteModels* psp, struct SiteModels* psd, int model_ps, int model_ds, int i, vector<double> &div)
{
  double tmp=0.0;
  double j;
  
  for(j=1.0;j<=species_n-1;j++){
    tmp=tmp+(1/j);
  }

  if (psp[i].sms[model_ps].p ==0.0)
  {
    return;
  }
  if (psd[i].sms[model_ds].p==0.0)
  {
    return;
  }
  double t=(psd[i].sms[model_ds].p / psp[i].sms[model_ps].p)*(tmp)-1-(1/(double)species_n);
  if (t<0) return;
  div[i] = t;
}
/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
void PRFCluster::DivergentTime(long species_n, long N, vector<double> &div){
  double tmp=0.0;
  double j;
  
  for(j=1.0;j<=species_n-1;j++){
    tmp=tmp+(1/j);
  }

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

  double t=(ds/ps)*(tmp)-1-(1/(double)species_n);
  for (int i = 0; i< N; i++)
    div.push_back(t);
}

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
//Calculate the 95% confidence interval for gamma stochastically, randomly pick dr and pr, calculate gamma using ModelSpecific_r_PRF, and get gamma weight; rank gamma, and sum weight, get the 95% weight boundary gamma.
int PRFCluster::r_stochastic(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N,long species_n){
	time_t time_start=time(NULL);
	if (verbose==1){
		cout<<"Site\tPolymorphismModelSize\tDivergenceModelSize\tAllUsedModels\tModelsWithSavedGammaValues\tr_ma\tr_ma_prob\tCIr_lower\tCIr_upper"<<endl;
	}

  vector<thread> stochastic_threads;
  int index=0;
  ostringstream* thread_outputs = new ostringstream[NUM_CPU];

  while(index<N)
  {
    int num_threads=0;
    for (int i=0;(i<NUM_CPU && index<N);i++)
    {
     stochastic_threads.push_back(std::thread(&PRFCluster::r_stochastic_threaded, 
      this,pr, dr, ps, ds, N, species_n, index,&thread_outputs[i]));
     index++;
     num_threads++;
   }
   for (long i=0;i<num_threads;i++){
    stochastic_threads[i].join();
    cout << thread_outputs[i].str();
	thread_outputs[i].str("");
    thread_outputs[i].clear();
  }
  stochastic_threads.clear();
}
return 1;
}
/***************************************************
* Function: use stochastic algorithm to calculate gamma, randomly pick one pr, dr, and ps, ds (if use -ssd)
* Input Parameter:
* Output:
***************************************************/
void PRFCluster::r_stochastic_threaded(struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n, long i,std::ostringstream* myout){
  int parameters=4; // p0,pc for dr and pr, rate-within-cluster, rate-without-cluster for polymorphism and divergence
  long modelCount_all=0;
  long modelCount_gamma=0;
  vector<rModels> vec_rModels2;

  vec_rModels2.clear();
  modelCount_all=0;
  modelCount_gamma=0;
    //Each time get one pr & dr, estimate r
  long end_num=Model_Number;
  double min_weight=1000000;
  if (verbose==1){
    *myout<<i+1<<"\t"<<pr[i].sms.size()<<"\t"<<dr[i].sms.size()<<"\t";
  }
    //Sort the models by weight (pr/dr)
  sort(dr[i].sms.begin(), dr[i].sms.end(), more_than_CI());
  sort(pr[i].sms.begin(), pr[i].sms.end(), more_than_CI());
  if (site_specific_flag)
  {
    sort(ds[i].sms.begin(), ds[i].sms.end(), more_than_CI());
    sort(ps[i].sms.begin(), ps[i].sms.end(), more_than_CI());
  }
  end_num=Model_Number;

    //Get Model_Number models from pr & dr,respectively
    //If all the selected models from pr & dr are zero, flag_mutation=0, won't calculate CIs for r.
    int flag_mutation=1; //if flag_mutation=0, there is no chance to mutate, no gamma will be calculated
	int crash_flag=0;
    int repeat_num=0;
    vector<double> PpWeightSums; //polymorphism weight sum
    vector<double> DpWeightSums; //divergence weight sum
    vector<double> PsWeightSums;
    vector<double> DsWeightSums;
    //vector with SummedWeight
    PpWeightSums=RandomModel_NumFastInit(pr,i);
    DpWeightSums=RandomModel_NumFastInit(dr,i);
    if (site_specific_flag)  
    {
      PsWeightSums=RandomModel_NumFastInit(ps,i);
      DsWeightSums=RandomModel_NumFastInit(ds,i);
    }
    int TotalSampleNum=Model_Number*1000;
    int TotalModel = dr[i].sms.size();
    if (TotalModel>pr[i].sms.size()) {
      TotalModel=pr[i].sms.size();//choose the smaller one of the numer of models of pr and dr
    }
    if (TotalSampleNum>TotalModel) {TotalSampleNum=TotalModel;}//choose the smaller one of Model_Number*1000 and Total Model to go through all models, as the cutoff model number when p=0.
    int jj=0; //count all the models sampled, including both p=0 and p!=0
    int ssd_jj=0;
    for(long j=1;j<=end_num;j++){
      long model_dr=RandomModel_NumFast(DpWeightSums); // create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
      long model_pr=RandomModel_NumFast(PpWeightSums); // create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
      
      if (site_specific_flag) 
      {
        long model_ps=RandomModel_NumFast(PsWeightSums);
        long model_ds=RandomModel_NumFast(DsWeightSums);
        parameters = 8;
		  if(ds[i].sms[model_ds].p==0.0 && ps[i].sms[model_ps].p==0.0){
			continue;
			ssd_jj++;
		  }else{
			SiteSpecificDivergentTime(species_n, ps, ds, model_ps, model_ds, i, divergent_time);
			divergent_time_sums[i] += divergent_time[i];
			divergent_time_weights[i] += 1;
		  }
		  if (ssd_jj==end_num){
		  cout<<"Error in calculating site specific divergence time from silent clusters! No silent clustering for the site."<<endl;
		  throw 1;
		  }
      }//end of if (site_specific_flag)
      jj++;//count all the models sampled, including both p=0 and p!=0
      //Re-generate random sample if both pr and dr=0
	  if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0){
		  j--;
		  if (jj == TotalSampleNum) //when Model_Number samples all have p=0, continue to find random samples till TotalSampleNum
		  {
			  if (repeat_num == TotalSampleNum) {
				  flag_mutation = 0; // only when all models sampled have p=0, the site has no chance to mutation, giving the N-INF for gamma
				  *myout << "For Site " << i
					  << ", all sampled models indicate no chance to mutate!"; // Every randomly extracted model (all Model_Number) gives p=0
				  *myout << "\tThe Number of Models with dr.p=0: "
					  << repeat_num
					  << ";\tThe total Number of Models sampled: "
					  << TotalSampleNum << endl;
			  } else {
				  flag_mutation = 1; //if there are models p!=0, calculate gamma.
			  }
			  break;
		  }
		  else{
			  repeat_num++;
		  }
		  continue;
	  }//end of if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0)
	  else{
		  //cout<<model_pr<<"\t"<<pr[i].sms[model_pr].weight<<"\t"<<pr[i].sms[model_pr].p<<endl;
		  double new_r=ModelSpecific_r_PRF(pr[i].sms[model_pr].p,dr[i].sms[model_dr].p,species_n, divergent_time[i]);
		  modelCount_all++;
		  if (new_r==50 or new_r==-50 or new_r==-99 or new_r==99 or new_r==-66) {
			  j--;
			  continue;
		  }
		  else {
			  rModels tmp_rm(1,new_r);
			  vec_rModels2.push_back(tmp_rm);
			  crash_flag = 1;
			  modelCount_gamma++;
		  }
	  }//end of else pr.p=0 and dr.p=0
	}//end of for(long j=1;j<=end_num;j++)	
if (verbose==1) {*myout<<modelCount_all<<"\t"<<modelCount_gamma<<"\t";}
//Find CI and 95% MA for r if wants CI, if CI is not wanted, just do ModelAverage_r
  if(vec_rModels2.size()==0){
   vec_r[i]=0;
   vec_lower_r[i]=0;
   vec_upper_r[i]=0;
   cout<< "The gamma models are empty. Put gamma and confidence intervals of gamma as 0 for SITE " << i << endl;
   return;
 }
  else if(flag_mutation==1 and ci_r==1 and crash_flag){
  CI_UpLow_r_thread(i,min_weight,vec_rModels2);
  if (verbose==1) {
	  double r_modelAveraged_p=ModelSpecific_r_PRF(vec_MA_rate_pr[i],vec_MA_rate_dr[i],species_n, divergent_time[i]);//get the gamma based on model averaged pr and dr
	  *myout<<vec_r[i]<<"\t"<<r_modelAveraged_p <<"\t"<<vec_lower_r[i]<<"\t"<<vec_upper_r[i]<<endl;
	}
  }
  else if(flag_mutation==1 and ci_r==0 and crash_flag){
	ModelAveraged_r_thread(i,min_weight,vec_rModels2);
  }
  else{
	*myout<<"Both pr and dr are 0. gamma was not calculated for model averaged gamma, 95% CI gamma!"<<endl;
	vec_r[i]=0;
	vec_lower_r[i]=0;
	vec_upper_r[i]=0;
  }
  vec_rModels2.clear();
}

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::r_exact (struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N, long species_n)
{
  time_t time_start=time(NULL);
  if (verbose==1){
	cout<<"Site\tPolymorphismModelSize\tDivergenceModelSize\tAllUsedModels\tModelsWithSavedGammaValues\tr_ma\tr_ma_prob\tCIr_lower\tCIr_upper"<<endl;
  }
  vector<thread> exact_threads;
  int index=0;
  ostringstream* thread_outputs = new ostringstream[NUM_CPU];

  while(index<N)
  {
	int num_threads=0;
	for (int i=0;(i<NUM_CPU && index<N);i++)
	{
	 exact_threads.push_back(std::thread(&PRFCluster::r_exact_threaded, 
	  this,pr, dr, ps, ds, N, species_n, index,&thread_outputs[i]));
	 index++;
	 num_threads++;
   }
   for (long i=0;i<num_threads;i++){
	exact_threads[i].join();
	cout << thread_outputs[i].str();
	thread_outputs[i].str("");
	thread_outputs[i].clear();
  }
  exact_threads.clear();
}
return 1;
}
/***************************************************
* Function: use exactly algorithm to calculate gamma, use four nested for loops if (-ssd)
* Input Parameter:
* Output:
***************************************************/
//Calculate the 95% confidence interval for gamma exactly, randomly pick dr and pr, calculate gamma using ModelSpecific_r_PRF, and get gamma weight; rank gamma, and sum weight, get the 95% weight boundary gamma.
void PRFCluster::r_exact_threaded (struct SiteModels *pr, struct SiteModels *dr, struct SiteModels *ps, struct SiteModels *ds, long N,long species_n, long i, std::ostringstream* myout){
	int parameters=4; // p0,pc for dr and pr, rate-within-cluster, rate-without-cluster for polymorphism and divergence
	double min_weight=1000000;
	long Dr_Model_Num=0;
	long Pr_Model_Num=0;
  long Ds_Model_Num=0;
  long Ps_Model_Num=0;
  time_t time_start=time(NULL);
  long modelCount_all=0;
  long modelCount_gamma=0;
		//Get all the models
  modelCount_all=0;
  modelCount_gamma=0;
  vector<rModels> vec_rModels2;
  vec_rModels2.clear();
	Dr_Model_Num=dr[i].sms.size();//The total number of divergence replacement models
	Pr_Model_Num=pr[i].sms.size();//The total number of polymorphism replacement models
	Ds_Model_Num=ds[i].sms.size();
    Ps_Model_Num=ps[i].sms.size();
    //If all the selected models from pr & dr are zero, flag_mutation=0, won't calculate CIs for r.
    if (verbose==1){
     cout<<i+1<<"\t"<<pr[i].sms.size()<<"\t"<<dr[i].sms.size()<<"\t";
   }
	int flag_mutation=1; //if flag_mutation=0, there is no chance to mutate, no gamma will be calculated
	for(long j=0;j<Dr_Model_Num;j++){
		long model_dr=j; //use the real model
		for(long jj=0;jj<Pr_Model_Num;jj++){
			long model_pr=jj;
			if (site_specific_flag)//use site specific divergence time calculated from silent clustering
			{
			  parameters=6;
			  for(long jjj=0;jjj<Ds_Model_Num;jjj++){
				long model_ds=jjj;
				for(long jjjj=0;jjjj<Ps_Model_Num;jjjj++)
				{
				  long model_ps=jjjj;			  
				  if(ds[i].sms[model_ds].p==0.0 && ps[i].sms[model_ps].p==0.0){
					continue;
				  }else{
					  SiteSpecificDivergentTime(species_n, ps, ds, model_ps, model_ds, i,divergent_time);
				  }
					//Re-generate random sample if both pr and dr=0
				  if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0){
					continue;
				  }
				  else{
				   //cout<<model_pr<<"\t"<<pr[i].sms[model_pr].weight<<"\t"<<pr[i].sms[model_pr].p<<endl;
					double new_r=ModelSpecific_r_PRF(pr[i].sms[model_pr].p,dr[i].sms[model_dr].p,species_n, divergent_time[i]);
					modelCount_all++;
					if (new_r==50 or new_r==-50 or new_r==-99 or new_r==99 or new_r==-66) {
					  continue;
					}
					else {
						double log_prob_tmp=pr[i].sms[model_pr].LogLikelihood+dr[i].sms[model_dr].LogLikelihood+ps[i].sms[model_ps].LogLikelihood+ds[i].sms[model_ds].LogLikelihood; //already log likelihood, no log any more
						double AIC_tmp  = -2*log_prob_tmp + 2*parameters;
						double AICc_tmp = AIC_tmp;
						if (N-parameters-1>0.0) AICc_tmp += 2*parameters*(parameters+1)/(N-parameters-1);
						else AICc_tmp = 2*AIC_tmp;
						double BIC_tmp = -2*log_prob_tmp + parameters*log(double(N));

						double weight_tmp = 0;
							  //BIC
						if (criterion_type==0){
						  weight_tmp = BIC_tmp;
						}
							  //AIC
						else if (criterion_type==1){
						  weight_tmp = AIC_tmp;
						}
							  //AICc
						else if (criterion_type==2){
						  weight_tmp = AICc_tmp;
						}

						if (site_specific_flag)
						{
						  divergent_time_weights[i] += weight_tmp;
						  divergent_time_sums[i] += divergent_time[i] * weight_tmp;
						} 
						if(weight_tmp<min_weight) min_weight=weight_tmp;
						rModels tmp_rm(weight_tmp,new_r);
						vec_rModels2.push_back(tmp_rm);
						modelCount_gamma++;
					  }
					}
				}//end of for loop for jjjj
			}//end of for loop for jjj
		}//end of if (site_specific_flag)
	  else{
				  //Re-generate random sample if both pr and dr=0
			if(pr[i].sms[model_pr].p==0.0 && dr[i].sms[model_dr].p==0.0){
			  continue;
			}else{
				  //cout<<model_pr<<"\t"<<pr[i].sms[model_pr].weight<<"\t"<<pr[i].sms[model_pr].p<<endl;
			  double new_r=ModelSpecific_r_PRF(pr[i].sms[model_pr].p,dr[i].sms[model_dr].p,species_n, divergent_time[i]);
			  modelCount_all++;
			  if (new_r==50 or new_r==-50 or new_r==-99 or new_r==99 or new_r==-66) {
				continue;
			  }
			  else {
					double log_prob_tmp=pr[i].sms[model_pr].LogLikelihood+dr[i].sms[model_dr].LogLikelihood; //already log likelihood, no log any more
					double weight_tmp=-2*log_prob_tmp+2*parameters;
					if(weight_tmp<min_weight) min_weight=weight_tmp;
					rModels tmp_rm(weight_tmp,new_r);
					vec_rModels2.push_back(tmp_rm);
					modelCount_gamma++;
				  }
				}//end of else dr.p and pr.p=0
		  }//end of else if (!site_specific_flag)
		}//end of for loop jj
	  }//end of for loop j
  if (verbose==1) {cout<<modelCount_all<<"\t"<<modelCount_gamma<<"\t";}
		//Calculate the CIs for r
  if(vec_rModels2.size()==0){
   vec_r[i]=0;
   vec_lower_r[i]=0;
   vec_upper_r[i]=0;
   cout<< "The gamma models are empty. Put gamma and confidence intervals of gamma as 0 for SITE " << i << endl;
   return;
 }
 //Find CI and 95% MA for r if wants CI, if CI is not wanted, just do ModelAverage_r
 if(ci_r==1){
   CI_UpLow_r_thread(i,min_weight,vec_rModels2);
   if (verbose==1) {
	  			double r_modelAveraged_p=ModelSpecific_r_PRF(vec_MA_rate_pr[i],vec_MA_rate_dr[i],species_n,divergent_time[i]);//get the gamma based on model averaged pr and dr
	  			cout<<vec_r[i]<<"\t"<<r_modelAveraged_p <<"\t"<<vec_lower_r[i]<<"\t"<<vec_upper_r[i]<<endl;
	  		}
      }
  if(ci_r==0){
   ModelAveraged_r_thread(i,min_weight,vec_rModels2);
   if (verbose==1) {cout<<i<<"\tr: "<<vec_r[i]<<endl;}
 }
 if (i%10==0){
   if (verbose==1) {Time(time_start);}
 }
 vec_rModels2.clear();
}

/***************************************************
* Function: Record time
* Input Parameter: start time
* Output: print time elapsed
***************************************************/
void PRFCluster::Time(int time_start) {
	int h, m, s;
	cout<<"Time elapsed: ";
	time_t t22 = time(NULL)-time_start;
	h=t22/3600, m=(t22%3600)/60, s=t22-(t22/60)*60;
	if(h)  cout<<h<<"h:"<<m<<"m:"<<s<<"s)"<<endl;
else   cout<<m<<"m:"<<s<<"s)"<<endl;
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
  vector<double> tempWeightSums;
  tempWeightSums.push_back(cip->weight);
  ++cip;
  for(long i = 1; i < pss; ++i, ++cip) {
    //if (cip->weight == 0.0) cerr << "null weight " << site << " " << i << endl;
    tempWeightSums.push_back(tempWeightSums[i-1] + cip->weight);
  }
  return (tempWeightSums);
}

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
  if(!(hi == 0 || (pWeightSums[hi-1] <= rand_num && rand_num <= pWeightSums[hi]))) {
    cerr << "blew fast model num " << lo << "  " << hi <<  " " << rand_num << ": " << pWeightSums[lo] << " " << pWeightSums[hi-1] << " " << pWeightSums[hi] << endl;
    exit(1);
  }
  return(hi);
}

 /***************************************************
* Function: Calculate gamma for the specific model from polymorphism and divergence probability/likelihood
* Input Parameter: polymorphism probability, divergence probability, the number of samples
* Output: gamma
***************************************************/
double PRFCluster::ModelSpecific_r_PRF (double p_pr, double p_dr,long species_n, double site_divergent_time){
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
  double r_optimal=-99;
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
       else if (-0.5<fraction < 0.0)
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
         f=site_divergent_time+gx1-tmp*gx2; //f(gamma)
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
   //if r_optimal is outside the -50 and 50, keep them as -50 and 50
    if (r_optimal>50){
      r_optimal=50;
    }
    if (r_optimal<-50){
      r_optimal=-50;
    }

    if (flag_root==false and min_dx<=MinDx){
      new_r=r_optimal;
    }

    if(flag_root==false and min_dx>MinDx){
     new_r=-99;
   }

   return (new_r);
 }

// model averaged gamma
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
    long Gamma_Model_Size=vec_rModels.size();
    if (ci_r_exact==0){ all_weight=Gamma_Model_Size;}
  //Calculate weight for each model and get a summed weight for gamma
    for(long i=0;i<vec_rModels.size();i++){
	  //For exact,use AIC weight
     if (ci_r_exact==1){
      vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
      all_weight+=vec_rModels[i].weight;
    }
	  //For stochastic, all weights are equal
    else{
      vec_rModels[i].weight=1;
    }
  }
  //Get the weight for the first model and put it into the list of rWeightSums
  vec_rModels[0].weight=vec_rModels[0].weight/all_weight;
  rWeightSums.push_back(vec_rModels[0].weight);

  	long last_item=vec_rModels.size()-1; //The real model number should be last_item+1, but inside the list ID, it should be last_item.
  	double lci=confidence_interval;
  	double uci=1-confidence_interval;
  	double averaged_gamma=0; // average only for weights between 95% CI

  	//consider the case that the first model already exceeding the weight sum of 0.975.
  	if (rWeightSums[0]>=uci) {
  		flag_lower=1;
  		lower_model=0;
  		upper_model=0;
  		double r1=vec_rModels[0].r;
  		vec_lower_r[site]=r1;
  		vec_upper_r[site]=r1;
  		vec_r[site]=r1;
  		return 1;
  	}

  	//lower_model, for the case the first item is already above the CI; record the lower model as 0;
  	if (rWeightSums[0]>=lci) {
      lower_model=0;
      flag_lower=1;
      double r1=vec_rModels[0].r;
      double w1=vec_rModels[0].weight;
      double r2=vec_rModels[1].r;
      double w2=vec_rModels[1].weight;
		vec_lower_r[site]= r1; //the Lower 95% CI gamma, not use interpolation, use the value directly
		averaged_gamma=r1*(w1-lci);//gamma, use interpolation, lci gamma*weight plus the first model gamma*weight
 }

  	//***Debug the model averaged gamma and lower CI gamma values
  	////cout<<"ModelID:\t"<<"ModelAveragedGamma:\t"<<"lowerCI:\t"<< "UpperCI:\t"<<"gamma[i]:\t"<< "weight[i]:\t"<<"gamma[i-1]:\t"<<"SummedWeight"<<endl;
    //starting from the second model for lower and upper CI_r and model_averaged_r
 for(long i=1;i<vec_rModels.size();i++){
  vec_rModels[i].weight=vec_rModels[i].weight/all_weight;
  double tmp_weight=vec_rModels[i].weight;
  double SummedWeight=rWeightSums[i-1] + tmp_weight;
  		//if (tmp_weight == 0.0) cerr << "null weight " << site << " " << i << endl;
  		rWeightSums.push_back(SummedWeight); // summed weight for Model i pushed into the list

  		//Get the lower CI gamma; record the model number for the rWeightSums at the two borders of the lower 95% CI lci 0.025, lower_model
  		if (rWeightSums[i-1]<lci && rWeightSums[i]>=lci) {
  			lower_model=i;
  			flag_lower=1;
  			vec_lower_r[site] = vec_rModels[i].r;
  			averaged_gamma=(rWeightSums[i]-lci)*vec_rModels[i].r;

  			//consider the case that the same model exceeding the weight sum of 0.975 and 0.025 at the same time, the all CI_r and gamma equal to the same model gamma.
  			if (rWeightSums[i]>=uci){
  				upper_model=i;
  				double r0=vec_rModels[i].r;
          vec_lower_r[site]=r0;
          vec_upper_r[site]=r0;
          vec_r[site]=r0;
          return 1;
        }
      }

  		//getting model averaged gamma between 95% CI
      if (flag_lower==1 && i>lower_model && rWeightSums[i]<uci){
       averaged_gamma+=vec_rModels[i].weight*vec_rModels[i].r;
     }
  		//Get the upper CI gamma; record the model number for the rWeightSums at the two sides of the upper 95% CI uci 0.975, upper_model
     if (flag_lower==1 && i>lower_model && rWeightSums[i-1]<uci && rWeightSums[i]>=uci) {
       upper_model=i;
       vec_upper_r[site] = vec_rModels[upper_model].r;
       averaged_gamma+=vec_rModels[i].r *(uci-rWeightSums[i-1]);
     }

  		//***Debug the model averaged gamma and lower CI gamma values
  		////cout<<i<<"\t"<< averaged_gamma<<"\t"<< vec_lower_r[site]<<"\t"<< vec_upper_r[site]<<"\t"<< vec_rModels[i].r<<"\t"<<tmp_weight<<"\t"<<vec_rModels[i-1].r<<"\t"<<rWeightSums[i]<<endl;
   }

   if (averaged_gamma!=0) {
		vec_r[site]=averaged_gamma/(uci-lci); // model averaged gamma only takes from lci to uci, so it should divide (uci-lci).
	}
	else {
		cout<<"Model averaged gamma: NULL\t"<<averaged_gamma<<endl;
    vec_r[site]=-66;
  }
  if (vec_r[site]!=-66 )
  {
  		//Make sure the LowerCI and UpperCI is two sides of value of model average; quit the program if model averaged gamma is not between 95% CI gamma

    if ( lower_model==-1 || upper_model==-1 || round(vec_lower_r[site])>round(vec_r[site]) || round(vec_upper_r[site]) < round(vec_r[site]) || flag_lower!=1 || upper_model<lower_model)
    {
     cout<<endl<<endl<<"Warning: Problem in getting proper values of model averaged gamma, lower and upper 95% CI gamma from all models!!!!"<<endl;
     cout<<"Site: "<<site<<"\tEstimatedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<"\nTotalModelNO: "<<last_item+1<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
     cout<<"1st Model r:\t"<<vec_rModels[0].r<<"\tSumWeight:\t"<<rWeightSums[0]<<"\t2rd Model r:\t"<<vec_rModels[1].r<<"\tSumWeight:\t"<<rWeightSums[1]<<"\tThe second last Model "<<last_item<<" r:\t"<<vec_rModels[last_item-1].r<<"\tSumWeight:\t"<<rWeightSums[last_item-1]<<"\tThe last Model "<<last_item+1<<"r:\t"<<vec_rModels[last_item].r<<"\tSumWeight:\t"<<rWeightSums[last_item]<<endl;
     cout<<"Lower95% CI r weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
     throw "Error in getting the right 95% confidence interval gamma in CI_UpLow_rc!";
   }
 }
 else
 {
  vec_lower_r[site]=-66;
  vec_upper_r[site]=-66;
  cout<<"Warning in getting proper gamma and 95% CI gamma: all -66"<<endl;
  cout<<"Site: "<<site<<"\tTotalModelNO: "<<last_item+1<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
  cout<<"Lower95% CI r weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
  cout<<"1st Model r:\t"<<vec_rModels[0].r<<"\tSumWeight:\t"<<rWeightSums[0]<<"\t2rd Model r:\t"<<vec_rModels[1].r<<"\tSumWeight:\t"<<rWeightSums[1]<<"\tThe second last Model "<<last_item<<" r:\t"<<vec_rModels[last_item-1].r<<"\tSumWeight:\t"<<rWeightSums[last_item-1]<<"\tThe last Model "<<last_item+1<<"r:\t"<<vec_rModels[last_item].r<<"\tSumWeight:\t"<<rWeightSums[last_item]<<endl;
  		//throw "Error in getting the right gamma in CI_UpLow_rc!";
}
return 1;
}

int PRFCluster::CI_UpLow_r_thread(long site,double min_weight, vector<rModels> vec_rModels) {
	double all_weight=0.0;
	long lower_model=-1;
	long upper_model=-1;
	vector<double> rWeightSums;

	rWeightSums.clear();
  int flag_lower=0;
  double lower=0.00, upper=0.00;
  	if (vec_rModels.size() == 0) return 1;

    sort(vec_rModels.begin(), vec_rModels.end()); // sort by gamma values incrementally
    long Gamma_Model_Size=vec_rModels.size();
    if (ci_r_exact==0){ all_weight=Gamma_Model_Size;}
  //Calculate weight for each model and get a summed weight for gamma
    for(long i=0;i<vec_rModels.size();i++){
	  //For exact,use AIC weight
     if (ci_r_exact==1){
      vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
      all_weight+=vec_rModels[i].weight;
    }
	  //For stochastic, all weights are equal
    else{
      vec_rModels[i].weight=1;
    }
  }
  //Get the weight for the first model and put it into the list of rWeightSums
  vec_rModels[0].weight=vec_rModels[0].weight/all_weight;
  rWeightSums.push_back(vec_rModels[0].weight);

  	long last_item=vec_rModels.size()-1; //The real model number should be last_item+1, but inside the list ID, it should be last_item.
  	double lci=confidence_interval;
  	double uci=1-confidence_interval;
  	double averaged_gamma=0; // average only for weights between 95% CI

  	//consider the case that the first model already exceeding the weight sum of 0.975.
  	if (rWeightSums[0]>=uci) {
  		flag_lower=1;
  		lower_model=0;
  		upper_model=0;
  		double r1=vec_rModels[0].r;
  		vec_lower_r[site]=r1;
  		vec_upper_r[site]=r1;
  		vec_r[site]=r1;
  		return 1;
  	}

  	//lower_model, for the case the first item is already above the CI; record the lower model as 0;
  	if (rWeightSums[0]>=lci) {
      lower_model=0;
      flag_lower=1;
      double r1=vec_rModels[0].r;
      double w1=vec_rModels[0].weight;
      double r2=vec_rModels[1].r;
      double w2=vec_rModels[1].weight;
		vec_lower_r[site]= r1; //the Lower 95% CI gamma, not use interpolation, use the value directly
		averaged_gamma=r1*(w1-lci);//gamma, use interpolation, lci gamma*weight plus the first model gamma*weight
 }

  	//***Debug the model averaged gamma and lower CI gamma values
  	////cout<<"ModelID:\t"<<"ModelAveragedGamma:\t"<<"lowerCI:\t"<< "UpperCI:\t"<<"gamma[i]:\t"<< "weight[i]:\t"<<"gamma[i-1]:\t"<<"SummedWeight"<<endl;
    //starting from the second model for lower and upper CI_r and model_averaged_r
 for(long i=1;i<vec_rModels.size();i++){
  vec_rModels[i].weight=vec_rModels[i].weight/all_weight;
  double tmp_weight=vec_rModels[i].weight;
  double SummedWeight=rWeightSums[i-1] + tmp_weight;
  		//if (tmp_weight == 0.0) cerr << "null weight " << site << " " << i << endl;
  		rWeightSums.push_back(SummedWeight); // summed weight for Model i pushed into the list

  		//Get the lower CI gamma; record the model number for the rWeightSums at the two borders of the lower 95% CI lci 0.025, lower_model
  		if (rWeightSums[i-1]<lci && rWeightSums[i]>=lci) {
  			lower_model=i;
  			flag_lower=1;
  			vec_lower_r[site] = vec_rModels[i].r;
  			averaged_gamma=(rWeightSums[i]-lci)*vec_rModels[i].r;

  			//consider the case that the same model exceeding the weight sum of 0.975 and 0.025 at the same time, the all CI_r and gamma equal to the same model gamma.
  			if (rWeightSums[i]>=uci){
  				upper_model=i;
  				double r0=vec_rModels[i].r;
          vec_lower_r[site]=r0;
          vec_upper_r[site]=r0;
          vec_r[site]=r0;
          return 1;
        }
      }

  		//getting model averaged gamma between 95% CI
      if (flag_lower==1 && i>lower_model && rWeightSums[i]<uci){
       averaged_gamma+=vec_rModels[i].weight*vec_rModels[i].r;
     }
  		//Get the upper CI gamma; record the model number for the rWeightSums at the two sides of the upper 95% CI uci 0.975, upper_model
     if (flag_lower==1 && i>lower_model && rWeightSums[i-1]<uci && rWeightSums[i]>=uci) {
       upper_model=i;
       vec_upper_r[site] = vec_rModels[upper_model].r;
       averaged_gamma+=vec_rModels[i].r *(uci-rWeightSums[i-1]);
     }

  		//***Debug the model averaged gamma and lower CI gamma values
  		////cout<<i<<"\t"<< averaged_gamma<<"\t"<< vec_lower_r[site]<<"\t"<< vec_upper_r[site]<<"\t"<< vec_rModels[i].r<<"\t"<<tmp_weight<<"\t"<<vec_rModels[i-1].r<<"\t"<<rWeightSums[i]<<endl;
   }

   if (averaged_gamma!=0) {
		vec_r[site]=averaged_gamma/(uci-lci); // model averaged gamma only takes from lci to uci, so it should divide (uci-lci).
	}
	else {
		cout<<"Model averaged gamma: NULL\t"<<averaged_gamma<<endl;
    vec_r[site]=-66;
  }
  if (vec_r[site]!=-66 )
  {
  		//Make sure the LowerCI and UpperCI is two sides of value of model average; quit the program if model averaged gamma is not between 95% CI gamma

    if ( lower_model==-1 || upper_model==-1 || round(vec_lower_r[site])>round(vec_r[site]) || round(vec_upper_r[site]) < round(vec_r[site]) || flag_lower!=1 || upper_model<lower_model)
    {
     cout<<endl<<endl<<"Warning: Problem in getting proper values of model averaged gamma, lower and upper 95% CI gamma from all models!!!!"<<endl;
     cout<<"Site: "<<site<<"\tEstimatedGamma: "<<vec_r[site]<<"\tLower 95% CI gamma: "<<vec_lower_r[site]<<"\tUpper CI gamma: "<<vec_upper_r[site]<<"\nTotalModelNO: "<<last_item+1<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
     cout<<"1st Model r:\t"<<vec_rModels[0].r<<"\tSumWeight:\t"<<rWeightSums[0]<<"\t2rd Model r:\t"<<vec_rModels[1].r<<"\tSumWeight:\t"<<rWeightSums[1]<<"\tThe second last Model "<<last_item<<" r:\t"<<vec_rModels[last_item-1].r<<"\tSumWeight:\t"<<rWeightSums[last_item-1]<<"\tThe last Model "<<last_item+1<<"r:\t"<<vec_rModels[last_item].r<<"\tSumWeight:\t"<<rWeightSums[last_item]<<endl;
     cout<<"Lower95% CI r weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
     throw "Error in getting the right 95% confidence interval gamma in CI_UpLow_rc!";
   }
 }
 else
 {
  vec_lower_r[site]=-66;
  vec_upper_r[site]=-66;
  cout<<"Warning in getting proper gamma and 95% CI gamma: all -66"<<endl;
  cout<<"Site: "<<site<<"\tTotalModelNO: "<<last_item+1<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
  cout<<"Lower95% CI r weight: "<<rWeightSums[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums[upper_model]<<endl;
  cout<<"1st Model r:\t"<<vec_rModels[0].r<<"\tSumWeight:\t"<<rWeightSums[0]<<"\t2rd Model r:\t"<<vec_rModels[1].r<<"\tSumWeight:\t"<<rWeightSums[1]<<"\tThe second last Model "<<last_item<<" r:\t"<<vec_rModels[last_item-1].r<<"\tSumWeight:\t"<<rWeightSums[last_item-1]<<"\tThe last Model "<<last_item+1<<"r:\t"<<vec_rModels[last_item].r<<"\tSumWeight:\t"<<rWeightSums[last_item]<<endl;
  		//throw "Error in getting the right gamma in CI_UpLow_rc!";
}
return 1;
}

// model averaged gamma
/***************************************************
* Function: Calculate  model averaged gamma
* Function: Use the same one for both Stochastic and Exact to find  model averaged gamma
* Input Parameter: SiteModels; site; min_weight_c
* Output: model averaged gamma
* Return Value: int 1
***************************************************/
int PRFCluster::ModelAveraged_r(long site,double min_weight) {
	double all_weight=0.0;
	rWeightSums.clear();
    sort(vec_rModels.begin(), vec_rModels.end()); // sort by gamma values incrementally
  //Calculate weight for each model and get a summed weight for gamma
    for(long i=0;i<vec_rModels.size();i++){
     vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
     all_weight+=vec_rModels[i].weight;
   }
   double averaged_gamma=0;
   double count_0s=0;
   double gamma=0;

   for(long i=0;i<vec_rModels.size();i++){
	  // Skip model averaging for models with special gamma values
     if (vec_rModels[i].r==50 or vec_rModels[i].r==-50 or vec_rModels[i].r==99 or vec_rModels[i].r==-99 or vec_rModels[i].r==-66 ){
      count_0s++;
      gamma+=vec_rModels[i].r;
      continue;
    }
    else {
      averaged_gamma+=vec_rModels[i].weight*vec_rModels[i].r;
    }
  }
  //If all gamma values are not valid for model averaging, keep that gamma value for the site.
  if (count_0s==vec_rModels.size()){
    vec_r[site]=gamma/count_0s;
    cout<<"For Site: "<<site<<" no valid gamma value averaged. Recorded averaged gamma value: "<<vec_r[site]<<endl;
  }

  if (averaged_gamma==0){
    vec_r[site]=0;
    vec_lower_r[site]=0;
    vec_upper_r[site]=0;
  }
  else{
   vec_r[site]=averaged_gamma;
 }
 return 1;
}
int PRFCluster::ModelAveraged_r_thread(long site,double min_weight,vector<rModels> vec_rModels) {
	double all_weight=0.0;
	rWeightSums.clear();
    sort(vec_rModels.begin(), vec_rModels.end()); // sort by gamma values incrementally
  //Calculate weight for each model and get a summed weight for gamma
    for(long i=0;i<vec_rModels.size();i++){
     vec_rModels[i].weight=exp(-0.5*(vec_rModels[i].weight-min_weight));
     all_weight+=vec_rModels[i].weight;
   }
   double averaged_gamma=0;
   double count_0s=0;
   double gamma=0;

   for(long i=0;i<vec_rModels.size();i++){
	  // Skip model averaging for models with special gamma values
     if (vec_rModels[i].r==50 or vec_rModels[i].r==-50 or vec_rModels[i].r==99 or vec_rModels[i].r==-99 or vec_rModels[i].r==-66 ){
      count_0s++;
      gamma+=vec_rModels[i].r;
      continue;
    }
    else {
      averaged_gamma+=vec_rModels[i].weight*vec_rModels[i].r;
    }
  }
  //If all gamma values are not valid for model averaging, keep that gamma value for the site.
  if (count_0s==vec_rModels.size()){
    vec_r[site]=gamma/count_0s;
    cout<<"For Site: "<<site<<" no valid gamma value averaged. Recorded averaged gamma value: "<<vec_r[site]<<endl;
  }

  if (averaged_gamma==0){
    vec_r[site]=0;
    vec_lower_r[site]=0;
    vec_upper_r[site]=0;
  }
  else{
   vec_r[site]=averaged_gamma;
 }
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
  site_specific_flag=1; //Use -SSD site-specific divergence time as a default option
  try {

    if (argc==2) {
      temp = stringtoUpper(argv[1]);
      if (temp=="-H") showHelpInfo();
      else {
       cout<<"Errors in argc!\n";
       throw 1;
     }
   }
   /*
   else if (argc!=5 && argc!=7 && argc!=9 && argc!=11 && argc!=13 && argc!=15 && argc!=17 && argc!=19 && argc!=21 && argc!=23 && argc!=25 && argc!=27 && argc!=29 && argc!=31 && argc!=33 && argc!=35 && argc!=37 && argc!=39 && argc!=41) {
    cout<<"argc number messed!\n";
    throw 1;
  }
  */
  else {
      //parse parameters
    int modelAveraged_p_gamma_flag=0, pol_num_flag=0, input_consensus_flag=0, output_flag=0, pol_flag=0, div_flag=0, code_flag=0, criterion_flag=0,ms_flag=0, synonymous_flag=0,ci_ma_flag=0,r_flag=0,ci_r_flag=0,ci_method_flag=0, nuc_replace_flag=0, ni_flag=0;
    model_num_flag=0;
    int verbose_flag=0;
    verbose=1;
    for (i=1; i<argc; i++) {
    	  //cout<<"argc: "<<argc<<" argv: "<<argv[i]<<endl;
    	 // cout<<"Parsing parameters i"<<i<<"\n";
     temp = stringtoUpper(argv[i]);
	//Polymorphic input fasta file
     if (temp=="-P" && (i+1)<argc && pol_flag==0) {
		//cout<<"Parsing parameters "<<temp<<"\n";
       pol_seqfile = argv[++i];
       pol_flag++;
     }
	//Divergent input fasta file
     else if (temp=="-D" && (i+1)<argc && div_flag==0) {
		//cout<<"Parsing parameters "<<temp<<"\n";
       div_seqfile = argv[++i];
       div_flag++;
     }
	//Use model averaged pr and dr calculate gamma
     else if (temp=="-RMAP" && (i+1)<argc && modelAveraged_p_gamma_flag==0) {
		//cout<<"Parsing parameters "<<temp<<"\n";
      int tmp = CONVERT<int>(argv[++i]);
		//cout<<"parameter for rMAp is "<<tmp<<endl;
      if (tmp==0) {
       modelAveraged_p_gamma_flag=0;
       modelAveraged_p_gamma=0;
     }
     else if (tmp==1)
     {
       modelAveraged_p_gamma_flag=1;
       modelAveraged_p_gamma=1;
     }
     else {
       cout<< "Error! The rMAp - parameter should be 0 or 1!\n";
       throw 1;
     }
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

	//model number for stochastic algorithm for gamma - the number of samples, as an option when stochastic method is used, default is 10000.
else if (temp=="-MN" && (i+1)<argc && model_num_flag==0) {
  Model_Number = CONVERT<int>(argv[++i]);
  model_num_flag++;
  if (Model_Number<10000){
   cout<< "Warning! The model number for stochastic algorithm should be bigger than 10000, 10000 is uded instead!\n";
   Model_Number=10000;
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
 if (num>=0 && num<=2) {
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
//Utilize site specific divergence time, calculated using silent clustering
else if (temp=="-SSD"){
  cout<<"Utilize site specific divergence time, calculated using silent clustering."<<endl;
  site_specific_flag=1;
  if (divtime_flag)
  {
    throw 1;
  }
}
//User input species divergence time
else if (temp=="-T" && (i+1)<argc && divtime_flag==0){
  double num=CONVERT<double>(argv[++i]);
  Div_time = num;
  divtime_flag=1;
  site_specific_flag=0;
  if (site_specific_flag)
  {
    cout<<"Error! User specified divergence time -T could not be used together with site-specific divergence time option -SSD."<<endl;
    throw 1;
  }
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

	//verbose output or concise
else if(temp=="-V" && (i+1)<argc && verbose_flag==0){
  int num=CONVERT<int>(argv[++i]);
  if(num==0 || num==1){
    verbose=num;
    verbose_flag=1;
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

/***************************************************
* Function:
* Input Parameter:
* Output:
***************************************************/
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
  cout<<"  -mn\tInput the number of models to be used in stochastic algorithm [integer, optional when -exact=0]"<<endl;

  cout<<"  -c\tCriterion used for clustering [integer, optional], {0:BIC || 1:AIC || 2:AICc}, default = 0"<<endl;
  cout<<"  -g\tGenetic code used for sequences [integer, optional], {1:standard}, default = 1"<<endl;
  cout<<"    \tMore information about the genetic codes can be found at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"<<endl;
  cout<<"  -m\tModel selection and model averaging  [integer, optional], {0: use both model selection and model averaging || 1: use only model selection}, default = 0"<<endl;
  cout<<"  -ci_m\tCalculate 95% confidence intervals for results of model averaging [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 0"<<endl;
  cout<<"  -s\tShow clustering results of synonymous sites from Polymorphism and Divergent sequences [integer, optional], {0: without clustering results of synonymous sites || 1: with clustering results of synonymous and replacement sites}, default = 0"<<endl;
  cout<<"  -ssd\tDefault option: use site-specific divergence time from silent site clustering (cannot be used in conjunction with -t flag)"<<endl;
  cout<<"  -r\tEstimate selection coefficient for each site [integer, optional], {0: NOT estimate selection coefficient || 1: estimate selection coefficient}, default=1"<<endl;
  cout<<"  -ci_r\tCalculate 95% confidence intervals for selection coefficient [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 1"<<endl;
  cout<<"  -exact\tAlgorithm for calculating 95% confidence intervals for selection coefficient [integer, optional], {0: use stochastic algorithm || 1: use exact algorithm}, default = 0"<<endl;
  cout<<"  -t\tUser specified species divergence time [Default: estimate site-specific divergence time using the option -ssd (cannot be used in conjunction with -ssd flag)]."<<endl;
  cout<<"  -n\tNucleotide is replaced or seen as gap when it is not A, T, G or C in the sequences [integer, optional], {0: see it as gap || 1: replace this nucleotide with the most frequently used nucleotide in other sequences}, default = 1"<<endl;
  cout<<"  -NI\tEstimate the Neutrality Index for each site [integer, optional], {0: NOT estimate Neutrality Index || 1: estimate Neutrality Index}, default=0"<<endl;
  cout<<"  -v\tVerbose output or not [integer, optional], {0: not verbose, concise output || 1: verbose output}, default=1"<<endl;

  cout<<"  -rMAp\tOutput gamma calculated using model averaged pr and dr [integer, optional], {0: not used model averaged pr and dr || 1: yes, default=0}"<<endl;

  cout<<"  -h\tShow help information"<<endl;
  cout<<endl;

  cout<<"COPYRIGHT & LICENSE:"<<endl;
  cout<<NAME<<"Source codes are under a Creative Commons CC BY-NC license."<<endl;
  cout<<"The codes can be attributed but not for commercial purpose."<<endl;
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
/***************************************************
* Function:use model averaged pr and dr to estimate gamma
* Input Parameter:
* Output:
***************************************************/
int PRFCluster::SitePRF(int species_n,long N){
  //Estimate r for each site one by one
  for(long i=0; i<N; i++){
   double new_r=ModelSpecific_r_PRF(vec_MA_rate_pr[i],vec_MA_rate_dr[i],species_n, divergent_time[i]);
   vec_r[i]=new_r;
 }
 return 1;
}