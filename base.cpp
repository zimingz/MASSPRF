/************************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: base.cpp
* Abstract: Definition of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

*************************************************************/

#include "base.h"


/******** Global variables ********/


/*						The Genetic Codes 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
	Last update of the Genetic Codes: October 05, 2000 */
int genetic_code=1; //from 1 to 23
/* Genetic standard codon table, !=stop codon */
const char* transl_table[] = {
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "1-Standard Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG", "2-Vertebrate Mitochondrial Code",
 "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "3-Yeast Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "4-Mold Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "5-Invertebrate Mitochondrial Code",
 "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "6-Ciliate, Dasycladacean and Hexamita Code",
 "", "7-",
 "", "8-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "9-Echinoderm and Flatworm Mitochondrial Code",
 "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "10-Euplotid Nuclear Code",
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "11-Bacterial and Plant Plastid Code",
 "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "12-Alternative Yeast Nuclear Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "13-Ascidian Mitochondrial Code",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "14-Alternative Flatworm Mitochondrial Code",
 "FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "15-Blepharisma Nuclear Code",
 "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "16-Chlorophycean Mitochondrial Code",
 "", "17-",
 "", "18-",
 "", "19-",
 "", "20-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "21-Trematode Mitochondrial Code",
 "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "22-Scenedesmus obliquus mitochondrial Code",
 "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "23-Thraustochytrium Mitochondrial Code"
};

//********End of Global variables**********

//Constructor function
Base::Base() {
	srand((unsigned)time(NULL));
}




/********************************************
* Function: addString
* Input Parameter: string, string, string
* Output: result = result + str + flag
* Return Value: void
* Note: flag = "\t" (default) or "\n"
*********************************************/
void Base::addString(string &result, string str, string flag) {
	result += str;
	result += flag;
}

/**********************************************************************
* Function: getAminoAcid
* Input Parameter: codon or codon's id
* Output: Calculate the amino acid according to codon or codon's id.
* Return Value: char 
***********************************************************************/
char Base::getAminoAcid(string codon) {
	return transl_table[2*(genetic_code-1)][getID(codon)];
}

char Base::getAminoAcid(int id) {
	return transl_table[2*(genetic_code-1)][id];
}

/**********************************
* Function: getNumNonsense
* Input Parameter: int
* Output: get the number of nonsense codons
* Return Value: int
***********************************/
int Base::getNumNonsense(int genetic_code) {

	int num, i;
	for(num=i=0; i<CODONCOUNT; i++) {
		if(getAminoAcid(i)=='!') num++;
	}

	return num;
}

/********************************************
* Function: getID
* Input Parameter: codon
* Output: Get codon's id in array of codon_table.
* Return Value: int
*********************************************/
int Base::getID(string codon) {
	return (convertChar(codon[0])*XSIZE + convertChar(codon[1])*DNASIZE + convertChar(codon[2]));
}

/********************************************
* Function: getCodon
* Input Parameter: int
* Output: Get the codon according to id;
		  a reverse funtion of getID.
* Return Value: string
*********************************************/
string Base::getCodon(int IDcodon) {
	
	string codon = "TTT";

	if (IDcodon>=0 && IDcodon<64) {
		codon[0]=convertInt(IDcodon/16); 
		codon[1]=convertInt((IDcodon%16)/4);
		codon[2]=convertInt(IDcodon%4);
	}

	return codon;
}

/*********************************************
* Function: convertChar
* Input Parameter: ch as char
* Output: Convert a char-T,C,A,G into a digit
*		  0,1,2,3, respectively.
* Return Value: int.
**********************************************/
int Base::convertChar(char ch) {
	int ret = -1;
	switch(ch) {
		case 'T':case 'U':
			ret = 0;
			break;
		case 'C':
			ret = 1;
			break;
		case 'A':
			ret = 2;
			break;
		case 'G':
			ret = 3;
			break;
	}
	return ret;
}

/********************************************
* Function: convertInt
* Input Parameter: int
* Output: Convert a digit- 0,1,2,3 into a 
*		  char-T,C,A,G, respectively.
* Return Value: char
*********************************************/
char Base::convertInt(int i) {
	char ch = '-';
	switch(i) {
		case 0:
			ch = 'T';
			break;
		case 1:
			ch = 'C';
			break;
		case 2:
			ch = 'A';
			break;
		case 3:
			ch = 'G';
			break;
	}
	return ch;
}

/********************************************
* Function: stringtoUpper
* Input Parameter: string
* Output: upper string
* Return Value: string
*********************************************/
string Base::stringtoUpper(string str) {
	int i;	
	for (i=0; i<str.length(); i++) str[i] = toupper(str[i]);
	return str;
}


string Base::Trim(string str) {
	int i, flag;
	//left trim
	flag=0;
	for (i=0; i<str.length() && flag==0;) {
		int c = (char)str[i];
		if (iscntrl(c) || c<33 || c>127) i++;
		else flag=1;
	}
	str = str.substr(i, str.length()-i);

	//right trim
	flag=0;
	for (i=str.length()-1; i>=0 && flag==0;) {
		int c = (char)str[i];
		if (iscntrl(c) || c<33 || c>127) i--;
		else flag=1;
	}
	str = str.substr(0, i+1);

	return str;
}

string Base::filterString(string str) {
	int i=0;
	string tmp="";
	for (;i<str.length();i++) {		
        int c = (char)(str[i]);
		if (c>32 && c<128) tmp += str[i];
	}
	return tmp;
}

string Base::replaceAll(string str, char from, char to) {
	int i=0;
	string tmp="";
	while (i<str.length()) {
		if (str[i]!=from) tmp += str[i];
		else tmp += to;
		i++;
	}

	return tmp;
}

string Base::replaceAll(string str, char from) {
	int i=0;
	string tmp="";
	while (i<str.length()) {
		if (str[i]!=from) tmp += str[i];		
		i++;
	}

	return tmp;
}
/********************************************
* Function: getRandom
* Input Parameter: void
* Output: Generate a radnom integer
* Return Value: int
*********************************************/
int Base::getRandom() {	
	return rand();
}

double Base::getRandDouble() {
	return (double)rand()/RAND_MAX;
}

/*******************

Note: begin<= getRandInt() < end

*******************/
int Base::getRandInt(int begin, int end) {
	return (int)((double)(end-begin)*getRandDouble()+begin);	
}

/********************************************
* Function: initArray
* Input Parameter: array of int/double, int, int/double(default=0)
* Output: Init the array x[0...n-1]=value
* Return Value: int
*********************************************/
int Base::initArray(double x[], int n, double value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

int Base::initArray(int x[], int n, int value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

/********************************************
* Function: sumArray
* Input Parameter: double/int, int, int(default=0)
* Output: Sum of array x[]
* Return Value: double/int
*********************************************/
double Base::sumArray(double x[], int end, int begin) { 
	int i;
	double sum=0.;	
	for(i=begin; i<end; sum += x[i], i++);    
	return sum;
}

int Base::sumArray(int x[], int end, int begin) { 
	int i, sum=0;	
	for(i=begin; i<end; sum += x[i], i++);    
	return sum;
}

/********************************************
* Function: norm
* Input Parameter: array of double, int
* Output: Sqrt of the sum of the elements' square 
           sqrt(x0*x0 + x1*x1 + ...)
* Return Value: double
*********************************************/
double Base::norm(double x[], int n) {
	int i; 
	double t=0; 

	for(i=0; i<n; t += square(x[i]), i++);

	return sqrt(t);
}

/********************************************
* Function: scaleArray
* Input Parameter: double, array of double, int
* Output: Elements in array are mutipled by scale 
* Return Value: int
*********************************************/
int Base::scaleArray(double scale, double x[], int n) {	
	int i; 	
	for (i=0; i<n; i++) x[i] *= scale;

	return 1; 
}

/********************************************
* Function: copyArray
* Input Parameter: array, array, int
* Output: Copy array's values one by one: to[] = from[]
* Return Value: int
*********************************************/
int Base::copyArray(double from[], double to[], int n) {	
	int i; 
	for(i=0; i<n; i++) to[i] = from[i];
	
	return 1; 
}

int Base::countChar(string str, char ch) {	
	int i, num; 
	for(i=num=0; i<str.length(); i++) if (str[i]==ch) num++;	
	return num; 
}


/********************************************
* Function: innerp
* Input Parameter: array, array, int
* Output: Sum of 'n' products multiplied by 
			two elements in x[], y[].
* Return Value: int
*********************************************/
double Base::innerp(double x[], double y[], int n) {
	
	int i; 
	double t=0;

	for(i=0; i<n; t += x[i]*y[i], i++); 

	return t; 
}

/********************************************
* Function: initIdentityMatrix
* Input Parameter: array of double, int
* Output: Set x[i,j]=0 when x!=j and 
			  x[i,j]=1 when x=j 
* Return Value: int
*********************************************/
int Base::initIdentityMatrix(double x[], int n) {
	
	int i,j;

	for (i=0; i<n; i++)  {
		for(j=0; j<n; x[i*n+j]=0, j++);
		x[i*n+i] = 1; 
	} 

	return 0; 
}


string Base::RGB(int r, int g, int b) {	
	char buffer[3];	

	sprintf(buffer, "%X", r);
	string rr = buffer;	
	if(rr.length()==1) rr = "0" + rr;

	sprintf(buffer, "%X", g);
	string gg = buffer;	
	if(gg.length()==1) gg = "0" + gg;

	sprintf(buffer, "%X", b);
	string bb = buffer;	
	if(bb.length()==1) bb = "0" + bb;

	return rr+gg+bb;	
}

/************************************************
* Function: writeFile
* Input Parameter: string, string
* Output: Write content into the given file.
* Return Value: True if succeed, otherwise false. 
*************************************************/
bool Base::writeFile(string output_filename, const char* result) {
	
	bool flag = true;
	try {
		//file name is ok
		if (output_filename!="" && output_filename.length()>0) {
			ofstream os(output_filename.c_str());
			if (!os.is_open()) throw 1;

			os<<result;
			os.close();					
		}
	}
	catch (...) {
		cout<<"Error in writing to file..."<<endl;
		flag = false;
	}	

	return flag;
}

bool Base::removeCntrl(string name, vector<string> &str) {
	bool flag = true;
	int i;

	try {
		//Check whether (sequence length)/3==0
		for (i=0; i<str.size(); i++) {
			str[i] = stringtoUpper(filterString(str[i]));
			if (str[0].length()!=str[i].length()) {
				cout<<endl<<"Error. The length of sequences "<<"'"<<name<<"' is not equal."<<endl;
				throw 1;
			}			
		}	
	}
	catch (...) {
		flag = false;
	}
	
	return flag;
}

bool Base::checkValid(string name, vector<string> &str) {

	bool flag = true;
	long i;

	try {

		if (removeCntrl(name, str)) {

			string str1 = str[0];
			string str2 = str[1];
			
			//Delete gap and stop codon
			bool found;
			int j;
			for(i=0; i<str1.length(); i+=3) {			
				for(found=false, j=0; j<3 && !found; j++) {
					if (str1[j+i]=='-' || str2[j+i]=='-') {
						found = true;
					}

					if (convertChar(str1[i+j])==-1 || convertChar(str2[i+j])==-1) {
						found = true;
					}
				}

				if ((getAminoAcid(str1.substr(i,3))=='!') || (getAminoAcid(str2.substr(i,3))=='!')) {
					found = true;
				}

				if (found) {
					str1 = str1.replace(i, 3, "");
					str2 = str2.replace(i, 3, "");
					i -= 3;
				}
			}

			str[0] = str1;
			str[1] = str2;

		}
		else {
			flag = true;
		}
	}
	catch (...) {
		flag = false;
	}
	
	return flag;
}

int Base::readAXT(string seqfile, vector<string> &seqname, vector<string> &seq) {

	int flag = 1;

	try	{

		ifstream is(seqfile.c_str());
		if (!is) {
			cout<<"Error in opening file..."<<endl;
			throw 1;
		}		
				
		string temp="", name="", str="";

		while (getline(is, temp, '\n'))	{
			
			name = temp;			
			
			getline(is, temp, '\n');
			cout<<"reading first seq"<<endl;
			while (temp!="") {				
				str += temp;
				temp = "";
				getline(is, temp, '\n');
			}
			cout<<"reading all seq"<<endl;
			cout<<str.length()<<endl;
			//Check str's validility
			vector<string> strr;
			strr.push_back(str.substr(0, str.length()/2));
			strr.push_back(str.substr(str.length()/2, str.length()/2));

			if (checkValid(name, strr)) {
				seq.push_back(generalInput(strr));
				seqname.push_back(name);
			}			
			name = str = "";
		}
		is.close();
		is.clear();	

	}
	catch (...) {
		cout<<"Error in reading AXT."<<endl;
		flag = 0;
	}

	return flag;
}

int Base::readFasta(string seqfile, vector<string> &seqname, vector<string> &seq) {
	int flag = 1;
	try	{
		ifstream is(seqfile.c_str());
		if (!is) throw "Error in opening file...";
				
		string temp="", name="";
		
		getline(is, temp, '\n');
		//cout<<"Debug - Each line: \n"<<temp<<endl;
		while ( (int)(temp.find('>')) > -1) {
			name = temp.substr(1, temp.length()-1);  //remove >; get the 'name' of the sequence
			string str="";
			getline(is, temp, '\n'); //get the line with the sequence 'str' for 'name'
			//cout<<"Debug - Seq line: \n"<<temp<<endl;
			while((int)(temp.find('>'))<0) {
				str += temp;
				//cout<<"Debug-Found seq: \n"<<temp<<endl;
				getline(is, temp, '\n'); //Get the sequence 
				if (is.eof()==true && temp=="") {
				cout<<"Error in readFasta. Check fasta file format."<<endl;
				cout<<"Standard fasta format should be used, '>' should be put before the taxa name and line break should be used for each sequence including the last sequence."<<endl;
				break;
				}
			} 

			//name = Trim(filterString(name));
			name = Trim(name);
			str = Trim(filterString(str));
			
			seqname.push_back(name);
			seq.push_back(str);
			//cout<<"Debug-readFasta: seqname: "<<name<<"; seq: "<<str<<endl;
		}		

		is.close();
		is.clear();
	}
	catch (const char* e) {
		cout<<e<<endl;
		flag = 0;
	}

	return flag;
}



int Base::readFastaConsensus(string seqfile, vector<string> &seqname, vector<string> &seq) {
	int flag = 1;
       
	try	{

		ifstream is(seqfile.c_str());
		if (!is) throw "Error in opening file in readFasta...";

		string temp="", name="";
		
		getline(is, temp, '\n');
		int c1=0;
		int c2=0;
		//cout<<"temp seq: "<<temp<<endl;
		while ( (int)(temp.find('>')) > -1) {
			c1=c1+1;
			//cout<<"The sequence count is "<<c1<<endl;

			name = temp.substr(1, temp.length()-1);  //remove >

			string str="";
			getline(is, temp, '\n');
			while((int)(temp.find('>'))<0) {
				c2=c2+1;
				//cout<<"The line count is "<<c2<<endl;
				str += temp;
				getline(is, temp, '\n');
				//cout<<"The line "<<c2<<": "<<temp<<endl;
				if (is.eof()==true) break;
				if (is.eof()==true && temp=="") {
				cout<<"Error in readFastaConsensus. Check fasta file format."<<endl;
				cout<<"Standard fasta format should be used, '>' should be put before the taxa name and line break should be used for each sequence including the last sequence."<<endl;
				break;}
				
				/*if (c2>10)
				{ cout<<"Error in sequence. Sequence line is greater than 1000, above the threshold setup, check for bugs."<<endl;
					break;}
					*/
			} 

			//name = Trim(filterString(name));
			name = Trim(name);
			cout<<"Sequence Name:"<<endl<<name<<endl;
			//str = Trim(filterString(str));
			cout<<"Sequence string:"<<endl<<str<<endl;
			
			seqname.push_back(name); //vector of seqname
			seq.push_back(str); //vector of seq
		}		

		is.close();
		is.clear();
	}
	catch (const char* e) {
		cout<<e<<endl;
		flag = 0;
	}

	cout<<"Finish reading for the fasta file "<<seqfile<<endl;
	return flag;
}


string Base::generalInput(vector<string> seq) {

	int i, j;
	string output = "";
	for (i=0; i<seq[0].length(); i++) {		
		int issame = 1;
		for (j=1; j<seq.size() && issame==1; j++) {
			if (seq[0][i]!=seq[j][i]) issame=0;
		}
		if (issame==1) output += "0";
		else output += "1";
	}

	return output;
}

/*

int matby (double a[], double b[], double c[], int n,int m,int k)
// a[n*m], b[m*k], c[n*k]  ......  c = a*b

{
   int i,j,i1;
   double t;
   FOR (i,n)  FOR(j,k) {
      for (i1=0,t=0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
      c[i*k+j] = t;
   }
   return (0);
}





void PMatrixTaylor(double P[], double n, double t) {

// Get approximate PMat using polynomial approximation
//   P(t) = I + Qt + (Qt)^2/2 + (Qt)^3/3!

   int nterms=1000, i,j,k, c[2],ndiff,pos=0,from[3],to[3];
   double *Q=space, *Q1=Q+n*n, *Q2=Q1+n*n, mr, div;

   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
      to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
      c[0]=GenetCode[com.icode][i];
      c[1]=GenetCode[com.icode][j];
      if (c[0]==-1 || c[1]==-1)  continue;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]*=kappa;
      if(c[0]!=c[1])  Q[i*n+j]*=omega;
      Q[j*n+i]=Q[i*n+j];
   }
   FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   for (i=0,mr=0;i<n;i++) { 
      Q[i*n+i]=-sum(Q+i*n,n); mr-=pi[i]*Q[i*n+i]; 
   }
   FOR(i,n*n) Q[i]*=t/mr;

   xtoy(Q,P,n*n);  FOR(i,n) P[i*n+i]++;   // I + Qt 
   xtoy(Q,Q1,n*n);
   for (i=2,div=2;i<nterms;i++,div*=i) {  // k is divisor 
      matby(Q, Q1, Q2, n, n, n);
      for(j=0,mr=0;j<n*n;j++) { P[j]+=Q2[j]/div; mr+=fabs(Q2[j]); }
      mr/=div;
      // if(debug) printf("Pmat term %d: norm=%.6e\n", i, mr); 
      if (mr<e) break;
      xtoy(Q2,Q1,n*n);
   }

   FOR(i,n*n) if(P[i]<0) P[i]=0;

}

*/
