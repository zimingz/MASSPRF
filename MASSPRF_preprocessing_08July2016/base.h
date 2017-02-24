/******************************************************************
* Copyright (c) 2007
* All rights reserved.
 
* Filename: base.h
* Abstract: Declaration of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang  (zhang.zhang@yale.edu)
* Date: Feb.2, 2005

******************************************************************/
#pragma warning(disable:4786)		//Disable warning when using vector

#if !defined(BASE_H)
#define  BASE_H

#define CODONSIZE 3			//Length of codon
#define DNASIZE 4				//A C G T
#define XSIZE DNASIZE*DNASIZE   //Size of one group AXX (X=A,C,G,T) 
#define CODONCOUNT 64				//Codon Count
#define NULL 0					//Zero
#define NA -1					//Not Available
#define NCODE	23				//Number of genetic codes
#define NNCODE  NCODE*2			//Double of the number genetic codes		
#define SMALLVALUE 1e-10			//Value near to zero
#define NUMBER_OF_RATES	6		//Number of substitution rates
#define MODELCOUNT	14			//Number of candidate models

#define gammap(x,alpha) (alpha*(1-pow(x,-1/alpha)))
#define square(a) ((a)*(a))

#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Stanard lib of C++ */
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<algorithm>
#include <iomanip>
#include <cstdlib>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
#include<cctype>
//#include<WCHAR.H>
//#include<WCTYPE.H>
#include<valarray>


using namespace std;

/*****Global variables*****/
extern const char* transl_table[];	//Genetic codon table
extern int genetic_code;			//ID of codon table from 1 to 23

//End of Global variables


/* Convert one type to any other type */
template<class out_type,class in_value>
	out_type CONVERT(const in_value & t) {
		stringstream stream;
		//Put the value 't' into the stream
		stream<<t;			
		out_type result;
		//Put the stream into the 'result'
		stream>>result;

		return result;
	}

class Base {

public:
	Base();

	/* Write the content into file */
	bool writeFile(string output_filename, const char* result);

	/* Parse results */
	string parseOutput();	
	/* Format string for outputing into file */

	void addString(string &result, string str, string flag="\t");

	/* Generate a radnom integer */
	int getRandom();
	double getRandDouble();
	int getRandInt(int begin, int end);
	int setRandSeed();

	int countChar(string str, char ch);
	
	/* Convert a char-T,C,A,G into a digit 0,1,2,3, respectively */
	int  convertChar(char ch);
	/* Convert a digit-0,1,2,3 into a char T,C,A,G, respectively */
	char convertInt(int ch);
	/* Convert a string to uppercase */
	string stringtoUpper(string str);
	string Trim(string str);
	string filterString(string str);
	string replaceAll(string str, char from, char to);
	string replaceAll(string str, char from);
	
	/* Calculate the amino acid of codon by codon */
	char getAminoAcid(string codon);
	/* Calculate the amino acid of codon by codon's id*/
	char getAminoAcid(int id);
	/* Get the number of stop codon in a given genetic code table */
	int getNumNonsense(int genetic_code);

	/* Return the codon's id from codon table */
	int getID(string codon);
	/* Return a codon according to the id */
	string getCodon(int IDcodon);

	/* Sum array's elements */
	double sumArray(double x[], int end, int begin=0);
	int sumArray(int x[], int end, int begin=0);

	/* Init value to array */
	int initArray(double x[], int n, double value=0.0);
	int initArray(int x[], int n, int value=0);

	/* Elements in array are mutipled by scale */
	int scaleArray(double scale, double x[], int n);
	/* Sqrt of the sum of the elements' square */
	double norm(double x[], int n);
	/* Copy array's values one by one: to[] = from[] */
	int copyArray(double from[], double to[], int n);
	/* Sum of 'n' products multiplied by two elements x[], y[] */
	double innerp(double x[], double y[], int n);
	/* Set x[i,j]=0 when x!=j and x[i,j]=1 when x=j */
	int initIdentityMatrix(double x[], int n);

	bool removeCntrl(string name, vector<string> &str);
	bool checkValid(string name, vector<string> &str);
	int readAXT(string seqfile, vector<string> &seqname, vector<string> &seq);
	int readFasta(string seqfile, vector<string> &seqname, vector<string> &seq);
	int readFastaConsensus(string seqfile, vector<string> &seqname, vector<string> &seq);
	string generalInput(vector<string> seq);
	string RGB(int r, int g, int b);

protected:
	/* Name of method for calculating ka/ks */
	//string name;
	/* Sysnonymous sites(S) and nonsysnonymous sites(N) */
	//double S, N;
	/* Number of sysnonymous differences(Sd) and nonsysnonymous differences(Nd), snp=Sd+Nd */
	//double Sd, Nd, snp;
	/* Number of sysnonymous substitutions(ks) and nonsysnonymous substitutions(ka) per site */
	//double Ka, Ks;
	/* Standard Errors for Ka and Ks */
	//double SEKa, SEKs;
	/* Transitional(Si) and transversional(Vi) substitutions */
	//double Si[5], Vi[5];
	/* 0-fold, 2-fold, 4-fold */
	//double L[5];
	/* Total Numbers of substitutions per i-th type site: K[i]=A[i]+B[i] */
private:
	unsigned long Seed;
	
};

#endif


