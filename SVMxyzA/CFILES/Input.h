#ifndef INPUT_H
#define INPUT_H
#include <vector>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>


const int ndpterms = 10;
const double pi = 4.0*atan(1.0);

//=============================================================================
struct potentialop {
	double vpot[ndpterms]={}, aquad[ndpterms]={}, alin[ndpterms]={};
	int    npow[ndpterms]={} , nterms = 0;
	std::string optype;
};
//=============================================================================
struct IsospinVector {
	IsospinVector(size_t size) :tz(size) {}
	IsospinVector() {}
	std::vector <int> tz;
};
//=============================================================================
struct SpinVector {
	SpinVector(size_t size) :sz(size) {}
	SpinVector() {}
	std::vector <int> sz;
};
//=============================================================================
struct IsospinState {
	int ncmp;
	std::vector<double> coef;
	std::vector<IsospinVector> cmp;
};
//=============================================================================
struct SpinState {
	int ncmp;
	std::vector<double> coef;
	std::vector<SpinVector> cmp;
};
//=============================================================================

class Input
{
private:
	int    rdai(std::istringstream& iss);
	double rdaf(std::istringstream& iss);
	void   rdpot(std::string& word, std::istringstream& iss);
	void   rdts(std::string& word, std::istringstream& iss);

public:
	Input(std::string const &jobname);
	~Input();
	void print();
	int  npar, irand, keycontinue, bosefermi;
	int  maxbasis, mm0, kk0;
	int dmax;
        double eB, momega, eBspin;
	double h2m, rndmin, rndmax;
	std::vector<double> mass;
	std::vector<double> charges;
	int nop, npt;
	std::vector<potentialop> potop;
	double vpot3b, apot3b;
	SpinState spin;
	IsospinState isospin;
	FILE *printfile;
	std::string jobname;
};

#endif 
