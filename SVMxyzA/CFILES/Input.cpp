#include "Input.h"
#include "BasisState.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
using namespace std;
  
Input::Input(string const &jobname_) : jobname(jobname_)
{
    size_t notfound = string::npos;
    //ofstream myfile;
    //myfile.open("example.txt");
    //myfile << "Writing this to a file.\n";
    
    ifstream inputfile;
    std::string buffer, line, word;
    
    /*                  open inputfile  */
    cout << "\t open input file " << jobname << endl;
    inputfile.open(jobname);
    
    
    /*    read inputfile and remove all comment lines  */
    buffer.clear();
    while (getline(inputfile, line))
	{
	    int ignore = line.find('!');
	    if (ignore != notfound)
		line = line.substr(0, ignore);
	    if (line.length() == 0)
		line = " ";
	    buffer.append(line);
	}
    //myfile << buffer << endl;
    //myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;

    /*                  remove all '=' signs from buffer         */
    std::replace(buffer.begin(), buffer.end(), '=', ' ');
    //myfile << buffer << endl;
    //myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    
    /*                  initialize some variables */
    dmax = 0;
    vpot3b = 0. ; apot3b = 0.;
    nop1b = 4;
    nop2b = 4;

    /*                  read data from buffer */
    istringstream iss(buffer);
    while (iss >> word){
	cout << word << endl;
	//myfile << "Input read: " << word << endl;
	if (word == "npar") {
	    npar = rdai(iss);
	    mass.resize(npar);
	    charges.resize(npar);
	}
	else if (word == "h2m") h2m = rdaf(iss);
	else if (word == "hh" ) h2m = rdaf(iss);
	else if (word == "eB" ) eB = rdaf(iss);  
	else if (word == "momega" ) momega = rdaf(iss); 
	else if (word == "eBspin") eBspin = rdaf(iss);
	else if (word == "dmax") dmax = rdai(iss);
	else if (word == "bmin") rndmin = rdaf(iss);
	else if (word == "bmax") rndmax = rdaf(iss);
	else if (word == "irand") irand = rdai(iss);
	else if (word == "ibf") bosefermi = rdai(iss);
	else if (word == "ico") keycontinue = rdai(iss);
	else if (word == "mnb") maxbasis = rdai(iss);
	else if (word == "mm0") mm0 = rdai(iss);
	else if (word == "kk0") kk0 = rdai(iss);
	else if (word == "vpot3b") vpot3b = rdaf(iss);
	else if (word == "apot3b") apot3b = rdaf(iss);
	else if (word == "nop") { nop = rdai(iss); potop.resize(nop);}
	else if (word == "npt") npt = rdai(iss);
	else if (word.find("xm") != notfound) {
	    for (int i = 0; i < npar; i++) mass[i] = rdaf(iss);
	}
	else if (word.find("z") != notfound) {
	    for (int i = 0; i < npar; i++) charges[i] = rdaf(iss);
	}
	else if (word.find("vpot") != notfound) rdpot(word, iss);
	else if (word.find("apot") != notfound) rdpot(word, iss);
	else if (word.find("bpot") != notfound) rdpot(word, iss);
	else if (word.find("npot") != notfound) rdpot(word, iss);
	else if (word == "nisc") {
	    int ncmp = rdai(iss);
	    isospin.ncmp = ncmp; isospin.coef.resize(ncmp); isospin.cmp.resize(ncmp, npar);
	}
	else if (word == "nspc") {
	    int ncmp = rdai(iss);
	    spin.ncmp = ncmp; spin.coef.resize(ncmp); spin.cmp.resize(ncmp, npar);
	}
	else if (word.find("cisc") != notfound) rdts(word, iss);
	else if (word.find("cspc") != notfound) rdts(word, iss);
	else if (word.find("iso(") != notfound) rdts(word, iss);
	else if (word.find("isp(") != notfound) rdts(word, iss);
	else if (word.find("nts_states") != notfound) nts_states = rdai(iss);
	else if (word.find("ts_state") != notfound) rdts_state(iss);
    }
    inputfile.close();
    if (nts_states != ts_states.size()) {
	cout << "\n\tInupt ERROR - number of ts_statees neq size of ts_states" 
	     << endl << endl; exit(0);
    }
    //myfile.close();
    //  printfile = fopen("out.txt", "w");
}

////=============================================================================
double Input::rdaf(istringstream& iss)
{
	double f;
	string c;

	if (!(iss >> f)) { iss.clear(); iss >> c; iss >> f; }
	return f;
}
////=============================================================================
int Input::rdai(istringstream& iss)
{
	int i;
	string c;

	if (!(iss >> i)) { iss.clear(); iss >> c; iss >> i; }
	return i;
}
////=============================================================================
void Input::rdpot(string& word, istringstream& iss)
{
	int iop, iterm, i1, i2, i3, key3b;
	string comp;
	double f;

	key3b = word.find("3b");
	if (key3b == string::npos) key3b = -1;

	if (key3b == -1) {
		i1 = word.find("("); i2 = word.find(","); i3 = word.find(")");
		iterm = stoi(word.substr(i1 + 1, i2 - i1 - 1));
		iop = stoi(word.substr(i2 + 1, i3 - i2 - 1));
		comp = word.substr(0, 4);
		if (comp == "vpot") {
			iss >> potop[iop - 1].vpot[iterm - 1];
			potop[iop - 1].nterms++;
		}
		if (comp == "apot") iss >> potop[iop - 1].aquad[iterm - 1];
		if (comp == "bpot") iss >> potop[iop - 1].alin[iterm - 1];
		if (comp == "npot") iss >> potop[iop - 1].npow[iterm - 1];
	}
}
////=============================================================================
void Input::rdts(string& word, istringstream& iss)
{
	int ipar, icmp, i1, i2, i3, keyts;
	size_t notfound = string::npos;
	string caux;
	double f;

	if (word.find("iso") != notfound) {
		i1 = word.find("("); i2 = word.find(","); i3 = word.find(")");
		ipar = stoi(word.substr(i1 + 1, i2 - i1 - 1));
		icmp = stoi(word.substr(i2 + 1, i3 - i1 - 1));
		iss >> isospin.cmp[icmp - 1].tz[ipar - 1];
	}
	if (word.find("isp") != notfound) {
		i1 = word.find("("); i2 = word.find(","); i3 = word.find(")");
		ipar = stoi(word.substr(i1 + 1, i2 - i1 - 1));
		icmp = stoi(word.substr(i2 + 1, i3 - i1 - 1));
		iss >> spin.cmp[icmp - 1].sz[ipar - 1];
	}
	if (word.find("cisc") != notfound) {
		i1 = word.find("("); i3 = word.find(")");
		icmp = stoi(word.substr(i1 + 1, i3 - i1 - 1));
		iss >> isospin.coef[icmp - 1];
	}
	if (word.find("cspc") != notfound) {
		i1 = word.find("("); i3 = word.find(")");
		icmp = stoi(word.substr(i1 + 1, i3 - i1 - 1));
		iss >> spin.coef[icmp - 1];
	}
}
//=============================================================================
void Input::rdts_state(istringstream& iss){
    int ncmp;
    size_t notfound = string::npos;
    string word,npud;
    MatrixXi sz,tz;
    VectorXd coef;
    
    iss >> word;
    if (word.find("ncmp") != notfound) {
	ncmp=rdai(iss);	
	coef.resize(ncmp); sz.resize(npar,ncmp); tz.resize(npar,ncmp);
    }
    else {cout << "Input error expecting ncmp (rdts_state)"; exit(0);}
    
    for (int icmp = 0; icmp < ncmp; icmp++) {
	iss >> word;
	if (word != "coef") {
	    cout << "\tInput ERROR expecting 'coef' (rdts_state)\n" ; exit(0);}
	coef(icmp)=rdaf(iss);
	iss >> word;
	if (word != "config") {
	    cout << "\tInput ERROR expecting 'config' (rdts_state)\n" ; exit(0);}
	for (int ipar = 0; ipar < npar; ipar++) {
	    iss >> npud;
	    cout << npud << " " << npud[1] << endl;
	    tz(ipar,icmp)=-9999; sz(ipar,icmp)=-9999;
	    if (npud.substr(0,1) == "p") tz(ipar,icmp)=1;
	    if (npud.substr(0,1) == "n") tz(ipar,icmp)=2;
	    if (npud.substr(1,2) == "1") sz(ipar,icmp)=1;
	    if (npud.substr(1,2) == "2") sz(ipar,icmp)=2;
	}
    }
    cout << ncmp << npar << sz.size() << coef.size() << endl;
    SpinIsospinState ts_state(ncmp,coef,sz,tz);
    ts_state.print();
    ts_states.push_back(ts_state);
}
//=============================================================================
void Input::print(void)
{
/*
	// print into output file
	fprintf(printfile, "Input Data: \n");
	fprintf(printfile, "\t input file = %s   \n", jobname.c_str());
	fprintf(printfile, "\t npar = %8d   \n", npar);
	fprintf(printfile, "\t h2m  = %8.3f \n", h2m);
	fprintf(printfile, "\t rndmin = %8.3f  rndmax = %8.3f \n", rndmin, rndmax);
	fprintf(printfile, "\t irand= %8d    bosefermi = %3d \n", irand, bosefermi);
	fprintf(printfile, "\t niters_out = %8d  niters_in = %8d \n", kk0, mm0);
	fprintf(printfile, "\t maxbasis = %8d\n", maxbasis);
	fprintf(printfile, "\t pi   = %16.12f \n", pi);
	for (int i = 0; i<npar; i++)
        fprintf(printfile, "\t mass[%d] = %8.3f", i, mass[i]);
	fprintf(printfile, "\n");

	// print isospin state
	fprintf(printfile, "\n\t spin-isospin states\n");
	fprintf(printfile, "\t isospin state: ncmp = %d  \n", isospin.ncmp);
	for (int icmp = 0; icmp < isospin.ncmp; icmp++)
	{
		fprintf(printfile, "\t\t cmp = %d   coef = %8.3f ", icmp, isospin.coef[icmp]);
		for (int ipar = 0; ipar<npar; ipar++)
			fprintf(printfile, "  tz[%d]=%d", ipar, isospin.cmp[icmp].tz[ipar]);
		fprintf(printfile, "\n");
	}
	// print spin state
	fprintf(printfile, "\t spin state: ncmp = %d  \n", spin.ncmp);
	for (int icmp = 0; icmp < spin.ncmp; icmp++)
	{
		fprintf(printfile, "\t\t cmp = %d   coef = %8.3f ", icmp, spin.coef[icmp]);
		for (int ipar = 0; ipar<npar; ipar++)
			fprintf(printfile, "  sz[%d]=%d", ipar, spin.cmp[icmp].sz[ipar]);
		fprintf(printfile, "\n");
	}
	// print potential params
	fprintf(printfile, "\n\t potential parameters\n");
	for (int iop = 0; iop<nop; iop++)
	{
	  potentialop &pop = potop[iop];
	  if (pop.nterms == 0) continue;
	  fprintf(printfile, "\t iop =%4d  nterms =%4d \n", iop, pop.nterms);
	  for (int jt = 0; jt < pop.nterms; jt++)
	    fprintf(printfile, "\t\t iterm = %4d   vpot = %8.3f   "
		    "aquad = %8.3f   alin = %8.3f   npow = %4d \n"
		    , jt, pop.vpot[jt], pop.aquad[jt], pop.alin[jt], pop.npow[jt]);
	}
	fprintf(printfile, "\t 3-body   vpot3b = %8.3f   "
		    "aquad3b = %8.3f   \n"
		    ,  vpot3b, apot3b);
*/

}
//=============================================================================

Input::~Input()
{
}
