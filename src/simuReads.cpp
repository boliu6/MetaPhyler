// simulate reads uniformly from genes for metaphyler training

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::istringstream;

#include <string>
using std::string;

#include <cstdlib>

typedef unsigned int Uint;

struct Cmdopts {
  string fastafile;
  Uint   length,
         stepsize;
};

void helpmsg();
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts);
void simulate(Uint length, Uint stepsize, string seqid, string &seq);


int main(int argc, char *argv[]) {

  // read in command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);


  // open fasta file
  ifstream ifs(cmdopts.fastafile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << cmdopts.fastafile << endl;
    exit(1);
  }


  // simulate reads from each fasta sequence
  string eachline, seqid, seq;
  istringstream iss;
  while (getline(ifs, eachline)) {

    if (eachline[0] == '>') {
      
      simulate(cmdopts.length, cmdopts.stepsize, seqid, seq); // process previous sequence

      iss.clear();
      iss.str(eachline);
      iss >> seqid;

      seqid.erase(0, 1); // remove '>'
      seq.clear();
    }
    else
      seq += eachline;   // store sequences
  }

  simulate(cmdopts.length, cmdopts.stepsize, seqid, seq);     // process last sequence

  return 0;
}


// simulate reads from a given sequence
void simulate(Uint length, Uint stepsize, string seqid, string &seq){

  if (seq.empty()) return;

  for (Uint i = 0; i+length-1 < seq.length(); i += stepsize) {
    cout << ">" << seqid << "_" << i/stepsize+1 << " " << i+1 << " " << i+length << endl;
    cout << seq.substr(i, length) << endl;
  }
}


// parse command line options
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 4) {
    helpmsg();
    exit(1);
  }

  cmdopts.length    = atoi(argv[1]);
  cmdopts.stepsize  = atoi(argv[2]);
  cmdopts.fastafile = argv[3];
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./simuReads <length> <step size> <FASTA file>" << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "        <length>     length of reads to be simulated." << endl << endl;
  cerr << "        <step size>  distance between two simulated reads." << endl << endl;
  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Bo Liu - boliu@umiacs.umd.edu" << endl;
  cerr << endl;
}
