#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::ios_base;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::istringstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <cstdlib>

struct pred {
  string id;
  float conf;
};

typedef unsigned int Uint;
typedef map<string, string> S2S;
typedef vector<string>      VS;
typedef map<string, Uint>   S2I;
typedef vector<S2I>         VS2I;
typedef map<string, Uint>   S2I;
typedef vector<pred>        VP;
typedef map<string, VP>     S2VP;

struct Cmdopts {
  VS files;
};

const Uint TLEV = 6;

void helpmsg();
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts);
void merge(string file, S2VP &confs);
void printconfs(const S2VP &confs);

int main(int argc, char *argv[]) {

  // read in command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);

  S2VP confs;
  for (Uint i = 0; i < cmdopts.files.size(); i++) {
    merge(cmdopts.files[i], confs);
  }

  printconfs(confs);
  return 0;
}

void printconfs(const S2VP &confs) {
  cout.setf(ios_base::fixed);
  cout.precision(3);
  for (S2VP::const_iterator citer = confs.begin(); citer != confs.end(); ++citer) {
    cout << citer->first << "\t";
    for (VP::const_iterator citer2 = citer->second.begin(); citer2 != citer->second.end(); ++citer2) {
      cout << citer2->id << "(" << citer2->conf << ")\t";
    }
    cout << endl;
  }
}

void merge(string file, S2VP &confs) {

  ifstream ifs(file.c_str());
  if (!ifs) {
    cerr << "Could not open file " << file << endl;
    exit(1);
  }

  string eachline, eachword, tid, seqid;
  istringstream iss;
  while (getline(ifs, eachline)) {
    
    iss.clear();
    iss.str(eachline);
    iss >> seqid;
    S2VP::iterator iter = confs.find(seqid);
    pred p;
    if (iter == confs.end()) {
      VP &vp = confs.insert(S2VP::value_type(seqid, VP())).first->second;
      while (iss >> eachword) {
	if (eachword == "NA") {
	  pred p = {"NA", 0.0};
	  vp.push_back(p);
	}
	else {
	  size_t pos = eachword.find('(');
	  pred p = {eachword.substr(0, pos), atof(eachword.substr(pos+1, 5).c_str())};
	  vp.push_back(p);
	}
      }
    }
    else {
      VP &vp = iter->second;
      Uint i = 0;
      while (iss >> eachword) {
	i++;
	if (eachword == "NA") continue;
	size_t pos = eachword.find('(');
	float conf = atof(eachword.substr(pos+1, 4).c_str());
	if (vp[i-1].conf < conf)
	  vp[i-1].conf = conf;
      }
    }
  }

}
// parse command line options
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc == 1) {
    helpmsg();
    exit(1);
  }

  for (int i = 1; i < argc; ++i) {
    cmdopts.files.push_back(argv[i]);
  }
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./combine <classification 1> <classification 2> ..." << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "        <classification> Result file from program metaphylerClassify." << endl;
  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Bo Liu - boliu@umiacs.umd.edu" << endl;
  cerr << endl;
}
