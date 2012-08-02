// BLAST file: the simulated reads are mapped to reference genes
// This program reports the scores in a histogram style

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::istringstream;

#include <map>
using std::map;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <algorithm>
#include <functional>
using std::greater;

#include <cstdlib>
#include <ctime>
#include <cmath>

#include <utility>
using std::pair;


typedef unsigned short int  Usint;
typedef unsigned int        Uint;
typedef vector<string>      VS;
typedef map<string, VS>     S2VS;
typedef vector<Usint>       VSI;
typedef vector<VSI>         VVSI;
typedef map<string, VVSI>   S2VVSI;
typedef map<string, Uint>   S2I;
typedef map<string, float>  S2F;

const unsigned int BITCUTOFF  = 1;


// stores command line options
struct Cmdopts{
  string taxfile,
         blastfile,
         blast,
         refseq;
  Usint  readlen;
};

void helpmsg();
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts);
void readTaxFile(string taxfile, S2VS &seq2tax, S2I &tid2num);
void train(string blastfile, S2VS &seq2tax, S2VVSI &seq2scores);
void printScores(S2VVSI &seq2scores, Cmdopts &cmdopts, const S2I &ref2len);
void printDistribution(S2VVSI &seq2scores, Cmdopts &cmdopts, const S2I &ref2len);
Usint findlca(S2VS::iterator riter, S2VS::iterator qiter);
void readRefseq(string refseqfile, S2I &ref2len, string blast);


int main(int argc, char *argv[]) {

  Cmdopts cmdopts;         // read in command line options
  getcmdopts(argc, argv, cmdopts);


  S2I ref2len;             // lengths of reference sequences
  readRefseq(cmdopts.refseq, ref2len, cmdopts.blast);


  S2VS seq2tax;            // taxonomic profile of each reference sequence
  S2I  tid2num;            // number of genes under each lowest taxonomic cluster
  readTaxFile(cmdopts.taxfile, seq2tax, tid2num);


  S2VVSI seq2scores;       // bit scores under each taxonomic level for each sequence
  train(cmdopts.blastfile, seq2tax, seq2scores);


  //printScores(seq2scores, cmdopts, ref2len); // summarize scores and print them out
  printDistribution(seq2scores, cmdopts, ref2len); // summarize scores and print them out

  return 0;
}


void readRefseq(string refseqfile, S2I &ref2len, string blast) {

  ifstream ifs(refseqfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << refseqfile << endl;
    exit(1);
  }


  // simulate reads from each fasta sequence
  string eachline, seqid("");
  Uint len = 0;
  istringstream iss;
  while (getline(ifs, eachline)) {

    if (eachline[0] == '>') {
      
      if (blast == "blastx")
	len *= 3;
      ref2len.insert(S2I::value_type(seqid, len));

      iss.clear();
      iss.str(eachline);
      iss >> seqid;

      seqid.erase(0, 1); // remove '>'
      len = 0;
    }
    else
      len += eachline.size();   // store sequences
  }

  if (blast == "blastx")
    len *= 3;
  ref2len.insert(S2I::value_type(seqid, len));
}


// load in taxonomy file
void readTaxFile(string taxfile, S2VS &seq2tax, S2I &tid2num) {

  ifstream ifs(taxfile.c_str());
  if (!ifs)
    cerr << "Could not open file " << taxfile << endl;

  string eachline, eachword, seqid;
  VS tlabs;
  istringstream iss;
  while (getline(ifs, eachline)) {

    tlabs.clear();
    iss.clear();
    iss.str(eachline);
    iss >> seqid >> eachword;
    
    // consider the lowest taxonomic level as clusters
    if (eachword != "NA") {
      if (tid2num.find(eachword) == tid2num.end())
	tid2num.insert(S2I::value_type(eachword, 1));
      
      else // store the abundance of each cluster
	++(tid2num.find(eachword)->second);
    }

    // store taxonomic labels of each gene
    tlabs.push_back(eachword);
    while(iss >> eachword)
      tlabs.push_back(eachword);
    seq2tax.insert( pair<string, VS> (seqid, tlabs) );
  }
}


// parse command line options
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 6) {
    helpmsg();
    exit(1);
  }


  cmdopts.taxfile   = argv[1];
  cmdopts.refseq    = argv[2];
  cmdopts.blastfile = argv[3];
  cmdopts.readlen   = atoi(argv[4]);
  cmdopts.blast     = argv[5];
}


// print sorted bit scores at each taxonomy level for each sequence
void printScores(S2VVSI &seq2scores, Cmdopts &cmdopts, const S2I &ref2len) {

  cout << "#Length " << cmdopts.readlen << endl;
  cout << "#BLAST " << cmdopts.blast << endl;
  
  // for each gene
  for (S2VVSI::iterator iter1 = seq2scores.begin(); iter1 != seq2scores.end(); ++iter1) {

    if (ref2len.find(iter1->first)->second < cmdopts.readlen) continue;
    
    cout << ">" << iter1->first << "\t" << iter1->second.size() - 1 << endl;

    // for each taxonomic level
    for (VVSI::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {

      if (iter2->empty()) {
	cout << endl;
	continue;
      }
      
      // sort the bit scores
      sort(iter2->begin(), iter2->end(), greater<int>());
      
      // print out each score, the number of scores that are >= this score
      Uint pre = *iter2->begin(), num = 1;
      for (VSI::iterator iter3 = iter2->begin()+1; iter3 != iter2->end(); ++iter3) {

	// this is a new socre
	if (*iter3 != pre) { 
	  cout << pre << " " << num << " "; // print out last record
	  pre = *iter3;
			num = 1;
	}
	num++;
      }
      cout << pre << " " << num << " " << endl; // print out last record
    }
  }
}

void printDistribution(S2VVSI &seq2scores, Cmdopts &cmdopts, const S2I &ref2len) {

  cout << "#Length " << cmdopts.readlen << endl;
  cout << "#BLAST " << cmdopts.blast << endl;
  
  // for each gene
  for (S2VVSI::iterator iter1 = seq2scores.begin(); iter1 != seq2scores.end(); ++iter1) {

    if (ref2len.find(iter1->first)->second < cmdopts.readlen) continue;
    
    cout << ">" << iter1->first << "\t" << iter1->second.size() - 1 << endl;

    // for each taxonomic level
    for (VVSI::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {

      if (iter2->empty()) {
	cout << endl;
	continue;
      }
      
      // print out each score, the number of scores that are >= this score
	Uint sum = 0, sqsum = 0;
      for (VSI::iterator iter3 = iter2->begin(); iter3 != iter2->end(); ++iter3) {
		sum += *iter3;
		sqsum += (*iter3)*(*iter3);
      }
	Uint num = iter2->size();
	cout << sum / num << "\t" << sqrt((sqsum/num) - (sum/num)*(sum/num)) << endl;
      //cout << pre << " " << num << " " << endl; // print out last record
    }
  }
}


// process blast bit scores, compare the tax labels between query and reference
// store them in corresponding tax level
void train(string blastfile, S2VS &seq2tax, S2VVSI &seq2scores) {
  
  ifstream blastfile_ifs(blastfile.c_str());
  if (!blastfile_ifs)
    cerr << "Could not open file: " << blastfile << endl;

  string eachline;
  while (getline(blastfile_ifs, eachline)) {

    size_t pos1 = eachline.find("\t");
    string qid  = eachline.substr(0, pos1);
    qid.resize(qid.rfind('_'));    // trim the end of the simulated ID
    
    size_t pos2 = eachline.find("\t", pos1+1);
    string rid  = eachline.substr(pos1+1, pos2-pos1-1);
    
    size_t pos3 = eachline.rfind(" ") != string::npos ? eachline.rfind(" ") : eachline.rfind("\t"); // sometimes an extract space before bit score
    string bit_str = eachline.substr(pos3+1, eachline.size()-pos3-1);
    Usint bit   = atoi(bit_str.c_str());

    if (bit < BITCUTOFF) continue; // ignore bad blast hit

    // taxonomic labels should be available for both sequences
    S2VS::iterator qiter = seq2tax.find(qid);
    S2VS::iterator riter = seq2tax.find(rid);
    if (qiter == seq2tax.end() || riter == seq2tax.end()) continue;

    // suppose sequence A has 5 tax labels, then we need 6 vectors to store bit scores
    // for each level plus an "other" level
    if (seq2scores.find(rid) == seq2scores.end())
      seq2scores.insert(S2VVSI::value_type(rid, VVSI(riter->second.size()+1, VSI())));

    Usint lca = findlca(riter, qiter);
    
    // their tax labels match at level lca of reference sequence
    S2VVSI::iterator scoreiter = seq2scores.find(rid);
    scoreiter->second[lca].push_back(bit);

  }
  
}


// Find lowest common ancester between query and reference w.r.t reference
// For example, reference: A B C D, query: E C F G,
// then this function returns 2, because reference[2] = C = query[1];
Usint findlca(S2VS::iterator riter, S2VS::iterator qiter) {

  Usint i = 0;
  bool tag = false;
  for (VS::const_iterator rtiter = riter->second.begin(); rtiter != riter->second.end(); ++rtiter, ++i) {

    if (*rtiter == "NA") continue;

    for (VS::const_iterator qtiter = qiter->second.begin(); qtiter != qiter->second.end(); ++qtiter) {
      if (*rtiter == *qtiter) {
	tag = true;
	break;
      }
    }

    if (tag) break; 
  }

  return i;
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./metaphylerTrain <taxonomy file> <ref seq> <BLAST file> <length> <BLAST program>" << endl;
  cerr << endl;

  cerr << "Options:" << endl;

  cerr << "        <taxonomy file> Taxonomy labels of reference sequences in the BLAST file." << endl << endl;

  cerr << "        <ref seq>       Reference sequence FASTA file." << endl << endl;

  cerr << "        <BLAST file>    BLAST alignment between simulated reads and reference sequences." << endl;
  cerr << "                        All simulated reads come from reference sequences, and are named as follows:" << endl;
  cerr << "                        If n reads come from A, then their IDs are A_0, A_1,..., A_n-1." << endl;

  cerr << "        <length>        Length of simulated reads." << endl << endl;

  cerr << "        <BLAST program> BLASTN, BLASTP, BLASTX or TBLASTX." << endl << endl;

  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Bo Liu - boliu@umiacs.umd.edu" << endl;
  cerr << endl;

}

/*
  todo list:
  1, automatically find bitcutoff
 */
