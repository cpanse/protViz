// Christian Panse <cp@fgcz.ethz.ch>
// 20171111

#include <Rcpp.h>
#include <numeric>
#include <string>
#include <fstream>
#include <iostream>

using namespace Rcpp;

/*
void read_FASTA(const String& filename, StringVector& fasta){
  std::string line;
  std::ifstream myfile (filename);
  
  if (myfile.is_open())
  {
    while ( std::getline (myfile, line) )
    {
      fasta.push_back(line);
    }
    
    myfile.close();
  }
}
*/
// [[Rcpp::export]]
StringVector fcat_FASTA(const StringVector& fasta){
  //std::vector<std::string> rv;
  StringVector rv;
  std::string aa_seq = "";
  
  for (auto s : fasta) {
    
    // if (s.at(0) == '>') {
    
    if (s[0] == '>') {
      if (aa_seq != ""){
        rv.push_back(aa_seq);
        aa_seq = "";
      }
    }else {
      aa_seq = aa_seq.append(s);
    }
  }
  rv.push_back(aa_seq);
  
  return rv;
}


// [[Rcpp::export]]
StringVector tryptic_digest_FASTA(const StringVector& fasta) {
  
  StringVector rv;
  
  char aa0 = '\0';
  std::string digest = "";
  
  for (auto sequence : fasta) {
    for (auto aa1 : sequence) {
      
      if (aa0 != '\0') {
        digest += aa0;
      }
      
      if ((aa1 != 'P' && aa0 == 'R') || aa0 == 'K') {
        if (digest != ""){
          rv.push_back(digest);
          digest = "";
        }
      }
      aa0 = aa1;
    }
  }
  
  return rv;
}


// [[Rcpp::export]]
int number_of_tryptic_peptides_FASTA(const StringVector& fasta0) {
  int n = 0;
  StringVector fasta = fcat_FASTA(fasta0);
  char aa0 = '\0';
  std::string digest = "";
  
  for (auto sequence : fasta) {
    for (auto aa1 : sequence) {
      
      if (aa0 != '\0') {
        digest += aa0;
      }
      
      if ((aa1 != 'P' && aa0 == 'R') || aa0 == 'K') {
        if (digest != ""){
          n++;
          digest = "";
        }
      }
      aa0 = aa1;
    }
  }
  
  return n;
}



class Fasta {
public:
  Fasta(std::string filename): filename_(filename) {
    read();
  }
  std::vector<std::string> getTrypticPeptides(){
    computeTrypticPeptides();
    return(Tryptic_);
  }
  
  int getNumberOfTrypticPeptides(){
    computeTrypticPeptides();
    return(Tryptic_.size());
  }
  
  int getDesc() {
    return Desc_.size(); }
  int getSeq() {
    return Seq_.size(); }
  void addValue(int y) { x_ += y; }
  void merge(const Fasta& rhs) { x_ += rhs.x_; }
  
private:
  int x_;
  std::string filename_;
  std::vector<std::string> Desc_;
  std::vector<std::string> Seq_;
  std::vector<std::string> Tryptic_;
  
  void computeTrypticPeptides() {
    if (Tryptic_.size() > 0) return;
    // int n = 0;
    char aa0 = '\0';
    std::string digest = "";
    
    for (auto sequence : Seq_) {
      for (auto aa1 : sequence) {
        
        if (aa0 != '\0') {
          digest += aa0;
        }
        
        if ((aa1 != 'P' && aa0 == 'R') || aa0 == 'K') {
          if (digest != ""){
            Tryptic_.push_back(digest);
            digest = "";
          }
        }
        aa0 = aa1;
      }
    }
  }
  
  
  void read(){
    //filename_ = "/Users/cp/p1875_db10_20170817.fasta";
    
    std::string line;
    std::ifstream myfile (filename_);
    std::string s = "";
    
    if (myfile.is_open())
    {
      while ( std::getline (myfile, line) )
      {
        if (line[0] == '>') {
          Desc_.push_back(line);
          
          if (s.length() > 0) {
            Seq_.push_back(s);
            s = "";
          }
        }
        else {
          s = s.append(line);
        }
      }
      
      myfile.close();
      if (s.length() > 0)
        Seq_.push_back(s);
    }
  }
};


// RCPP_EXPOSED_CLASS(Fasta)
  RCPP_MODULE(mod_fasta) {
	using namespace Rcpp;
    
    class_<Fasta>("Fasta")
    .default_constructor("Default constructor") // This exposes the default constructor
    .constructor<std::string>("sets initial value")
    .method("getDesc", &Fasta::getDesc, "Returns the value")
    .method("getSeq", &Fasta::getSeq, "Returns the value")
    .method("getNumberOfTrypticPeptides", &Fasta::getNumberOfTrypticPeptides, "Returns the value")
    .method("getTrypticPeptides", &Fasta::getTrypticPeptides, "Returns the value")
    .method("addValue", &Fasta::addValue, "Adds a value")
    .method("merge", &Fasta::merge, "Merges another Test into this object")
    ;
  }
