// Christian Panse <cp@fgcz.ethz.ch>
// 20171114 CNL470 ZH -> HAL ; 20171115 ICE77 F -> ZH


#include <numeric>
#include <string>
#include <fstream>
#include <iostream>

#include <Rcpp.h>

using namespace Rcpp;

class Fasta {
public:
  Fasta(){}
  Fasta(const String filename): filename_(filename) {
    read();
  }
  
  StringVector getTrypticPeptides(){
    computeTrypticPeptides();
    return(Tryptic_);
  }
  StringVector getSequences() {
    return(Seq_);
  }
  int getNumberOfTrypticPeptides(){
    
    if (Tryptic_.size() > 0) {
      return (Tryptic_.size());}
   
    int n = 0;
    char aa0 = '\0';
    std::string digest = "";
    
    for (auto sequence : Seq_) {
      for (auto aa1 : sequence) {
        
        if (aa0 != '\0') {
          digest += aa0;
        }
        
        if ((aa1 != 'P' && aa0 == 'R') || aa0 == 'K') {
          if (digest != ""){
            //Tryptic_.push_back(digest);
            n++;
            digest = "";
          }
        }
        aa0 = aa1;
      }
    }
    return(n);
  }
  
  int getNumberOfAminoAcids() {
    int sumAA = 0;
    
    for (auto a : Seq_)
      sumAA += a.size();
      
      return sumAA;
  }
  int getNumberOfDescriptions() {
    return Desc_.size(); }
  
  int getNumberOfSequences() {
    return Seq_.size(); }
  
  void addValue(int y) { x_ += y; }
  
  void merge(const Fasta& rhs) { x_ += rhs.x_; }
  
  List summary(){
    return List::create(Named("filename") = filename_,
                        Named("number of amino acids") = getNumberOfAminoAcids(),
                        Named("number of proteins") = getNumberOfDescriptions(),
                        Named("number of tryptic peptides") = getNumberOfTrypticPeptides()
                        );
  }
  
private:
  int x_;
  std::string filename_;
  StringVector Desc_;
  StringVector Seq_;
  //std::vector<std::string> Tryptic_;
  StringVector Tryptic_;
  
  void computeTrypticPeptides() {
    if (Tryptic_.size() > 0) return;
    // int n = 0;
    char aa0 = '\0';
    std::string digest;
    
    for (auto sequence : Seq_) {
      digest = "";
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
      // this is the end of the seq
      if (aa0 != '\0'){
        digest += aa0;
      }
      Tryptic_.push_back(digest);
    }
   
  }
  
  
  void read(){
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


// Expose the classes
RCPP_MODULE(FastaMod) {
  class_<Fasta>("Fasta")
  .default_constructor("Default constructor") // This exposes the default constructor
  .constructor<std::string>("FASTA filename")
  .method("getNumberOfDescriptions", &Fasta::getNumberOfDescriptions, "Returns the value.")
  .method("getNumberOfSequences", &Fasta::getNumberOfSequences, "Returns the value.")
  .method("getNumberOfTrypticPeptides", &Fasta::getNumberOfTrypticPeptides, "Returns the number of tryptic peptides.")
  .method("getNumberOfAminoAcids", &Fasta::getNumberOfAminoAcids, "Returns the number of AAs")
  .method("getSequences", &Fasta::getSequences, "Returns a vector of amino acid sequences.")
 // .method("addValue", &Fasta::addValue, "Adds a value")
  //.method("merge", &Fasta::merge, "Merges another Test into this object")
  .method("getTrypticPeptides", &Fasta::getTrypticPeptides, "Returns tryptic peptides.")
  .method("summary", &Fasta::summary, "computes a summary of the FASTA object.")
  ;
}


