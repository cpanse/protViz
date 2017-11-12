// Christian Panse <cp@fgcz.ethz.ch>
// 20171111

#include <Rcpp.h>
#include <numeric>
#include <string>

using namespace Rcpp;

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

