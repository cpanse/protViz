// Christian Panse <cp@fgcz.ethz.ch>
// 20171111

#include <Rcpp.h>
#include <numeric>
#include <string>


namespace protVizFasta{

  std::vector<std::string> fcat(const std::vector<std::string>& fasta) {
    //StringVector rv;
    std::vector<std::string> rv;
      std::string aa_seq = "";
  
      for (auto s : fasta) {
      	if (s.at(0) == '>') {
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

}

using namespace Rcpp;

// [[Rcpp::export]]
StringVector fastacat(const StringVector& fasta){
  StringVector rv;
  
  return rv;
}
