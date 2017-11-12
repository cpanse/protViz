/// Christian Panse <cp@fgcz.ethz.ch>
/// 2006-2017

#include <Rcpp.h>
#include <numeric>
#include <string>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector findNN_(const NumericVector &q, const NumericVector & vec, bool check = false){
  IntegerVector NN(q.size(), -1);
  int n = vec.size();
  
  if (check){
    if (!std::is_sorted(vec.begin(), vec.end())){
      // std::sort(vec.begin(), vec.end());
      return (NN);
    }
  }
  
  size_t dist;
  double d;
  
  for(int i = 0; i < q.size(); i++) {
    
    dist = std::distance (vec.begin(), std::lower_bound(vec.begin(), vec.end(), q[i]));
    
    NN[i] = dist;
    
    if (dist > 0){
      d = std::fabs(q[i] - vec[dist - 1]);
      
      if (dist <  n){
        if (d < std::fabs(q[i] - vec[dist]))
          NN[i] = dist - 1;
      }else if (NN[i] >= n){
        NN[i]--;
      }
    }
  }
  
  std::for_each(NN.begin(), NN.end(), [&](int &c){c++;});
  
  return(NN);
}
