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
  
  if (check){
    if (!std::is_sorted(vec.begin(), vec.end())){
      // std::sort(vec.begin(), vec.end());
      return (NN);
    }
  }
  
  double d;
  size_t dist;
  size_t n = vec.size();
  
  for(int i = 0; i < static_cast<int> (q.size()); i++) {
    
    dist = std::distance (vec.begin(), std::lower_bound(vec.begin(), vec.end(), q[i]));
    
    NN[i] = dist;
    
    if (dist > 0){
      d = std::fabs(q[i] - vec[dist - 1]);
      
      if (static_cast<size_t>(dist) <  n){
        if (d < std::fabs(q[i] - vec[dist]))
          NN[i] = dist - 1;
      }else if (static_cast<size_t>(NN[i]) >= n){
        NN[i]--;
      }
    }
  }
  
  std::for_each(NN.begin(), NN.end(), [&](int &c){c++;});
  
  return(NN);
}


// [[Rcpp::export]]
IntegerVector lower_bound__(const NumericVector xq, const NumericVector xvec)
{
  IntegerVector idxvec;
  try {
    
    for (int x:xq){
      
      idxvec.push_back(std::distance(xvec.begin(),
                                     std::lower_bound(xvec.begin(), xvec.end(), x, [](const double& a, const double& b){ return a <= b; })));
    }
    return (idxvec);
  }
  catch( ...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  
  // TODO(cp): handle possible exceptions
}
