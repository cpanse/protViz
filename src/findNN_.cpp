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


class MyComparator : public std::binary_function<double, double, bool> {
public : 
  bool operator()(double a, double b){
    return (a <= b);
  }
};

// [[Rcpp::export]]
IntegerVector lower_bound__(const NumericVector xq, const NumericVector xvec)
{
  IntegerVector idxvec;
  MyComparator comparator;
  try {
    
    for (int x:xq){
      
      idxvec.push_back(std::distance(xvec.begin(),
                                     std::lower_bound(xvec.begin(), xvec.end(), x, comparator)));
    }
    return (idxvec);
  }
  catch( ...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  
  // TODO(cp): handle possible exceptions
}