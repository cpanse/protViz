
#ifndef RCPPDEISOTOPERENVELOPE_H
#define RCPPDEISOTOPERENVELOPE_H


/*

Authors   : Witold Eryk Wolski <wewolski@gmail.com> and Christian Panse <cp@fgcz.ethz.ch>

This code belongs to two projects:
o http://cran.r-project.org/web/packages/protViz/index.html under GPLv3 
o https://github.com/wolski/rl under three-clause BSD license

Copyright : Functional Genomics Center Zurich | UZH | ETHZ and  ETH Zurich 2013


$HeadURL: $
$Id: $


*/

#include <math.h>
#include <stdint.h>

#include <iostream>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
 #include <sstream>

#include <Rcpp.h>

#include <base/chemistry/iisotopeenvelope.h>
#include <rcppdeisotoperenvelope.h>


namespace ralab{
  namespace base{
    namespace ms{

      class RcppIsotopeenvelope : public  ralab::base::chemistry::IIsotopeEnvelope {
            Rcpp::DataFrame df_;

            public:


            RcppIsotopeenvelope(Rcpp::DataFrame & df):df_(df){
                
            }

        template < class Iterator, class T > Iterator findNearestNeighbor(Iterator first,
                                                                          Iterator last,
                                                                          const T & val) {
          Iterator it = std::upper_bound(first, last, val);

          if (it != first) {
              Iterator it_pre = it - 1;
              if (fabs(val - *it_pre) < fabs(val - *it) || it == last) {
                  return (it_pre);
                }
            }

          return (it);
        }

            double string_to_double( const std::string& s ) { std::istringstream i(s); double x; if (!(i >> x)) return 0; return x; } 

            std::vector<double> isotopeEnvelope(double mass ) 
            {
               std::vector<std::string> sNames = df_.attr("names");
               std::vector< double > dNames;

               for(std::vector<std::string>::iterator it = sNames.begin(); it != sNames.end(); ++it){
                    dNames.push_back(string_to_double(*it));
               }

               std::vector<double>::iterator idx = findNearestNeighbor(dNames.begin(), dNames.end(), mass);
               std::vector<double> L =  df_[ (int)std::distance(dNames.begin(), idx) ];
               return L;
            }

      };
    }
  }
}
#endif				// RCPPDEISOTOPERENVELOPE_H
