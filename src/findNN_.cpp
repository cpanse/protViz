#include <algorithm>
#include <cmath>


/*

$HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/findNN_.cpp $
$Id: findNN_.cpp 6603 2014-08-12 11:17:31Z cpanse $

*/

extern "C" {
    void findNN_ (int *m_, int *n_, double *q_, double *vec_, int *NN_) {

        size_t dist;
        double d;

        for (int i = 0; i < *m_; i++){

            
            dist =  std::distance (vec_, std::lower_bound(vec_, vec_ + *n_, q_[i]) );
            NN_[i] = dist ;

            if (dist > 0){
                d = std::fabs(q_[i] - vec_[dist - 1]);

                if (dist <  *n_){
                    if (d < std::fabs(q_[i] - vec_[dist]))
                        NN_[i] = dist - 1;
                }else if (NN_[i] >= *n_){
                        NN_[i]--;
                }

            }
        }
    }
}
