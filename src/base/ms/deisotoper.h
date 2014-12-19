#ifndef DEISOTOPER_H
#define DEISOTOPER_H

#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stack>
#include <cmath>

#include <iostream>

#include <base/chemistry/iisotopeenvelope.h>

/*
    Authors : Witold Eryk Wolski <wewolski@gmail.com> 
              Christian Panse <cp@fgcz.ethz.ch>

This code belongs to two projects:
    o http://cran.r-project.org/web/packages/protViz/index.html under GPLv3 
    o https://github.com/wolski/rl under three-clause BSD license

    Copyright : Functional Genomics Center Zurich | UZH | ETHZ 2013

    $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/base/ms/deisotoper.h $
    $Id: deisotoper.h 6678 2014-09-17 12:01:18Z cpanse $


*/

namespace ralab{
  namespace base{
    namespace ms{

      class Deisotoper {

        std::vector<std::string> warnings;
        std::vector < int > Z_; // charge states
        std::vector<std::vector < double > > isotopPatternMap;

        double massError_;
        double hydrogenMass_;
        ralab::base::chemistry::IIsotopeEnvelope * ie_;

        struct Peaks{
            std::vector < double > m;
            std::vector < double > aggregatedIntensity;
            std::vector < int > z;
        };

      public:

        typedef std::vector< std::vector< std::vector < int > > > ListOfVectors;
        typedef std::vector< std::vector < double > > ListOfDoubleVectors;


        struct isotopResult {
          ListOfVectors isotopChains;
          ListOfDoubleVectors isotopInnerProducts;
          ListOfDoubleVectors isotopInnerProducts1;

          std::vector < std::vector < int > >isotopPatternIdx;
          std::vector < std::vector < double > >isotopPeakScores;
          std::vector < double >isotopScores;

          std::vector < std::vector < int > >isotopGroups;

          Peaks peaks;


          const std::vector < int > & getIsotopPatten(int i) {
            if (isotopPatternIdx[i].size() > 0)
              return isotopPatternIdx[i];
            return isotopPatternIdx[i];
          }
        };

        isotopResult result_;

	const std::vector<std::string> & getWarnings(){
		return warnings;
	}

        const std::vector < std::vector < int > >  & getIsotopGroups() {
          return result_.isotopGroups;
        }

        const std::vector < double > & getPeaksMass() {
          return result_.peaks.m;
        }

        const std::vector < double > & getPeaksIntensity() {
          return result_.peaks.aggregatedIntensity;
        }

        const std::vector < int > & getPeaksCharge() {
          return result_.peaks.z;
        }

        const ListOfVectors & getIsotopChainResults() {
          return result_.isotopChains;
        }

        const ListOfDoubleVectors & getIsotopInnerProductsResults() {
          return result_.isotopInnerProducts;
        }

        const ListOfDoubleVectors & getIsotopInnerProducts1Results() {
          return result_.isotopInnerProducts1;
        }

        // constructor
        Deisotoper(double massError = 0.1 ):massError_(massError), hydrogenMass_(1.008){}

        template<class TintI>
        Deisotoper(TintI beginz , TintI endz, double massError = 0.1):
          Z_( beginz , endz ),
          massError_(massError),hydrogenMass_(1.008)
        {}

        void setIsotopPatternMap(ralab::base::chemistry::IIsotopeEnvelope * isotopeEnvlope){
          ie_ = isotopeEnvlope;
        }


        template < typename TmassI, typename TinteI >
        void 
        computeIsotopChains(TmassI beginMass, TmassI endMass, TinteI beginIntensity) {

          std::vector < int > isoptopGroupsIdx;

          std::vector <std::vector < int > >isoptopGroupsIdxChargeBegin(Z_.size());
          std::vector <std::vector < int > >isoptopGroupsIdxChargeEnd(Z_.size());

          std::vector<std::vector< int > > results ;
          for (unsigned int zi = 0; zi < Z_.size(); zi++) {


              computeDeisotopCandidates(beginMass, endMass, Z_[zi], results);

              result_.isotopChains.push_back(results);

              // preprocessing for ISOTOP GROUPING
              // step 0: determine all isotop groups for all charge states
              for (std::vector<std::vector< int > >::iterator it = results.begin(); it !=  results.end(); ++it){

                isoptopGroupsIdx.push_back((*it).front());
                isoptopGroupsIdxChargeBegin[zi].push_back((*it).front());
                isoptopGroupsIdxChargeEnd[zi].push_back((*it).back());

                (*it).clear();
              }
              results.clear();
          }

            // ISOTOP GROUPING 
            std::sort (isoptopGroupsIdx.begin(), isoptopGroupsIdx.end());
            std::unique (isoptopGroupsIdx.begin(), isoptopGroupsIdx.end());

            std::vector<int>::iterator it_NN_begin;
            std::vector<int>::iterator it_NN_end;

            int oldGroupIdx=-1;
            // for all isotop groups
            for(std::vector< int >::iterator it = isoptopGroupsIdx.begin(); it !=isoptopGroupsIdx.end(); ++it){

                // std::unique cannot alter the properties of the object containing the range of elements
                if (*it <= oldGroupIdx) break;
                oldGroupIdx = *it;
                
                // default -1 means no isotop for this charge state
                std::vector< int > group(Z_.size(), -1);

                // for all charge states
                for (unsigned int zi = 0; zi < Z_.size(); zi++) {

                    // if we have a isotop group for charge state zi
                    if (std::distance(isoptopGroupsIdxChargeBegin[zi].begin(), isoptopGroupsIdxChargeBegin[zi].end()) > 0){

                        // compute lower and upper bound for a isotop group intersection
                        it_NN_begin = std::lower_bound(isoptopGroupsIdxChargeBegin[zi].begin(), isoptopGroupsIdxChargeBegin[zi].end(), *it);
                        it_NN_end = std::upper_bound(isoptopGroupsIdxChargeEnd[zi].begin(), isoptopGroupsIdxChargeEnd[zi].end(), *it);

                        // the end of the list is not the last entry!
                        if (it_NN_begin == isoptopGroupsIdxChargeBegin[zi].end())
                            --it_NN_begin;
                        if (it_NN_end == isoptopGroupsIdxChargeEnd[zi].end())
                            --it_NN_end;


                        // determine if thhere is a isotop group intersection;  
                        if ( (*it_NN_begin) <=  (*it) && (*it) <=  (*it_NN_end) )
                            group[zi] = std::distance(isoptopGroupsIdxChargeBegin[zi].begin(),it_NN_begin);

                    } 
                }
                // save the result
                result_.isotopGroups.push_back(group);
            }
        }

        template < typename TmassI, typename TinteI >
        void 
        assignIsotopIntensities(TmassI beginMass, TmassI endMass, TinteI beginIntensity){

            double dIntensityL2,  dIntensityInnerProduct0, dIntensityInnerProduct1;

            // for all charge states
            for (std::vector< std::vector< std::vector< int > > >::iterator it_z = result_.isotopChains.begin(); it_z !=  result_.isotopChains.end(); ++it_z){
                double z = Z_[std::distance(result_.isotopChains.begin(), it_z)];

                std::vector < double > Lscore, Lscore1;

                // for all chains
                for (std::vector< std::vector< int > >::iterator it_c = (*it_z).begin(); it_c !=  (*it_z).end(); ++it_c){
                    // determine queryMass for insilico isotope envelope query
                    double queryMass = *(beginMass + *(*it_c).begin());
                    std::vector<double> LIsotopeEnvelope = ie_ -> isotopeEnvelope(z * queryMass);

                    // compute L2 norm for selected isotop chain
                    dIntensityL2 = 0.0;
                    dIntensityL2 = sqrt (computeL2norm(beginIntensity, (*it_c).begin(), (*it_c).end(), dIntensityL2));
                    
                    // TODO here we have to find the mono iso peak!!!  
                    // the winner is tbd
                    dIntensityInnerProduct0 = 0.0;
                    dIntensityInnerProduct0 = inner_product(beginIntensity, (*it_c).begin(), (*it_c).end(), LIsotopeEnvelope.begin(), LIsotopeEnvelope.end(), dIntensityInnerProduct0, dIntensityL2);
                    Lscore.push_back(dIntensityInnerProduct0);

                    dIntensityInnerProduct1 = 0.0;
                    std::vector< int > ::iterator it_c_begin_succ = (*it_c).begin();
                    ++it_c_begin_succ;
                    dIntensityInnerProduct1 = inner_product(beginIntensity, it_c_begin_succ, (*it_c).end(), LIsotopeEnvelope.begin(), LIsotopeEnvelope.end(), dIntensityInnerProduct1, dIntensityL2);
                    Lscore1.push_back(dIntensityInnerProduct1);

                    /*
                    int l0 = std::distance((*it_c).begin(), (*it_c).end());
                    int l1 = std::distance(it_c_begin_succ, (*it_c).end());

                    if (dIntensityInnerProduct0 > 0.95 || (l1 > 2 && dIntensityInnerProduct1 > 0.95)){
                        std::cout << z << "\t" << l0 << "\tscore0=" << l0 * dIntensityInnerProduct0 << "\tscore1=" << l1 * dIntensityInnerProduct1 << std::endl;
                    }
                    */

                } //  END for all chains
                result_.isotopInnerProducts.push_back(Lscore);
                result_.isotopInnerProducts1.push_back(Lscore1);
            } // END for all charge states z
         }

        template < typename TmassI, typename TinteI >
        void 
        deisotopPeaklist(TmassI beginMass, TmassI endMass, TinteI beginIntensity){
            int nMass = std::distance(beginMass, endMass);

            for (int i = 0; i < nMass; i++){   
                result_.peaks.m.push_back(*(beginMass + i));
                result_.peaks.aggregatedIntensity.push_back(*(beginIntensity + i));
                result_.peaks.z.push_back(0);
            }
            //for (int i=0; i<nMass;i++)   
            //   std::cout << result_.peaks.m[i] << "\t" << mark[i] << "\t" << result_.peaks.aggregatedIntensity[i] << std::endl;

             std::vector < int > mark(nMass, 0);

             int groupIdx, chainIdx, zIdx;

            double monoIsotopicMass, monoIsotopicIntensity;  
            int nPeakCand = -1; 
            double eps = 0.90;

            // assume that the first peak is the monoisotopic peak
            int max_zIdx = -1;
            int max_chainIdx  = -1;
            int max_peakShift = 0;
            double max_score = -1.0;

            int peakShift;
            double score;

            for(std::vector< std::vector< int > >::iterator itGroup=result_.isotopGroups.begin(); itGroup != result_.isotopGroups.end(); ++itGroup){
                max_score = -1.0; max_zIdx = -1; max_peakShift = 0; 
                max_chainIdx = -1;
                groupIdx = std::distance(result_.isotopGroups.begin(), itGroup);
                
                // forall charge idx
                for (std::vector< int >::iterator it_z = (*itGroup).begin(); it_z !=  (*itGroup).end(); ++it_z){

                    chainIdx = *it_z;
                    zIdx = std::distance((*itGroup).begin(), it_z);

                    if (chainIdx > -1 ){
                        nPeakCand = std::distance(result_.isotopChains[zIdx][chainIdx].begin(), result_.isotopChains[zIdx][chainIdx].end());

                        if (result_.isotopInnerProducts[zIdx][chainIdx] > result_.isotopInnerProducts1[zIdx][chainIdx] || nPeakCand < 3){
                            score = result_.isotopInnerProducts[zIdx][chainIdx] * nPeakCand;
                            peakShift = 0;
                        } else {
                            score = result_.isotopInnerProducts1[zIdx][chainIdx] * (nPeakCand - 1); 
                            peakShift = 1;
                        }

                        if (score > max_score){
                            max_score = score;
                            max_chainIdx = chainIdx;
                            max_zIdx  = zIdx; 
                            max_peakShift = peakShift;
                        }
                    }
                }

// monoisotop peak is determined prepare the output
                nPeakCand = std::distance(result_.isotopChains[max_zIdx][max_chainIdx].begin(), result_.isotopChains[max_zIdx][max_chainIdx].end());
                if (result_.isotopInnerProducts.size() <= max_zIdx && result_.isotopInnerProducts[max_zIdx].size() <= max_chainIdx){
                    if (result_.isotopInnerProducts[max_zIdx][max_chainIdx] > eps || (nPeakCand > 3 && result_.isotopInnerProducts1[max_zIdx][chainIdx] > eps)){
                    // chose mono isotopic peak

                    int monoIsotopicPeakIdx = *(result_.isotopChains[max_zIdx][max_chainIdx].begin()) + max_peakShift;
                    monoIsotopicMass=*(beginMass + monoIsotopicPeakIdx);
                    monoIsotopicIntensity = *(beginIntensity  + monoIsotopicPeakIdx);

		
		    if (std::distance(result_.isotopChains[max_zIdx][max_chainIdx].begin() , result_.isotopChains[max_zIdx][max_chainIdx].end())< max_peakShift + 2){

				std::ostringstream sWarning;
                                sWarning << "isotop peak group to small."; 
				warnings.push_back(sWarning.str());

			} else{
                    std::vector<int>::iterator it_BeginIsotopPeakGroup = (result_.isotopChains[max_zIdx][max_chainIdx].begin() +  max_peakShift + 1);
                    std::vector<int>::iterator it_EndIsotopPeakGroup = result_.isotopChains[max_zIdx][max_chainIdx].end();
		    std::vector<int>::iterator it_chain;
			
                    for (it_chain = it_BeginIsotopPeakGroup; it_chain !=  it_EndIsotopPeakGroup; ++it_chain){

                            // peak can only be used once
                            if (mark[*it_chain] == 0){
                                monoIsotopicIntensity += *(beginIntensity + *it_chain);
                                mark[*it_chain] = monoIsotopicPeakIdx;
                                result_.peaks.aggregatedIntensity[*it_chain] = 0.0;
                            }
                            else{
				std::ostringstream sWarning;

                                sWarning << "Peak " << result_.peaks.m[*it_chain] << "[" \
<< *it_chain << "] / " << Z_[max_zIdx] << " can only be used once. Already in use with (idx, m, Z, I) = (" \
<< *it_chain << ", " << result_.peaks.m[mark[*it_chain]] << ", "  << result_.peaks.z[mark[*it_chain]] \
<< ", " << result_.peaks.aggregatedIntensity[mark[*it_chain]] << ").";

				warnings.push_back(sWarning.str());
                            }
                    }

                    if (result_.peaks.m[monoIsotopicPeakIdx] != monoIsotopicMass){
			warnings.push_back("ERROR: contradiction in monoIsotopicMass");
		    }
                    result_.peaks.aggregatedIntensity[monoIsotopicPeakIdx] = monoIsotopicIntensity;
                    result_.peaks.z[monoIsotopicPeakIdx] = Z_[max_zIdx];

                    it_chain = result_.isotopChains[max_zIdx][max_chainIdx].begin();
                }
		}}else{
			//warnings.push_back("prevent use of uninitialised memory in methode 'result_.isotopInnerProducts'.");
            }
            }//foreach group (isotop cluster)
         }

      private:
        template <typename TinteI, class InputIterator1, class T>
        T computeL2norm (TinteI beginIntensity, InputIterator1 first1, InputIterator1 last1, T init)
        {
            while (first1!=last1) {
                init = init + ((*(beginIntensity + *first1)) * (*(beginIntensity + *first1)));
                ++first1; 
            }
            return init;
        }

        // modified code taken from http://www.cplusplus.com/reference/numeric/inner_product/
        // first2 is reqiered because both vectors can have different length
        template <typename TinteI, class InputIterator1, class InputIterator2, class T>
        T inner_product (TinteI beginIntensity, InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, T init, T norm)
        {
            while (first1!=last1 && first2!=last2) {
                init = init + ((*(beginIntensity + *first1)) * (*first2) / norm);
                ++first1; ++first2;
            }
            return init;
        }

        // non recursive depth first search for special case 
        int explore(int v, std::vector < int >&G,
                    std::vector < bool > & VISITED,
                    std::vector< int > & result
                    ){

          std::stack < int > S;
          S.push(v);

          while (!S.empty()) {
              v = S.top();
              S.pop();
              result.push_back(v);

              if (!VISITED[v]) {
                  VISITED[v] = true;
                  if (G[v] > -1)
                    S.push(G[v]);
                }
            }

          return 0;
        }

        template < class Iterator, class T > Iterator findNearestNeighbor(Iterator first,
                                                                          Iterator last,
                                                                          const T & val) {
          Iterator it = std::upper_bound(first, last, val);

          if (it == last)
              return (it--);
          if (it != first && it - 1 != first) {
              if (std::fabs(val - *(it - 1)) < std::fabs(val - *it) || it == last) {
                 // return (it - 1);
                 it--;
                }
            }

          return (it);
        }

        /*
       graph G(V,E);
       each peak represents one node of a G so V = {1, ..., n}
       the edges are defined by (v, node_array_G[v]) iff node_array_G[v] > -1;

       simply speaking an edge is defined between two peaks if the mZ difference is
       below the given massError.

       TODO: for efficiency we can have a list containg only the candidates within the massError
       having those list we do not have to iterate over the whole graph again.
     */
        template < class Iterator >
        int generateGraph(std::vector <int >&G,
                          int z,
                          Iterator firstMz,
                          Iterator lastMz) {
          int idx=-1;
          for (Iterator itMz = firstMz; itMz != lastMz; ++itMz)
            {

              double dQmZ = (*itMz + (hydrogenMass_ / z));
              Iterator itNNmZ = findNearestNeighbor(firstMz, lastMz, dQmZ);

              idx=std::distance(firstMz, itMz);

              if (itNNmZ != lastMz && std::fabs(*itNNmZ - dQmZ) < massError_ && idx <= G.size()){
                G[idx] = std::distance(firstMz, itNNmZ);
              }
            }
          return 0;
        }

        template < typename TmassI > int computeDeisotopCandidates(TmassI beginMass,
                                                                   TmassI endMass,
                                                                   int z,
                                                                   std::vector< std::vector<int> > & candidates) {
          int n = std::distance(beginMass, endMass);
          std::vector < int >G(n, -1);
          generateGraph(G, z, beginMass, endMass);

          std::vector < bool > VISITED(n, false);

          for (unsigned int i = 0; i < G.size(); i++) {
              if (G[i] > -1) {
                  std::vector<int> result;
                  explore(i, G, VISITED,result);
                  if(result.size() > 1){
                      candidates.push_back(result);
                    }
                }
            }
          return 0;
        }

      };
    }
  }
}
#endif				// DEISOTOPER_H
