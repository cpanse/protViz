/*
 * aa2mass.cpp R package protViz parts taken from msmsseqass.c deltatectra
 * 
 * Copyright  2013 Christian Panse <cp@fgcz.ethz.ch>
 * 
 * This file is part of deltatectra.
 * 
 * deltatectra is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version. This program is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details. You should have received a copy of the
 * GNU General Public License along with this program; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301  USA *
 * 
 * 
 * $HeadURL:
 * http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/aa2mas
 * s.cpp $ $Id: aa2mass.cpp 6179 2014-02-27 09:34:04Z cpanse $
 * 
 */

#include <string>
#include <iterator>
#include <vector>

#include <Rcpp.h>


namespace protViz{
class Aa2Mass {
    std::vector < std::string > peptides_;
    std::vector < double >AA_mass_;
    std::vector < std::string > AA_names_;


    /*
        dAA_mass_ : {Alanin, ..., Valine} -> {A, ..., Z}

        |{Alanin, ..., Valine}| = 20  
        |{A, ..., Z}| = 26
    */
    double dAA_mass_[26];

  public:
    typedef std::vector< std::vector < double  > > ListOfVectors;

     template < class TaaD, class TaaC >
        Aa2Mass(TaaD beginMass,
                TaaD endMass, TaaC beginAA,
                TaaC endAA): AA_mass_(beginMass, endMass), AA_names_(beginAA, endAA) {
            init();
    }

    /*

    initialze the dAA_mass_[20] array
    */
    void init(void) {

        int letter;

    long unsigned int i;
    for (i = 0; i < 26; i++){
        dAA_mass_[i] = 0;
    }

        for (i = 0; i < AA_names_.size(); i++) {
            letter = *(AA_names_[i].c_str());

            if (64 < letter && letter < 92)
                    dAA_mass_[letter - 65] = AA_mass_[i];
        }
    } // init

    template < class TaaS >
    ListOfVectors process(TaaS beginPeptides, TaaS endPeptides) {

        int letter;

    ListOfVectors result;

    std::string pep;
        std::string::iterator s_it;

        for (TaaS  vs_it = beginPeptides; vs_it < endPeptides; ++vs_it) {
        std::vector<double> L;
            pep = *vs_it;

        for ( s_it = pep.begin(); s_it < pep.end(); s_it++){
            letter = *s_it;

                if (64 < letter && letter < 92)
                L.push_back(dAA_mass_[letter - 65]);
            else
            {
                L.push_back(0.0);
                Rcpp::stop("Inadmissible value");
            }
        }

        result.push_back(L);
    }
        return result;
    } // process
};
}


//' Determines the weight of each amino acid sequence given as input.
//'
//' @details For the computation no C-Term and N-Term is considered.
//' \code{aa2mass} is useful if you have AA modifications. 
//'
//' @author Christian Panse 2014, 2023
//' @param peptides_ a amino acid sequence (peptide) 
//' @param aa_ amino acid (letter 1 code) in the same size and order as the \code{mass_} attribute.
//' @param mass_ vector of size 20 containing the weight of the AA.
//' @importFrom Rcpp evalCpp
//' @useDynLib protViz, .registration=TRUE
//' @export
//' @examples
//' mass <- c(71.03711, 156.10111, 114.04293, 115.02694, 103.00919,
//' 129.04259, 128.05858, 57.02146, 137.05891, 113.08406, 113.08406,
//' 128.09496, 131.04049, 147.06841, 97.05276, 87.03203, 101.04768,
//' 186.07931, 163.06333, 99.06841)
//'
//' letter1 <- c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
//' 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
//'
//' peptide <- "HTLNQIDSVK"
//'
//' aa2mass(peptide, mass, letter1)
//' peptides <- c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')
//' C_term <- 17.002740
//' N_term <- 1.007825
//' H_ <- 1.008
//' unlist(lapply(aa2mass(peptides), sum)) + C_term + N_term + H_ - parentIonMass(peptides)
//'
//' ## determine the fragment ions
//' lapply(aa2mass(peptides), function(x){protViz::fragmentIon(x)[[1]]})
//'
//' ## an example with [STY] AA modification
//' peptide <- ' HTLNQIDSVK '
//' mod <- rep(0.0, nchar(peptide)); mod[8] <- 79.998;
//' aa2mass(peptide)[[1]] + mod
//[[Rcpp::export]]
Rcpp::List aa2mass(Rcpp::StringVector peptides_ , Rcpp::Nullable<Rcpp::NumericVector> mass_ = R_NilValue, Rcpp::Nullable<Rcpp::StringVector> aa_=R_NilValue)
{
    	std::vector<double> mass({71.03711, 156.10111, 114.04293, 115.02694, 103.00919,
129.04259, 128.05858, 57.02146, 137.05891, 113.08406, 113.08406,
128.09496, 131.04049, 147.06841, 97.05276, 87.03203, 101.04768,
186.07931, 163.06333, 99.06841});
  	std::vector<std::string> aastd({"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"});
	Rcpp::StringVector aa(aastd.size());
	aa = aastd;

     if (mass_.isNotNull()) {
        Rcpp::NumericVector mmm (mass_);
	mass =Rcpp::as<std::vector<double> >(mmm);
     }else{
  }
     if (aa_.isNotNull()) {
	     Rcpp::StringVector aaa (aa_);
	     aa = aaa;
     }else{

     }

	try{

        protViz::Aa2Mass myAa2Mass(mass.begin(), mass.end(), aa.begin(), aa.end());
        //myAa2Mass.init();

        Rcpp::List NDF = Rcpp::List::create(Rcpp::Named("input") = peptides_, Rcpp::Named("output") = myAa2Mass.process(peptides_.begin(), peptides_.end()));

        return (NDF);

    }
    catch(std::exception & ex) {        // or use END_RCPP macro
        forward_exception_to_r(ex);
    }
    catch( ...) {
        ::Rf_error("C++ exception (unknown reason)");
    }
    return R_NilValue;          // -Wall
}

