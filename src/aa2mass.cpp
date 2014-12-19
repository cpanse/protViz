/*
 * aa2mass.cpp R package protViz
 * parts taken from msmsseqass.c deltatectra
 *
 * Copyright  2013
 * Christian Panse <cp@fgcz.ethz.ch>
 *
 * This file is part of deltatectra.
 *
 * deltatectra is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA *


 $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/aa2mass.cpp $
 $Id: aa2mass.cpp 6179 2014-02-27 09:34:04Z cpanse $

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
    void init() {

	int letter;

    int i;
    for (i = 0; i < 26; i++)
        dAA_mass_[i] = 0;

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

RcppExport SEXP aa2mass_main(SEXP Dsexp_peptides_, SEXP Dsexp_mass_,
			     SEXP Dsexp_aa_)
{
    try {			// or use BEGIN_RCPP macro

	Rcpp::StringVector peptides_ = Rcpp::StringVector(Dsexp_peptides_);

	Rcpp::NumericVector mass_ = Rcpp::NumericVector(Dsexp_mass_);

	Rcpp::StringVector aa_ = Rcpp::StringVector(Dsexp_aa_);

	protViz::Aa2Mass myAa2Mass(mass_.begin(), mass_.end(), aa_.begin(), aa_.end());
	//myAa2Mass.init();

	Rcpp::List NDF = Rcpp::List::create(Rcpp::Named("input") = peptides_, Rcpp::Named("output") = myAa2Mass.process(peptides_.begin(), peptides_.end()));
    
	return (NDF);

    }
    catch(std::exception & ex) {	// or use END_RCPP macro
	forward_exception_to_r(ex);
    }
    catch( ...) {
	::Rf_error("C++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
