#include <Rcpp.h>
#include <base/ms/deisotoper.h>
#include <rcppdeisotoperenvelope.h>

/*

$HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/deisotoper.cpp $
$Id: deisotoper.cpp 6179 2014-02-27 09:34:04Z cpanse $

*/


void printWarning(std::string s){
	::Rf_warning(s.c_str());
}

RcppExport SEXP deisotoper_main(SEXP Dsexp_mZ_, SEXP Dsexp_intensity_, SEXP Dsexp_Z_, SEXP Dsexp_DF_, SEXP Dsexp_massError_)
{
    try {			
	
    Rcpp::DataFrame DF = Rcpp::DataFrame(Dsexp_DF_);

	Rcpp::NumericVector mZ_ = Rcpp::NumericVector(Dsexp_mZ_);
	Rcpp::NumericVector intensity_ =
	    Rcpp::NumericVector(Dsexp_intensity_);

	Rcpp::IntegerVector Z = Rcpp::IntegerVector(Dsexp_Z_);

    double dMassError = Rcpp::as<double>(Dsexp_massError_);

    double dChargeCutOffMass;
    for (unsigned int zi = 0; zi < Z.size(); zi++){
	dChargeCutOffMass = (1/(double)Z[zi]);
	if (dMassError > dChargeCutOffMass){
		//dMassError = dChargeCutOffMass;
		std::ostringstream sWarning;
		sWarning << "deisotoper_main: mass error schould be smaller than 1 / Z. set massError=" \
		<< dMassError << ".";
		printWarning(sWarning.str());
	}
    }

    ralab::base::ms::RcppIsotopeenvelope ril(DF);
    ralab::base::ms::Deisotoper x(Z.begin(), Z.end(), dMassError);
    x.setIsotopPatternMap(& ril);

    x.computeIsotopChains(mZ_.begin() , mZ_.end(), intensity_.begin());
    x.assignIsotopIntensities(mZ_.begin() , mZ_.end(), intensity_.begin());
    x.deisotopPeaklist(mZ_.begin() , mZ_.end(), intensity_.begin());

	Rcpp::List NDF =
	    Rcpp::List::create(Rcpp::Named("result") = x.getIsotopChainResults(), 
        Rcpp::Named("score") = x.getIsotopInnerProductsResults(),
        Rcpp::Named("score1") = x.getIsotopInnerProducts1Results(),
        Rcpp::Named("group") = x.getIsotopGroups(),
        Rcpp::Named("peaksMass") = x.getPeaksMass(),
        Rcpp::Named("peaksIntensity") = x.getPeaksIntensity(),
        Rcpp::Named("peaksCharge") = x.getPeaksCharge()
        );

	
	std::vector<std::string> warnings = x.getWarnings();
	std::for_each(warnings.begin(), warnings.end(), printWarning);

	return (NDF);
    }
    catch(std::exception & ex) {	// or use END_RCPP macro
	    forward_exception_to_r(ex);
    }
    catch( ...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
