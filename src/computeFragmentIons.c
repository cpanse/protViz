/*
 * computeFragmentIons.c R package protViz
 * parts taken from msmsseqass.c deltatectra
 *
 * Copyright 2007, 2008, 2009, 2010, 2012
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
*/

#include <stdio.h>
#include <stdlib.h>

/*

NOTE:

-- setting fixed modification
    ->Carbamidomethyl (C)
        ->C 160.030649


    M[2] = 160.030649;

    Cysteine mono mass is 103.00919


$HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/computeFragmentIons.c $
$Id: computeFragmentIons.c 6535 2014-06-18 15:21:13Z cpanse $


*/

double *initAminoAcidMass(void)
{

    double *M;
    if ((M = (double *) malloc(26 * sizeof(double))) == NULL) {
	return (NULL);
    }

    /* A Alanine Ala */
    M[0] = 71.037110;
    /* B ***************** */
    M[1] = 114.534940;
    /* C Cysteine Cys */
    M[2] = 160.030649;
    /* D Aspartic Asp */
    M[3] = 115.026940;
    /* E Glu Gluteami Acid */
    M[4] = 129.042590;
    /* F Phe Phenyalanine */
    M[5] = 147.068410;
    /* G Gly Glycine */
    M[6] = 57.021460;
    /* H His Histidine  */
    M[7] = 137.058910;
    /* I Isoleucince !!! same mass as L */
    M[8] = 113.084060;
    /* J ****************** */
    M[9] = 0.000000;
    /* K Lys Lysine */
    M[10] = 128.094960;
    /* L Leu Leucine */
    M[11] = 113.084060;
    /* M Met Methionine */
    M[12] = 131.040480;
    /* N Asparagine Asn */
    M[13] = 114.042930;
    /* O ***************** */
    M[14] = 0.000000;
    /* P Pro Proline */
    M[15] = 97.052760;
    /* Q Gln Glutamine */
    M[16] = 128.058580;
    /* R Arg Arginine */
    M[17] = 156.101110;
    /* S Ser Serine */
    M[18] = 87.032030;
    /* Thr Threonine */
    M[19] = 101.047680;
    /* U ****************** */
    M[20] = 150.953630;
    /* V Val Valine  */
    M[21] = 99.068410;
    /* W Trp Tryptophan */
    M[22] = 186.079310;
    /* X ******************* */
    M[23] = 111.000000;
    /* Y Tyr Tyrosine */
    M[24] = 163.063330;
    /* Z ******************* */
    M[25] = 128.550590;

    return (M);
}

void computeFragmentIons(int *n_, char **seq_, double *pim_, double *b_,
			 double *y_)
{

    int i;
    double b, y;
    double C_term;
    double N_term;
    double Electron;
    double Hydrogen;
    int letter;

    double *M;
    M = initAminoAcidMass();

    if (M != NULL) {

	C_term = 17.002740;
	N_term = 1.007825;
	Electron = 0.000549;
	Hydrogen = 1.007825;

	b = N_term - Electron;
	y = *pim_;

	for (i = 0; i < *n_; i++) {
	    letter = seq_[0][i];
	    if (64 < letter && letter < 92) {
		b += M[letter - 65];
		b_[i] = b;

		y_[*n_ - i - 1] = y;
		y -= M[letter - 65];
	    }
	}
    }
    free(M);
}



void computeFragmentIonsModification(int *n_, char **seq_, double *pim_,
				     double *b_, double *y_,
				     int *modified_, double *modification_)
{

    int i;
    double b, y;
    double C_term;
    double N_term;
    double Electron;
    double Hydrogen;
    int letter;

    double *M;
    M = initAminoAcidMass();
    if (M != NULL) {

	C_term = 17.002740;
	N_term = 1.007825;
	Electron = 0.000549;
	Hydrogen = 1.007825;

	b = N_term - Electron;
	y = *pim_;

	for (i = 0; i < *n_; i++) {

	    letter = seq_[0][i];

	    if (64 < letter && letter < 92) {
		b += M[letter - 65] + modification_[modified_[i]];
		b_[i] = b;

		y_[*n_ - i - 1] = y;
		y -= (M[letter - 65] + modification_[modified_[i]]);
	    }
	}
    }
    free(M);
}

void computeFragmentIonsFixedVariableModification(int *n_, char **seq_,
						  double *pim_, double *b_,
						  double *y_,
						  int *modified_,
						  double *modification_,
						  double
						  *fixedMmodification)
{

    int i;
    double b, y;
    double C_term;
    double N_term;
    double Electron;
    double Hydrogen;
    int letter;

    double *M;
    M = initAminoAcidMass();

    for (i = 0; i < 26; i++) {
	if (fixedMmodification[i] > 0.0) {
	    M[i] = fixedMmodification[i];
	} else {
	    fixedMmodification[i] = M[i];
	}
    }

    C_term = 17.002740;
    N_term = 1.007825;
    Electron = 0.000549;
    Hydrogen = 1.007825;

    b = N_term - Electron;
    y = *pim_;

    for (i = 0; i < *n_; i++) {

	letter = seq_[0][i];

	if (64 < letter && letter < 92) {
	    b += M[letter - 65] + modification_[modified_[i]];
	    b_[i] = b;

	    y_[*n_ - i - 1] = y;
	    y -= (M[letter - 65] + modification_[modified_[i]]);
	}
    }
    free(M);
}



/*
debug:

library(protViz)
peptide.AA<-"KINHSFLR"

peptide.AA.weights<-c(128.09496,113.08406,114.04293,137.05891,87.03203,147.06841,113.08406,156.10111)

n<-as.integer(length(peptide.AA.weights))

fragmentIons(peptide.AA)
.C("_computeFragmentIons", n=n, W_=as.double(peptide.AA.weights), b_=as.double(rep(0,n)), y_=as.double(rep(0,n)))

*/
void _computeFragmentIons(int *n_, double *W_, 
				     double *b_, double *y_)
{
    int i;
    double b;
    double y;
    double C_term;
    double N_term;
    double Electron;
    double Hydrogen;

    double *M;
    M = initAminoAcidMass();

	C_term = 17.002740;
	N_term = 1.007825;
	Electron = 0.000549;
	Hydrogen = 1.007825;

	b = N_term - Electron;

    /* compute parent ion mass */
	y = C_term + N_term + Hydrogen - Electron;

	for (i = 0; i < *n_; i++) {
        y += W_[i];
    }

	for (i = 0; i < *n_; i++) {

		b += W_[i]; 
		b_[i] = b;

		y_[*n_ - i - 1] = y;
		y -= W_[i];
	}
    free(M);
}

