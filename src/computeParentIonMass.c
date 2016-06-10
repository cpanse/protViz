#include <stdio.h>

/*
 *  computeParentIonMass.c
 *
 * Copyright 2006
 * Christian Panse <cp@fgcz.ethz.ch>
 *
 * This file is part of the R-package protViz.
 *
 * protViz is free software; you can redistribute it and/or modify
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


/*
 * Christian Panse <cp@fgcz.ethz.ch> 
 * Mit Nov  1 17:50:33 CET 2006
 *
 * DESCRIPTION:
 * digest conditions:
 * cleave AFTER arginine (R) and lysine (K) except followed by proline(P)
 * xxxRxx --> xxxR and xx
 * xxxKxx --> xxxK and xx
 * xxxRPx --> xxxRPx
 *
 * COMPILER OPTION:
 * gcc -opim pim.c -O3 -pedantic -Wall
 *
 * EXAMPLE:
 dyn.load("parentIonMass.so")
 r<-.C("computeParentIonMass", "TEST", n=as.integer(1), pim_=as.double(1))
#include <stdio.h>
 *

 $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/computeParentIonMass.c $
 $Id: computeParentIonMass.c 6179 2014-02-27 09:34:04Z cpanse $


 */

static double M[26];
static double C_term;
static double N_term;
static double Electron;
static double Hydrogen;
static double pepmass0;


int init()
{
    /* taken from mascot server
     * A=71.037114
     * B=114.534940
     * C=160.030649
     * D=115.026943
     * E=129.042593
     * F=147.068414
     * G=57.021464
     * H=137.058912
     * I=113.084064
     * J=0.000000
     * K=128.094963
     * L=113.084064
     * M=131.040485
     * N=114.042927
     * O=0.000000
     * P=97.052764
     * Q=128.058578
     * R=156.101111
     * S=87.032028
     * T=101.047679
     * U=150.953630
     * V=99.068414
     * W=186.079313
     * X=111.000000
     * Y=163.063329
     * Z=128.550590
     * Hydrogen=1.007825
     * Carbon=12.000000
     * Nitrogen=14.003074
     * Oxygen=15.994915
     * Electron=0.000549
     * C_term=17.002740
     * N_term=1.007825
     * delta1=0.984009,Deamidated (NQ)
     */
    M[0] = 71.037110;
    M[1] = 114.534940;
    M[2] = 160.030649;
    M[3] = 115.026940;
    M[4] = 129.042590;
    M[5] = 147.068410;
    M[6] = 57.021460;
    M[7] = 137.058910;
    M[8] = 113.084060;
    M[9] = 0.000000;
    M[10] = 128.094960;
    M[11] = 113.084060;
    M[12] = 131.040480;
    M[13] = 114.042930;
    M[14] = 0.000000;
    M[15] = 97.052760;
    M[16] = 128.058580;
    M[17] = 156.101110;
    M[18] = 87.032030;
    M[19] = 101.047680;
    M[20] = 150.953630;
    M[21] = 99.068410;
    M[22] = 186.079310;
    M[23] = 111.000000;
    M[24] = 163.063330;
    M[25] = 128.550590;

    C_term = 17.002740;
    N_term = 1.007825;
    Electron = 0.000549;
    Hydrogen = 1.007825;

    pepmass0 = C_term + N_term + Hydrogen - Electron;

    return 0;
}

/*
void computeParentIonMass(char *input, int *n, double *pim_)

fetuinPeptides<-c('LCPGR', 'GSVIQK', 'QYGFCK', 'AHYDLR', 'EVVDPTK', 'CNLLAEK', 'ALGGEDVR', 'HTLNQIDSVK', 'QDGQFSVLFTK', 'CDSSPDSAEDVR', 'TPIVGQPSIPGGPVR', 'LCPDCPLLAPLNDSR', 'QQTQHAVEGDCDIHVLK', 'HTFSGVASVESSSGEAFHVGK', 'EPACDDPDTEQAALAAVDYINK', 'AQFVPLPVSVSVEFAVAATDCIAK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK', 'VVHAVEVALATFNAESNGSYLQLVEISR', 'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR', 'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR')

dyn.load("parentIonMass.so")
source("computeParentIonMass.R")
computeParentIonMass(fetuinPeptides)


*/
void computeParentIonMass(int *n_, char **seq_,  double *pim_)
{

    int i;
    int j;
    int letter;
    init();

    /* for each sequence */
    for (i=0; i < n_[0]; i++){ 
        pim_[i] = pepmass0;
        
        /* for each letter */
        for (j = 0; seq_[i][j] != '\0'; j++){
            letter=seq_[i][j];
	        if (64 < letter && letter < 92) {
		        pim_[i] += M[letter - 65];
            }
        }
    }
}

void computeParentIonMass2(int *n_, char **seq_,  double *pim_, double *M_, double *N_term_)
{

    int i;
    int j;
    int letter;
    double C_term, Electron, Hydrogen;

    C_term = 17.002740;
    Electron = 0.000549;
    Hydrogen = 1.007825;

    /* for each sequence */
    for (i=0; i < n_[0]; i++){ 
        pim_[i] = C_term + N_term_[0] + Hydrogen - Electron;
        
        /* for each letter */
        for (j = 0; seq_[i][j] != '\0'; j++){
            letter=seq_[i][j];
	        if (64 < letter && letter < 92) {
		        pim_[i] += M_[letter - 65];
            }
        }
    }
}
