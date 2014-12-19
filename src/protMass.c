#include <stdio.h>

/*
 * Christian Panse <cp@fgcz.ethz.ch> 
 * Mit Nov  1 17:50:33 CET 2006
 *
 * DESCRIPTION:
 * The program reads AA on stdin and computes the weight.
 * It adds C and N term.
 * AA weight is defined in the code and can not be changed 
 * on the program execution.
 *
 * COMPILER OPTION:
 * gcc -oprotMass protMass.c -O3 -pedantic -Wall
 *
 * EXAMPLE:
 * cat SWISPPROT.fasta | fcat | protMass | head
 *
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

    /*
    M[2] = 160.030649;
     Cys   C3H5ONS    103.00919 103.1388
     */

    M[2] = 103.00919;

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

int main()
{

    int i;
    double protMass;

    init();

	protMass = 0.0;

    while ((i = getchar()) != EOF) {
	    if (64 < i && i < 92) {
         protMass += M[i - 65];
	    }else{
            printf("%f\n", protMass + C_term + N_term);
            protMass = 0.0;
        }

	}
    return 0;
}
