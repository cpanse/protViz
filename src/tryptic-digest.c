#include <stdio.h>


/*
 *  tryptic-digest.c
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

 $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/tryptic-digest.c $
 $Id: tryptic-digest.c 6179 2014-02-27 09:34:04Z cpanse $

	gcc -o tryptic-digest tryptic-digest.c -O3 -pedantic -Wall

*/

int trypticdigest ()
{

    int i;
    int j;

    j = 0;


    while ((i = getchar()) != EOF) {
	if (i != 'P' && j == 'R') {
	    putchar(j);
	    printf("\n");
	    j = i;
	} else if (j == 'K') {
	    putchar(j);
	    printf("\n");
	    j = i;
	} else if (j != 0) {
	    putchar(j);
	    j = i;
	} else {
	    j = i;
	}

    }

    printf("\n");

    return 0;
}

int main(int argc, char *argv[]) { 
return (trypticdigest());
}
