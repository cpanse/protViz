#include <stdio.h>
#include <float.h>

/*
 *  retrival.c
 *  deltatectra
 *
 * Copyright 2007, 2008, 2009, 2010
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


/*
    NOTES by <cp@fgcz.ethz.ch>
    
    Q: Why did we reinvent the wheel? Why do we not use bsearch instead?
    A:
    - ansi-c stdlib bsearch does an exact match otherwise it returns NULL.
    - using core R findInterval gives by definition only positive values.
    - NNQuery returns the index so that there is no other index with a closer distance. 


    $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/findNN.c $
    $Id: findNN.c 6179 2014-02-27 09:34:04Z cpanse $



*/
int NNQuery(void *key, void *base, size_t n, size_t size,
	    double (*dst) (const void *, const void *))
{
    int mid = -1;
    double d = -1;
    int low = 0;
    int high = n;
    double minDist = DBL_MAX;
    int minIdx = -1;
    char *cbase = base;
    char *ckey = key;

    if (high == 0)
       return (-1);

    while (low <= high) {
	    mid = ((low + high) / 2);
        if (mid >= n) break;

	    d = (*dst) (&cbase[mid * size], ckey);

	    if (d < 0) {
	        d = -d;
	    }
	    if (d < minDist) {
	        minDist = d;
	        minIdx = mid;
	    }

	    d = (*dst) (&cbase[mid * size], ckey);
	    if (d < 0) {
	        high = mid - 1;
	    } else if (d > 0) {
	        low = mid + 1;
	    } else {
	        return (minIdx);
	    }
    }

    return (minIdx);
}

double distd (const void *arg_a, const void *arg_b) {
    double *b = (double *) arg_b;
    double *a = (double *) arg_a;

    return (double) (*b - *a);
}



void __findNN (int *m_, int *n_, double *q_,  double *vec_, int *NN_)
{
    int i;

    for (i = 0; i < m_[0]; i++)
        NN_[i] = NNQuery(&q_[i], vec_, n_[0], sizeof(vec_[0]), distd);
}
