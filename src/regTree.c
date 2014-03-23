/*******************************************************************
   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/

/******************************************************************
 * regression tree is modifided for learning time series patterns
 ******************************************************************/
 
#include <Rmath.h>
#include <R.h>
#include "rf.h"

/*--------------------------------------------------------------*/
void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
		   int ndstart, int ndend, int *msplit, double *decsplit,
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat) {
    int last, lc, nl, nr, npopl, npopr;
    int i, j, kv, *mind, *ncase;
    double *xt, *ut, *v, *yl, avcat[32], tavcat[32], ubestt;
    double crit, critmax, critvar, suml, sumr, d, critParent;

    ut = (double *) Calloc(nsample, double);
    xt = (double *) Calloc(nsample, double);
    v  = (double *) Calloc(nsample, double);
    yl = (double *) Calloc(nsample, double);
    mind  = (int *) Calloc(mdim, int);
    ncase = (int *) Calloc(nsample, int);
    zeroDouble(avcat, 32);
    zeroDouble(tavcat, 32);

    /* START BIG LOOP */
    *msplit = -1;
    *decsplit = 0.0;
    critmax = 0.0;
    ubestt = 0.0;
    for (i=0; i < mdim; ++i) mind[i] = i;

    last = mdim - 1;
    for (i = 0; i < mtry; ++i) {
		critvar = 0.0;
		j = (int) (unif_rand() * (last+1));
		kv = mind[j];
        swapInt(mind[j], mind[last]);
		last--;
		
		lc = cat[kv];
		if (lc == 1) {
			/* numeric variable */
			for (j = ndstart; j <= ndend; ++j) {
				xt[j] = x[kv + (jdex[j] - 1) * mdim];
				yl[j] = y[jdex[j] - 1];
			}
		}
		
        /* copy the x data in this node. */
		for (j = ndstart; j <= ndend; ++j) v[j] = xt[j];
		for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
		R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
		if (v[ndstart] >= v[ndend]) continue;
		/* ncase(n)=case number of v nth from bottom */
		/* Start from the right and search to the left. */
		critParent = sumnode * sumnode / nodecnt;
		suml = 0.0;
		sumr = sumnode;
		npopl = 0;
		npopr = nodecnt;
		crit = 0.0;
		/* Search through the "gaps" in the x-variable. */
		for (j = ndstart; j <= ndend - 1; ++j) {
			d = yl[ncase[j] - 1];
			suml += d;
			sumr -= d;
			npopl++;
			npopr--;
			if (v[j] < v[j+1]) {
				crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
					critParent;
				if (crit > critvar) {
					ubestt = (v[j] + v[j+1]) / 2.0;
					critvar = crit;
				}
			}
		}
		if (critvar > critmax) {
			*ubest = ubestt;
			*msplit = kv + 1;
			critmax = critvar;
			for (j = ndstart; j <= ndend; ++j) {
				ut[j] = xt[j];
			}
			if (cat[kv] > 1) {
				for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
			}
		}
    }
    *decsplit = critmax;

    /* If best split can not be found, set to terminal node and return. */
    if (*msplit != -1) {
        nl = ndstart;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] <= *ubest) {
                nl++;
                ncase[nl-1] = jdex[j];
            }
        }
        *ndendl = imax2(nl - 1, ndstart);
        nr = *ndendl + 1;
        for (j = ndstart; j <= ndend; ++j) {
            if (ut[j] > *ubest) {
                if (nr >= nsample) break;
                nr++;
                ncase[nr - 1] = jdex[j];
            }
	    }
        if (*ndendl >= ndend) *ndendl = ndend - 1;
        for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];
		
    } else *jstat = 1;
	
    Free(ncase);
    Free(mind);
    Free(v);
    Free(yl);
    Free(xt);
    Free(ut);
}
                
void regTree_time_series(double *x, double *segfactor, int targetdiff, int segmentdiff, int maxdepth, 
			int mdim, int nsample, int *lDaughter, int *rDaughter, double *upper, double *avnode, 
			int *nodedepth, int *nodestatus, int *splitType, int nrnodes, int *treeSize, int nthsize, 
			int mtry, int *mbest, int *currentTarget, int *currentTargetType, int *cat) {
			 
    int i, j, k, m, ncur, *jdex, *nodestart, *nodepop, nofsampleobs, segmentlen, cur_target, cur_x;
    int ndstart, ndend, ndendl, nodecnt, jstat, msplit, cnt, totalnodes, sel, add, rotate;
    double d, ss, av, decsplit, ubest, sumnode, *y, *obsx, rndm;

	segmentlen = (int) (*segfactor * mdim);
	nofsampleobs = segmentlen * nsample;
	
	rotate=0; // experimental 1 if we choose any time point 
			  // as the starting point of the segment 

    nodestart = (int *) Calloc(nrnodes, int);
    nodepop   = (int *) Calloc(nrnodes, int);
     
    obsx = (double *) Calloc(nofsampleobs, double);
    y   = (double *) Calloc(nofsampleobs, double);
      
   // decide if target is difference series
   rndm=unif_rand();
   if(targetdiff==1 && rndm<0.5) {
	    if(rotate<1){
			cur_target=(int) (unif_rand()*(mdim-segmentlen)); //unequal observation probabilities
		} else {
			cur_target=(int) (unif_rand()*(mdim-1)); //unequal observation probabilities
		}
		*currentTargetType=DIFF_SERIES;
   } else {	   
	    if(rotate<1){
			cur_target=(int) (unif_rand()*(mdim-segmentlen+1));
		} else {
			cur_target=(int) (unif_rand()*mdim); //unequal observation probabilities
		}
		*currentTargetType=OBS_SERIES;
   }

   cnt=0;
   for (i = 0; i < nsample; i++) {
		for (j = 0; j < segmentlen; j++) {
			 if(*currentTargetType==OBS_SERIES){
				if(cur_target+j>mdim-1)
					y[cnt++]=x[cur_target+j-mdim+i*mdim];
				else
					y[cnt++]=x[cur_target+j+i*mdim];
						
			} else if(*currentTargetType==DIFF_SERIES){
				if(cur_target+j>mdim-2)
					y[cnt++]=x[cur_target+(j+2)-mdim+i*mdim]-x[cur_target+j+1-mdim+i*mdim];
				else
					y[cnt++]=x[cur_target+(j+1)+i*mdim]-x[cur_target+j+i*mdim];	
			}
		}
	}

    *currentTarget=cur_target+1;
    		
    /* initialize some arrays for the tree */
    zeroInt(nodedepth, nrnodes);
    zeroInt(nodestatus, nrnodes);
    zeroInt(nodestart, nrnodes);
    zeroInt(nodepop, nrnodes);
    zeroDouble(avnode, nrnodes);

    jdex = (int *) Calloc(nofsampleobs, int);
    for (i = 1; i <= nofsampleobs; ++i) jdex[i-1] = i;

    ncur = 0;
    nodestart[0] = 0;
    nodepop[0] = nofsampleobs;
    nodestatus[0] = NODE_TOSPLIT;

    /* compute mean and sum of squares for Y */
    av = 0.0;
    ss = 0.0;
    for (i = 0; i < nofsampleobs; ++i) {
		d = y[jdex[i] - 1];
		ss += i * (av - d) * (av - d) / (i + 1);
		av = (i * av + d) / (i + 1);
    }
    avnode[0] = av;

    /* start main loop */
    for (k = 0; k < nrnodes - 2; ++k) {	
		if (k > ncur || ncur >= nrnodes - 2) break;
		/* skip if the node is not to be split */
		if (nodestatus[k] != NODE_TOSPLIT) continue;

		if (nodedepth[k]==maxdepth && nodestatus[k] == NODE_TOSPLIT) {
			nodestatus[k] = NODE_TERMINAL;
			continue;
		} 
		
		/* initialize for next call to findbestsplit */
		ndstart = nodestart[k];
		ndend = ndstart + nodepop[k] - 1;
		nodecnt = nodepop[k];
		sumnode = nodecnt * avnode[k];
		jstat = 0;
		decsplit = 0.0;

	    // select predictor segment (should be different than target segment)
	    cur_x=cur_target;
	    rndm=unif_rand();
	    splitType[k]=*currentTargetType;
	    while(cur_x == cur_target && *currentTargetType==splitType[k]) {
			if(segmentdiff==1 && rndm<0.5) {
				if(rotate<1){
					cur_x=(int) (unif_rand()*(mdim-segmentlen)); //unequal observation probabilities
				} else {
					cur_x=(int) (unif_rand()*(mdim-1)); //unequal observation probabilities
				}
				splitType[k]=DIFF_SERIES;
			} else {	   
				if(rotate<1){
					cur_x=(int) (unif_rand()*(mdim-segmentlen+1)); //unequal observation probabilities
				} else {
					cur_x=(int) (unif_rand()*mdim); //unequal observation probabilities
				}
				splitType[k]=OBS_SERIES;
			}
		}

		cnt=0;
		for (i = 0; i < nsample; i++) {
			for (j = 0; j < segmentlen ; j++) {
				if(splitType[k]==OBS_SERIES){
					if(cur_x+j>mdim-1)
						obsx[cnt++]=x[cur_x+j-mdim+i*mdim];
					else
						obsx[cnt++]=x[cur_x+j+i*mdim];
				} else if(splitType[k]==DIFF_SERIES){
					if(cur_x+j>mdim-2)
						obsx[cnt++]=x[cur_x+(j+2)-mdim+i*mdim]-x[cur_x+j+1-mdim+i*mdim];
					else
						obsx[cnt++]=x[cur_x+(j+1)+i*mdim]-x[cur_x+j+i*mdim];			
				}
			}
		}

		findBestSplit(obsx, jdex, y, 1, nofsampleobs, ndstart, ndend, &msplit,
                      &decsplit, &ubest, &ndendl, &jstat, mtry, sumnode,
                      nodecnt, cat);
                 	
		if (jstat == 1) {
			/* Node is terminal: Mark it as such and move on to the next. */
			nodestatus[k] = NODE_TERMINAL;
			continue;
		}
		
        /* Found the best split. */
		mbest[k] = cur_x + 1;
		upper[k] = ubest;
		nodestatus[k] = NODE_INTERIOR;
		
		/* leftnode no.= ncur+1, rightnode no. = ncur+2. */
		nodepop[ncur + 1] = ndendl - ndstart + 1;
		nodepop[ncur + 2] = ndend - ndendl;
		nodestart[ncur + 1] = ndstart;
		nodestart[ncur + 2] = ndendl + 1;
		
		nodedepth[ncur + 1] = nodedepth[k] + 1;
		nodedepth[ncur + 2] = nodedepth[k] + 1;
		
		/* compute mean and sum of squares for the left daughter node */
		av = 0.0;
		ss = 0.0;
		for (j = ndstart; j <= ndendl; ++j) {
			d = y[jdex[j]-1];
			m = j - ndstart;
			ss += m * (av - d) * (av - d) / (m + 1);
			av = (m * av + d) / (m+1);
		}
		avnode[ncur+1] = av;
		nodestatus[ncur+1] = NODE_TOSPLIT;
		if (nodepop[ncur + 1] <= nthsize) {
			nodestatus[ncur + 1] = NODE_TERMINAL;
		}

		/* compute mean and sum of squares for the right daughter node */
		av = 0.0;
		ss = 0.0;
		for (j = ndendl + 1; j <= ndend; ++j) {
			d = y[jdex[j]-1];
			m = j - (ndendl + 1);
			ss += m * (av - d) * (av - d) / (m + 1);
			av = (m * av + d) / (m + 1);
		}
		avnode[ncur + 2] = av;
		nodestatus[ncur + 2] = NODE_TOSPLIT;
		if (nodepop[ncur + 2] <= nthsize) {
			nodestatus[ncur + 2] = NODE_TERMINAL;
		}
		
		/* map the daughter nodes */
		lDaughter[k] = ncur + 1 + 1;
		rDaughter[k] = ncur + 2 + 1;
		
		/* Augment the tree by two nodes. */			
		ncur += 2;
    }
    totalnodes = nrnodes;
    
    for (k = nrnodes - 1; k >= 0; --k) {
        if (nodestatus[k] == 0) totalnodes--;
        if (nodestatus[k] == NODE_TOSPLIT) {
            nodestatus[k] = NODE_TERMINAL;
        }
    }
    *treeSize=totalnodes;
    
    Free(nodestart);
    Free(jdex);
    Free(nodepop);
    Free(obsx);
    Free(y);
}
	
						
void predictRepresentation_time_series(double *x, int segmentlen, int nsample, int mdim, 
		int *lDaughter, int *rDaughter, int *nodedepth, int *nodestatus,
		double *split, int *splitVar, int *splitType, int *nodex, int maxdepth) {
						
    int i, j, k, m, rotate=0;
	
    for (i = 0; i < nsample; i++) {
		for (j = 0; j < segmentlen ; j++) {
			k = 0;
			while (nodestatus[k] != NODE_TERMINAL && nodedepth[k] < maxdepth) { /* go down the tree */
				m = splitVar[k] - 1;
				if(splitType[k]==OBS_SERIES){
					if(m+j>mdim-1){
						k = (x[m+j-mdim+i*mdim] <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;
					} else {
						k = (x[m+j+i*mdim] <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;						
					}
				}  else if(splitType[k]==DIFF_SERIES){
					if(m+j>mdim-2){
						k = ((x[m+(j+2)-mdim+i*mdim]-x[m+j-mdim+1+i*mdim]) <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;
					} else {
						k = ((x[m+(j+1)+i*mdim]-x[m+j+i*mdim]) <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;					
					}				
				}
			}
			nodex[k*nsample+i]++;
		}
	}
}
							
void predict_time_series(double *x, int segmentlen, int nsample, int mdim, 
		int *lDaughter, int *rDaughter, int *nodedepth, int *nodestatus,
		double *split, int *splitVar, int *splitType, double *nodepred, 
		int maxdepth, int target, double *prediction, int *targetcount) {
						
    int i, j, k, m, t;
    
	t = target - 1;
	for (j = 0; j < segmentlen ; j++) {
		if(t+j>mdim-1){
			targetcount[t+j-mdim]++;
		} else {
			targetcount[t+j]++;	
		}			
		for (i = 0; i < nsample; i++) {
			k = 0;
			while (nodestatus[k] != NODE_TERMINAL && nodedepth[k] < maxdepth) { /* go down the tree */
				m = splitVar[k] - 1;
				if(splitType[k]==OBS_SERIES){
					if(m+j>mdim-1){
						k = (x[m+j-mdim+i*mdim] <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;
					} else {
						k = (x[m+j+i*mdim] <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;						
					}
				}  else if(splitType[k]==DIFF_SERIES){
					if(m+j>mdim-2){
						k = ((x[m+(j+2)-mdim+i*mdim]-x[m+j-mdim+1+i*mdim]) <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;
					} else {
						k = ((x[m+(j+1)+i*mdim]-x[m+j+i*mdim]) <= split[k]) ?
						lDaughter[k] - 1 : rDaughter[k] - 1;					
					}				
				}
			}
			if(t+j>mdim-1){
				prediction[t+j-mdim+i*mdim] += nodepred[k];
			} else {
				prediction[t+j+i*mdim] += nodepred[k];
			}				
		}
	}
}


    /*************************************
     * Codes to use for future development
     *************************************/
     

/*
void predictRegTree_time_series(double *x, double segfactor, int nsample, int mdim,
		    int *lDaughter, int *rDaughter, int *nodestatus, double *split,
                    int *splitVar, int treeSize,
                    int *nodex, int max_depth) {
    int i, j, k, m,cnt,segmentlen,isdiff,cur_depth;
	unsigned int npack;
	cnt=0;
	
	segmentlen=(int) (segfactor*mdim);

    for (i = 0; i < nsample; i++) {
		for (j = 0; j < segmentlen ; j++) {
			k = 0;
			cur_depth=0;
			while (nodestatus[k] != NODE_TERMINAL && cur_depth<max_depth) { 
				m = splitVar[k] - 1;
				if(m>=0){
					k = (x[m+j+i*mdim]<=split[k]) ?
					lDaughter[k] - 1 : rDaughter[k] - 1;
				} else {
					k = ((x[-m+j+1+i*mdim]-x[-m+j+i*mdim]) <= split[k]) ?
					lDaughter[k] - 1 : rDaughter[k] - 1;					
				}
				cur_depth++;
			}
			nodex[cnt] = k + 1;
			cnt++;
		}
	}
}

void predictRegTree_time_series_predict(double *x, double segfactor, int nsample, int mdim,
		    int *lDaughter, int *rDaughter, int *nodestatus, double *split,
                    int *splitVar, double *predict, int treeSize,
                    double *avnodes, int max_depth) {
    int i, j, k, m,cnt,segmentlen,isdiff,cur_depth;
	unsigned int npack;
	cnt=0;
	
	segmentlen=(int) (segfactor*mdim);

    for (i = 0; i < nsample; i++) {
		for (j = 0; j < segmentlen ; j++) {
			k = 0;
			cur_depth=0;
			while (nodestatus[k] != NODE_TERMINAL && cur_depth<max_depth) { 
				m = splitVar[k] - 1;
				if(m>=0){
					k = (x[m+j+i*mdim]<= split[k]) ?
					lDaughter[k] - 1 : rDaughter[k] - 1;
				} else {
					k = ((x[-m+j+1+i*mdim]-x[-m+j+i*mdim]) <= split[k]) ?
					lDaughter[k] - 1 : rDaughter[k] - 1;					
				}
				cur_depth++;
			}
			avnodes[cnt]=predict[k];
			cnt++;
		}
	}
}
*/

