/*******************************************************************
   Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

#include <R.h>
#include "rf.h"

void regRF_time_series(double *x, double *seglength, int *tardiff, int *segdiff, 
		int *xdim, int *sampsize, int *nthsize, int *nrnodes, int *nTree, int *mtry, 
		int *cat, int *jprint, int *target, int *targetType, double *yptr, int *treeSize, 
		int *nodedepth, int *nodestatus, int *splitType, int *lDaughter, int *rDaughter, 
		double *avnode, int *mbest, double *upper, int *keepf, int *replace, int *inbag) {

    double  *xb, xrand;
    int *in, k, m, n, j, nsample, idx, mdim, keepF, keepInbag, maxdepth;

    nsample = xdim[0];    		//number of series
    mdim = xdim[1];	      		//length series
    keepF = keepf[0];     		//keep forest
    keepInbag = keepf[1]; 		//keep inbag data
    maxdepth=log2(*nrnodes);	//maximum depth
    
    GetRNGstate();

    /*************************************
     * Start the loop over trees.
     *************************************/
    for (j = 0; j < *nTree; ++j) {     
		idx = keepF ? j * *nrnodes : 0;		
		/* Draw a random sample for growing a tree. */
		if (*replace) { /* sample */		
			xb = (double *) S_alloc(mdim * *sampsize, sizeof(double));  //inbag x info
			in = (int *) S_alloc(nsample, sizeof(int));
			zeroInt(in, nsample);
			for (n = 0; n < *sampsize; ++n) {
				xrand = unif_rand();
				k = xrand * nsample;
				in[k] = 1;
				for(m = 0; m < mdim; ++m) {
					xb[m + n * mdim] = x[m + k * mdim];					
				}
			}
			if (keepInbag) {
				for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
			}
			/* grow the regression tree */
			regTree_time_series(xb, seglength[j], *tardiff, *segdiff, maxdepth, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, nthsize[j], *mtry, mbest + idx, target + j, targetType + j, cat);
		}  else { //use all data
			regTree_time_series(x, seglength[j], *tardiff, *segdiff, maxdepth, mdim, nsample, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, nthsize[j], *mtry, mbest + idx, target + j, targetType + j, cat);
		}                
        
        if(*jprint>0&&(j+1)%(*jprint)==0)
			Rprintf("Tree %d over\n",j+1);
    }
    PutRNGstate();
    /* end of tree iterations=======================================*/

}
		
void regForest_similarity(double *x, double *y, int *n, int *ny,
			double *seglength, int *mdim, int *ntree, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			int *treeSize, int *maxdepth, int *similarity) {
				   		                    
    int i, j, k, m, idx1, segmentlength;
    int *noderef, *nodetest, *tempnodestatus, totx, totxtst;
	
	totx=(*n)*(*nrnodes);
	totxtst=(*ny)*(*nrnodes);

    noderef = (int *) Calloc(totx, int);	
    nodetest = (int *) Calloc(totxtst, int);
    
    tempnodestatus = (int *) Calloc(*nrnodes, int);

    zeroInt(similarity, (*ny) * (*n));
	idx1 = 0;
	for (i = 0; i < *ntree; i++) {
		segmentlength=(int) (*mdim*seglength[i]);
		
		zeroInt(noderef, totx);
		zeroInt(nodetest, totxtst);
		zeroInt(tempnodestatus, *nrnodes);
		
		// based on the maxdepth setting identify terminal nodes
		for (k = 0; k < *nrnodes; k++) {
			if(nodedepth[idx1+k]==*maxdepth || nodestatus[idx1+k]==NODE_TERMINAL){
				tempnodestatus[k]=NODE_TERMINAL;
			}
		}
		
		predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
				nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1, noderef, *maxdepth);
		predictRepresentation_time_series(y, segmentlength, *ny, *mdim, lDaughter + idx1, rDaughter + idx1,
				nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1, nodetest, *maxdepth); 


		for (k = 0; k < *nrnodes; k++) {
			if(tempnodestatus[k]==NODE_TERMINAL){ 
				for(j = 0; j < (*ny); j++){	
					for(m = 0; m < (*n); m++){
						similarity[j+(*ny)*m]=similarity[j+(*ny)*m]+abs(noderef[(*n)*k+m]-nodetest[(*ny)*k+j]);
					}
				}
			}
	   }
	   idx1 += *nrnodes; 
	}	

	Free(noderef);
	Free(nodetest);
	Free(tempnodestatus);
}

void regForest_represent(double *x, int *n, int *whichtree,
			double *seglength, int *mdim, int *ntree, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			int *treeSize, int *maxdepth, int *representation, int *repLength) {
				
	int i, k, m, idx1, startind, endind,nodecount, termnodecount, tempidx;
    int *noderef, *tempnodestatus, segmentlength, totx;

	totx=(*n)*(*nrnodes);
	
	if(*whichtree>0){ //single tree
		startind=*whichtree-1;
		endind=startind+1;
		idx1=*nrnodes*(*whichtree-1);
	}
	else {	//ensemble
		startind=0;
		endind=*ntree;
		idx1 = 0;	
	}
	
	// compute the size of the representation  by computing 
	// the total number of terminal nodes (termnodecount) over trees
	tempidx = idx1;
	termnodecount = 0;
	for (i = startind; i < endind; i++) {
		for (k = 0; k < *nrnodes; k++) {
			if(nodedepth[tempidx+k]==*maxdepth || nodestatus[tempidx+k]==NODE_TERMINAL)
				termnodecount++;
		}	
		tempidx += *nrnodes; 
	}
	*repLength = termnodecount;
	
    noderef = (int *) Calloc(totx, int);	
    tempnodestatus = (int *) Calloc(*nrnodes, int);
    
	nodecount=0;
	for (i = startind; i < endind; i++) {
		segmentlength=(int) (*mdim*seglength[i]);
		zeroInt(noderef, totx);
		zeroInt(tempnodestatus, *nrnodes);
		
		// based on the maxdepth setting identify terminal nodes
		for (k = 0; k < *nrnodes; k++) {
			if(nodedepth[idx1+k]==*maxdepth || nodestatus[idx1+k]==NODE_TERMINAL){
				tempnodestatus[k]=NODE_TERMINAL;
			}
		}
		predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
				nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1,
				noderef, *maxdepth);
		

		for (k = 0; k < *nrnodes; k++) {
			if(tempnodestatus[k]==NODE_TERMINAL){ 
				for(m = 0; m < (*n); m++){
					representation[termnodecount*m+nodecount]=noderef[(*n)*k+m];
				}
				nodecount++;
			}
	   }
	   idx1 += *nrnodes; 
	}	

	Free(noderef);
	Free(tempnodestatus);				
}

/*************************************
 * Codes to use for future development
*************************************/
     
/*
void regForest_represent(double *x, double *seglength, int *mdim, int *n,
               int *ntree, int *lDaughter, int *rDaughter,
               int *nodestatus, int *nrnodes, int *nrterminalnodes, double *xsplit, 
               int *mbest, int *treeSize, int *max_depth, int *representation) {
				   		                    
    int i, j, k, m, idx1, cnt;
    int *noderef, *tempnodestatus, totx;

	totx=(*n)*(*nrnodes);
    noderef = (int *) Calloc(totx, int);	
    tempnodestatus = (int *) Calloc(*nrnodes, int);

	if(*ntree<0){
		zeroInt(representation, (*n));
		idx1 = 0;
		for (i = 0; i < -(*ntree); i++) {
			zeroInt(noderef, totx);
			zeroInt(tempnodestatus, *nrnodes);
			predictRepresentation_time_series(x, seglength[i], *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodestatus + idx1, xsplit + idx1, mbest + idx1, treeSize[i],
					noderef, *max_depth, tempnodestatus);		
			cnt=0;
			for(k=0;k<*nrnodes;k++){
				if(tempnodestatus[k]==NODE_TERMINAL){ 
					for(m=0;m<*n;m++){
						representation[cnt++]=noderef[(*n)*k+m];
					}
				}
			}
			idx1 += *nrnodes; 
		}	
	} else {
		zeroInt(representation, (*n));
		i=*ntree-1;
		idx1 = *nrnodes*i;
		zeroInt(noderef, totx);
		zeroInt(tempnodestatus, *nrnodes);
		predictRepresentation_time_series(x, seglength[i], *n, *mdim, lDaughter + idx1, rDaughter + idx1,
				nodestatus + idx1, xsplit + idx1, mbest + idx1, treeSize[i],
				noderef, *max_depth, tempnodestatus);		
		cnt=0;
		for(k=0;k<*nrnodes;k++){
			if(tempnodestatus[k]==NODE_TERMINAL){ 
				for(m=0;m<*n;m++){
					representation[cnt++]=noderef[(*n)*k+m];
				}
			}
		}
		idx1 += *nrnodes; 
	}

	Free(noderef);
	Free(tempnodestatus);
}


void regForest_time_series(double *x, double *seglength, int *mdim, int *n,
               int *ntree, int *lDaughter, int *rDaughter,
               int *nodestatus, int *nrnodes, double *xsplit, int *mbest, int *treeSize, int *max_depth, int *nodes, int *nodex) {
				   		                    
    int i, j, idx1, idx2, newmdim, noftree;
  
    if(*ntree<0){
		noftree=-(*ntree);
		zeroInt(nodex, *n * noftree * newmdim);
		idx1 = 0;
		idx2 = 0;
		for (i = 0; i < noftree; i++) {
			newmdim=(int) (seglength[i] * (*mdim));
			predictRegTree_time_series(x, seglength[i], *n, *mdim, lDaughter + idx1, rDaughter + idx1,
						   nodestatus + idx1, xsplit + idx1, mbest + idx1, treeSize[i],
						   nodex + idx2, *max_depth);
			idx1 += *nrnodes;
			if (*nodes) idx2 += (*n * newmdim);
		}
	} else {
		newmdim=(int) (seglength[*ntree-1] * (*mdim));
		zeroInt(nodex, *n * newmdim);
		idx1=(*ntree-1)* (*nrnodes);
		predictRegTree_time_series(x, seglength[*ntree-1], *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodestatus + idx1, xsplit + idx1, mbest + idx1, treeSize[(*ntree-1)],
					nodex, *max_depth);		
	}
}

void regForest_time_series_predict(double *x,double *segfactor, int *mdim, int *n,
               int *ntree, int *lDaughter, int *rDaughter,
               int *nodestatus, int *nrnodes, double *xsplit, int *mbest, 
               double *predicted, int *treeSize, int *max_depth, double *avnodes) {
				   		                    
    int i, j, idx1, idx2,newmdim, noftree;
    double seg;
    
    seg=*segfactor;
    newmdim=(int) (seg * (*mdim));

    if(*ntree<0){    
		noftree=-(*ntree);
		zeroDouble(avnodes, *n * noftree * newmdim);
		idx1 = 0;
		idx2 = 0;
		for (i = 0; i < noftree; i++) {
			predictRegTree_time_series_predict(x, seg, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
						   nodestatus + idx1, xsplit + idx1, mbest + idx1, predicted + idx1, treeSize[i], 
						   avnodes + idx2, *max_depth);
			idx1 += *nrnodes; // increment the offset 
			idx2 += (*n * newmdim);
		}
	} else {
		zeroDouble(avnodes, *n * newmdim);
		idx1=(*ntree-1)* (*nrnodes);
		predictRegTree_time_series_predict(x, seg, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodestatus + idx1, xsplit + idx1, mbest + idx1, predicted + idx1, treeSize[(*ntree-1)],
					avnodes, *max_depth);			
	}
    
} 

*/
