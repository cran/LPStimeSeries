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

void regRF_time_series(double *x, double *seglength, int *isRand, int *tardiff, int *segdiff, 
		int *xdim, int *sampsize, int *nthsize, int *nrnodes, int *nTree, int *mtry, 
		int *cat, int *jprint, int *oobpred, int *treeErrors, int *target, int *targetType, int *treeSize, 
		int *nodedepth, int *nodestatus, int *splitType, int *lDaughter, int *rDaughter, 
		double *avnode, int *mbest, double *upper, int *keepf, int *replace, double *oobpredictions,
		double *ooberrors, int *inbag, double *errorTree) {

    double  *xb, *predictions, xrand, temp;
    int *in, *targetcount, k, m, n, j, nsample, idx, mdim;
    int segmentlength, keepF, keepInbag, maxdepth, oobcount=0,oobcount2=0, oobcountperTree;

    nsample = xdim[0];    		//number of series
    mdim = xdim[1];	      		//length series
    keepF = keepf[0];     		//keep forest
    keepInbag = keepf[1]; 		//keep inbag data
    maxdepth=log2(*nrnodes+1);	//maximum depth
    
    GetRNGstate();
	
	if (*replace) { /* sample */		
		xb = (double *) Calloc(mdim*(*sampsize), double); //inbag x info
		in = (int *) Calloc(nsample, int);
		targetcount = (int *) Calloc(mdim*nsample, int);
		zeroInt(targetcount, mdim*nsample);
		predictions = (double *) Calloc(mdim*nsample, double);
	} else if (*treeErrors) {
		targetcount = (int *) Calloc(mdim, int);
		predictions = (double *) Calloc(mdim*nsample, double); 
	}
    /*************************************
     * Start the loop over trees.
     *************************************/
   
    for (j = 0; j < *nTree; ++j) {     
		idx = keepF ? j * *nrnodes : 0;		
		/* Draw a random sample for growing a tree. */
		if (*replace) { /* sample */		
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
			regTree_time_series(xb, seglength + j, *tardiff, *segdiff, maxdepth, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, target + j, targetType + j, cat, *isRand);
			
			if (*oobpred) {
			//	Rprintf("Yalan\n");
				segmentlength=(int) (mdim*seglength[j]);
				oobcount=0;    
				errorTree[j]=0;  
				ooberrors[j]=0;   
				zeroDouble(predictions, mdim*nsample);
				oobcount2=0; 
				for (n = 0; n < nsample; ++n) {
					if(in[n]==0){
						predict_time_series(x + n*mdim, segmentlength, 1, mdim, lDaughter + idx, rDaughter + idx,
							nodedepth + idx, nodestatus + idx, upper + idx, mbest + idx, splitType + idx,
							avnode + idx, maxdepth, target[j], predictions + n*mdim, targetcount + n*mdim, 0); 			
						 	 	
						for (m = 0; m < mdim; ++m) {
							if(m>=(target[j]-1)&&m<(target[j]+segmentlength-1)){		
								errorTree[j]=errorTree[j]+pow(x[m+n*mdim]-predictions[m+n*mdim],2);
								oobpredictions[m+n*mdim]=oobpredictions[m+n*mdim]+predictions[m+n*mdim];
								oobcount2++;
							}
						}	
					}
					for (m = 0; m < mdim; ++m) {
						if(targetcount[m+n*mdim]>0){
							temp=oobpredictions[m+n*mdim]/targetcount[m+n*mdim];
							ooberrors[j]=ooberrors[j]+pow(x[m+n*mdim]-temp,2);
							oobcount++;
						}
					}	
				}  
				ooberrors[j]=ooberrors[j]/oobcount;   
				errorTree[j]=errorTree[j]/oobcount2;
			}
					       
		}  else { //use all data
			
			regTree_time_series(x, seglength + j, *tardiff, *segdiff, maxdepth, mdim, nsample, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, target + j, targetType + j, cat, *isRand);
           
            if(*isRand==0&&*treeErrors) {
				zeroInt(targetcount, mdim);	
				zeroDouble(predictions, mdim*nsample);	
				segmentlength=(int) (mdim*seglength[j]);
				predict_time_series(x, segmentlength, nsample, mdim, lDaughter + idx, rDaughter + idx,
						nodedepth + idx, nodestatus + idx, upper + idx, mbest + idx, splitType + idx,
						avnode + idx, maxdepth, target[j], predictions, targetcount, 0);
						
				oobcount=0;  	 	
				for (m = 0; m < mdim; ++m) {
					if(targetcount[m]>0){
						for (n = 0; n < nsample; ++n) {		
							errorTree[j]=errorTree[j]+pow(x[m+n*mdim]-predictions[m+n*mdim],2);
							oobcount++;
						}
					}
				}	
				errorTree[j]=errorTree[j]/oobcount; 
			}

		}                
		
        if(*jprint>0&&(j+1)%(*jprint)==0)
			Rprintf("Tree %d over\n",j+1);
    }
   //  Rprintf("Tree construction Elapsed: %f seconds\n", temp / CLOCKS_PER_SEC);
    PutRNGstate();
    /* end of tree iterations=======================================*/

	if (*replace) { /* free memory */
	   if (*oobpred) {
		   for (n = 0; n < nsample; ++n) {
			   for (m = 0; m < mdim; ++m) {
					if(targetcount[m+n*mdim]>0){
						oobpredictions[m+n*mdim]=oobpredictions[m+n*mdim]/targetcount[m+n*mdim];
					} else {
						oobpredictions[m+n*mdim]=-999;
					}  
			   }		
		   }
	   }
       Free(xb);
       Free(in);
       Free(targetcount);
       Free(predictions);
	}
	else if (*treeErrors) {
	   Free(predictions);
       Free(targetcount);
	}
}
		
void regForest_similarity(double *x, double *y, int *n, int *ny,
			double *seglength, int *mdim, int *ntree, int *usedtrees, 
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
		if(usedtrees[i]==1) {
			segmentlength=(int) (*mdim*seglength[i]);
					
			zeroInt(noderef, totx);
			zeroInt(nodetest, totxtst);
			zeroInt(tempnodestatus, *nrnodes);
			
			// based on the maxdepth setting identify terminal nodes
			for (k = 0; k < *nrnodes; k++) {
				if(nodedepth[idx1+k]==*maxdepth) {
					tempnodestatus[k]=NODE_TERMINAL;
				}
				else if(nodestatus[idx1+k]==NODE_TERMINAL){
					tempnodestatus[k]=NODE_TERMINAL;
				}
			}
						
			predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, tempnodestatus, xsplit + idx1, mbest + idx1, splitType + idx1, noderef, *maxdepth);
			predictRepresentation_time_series(y, segmentlength, *ny, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, tempnodestatus, xsplit + idx1, mbest + idx1, splitType + idx1, nodetest, *maxdepth); 

			for (k = 0; k < *nrnodes; k++) {
				if(tempnodestatus[k]==NODE_TERMINAL){ 
					for(j = 0; j < (*ny); j++){	
						for(m = 0; m < (*n); m++){
							similarity[j+(*ny)*m]=similarity[j+(*ny)*m]+abs(noderef[(*n)*k+m]-nodetest[(*ny)*k+j]);
						}
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
			double *seglength, int *mdim, int *ntree, int *usedtrees, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			int *treeSize, int *maxdepth, int *representation, int *repLength) {
				
	int i, k, m, idx1, startind, endind,nodecount, termnodecount, tempidx;
    int *noderef, *tempnodestatus, segmentlength, totx;

	totx=(*n)*(*nrnodes);
	idx1 = 0;
	// compute the size of the representation  by computing 
	// the total number of terminal nodes (termnodecount) over trees
	tempidx = idx1;
	termnodecount = 0;
	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			for (k = 0; k < *nrnodes; k++) {
				if(nodedepth[tempidx+k]==*maxdepth) {
					termnodecount++;
				}
				else if(nodedepth[tempidx+k]<=*maxdepth&&nodestatus[tempidx+k]==NODE_TERMINAL){
					termnodecount++;
				}
			}	
		//	Rprintf("Tree %d, term node count %d\n",i+1,termnodecount);
		}
		tempidx += *nrnodes; 
	}
	*repLength = termnodecount;
	
	
    noderef = (int *) Calloc(totx, int);	
    tempnodestatus = (int *) Calloc(*nrnodes, int);
    
	nodecount=0;
	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			segmentlength=(int) (*mdim*seglength[i]);
			zeroInt(noderef, totx);
			zeroInt(tempnodestatus, *nrnodes);
			
			// based on the maxdepth setting identify terminal nodes
			for (k = 0; k < *nrnodes; k++) {
				if(nodedepth[idx1+k]==*maxdepth) {
					tempnodestatus[k]=NODE_TERMINAL;
				}
				else if(nodedepth[idx1+k]<=*maxdepth&&nodestatus[idx1+k]==NODE_TERMINAL){
					tempnodestatus[k]=NODE_TERMINAL;
				}
			}
			predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, tempnodestatus, xsplit + idx1, mbest + idx1, splitType + idx1,
					noderef, *maxdepth);
			

			for (k = 0; k < *nrnodes; k++) {
				if(tempnodestatus[k]==NODE_TERMINAL){ 
					for(m = 0; m < (*n); m++){
						representation[termnodecount*m+nodecount]=noderef[(*n)*k+m];
					}
					nodecount++;
				}
		   }
	   }
	   idx1 += *nrnodes; 
	}	

	Free(noderef);
	Free(tempnodestatus);				
}

void regForest_predict(double *x, int *n, int *whichtree,
			double *seglength, int *mdim, int *ntree, int *usedtrees, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			double *nodepred, int *treeSize, int *target, int *maxdepth, 
			double *prediction, int *targetcount) {
				
	int i, j, idx1, startind, endind;
    int segmentlength;


	idx1 = 0;	
	zeroDouble(prediction,(*n)*(*mdim));
	zeroInt(targetcount,(*mdim));
	
	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			segmentlength=(int) (*mdim*seglength[i]);
			predict_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1,
					nodepred + idx1, *maxdepth, target[i], prediction, targetcount, 1);
	   }		
	   idx1 += *nrnodes; 
	}	

	for (i = 0; i < *n; i++) {
		for (j = 0; j < *mdim; j++) {
			if(targetcount[j]>0){
				prediction[j+i*(*mdim)]=prediction[j+i*(*mdim)]/targetcount[j];
			} else {
				prediction[j+i*(*mdim)]=NA_indicator;
			}
		}
	}
		
}


void regForest_pattern(double *x, int *n, int *whichtree,
			int *whichterminal, double *seglength, int *mdim, int *ntree, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			int *treeSize, int *maxdepth, int *target, int *targetType, 
			double *predictpattern, double *targetpattern) {
				
	int i, j, k, m, idx1, startind, endind, nodecount, termnodecount;
    int segmentlength, totx, termid, curtree, lastk;

	curtree=*whichtree-1;
	idx1=*nrnodes*curtree;

	for (i = 0; i < *n; i++) {
		for (j = 0; j < *mdim; j++) {
			predictpattern[j+i*(*mdim)]=NA_indicator;
			targetpattern[j+i*(*mdim)]=NA_indicator;
		}
	}

	termnodecount = 0;
	for (k = 0; k < *nrnodes; k++) {
		if(nodedepth[idx1+k]==*maxdepth || nodestatus[idx1+k]==NODE_TERMINAL)
			termnodecount++;
				
		if(termnodecount==*whichterminal)
			break;
	}	
	termid=k;
	segmentlength=(int) (*mdim*seglength[curtree]);
	for (j = 0; j < segmentlength ; j++) {
		for (i = 0; i < *n; i++) {
			k = 0;
			while (nodestatus[idx1+k] != NODE_TERMINAL && nodedepth[idx1+k] < *maxdepth) {
				m = mbest[idx1+k] - 1;
				lastk = k;
				if(splitType[idx1+k]==OBS_SERIES){
					if(m+j>(*mdim)-1){
						k = (x[m+j-(*mdim)+i*(*mdim)] <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					} else {
						k = (x[m+j+i*(*mdim)] <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;						
					}
				}  else if(splitType[idx1+k]==DIFF_SERIES){
					if(m+j>(*mdim)-2){
						k = ((x[m+(j+2)-(*mdim)+i*(*mdim)]-x[m+j-(*mdim)+1+i*(*mdim)]) <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					} else {
						k = ((x[m+(j+1)+i*(*mdim)]-x[m+j+i*(*mdim)]) <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;					
					}				
				}		
			} 
			
			if(termid==k){
				m = target[curtree] - 1;
				if(m+j>(*mdim)-1){
					targetpattern[m+j-(*mdim)+i*(*mdim)] = x[m+j-(*mdim)+i*(*mdim)];
				} else {
					targetpattern[m+j+i*(*mdim)] = x[m+j+i*(*mdim)];			
				}

				m = mbest[idx1] - 1;
				if(m+j>(*mdim)-1){
					predictpattern[m+j-(*mdim)+i*(*mdim)] = x[m+j-(*mdim)+i*(*mdim)];
				} else {
					predictpattern[m+j+i*(*mdim)] = x[m+j+i*(*mdim)];			
				}
		 			
			}
		}
	}
		
}

/*************************************
 * Codes to use for future development
*************************************/
/*
void regForest_findNN(double *x, double *y, int *n, int *ny,
			double *seglength, int *mdim, int *ntree, 
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType, 
			int *treeSize, int *maxdepth, int *K, int *result) {
				   		                    
    int i, j, k, m, idx1, segmentlength;
    int *noderef, *nodetest, *tempnodestatus, dist, bestsofar;
	

    noderef = (int *) Calloc((*nrnodes), int);	
    nodetest = (int *) Calloc((*nrnodes), int);
    tempnodestatus = (int *) Calloc((*ntree)*(*nrnodes), int);

	zeroInt(tempnodestatus, (*ntree)*(*nrnodes));
		
	// based on the maxdepth setting identify terminal nodes
	for (i = 0; i < *ntree; i++) {
		idx1 = 0;
		for (k = 0; k < *nrnodes; k++) {
			if(nodedepth[idx1+k]==*maxdepth || nodestatus[idx1+k]==NODE_TERMINAL){
				tempnodestatus[idx1+k]=NODE_TERMINAL;
			}
		}
		idx1 += *nrnodes; 
	}
		
	for(j = 0; j < (*ny); j++){	
		bestsofar=30000;
		for(m = 0; m < (*n); m++){
			dist=0;
			for (i = 0; i < *ntree; i++) {
				idx1 = 0;
				segmentlength=(int) (*mdim*seglength[i]);			
				zeroInt(noderef, (*nrnodes));
				zeroInt(nodetest, (*nrnodes));
				predictRepresentation_time_series(x + m*(*mdim), segmentlength, 1, *mdim, lDaughter + idx1, rDaughter + idx1,
						nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1, noderef, *maxdepth);
				predictRepresentation_time_series(y+ j*(*mdim), segmentlength, 1, *mdim, lDaughter + idx1, rDaughter + idx1,
						nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1, nodetest, *maxdepth); 		
				for (k = 0; k < *nrnodes; k++) {
					if(tempnodestatus[idx1+k]==NODE_TERMINAL){ 
						dist=dist+abs(noderef[k]-nodetest[k]);
					}
				}
				idx1 += *nrnodes; 	
				if(bestsofar<dist)
					break;
			}

			if(dist<bestsofar){
				bestsofar=dist;
				result[j]=m+1;
			}
		}
	}	
	
	Free(noderef);
	Free(nodetest);
	Free(tempnodestatus);
}
* 
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
