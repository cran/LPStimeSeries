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

void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}

void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}

void compute_similarity(int *testterminal, int *noftest, int *trainterminal, int *noftrain, int *nofterminal, int *result){
    int i,k,j;
    zeroInt(result, (*noftest) * (*noftrain));
    
    for(j=0;j<*noftest;j++){	
        for(i=0;i<*noftrain;i++){
			for(k=0;k<*nofterminal;k++){
				result[j+(*noftest)*i]=result[j+(*noftest)*i]+abs(trainterminal[*noftrain*k+i]-testterminal[*noftest*k+j]);
			}
		}
   }
}
    /*************************************
     * Codes to use for future development
     *************************************/
  
/*   
void minusDouble(double *x, int length) {
	int i;
	for(i=0;i<length;i++)
		x[i]=-999;
}
      
void generate_codebook(int *nodestatus, int *nofnode, int *noftree, int *terminal,  
	int *nofterminal, int *nofobservations, int *total, int *nofseries, int *result){
     int i,j,k,index,temp,ind;
     double tmp;
     for(i=0;i<*noftree;i++){
		temp=0;
		for(k=0;k<*nofseries;k++){
			for(j=0;j<*nofterminal;j++){
				result[(*nofseries)*(i*(*nofterminal)+j)+k]=0;
			}
			for(j=0;j<nofobservations[k];j++){
				ind=terminal[(*total*i)+temp+j]-1;
				index=nodestatus[i*(*nofnode)+ind];
				result[(*nofseries)*(i*(*nofterminal)+index)+k]=result[(*nofseries)*(i*(*nofterminal)+index)+k]+1;
			}			
			temp=temp+nofobservations[k];
		}
     }
}

void NNclassify_series(int *testterminal, int *noftest, int *trainterminal, int *noftrain, int *nofterminal, int *result){
    int i,k,j,bestsofar,bestid,similarity;
    
    for(j=0;j<*noftest;j++){	
		bestsofar=1000000;
        for(i=0;i<*noftrain;i++){
			similarity=0;
			for(k=0;k<*nofterminal;k++){
				similarity=similarity+abs(trainterminal[*noftrain*k+i]-testterminal[*noftest*k+j]);			
				if(similarity>bestsofar)
					break;
			}
			if(similarity<bestsofar){
				bestsofar=similarity;
				bestid=i+1;
			}
		}
		result[j]=bestid;
   }
}

void generate_segment(double *series, int *nofseries, int *intlen, int *nofint, double *result){
    int i,k,j,cnt;
    cnt=0;
    for(j=0;j<*nofseries;j++){	
        for(i=0;i<*nofint;i++){
			for(k=0;k<*intlen;k++){
				result[cnt++]=series[j+(*nofseries)*(i+k)];
			}
		}
   }
}

double findmean(double *input, int id, int len, int start, int end){
    int k,stin;
    double avy;
    
    stin=id*len;
    avy=0;
    for(k=(start+stin);k<(end+stin);k++){                        
       avy=avy+input[k];                         
    }
    avy=avy/(end-start);            
    return avy;   
}

double findvariance(double *input, int id, int len, double mean, int start, int end){
    int k,stin;
    double avy;
    
    stin=id*len;
    avy=0;
    for(k=(start+stin);k<(end+stin);k++){                         
       avy=avy+pow((input[k]-mean),2);                         
    }
    avy=avy/(end-start-1);
            
    return avy;   
}

double findslope(double *input, int id, int len, int start, int end){
    int k,stin;
    double sx,sxy,avx,avy;
    
    stin=id*len;
    avx=0;
    avy=0;
    for(k=(start+stin);k<(end+stin);k++){    
       avx=avx+k;                         
       avy=avy+input[k];                         
    }
    avx=avx/(end-start);
    avy=avy/(end-start);
    
    sx=0;
    sxy=0;
    for(k=(start+stin);k<(end+stin);k++){  
       sx=sx+pow((k-avx),2);
       sxy=sxy+((k-avx)*(input[k]-avy));                        
    }

    return sxy/sx;    
}

void generate_interval_features(double *series, int *nofseries, int *lenseries, 
									int *intlen, int *slidestep, int *ismean, int *isvar, int *isslope, double *result){
	int st,cnt,i;
	double mean;
	
	for(i=0;i<*nofseries;i++){
		st=0;
		cnt=0;
		while(*intlen+st<*lenseries){	    
			if(*isslope>0)
				result[*nofseries*(cnt++)+i]=findslope(series,i,*lenseries,st,st+*intlen);
			if(*ismean>0){
				mean=findmean(series,i,*lenseries,st,st+*intlen);
				result[*nofseries*(cnt++)+i]=mean;
			}
			if(*isvar>0){
				if(*ismean==0) mean=findmean(series,i,*lenseries,st,st+*intlen);
				result[*nofseries*(cnt++)+i]=findvariance(series,i,*lenseries,mean,st,st+*intlen);  
			}
			st=st+(*slidestep); 
		}
		
		if(st<*lenseries){	
			if(*isslope>0)
				result[*nofseries*(cnt++)+i]=findslope(series,i,*lenseries,st,*lenseries);
			if(*ismean>0){
				mean=findmean(series,i,*lenseries,st,*lenseries);
				result[*nofseries*(cnt++)+i]=mean;
			}
			if(*isvar>0){
				if(*ismean==0) mean=findmean(series,i,*lenseries,st,*lenseries);
				result[*nofseries*(cnt++)+i]=findvariance(series,i,*lenseries,mean,st,*lenseries); 
			} 
		}
	}
}

void compute_shapelet_distance(double *series, int *nseries, int *lenseries, 
						int *nshapelet, int *isRandom, int *impOrder, 
						int *intLen, int *slideLen, int *intCount,
						int *selected, int *level, double *result)
{

	int i, j, k, m, l, curCnt, curTS, curInterval, st, en, nofbreaks, ShStart, ShEnd, ShLen, cnt;
	double *shapelet, minsim, tempsim;
	
	shapelet = (double *) Calloc(*lenseries, double);	
	
	for(i=0;i<*nshapelet;i++){
		ShStart=10000;
		ShEnd=-10000;

		curTS = (int) (unif_rand()* (*nseries));
		curCnt = (int) (unif_rand()* (*intCount))+1;
		selected[i]=curTS;
		level[i]=curCnt;		
		minusDouble(shapelet,*lenseries);

		for(j=0;j<curCnt;j++){
			curInterval=impOrder[*nseries*j+curTS];	
			st=(curInterval-1)*(*slideLen);
			if(st+(*intLen)>=*lenseries){
				en=*lenseries;
			} else {
				en=st+(*intLen);			
			}
			if(st<ShStart) ShStart=st;
			if(en>ShEnd) ShEnd=en;
						
			for(m=st;m<en;m++)
				shapelet[m]=series[*lenseries*curTS+m];
		}
		ShLen=ShEnd-ShStart;

		for(k=0;k<*nseries;k++){
			minsim=1000;
			for(l=0;l<(*lenseries-ShLen);l++){
				tempsim=0;
				cnt=0;
				for(m=ShStart;m<ShEnd;m++){
					if(shapelet[m]!=-999){
					   tempsim=tempsim+pow(shapelet[m]-series[*lenseries*k+l+cnt],2);   				   
					   if(tempsim>=minsim){
						   break;
					   }
					}
					cnt++;
				}
				if(tempsim<minsim)
					minsim=tempsim;				
			}
			result[*nseries*i+k]=minsim;
		}		
	}
	
	Free(shapelet);
}

void compute_shapelet_distance_test(double *series, double *refseries, int *nseries, int *nrefseries,  int *lenseries, 
						int *nshapelet, int *impOrder, 
						int *intLen, int *slideLen, int *selected, int *level, double *result)
{

	int i, j, k, m, l, curCnt, curTS, curInterval, st, en, nofbreaks, ShStart, ShEnd, ShLen, cnt;
	double *shapelet, minsim, tempsim;
	
	shapelet = (double *) Calloc(*lenseries, double);	
	
	for(i=0;i<*nshapelet;i++){
		ShStart=10000;
		ShEnd=-10000;
		minusDouble(shapelet,*lenseries);	
		
		curTS = selected[i];
		curCnt = level[i];			
		for(j=0;j<curCnt;j++){
			curInterval=impOrder[*nrefseries*j+curTS];	
			st=(curInterval-1)*(*slideLen);
			if(st+(*intLen)>=*lenseries){
				en=*lenseries;
			} else {
				en=st+(*intLen);			
			}
			if(st<ShStart) ShStart=st;
			if(en>ShEnd) ShEnd=en;
						
			for(m=st;m<en;m++)
				shapelet[m]=refseries[*lenseries*curTS+m];
		}
		ShLen=ShEnd-ShStart;
		
		for(k=0;k<*nseries;k++){
			minsim=1000;
			for(l=0;l<(*lenseries-ShLen);l++){
				tempsim=0;
				cnt=0;
				for(m=ShStart;m<ShEnd;m++){
					if(shapelet[m]!=-999){
					   tempsim=tempsim+pow(shapelet[m]-series[*lenseries*k+l+cnt],2);
					   if(tempsim>=minsim){
						   break;
					   }  
					}
					cnt++;
				}
				if(tempsim<minsim) 
					minsim=tempsim;				
			}
			result[*nseries*i+k]=minsim;
		}		
	}
	
	Free(shapelet);
}
  
void createClass(double *x, int realN, int totalN, int mdim) {
 // Create the second class by bootstrapping each variable independently. 
    int i, j, k;
    for (i = realN; i < totalN; ++i) {
        for (j = 0; j < mdim; ++j) {
            k = (int) (unif_rand() * realN);
            x[j + i * mdim] = x[j + k * mdim];
        }
    }
}

void normClassWt(int *cl, const int nsample, const int nclass,
                 const int useWt, double *classwt, int *classFreq) {
    int i;
    double sumwt = 0.0;

    if (useWt) {
        // Normalize user-supplied weights so they sum to one. 
        for (i = 0; i < nclass; ++i) sumwt += classwt[i];
        for (i = 0; i < nclass; ++i) classwt[i] /= sumwt;
    } else {
        for (i = 0; i < nclass; ++i) {
            classwt[i] = ((double) classFreq[i]) / nsample;
        }
    }
    for (i = 0; i < nclass; ++i) {
        classwt[i] = classFreq[i] ? classwt[i] * nsample / classFreq[i] : 0.0;
    }
}

void makeA(double *x, const int mdim, const int nsample, int *cat, int *a,
           int *b) {
    // makeA() constructs the mdim by nsample integer array a.  For each
    // numerical variable with values x(m, n), n=1, ...,nsample, the x-values
    // are sorted from lowest to highest.  Denote these by xs(m, n).  Then
    // a(m,n) is the case number in which xs(m, n) occurs. The b matrix is
    // also contructed here.  If the mth variable is categorical, then
    // a(m, n) is the category of the nth case number. 
    int i, j, n1, n2, *index;
    double *v;

    v     = (double *) Calloc(nsample, double);
    index = (int *) Calloc(nsample, int);

    for (i = 0; i < mdim; ++i) {
        if (cat[i] == 1) { // numerical predictor 
            for (j = 0; j < nsample; ++j) {
                v[j] = x[i + j * mdim];
                index[j] = j + 1;
            }
            R_qsort_I(v, index, 1, nsample);

           //  this sorts the v(n) in ascending order. index(n) is the case
           //     number of that v(n) nth from the lowest (assume the original
          //      case numbers are 1,2,...).  
            for (j = 0; j < nsample-1; ++j) {
                n1 = index[j];
                n2 = index[j + 1];
                a[i + j * mdim] = n1;
                if (j == 0) b[i + (n1-1) * mdim] = 1;
                b[i + (n2-1) * mdim] =  (v[j] < v[j + 1]) ?
                    b[i + (n1-1) * mdim] + 1 : b[i + (n1-1) * mdim];
            }
            a[i + (nsample-1) * mdim] = index[nsample-1];
        } else { // categorical predictor 
            for (j = 0; j < nsample; ++j)
                a[i + j*mdim] = (int) x[i + j * mdim];
        }
    }
    Free(index);
    Free(v);
}


void modA(int *a, int *nuse, const int nsample, const int mdim,
	  int *cat, const int maxcat, int *ncase, int *jin) {
    int i, j, k, m, nt;

    *nuse = 0;
    for (i = 0; i < nsample; ++i) if (jin[i]) (*nuse)++;

    for (i = 0; i < mdim; ++i) {
      k = 0;
      nt = 0;
      if (cat[i] == 1) {
          for (j = 0; j < nsample; ++j) {
              if (jin[a[i + k * mdim] - 1]) {
                  a[i + nt * mdim] = a[i + k * mdim];
                  k++;
              } else {
                  for (m = 0; m < nsample - k; ++m) {
                      if (jin[a[i + (k + m) * mdim] - 1]) {
                          a[i + nt * mdim] = a[i + (k + m) * mdim];
                          k += m + 1;
                          break;
                      }
                  }
              }
              nt++;
              if (nt >= *nuse) break;
          }
      }
    }
    if (maxcat > 1) {
        k = 0;
        nt = 0;
        for (i = 0; i < nsample; ++i) {
            if (jin[k]) {
                k++;
                ncase[nt] = k;
            } else {
                for (j = 0; j < nsample - k; ++j) {
                    if (jin[k + j]) {
                        ncase[nt] = k + j + 1;
                        k += j + 1;
                        break;
                    }
                }
            }
            nt++;
            if (nt >= *nuse) break;
        }
    }
}

void Xtranslate(double *x, int mdim, int nrnodes, int nsample,
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize) {

 // this subroutine takes the splits on numerical variables and translates them
 //back into x-values.  It also unpacks each categorical split into a
 //32-dimensional vector with components of zero or one--a one indicates that
 //the corresponding category goes left in the split.


    int i, m;

    for (i = 0; i < treeSize; ++i) {
	if (nodestatus[i] == 1) {
	    m = bestvar[i] - 1;
	    if (cat[m] == 1) {
		xbestsplit[i] = 0.5 * (x[m + (bestsplit[i] - 1) * mdim] +
				       x[m + (bestsplitnext[i] - 1) * mdim]);
	    } else {
		xbestsplit[i] = (double) bestsplit[i];
	    }
	}
    }
}

void permuteOOB(int m, double *x, int *in, int nsample, int mdim) {
 // Permute the OOB part of a variable in x.
 // Argument:
 //   m: the variable to be permuted
 //   x: the data matrix (variables in rows)
 //  in: vector indicating which case is OOB
 //   nsample: number of cases in the data
 //    mdim: number of variables in the data
 
    double *tp, tmp;
    int i, last, k, nOOB = 0;

    tp = (double *) Calloc(nsample, double);

    for (i = 0; i < nsample; ++i) {
		// make a copy of the OOB part of the data into tp (for permuting) 
		if (in[i] == 0) {
            tp[nOOB] = x[m + i*mdim];
            nOOB++;
        }
    }

    last = nOOB;
    for (i = 0; i < nOOB; ++i) {
		k = (int) last * unif_rand();
		tmp = tp[last - 1];
		tp[last - 1] = tp[k];
		tp[k] = tmp;
		last--;
    }

    // Copy the permuted OOB data back into x.
    nOOB = 0;
    for (i = 0; i < nsample; ++i) {
		if (in[i] == 0) {
            x[m + i*mdim] = tp[nOOB];
            nOOB++;
		}
    }
    Free(tp);
}

//Compute proximity. 
void computeProximity(double *prox, int oobprox, int *node, int *inbag,
                      int *oobpair, int n) {
 //  Accumulate the number of times a pair of points fall in the same node.
 //  prox:    n x n proximity matrix
 // oobprox: should the accumulation only count OOB cases? (0=no, 1=yes)
 //  node:    vector of terminal node labels
 //  inbag:   indicator of whether a case is in-bag
 //  oobpair: matrix to accumulate the number of times a pair is OOB together
 //  n:       total number of cases

    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = i+1; j < n; ++j) {
            if (oobprox) {
                if (! (inbag[i] || inbag[j]) ) {
                    oobpair[j*n + i] ++;
                    oobpair[i*n + j] ++;
                    if (node[i] == node[j]) {
                        prox[j*n + i] += 1.0;
                        prox[i*n + j] += 1.0;
                    }
                }
            } else {
                if (node[i] == node[j]) {
                    prox[j*n + i] += 1.0;
                    prox[i*n + j] += 1.0;
                }
            }
        }
    }
}

unsigned int pack(int nBits, int *bits) {
    int i = nBits;
	unsigned int pack = 0;
    while (--i >= 0) pack += bits[i] << i;
    return(pack);
}

void unpack(int nBits, unsigned int pack, int *bits) {
// pack is a 4-byte integer.  The sub. returns icat, an integer 
//   array of zeroes and ones corresponding to the coefficients 
//   in the binary expansion of pack. 
    int i;
    for (i = 0; i < nBits; pack >>= 1, ++i) bits[i] = pack & 1;
}

void F77_NAME(unpack)(int *nBits, unsigned int *pack, int *bits) {
	unpack(*nBits, *pack, bits);
}

#ifdef OLD

double oldpack(int l, int *icat) {
    // icat is a binary integer with ones for categories going left
    // and zeroes for those going right.  The sub returns npack- the integer
    int k;
    double pack = 0.0;

    for (k = 0; k < l; ++k) {
	if (icat[k]) pack += R_pow_di(2.0, k);
    }
    return(pack);
}


void oldunpack(int l, int npack, int *icat) {
//
// npack is a long integer.  The sub. returns icat, an integer of zeroes and
// ones corresponding to the coefficients in the binary expansion of npack.
 
    int i;
    zeroInt(icat, 32);
    icat[0] = npack % 2;
    for (i = 1; i < l; ++i) {
	npack = (npack - icat[i-1]) / 2;
	icat[i] = npack % 2;
    }
}

#endif // OLD */
