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
#ifndef RF_H
#define RF_H

// swap two integers 
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b));
           
void regTree_time_series(double *x, double *segfactor, int targetdiff, int segmentdiff, int maxdepth, 
			int mdim, int nsample, int *lDaughter, int *rDaughter, double *upper, double *avnode, 
			int *nodedepth, int *nodestatus, int *splitType, int nrnodes, int *treeSize, int nthsize, 
			int mtry, int *mbest, int *currentTarget, int *currentTargetType, int *cat);
			
void predictRepresentation_time_series(double *x, int segmentlen, int nsample, int mdim, 
		int *lDaughter, int *rDaughter, int *nodedepth, int *nodestatus,
		double *split, int *splitVar, int *splitType, int *nodex, int maxdepth);
		
void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample, 
		   int ndstart, int ndend, int *msplit, double *decsplit, 
		   double *ubest, int *ndendl, int *jstat, int mtry,
		   double sumnode, int nodecnt, int *cat);

void zeroInt(int *x, int length);
void zeroDouble(double *x, int length);



/* 
// test if the bit at position pos is turned on 
#define isBitOn(x,pos) (((x) & (1 << (pos))) > 0)

void minusDouble(double *x, int length) ;
unsigned int pack(int l, int *icat);
void unpack(int nBits, unsigned int npack, int *icat);
void createClass(double *x, int realN, int totalN, int mdim);
void prepare(int *cl, const int nsample, const int nclass, const int ipi, 
	     double *pi, double *pid, int *nc, double *wtt);
void makeA(double *x, const int mdim, const int nsample, int *cat, int *a, 
           int *b);
void modA(int *a, int *nuse, const int nsample, const int mdim, int *cat, 
          const int maxcat, int *ncase, int *jin);
void Xtranslate(double *x, int mdim, int nrnodes, int nsample, 
		int *bestvar, int *bestsplit, int *bestsplitnext,
		double *xbestsplit, int *nodestatus, int *cat, int treeSize);
void permuteOOB(int m, double *x, int *in, int nsample, int mdim);
void computeProximity(double *prox, int oobprox, int *node, int *inbag, 
                      int *oobpair, int n); */

/* Template of Fortran subroutines to be called from the C wrapper */
extern void F77_NAME(buildtree)(int *a, int *b, int *cl, int *cat, 
				int *maxcat, int *mdim, int *nsample, 
				int *nclass, int *treemap, int *bestvar, 
				int *bestsplit, int *bestsplitnext, 
				double *tgini, int *nodestatus, int *nodepop, 
				int *nodestart, double *classpop, 
				double *tclasspop, double *tclasscat, 
				int *ta, int *nrnodes, int *, 
				int *, int *, int *, int *, int *, int *, 
				double *, double *, double *,
				int *, int *, int *); 

/* Node status */
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

/* Input Type */
#define OBS_SERIES 1
#define DIFF_SERIES 2

#endif /* RF_H */
