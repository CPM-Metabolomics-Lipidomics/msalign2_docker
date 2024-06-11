/*                                                                                                 */                     
/*  msalign_uTOF - automatic alignment of 2 LC-MS datasets using feature extraction and GA         */
/*                                                                                                 */
/* Copyright (c) Magnus Palmblad 2007, 2008                                                        */ 
/*                                                                                                 */
/* This program is free software; you can redistribute it and/or modify it under the terms of the  */
/* GNU General Public License as published by the Free Software Foundation; either version 2 of    */
/* the License, or (at your option) any later version.                                             */
/*                                                                                                 */
/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;       */
/* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.       */
/* See the GNU General Public License for more details.                                            */
/*                                                                                                 */
/* You should have received a copy of the GNU General Public License along with this program; if   */
/* not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA       */
/* 02111-1307, USA.                                                                                */
/*                                                                                                 */
/* Contact information: magnus.palmblad@gmail.com                                                  */
/*                                                                                                 */
/*                                                                                                 */
/* compile with e.g. gcc -o msalign base64.c ramp.c msalign_uTOF.c -I. -lgd -lm -lz -std=gnu99     */
/*                                                                                                 */
/* TODO list (suggested improvements)                                                              */
/* 1) clean up code (change variable names, move code to main(), ...                               */
/* 2) interpolation to get better mass accuracy in MS (presently "max" peak picking)               */
/*                                                                                                 */
/*                                                                                                 */


#include <stdio.h>
#include <stdlib.h>  
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "ramp.h"
#include "base64.h"

#define HPLUS_MASS 1.00727646688

#define MAX_DATAPOINTS 50000
#define N_CANDIDATES 300
#define MAX_BREAKPOINTS 24      // original value was 12
#define N_GENERATIONS 10000    // original value was 10000 (before even 5000)
#define FRACTION_KEPT 0.5 
#define MAX_FEATURES 100000

/* mzXML read functions adopted from Pep3D */

typedef struct {
  int size;
  double * xval;
  double * yval;
} spectStrct;

void freeSpectStrct(spectStrct spectrum)
{
  free(spectrum.xval);
  free(spectrum.yval);

  return;
}

ramp_fileoffset_t *pScanIndex;
int iLastScan;
struct ScanHeaderStruct scanHeader;
struct RunHeaderStruct runHeader;
ramp_fileoffset_t indexOffset;

spectStrct *spects;
double *wghs;
int spectNum;
long i;
RAMPREAL *pPeaks;
int n, n_MS_spectra=0;
RAMPFILE *pFI;

int getMsSpect(spectStrct *msSpect, RAMPFILE *pFI, int scanNum[2])
{
  void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs);

  msSpect->size = -1;

  // initialize
  scanNum[0] = scanNum[0] > 1 ? scanNum[0] : 1;
  scanNum[1] = scanNum[1] < iLastScan ? scanNum[1] : iLastScan;
  spectNum = scanNum[1] - scanNum[0] + 1;
  if(spectNum < 1){
    printf("invalid scan number: %d-%d (full scan range: 1-%d)\n",scanNum[0], scanNum[1], iLastScan);
    fflush(stdout);
    free (pScanIndex);
    return -1;    
  }

  spects = (spectStrct *) calloc(spectNum, sizeof(spectStrct));
  spectNum = 0;
  for (i = scanNum[0]; i <= scanNum[1]; ++i) 
    {
      if((scanHeader.msLevel==1)&&(scanHeader.peaksCount>0)) /* MS ? */
	{                         
	  spects[spectNum].size = scanHeader.peaksCount;
	  spects[spectNum].xval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  spects[spectNum].yval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  
	  pPeaks = readPeaks (pFI, pScanIndex[i]);
	  
	  spects[spectNum].size = 0;
	  n = 0;
	  while (pPeaks[n] != -1)
	    {
	      spects[spectNum].xval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      spects[spectNum].yval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      ++(spects[spectNum].size);
	    }
	  free (pPeaks);
	  n_MS_spectra++;
	  if(spects[spectNum].size > 0) ++spectNum; 
	  else freeSpectStrct(spects[spectNum]);
	}
      
    } 
  
  if(spectNum > 0) {
    wghs = (double *) calloc(spectNum, sizeof(double));
    for (i = 0; i < spectNum; ++i)
      wghs[i] = 1.;
    getCmbSpect(msSpect, spectNum, spects, wghs);
    free(wghs);
  }
  else 
    {
      printf("cannot find an MS spectrum\n"); fflush(stdout);
      for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
      free(spects);
      return -1;
    }
  
  for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
  free(spects);
  
  return 0;
}


int getCidSpect(double *mz, double *et, spectStrct *cidSpect, RAMPFILE *pFI, int scanNum[2])
{
  void getCmbSpect(spectStrct *cmbSpect,int spectNum, spectStrct *spects, double *wghs);

  cidSpect->size = -1;

  scanNum[0] = scanNum[0] > 1 ? scanNum[0] : 1;
  scanNum[1] = scanNum[1] < iLastScan ? scanNum[1] : iLastScan;
  spectNum = scanNum[1] - scanNum[0] + 1;
  if(spectNum < 1){
    printf("invalid scan number: %d-%d (full scan range: 1-%d)\n",scanNum[0], scanNum[1], iLastScan);
    fflush(stdout);
    free (pScanIndex);
    return -1;    
  }
  spects = (spectStrct *) calloc(spectNum, sizeof(spectStrct));
  *mz = 0.;
  *et = 0.;
  spectNum = 0;
  for (i = scanNum[0]; i <= scanNum[1]; ++i) 
    {
      if((scanHeader.msLevel==2)&&(scanHeader.peaksCount>0)) 
	{                         
	  *mz += scanHeader.precursorMZ;
	  *et += scanHeader.retentionTime/60;
	  
	  spects[spectNum].size = scanHeader.peaksCount;
	  spects[spectNum].xval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  spects[spectNum].yval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  
	  pPeaks = readPeaks (pFI, pScanIndex[i]);
	  
	  spects[spectNum].size = 0;
	  n = 0;
	  while (pPeaks[n] != -1)
	    {
	      spects[spectNum].xval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      spects[spectNum].yval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      ++(spects[spectNum].size);
	    }
	  free (pPeaks);
	  
	  if(spects[spectNum].size > 0) ++spectNum;
	  else freeSpectStrct(spects[spectNum]);
	} 
      
    } 
  
  if(spectNum > 0) {
    *mz /= spectNum;
    *et /= spectNum;
    wghs = (double *) calloc(spectNum, sizeof(double));
    for (i = 0; i < spectNum; ++i)
      wghs[i] = 1.;
    getCmbSpect(cidSpect, spectNum, spects, wghs);
    free(wghs);
  }
  else 
    {
      printf("cannot find an MS/MS spectrum\n"); fflush(stdout); 
      for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
      free(spects);
      return -1;
    }
  
  for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
  free(spects);
  
  return 0;
}


void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs)
{
  void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect);
  spectStrct tmpSpect[2];
  int indx, indx1, indx2;
  double tmpWghs[2] = {1., 1.};
  int i;
  
  if (spectNum < 1)
    return;
  
  // single spectrum
  if(spectNum == 1) {
    copySpectStrct(cmbSpect, spects[0]);
    if(wghs[0] != 1.) {
      for(i = 0; i < cmbSpect->size; ++i)
	cmbSpect->yval[i] *= wghs[0];
    }
    return;
  } 
  
  // 2 spectra
  if(spectNum == 2) {
    tmpSpect[0].size = spects[0].size + spects[1].size;
    tmpSpect[0].xval = (double *) calloc(tmpSpect[0].size, sizeof(double));
    tmpSpect[0].yval = (double *) calloc(tmpSpect[0].size, sizeof(double));

    indx1 = 0;
    indx2 = 0;
    indx = 0;
    while(indx1 < spects[0].size || indx2 < spects[1].size) {
      
      if(indx1 >= spects[0].size){
	tmpSpect[0].xval[indx] = spects[1].xval[indx2];
	tmpSpect[0].yval[indx] = spects[1].yval[indx2]*wghs[1];      
	++indx2;
	++indx;
      }
      else if (indx2 >= spects[1].size) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0];      
	++indx1;
	++indx;
      }
      else if(spects[0].xval[indx1] == spects[1].xval[indx2]) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0] 
	  + spects[1].yval[indx2]*wghs[1];      
	++indx1;
	++indx2;
	++indx;
      }
      else if(spects[0].xval[indx1] < spects[1].xval[indx2]) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0];      
	++indx1;
	++indx;
      }
      else {
	tmpSpect[0].xval[indx] = spects[1].xval[indx2];
	tmpSpect[0].yval[indx] = spects[1].yval[indx2]*wghs[1];      
	++indx2;
	++indx;
      }
    } 
    tmpSpect[0].size = indx;
    
    copySpectStrct(cmbSpect, tmpSpect[0]);
    freeSpectStrct(tmpSpect[0]);

    return;
  }

  // at least three spectra
  indx1 = spectNum/2;
  indx2 = spectNum - spectNum/2;
  getCmbSpect(&(tmpSpect[0]), indx1, spects, wghs);
  getCmbSpect(&(tmpSpect[1]), indx2, spects+indx1, wghs+indx1);
  getCmbSpect(cmbSpect, 2, tmpSpect, tmpWghs);
  
  freeSpectStrct(tmpSpect[0]);
  freeSpectStrct(tmpSpect[1]);

  return;
}

void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect)
{
  int i;

  tgtSpect->size = srcSpect.size;
  tgtSpect->xval = (double *) calloc(tgtSpect->size, sizeof(double));
  tgtSpect->yval = (double *) calloc(tgtSpect->size, sizeof(double));

  for (i = 0; i < tgtSpect->size; ++i) {
    tgtSpect->xval[i] = srcSpect.xval[i];
    tgtSpect->yval[i] = srcSpect.yval[i];
  }

  return;
}




/* main program starts here */

int main(int argc, char *argv[]) 
{
  FILE *inp, *outp;

  typedef struct {
    double B[MAX_BREAKPOINTS][2];
    int nB;
  } PLF_type; /* piecewise linear function */
 
  PLF_type best_C;
  spectStrct mzXML_spectrum;
  RAMPFILE *mzXML_file;
  char pepXML_filename[100], mzXML_filename_1[100], mzXML_filename_2[100], output_filename[100], scan_range[100], line[4000], *p, temp[30], new_mass;
  ramp_fileoffset_t offset, *scan_index;
  struct ScanHeaderStruct scan_header;
  struct RunHeaderStruct run_header;
  int nographics, MS_start_scan, MS_end_scan, SetFeatures, CountFeatures, assumed_charge, n_masses, n_unique_masses, n_unique_masses1, n_unique_masses2, *SIC_1_argmax, *SIC_2_argmax, *SIC_2a_argmax, range[2], seen_before, LC_sigma_set, Xmax_set, Ymax_set, non_matched;
  long i, j, k, n;
  double max_mass_measurement_error, start_scan, end_scan, intensity_sum, *SIC_1_max, *SIC_2_max, *SIC_2a_max, A[MAX_DATAPOINTS][2], scan_diff[MAX_DATAPOINTS], scan_diff_median, intensity_max, scan_intensity_max, *SIC_1_max_mz, *SIC_2_max_mz, *SIC_2a_max_mz, error[MAX_DATAPOINTS], swap, sum_squared_deviation, LC_sigma, Xmax, Xmax_temp, Ymax, background, oldbackground1, background1, background2, oldbackground2, costs, argFitness;
  int go, r, SumFeatures, SetBackground, feature_index, CountAttempts; 
  
  double mimass[5]={1.0078250321,12,14.0030740052,15.9949146221,31.97207069};
  int fragment[5]; 

 
  double fitness(double A[MAX_DATAPOINTS][2], int nA, PLF_type S)
  {
    double f=0,k,m;
    int i,j,bp1,bp2;
      
    //printf("calculating fitness... nA=%d, nB=%d, fitness=",nA,S.nB); fflush(stdout);
    for(i=0;i<nA;i++)
      {
	j=0;
	while(S.B[j][1]<0) j=j+1;
	bp1=j; bp2=j+1;
	for(j=0;j<MAX_BREAKPOINTS-1;j++)
	  {
	    if((A[i][1]>S.B[j][1])&&(A[i][1]<S.B[j+1][1]))
	      { 
		bp1=j;
		bp2=j+1; 
	      }
	  }

	/* define y=k*x+m for line segment between breakpoints */
	if((S.B[bp2][1]-S.B[bp1][1])>0)
	  {
	    k=(S.B[bp2][2]-S.B[bp1][2])/(S.B[bp2][1]-S.B[bp1][1]);
	    m=S.B[bp1][2]-k*S.B[bp1][1];
	    if (A[i][2]>0) {f+=argFitness*exp(-((A[i][2]-(k*A[i][1]+m))*(A[i][2]-(k*A[i][1]+m)))/(LC_sigma*LC_sigma));} /* fitness (are we near point?) */
	  }
      }
    
    f=f-costs*S.nB; /* cost per breakpoint !!!original value was 0.5!!!! changed to costs variable*/
    // for(j=0;j<S.nB-1;j++) if(S.B[j][1]>0) f=f-0.1; /* cost per breakpoint */
    //printf("%f\n",f); fflush(stdout);
    return f;
  }

  PLF_type sort_breakpoints(PLF_type S) 
  {
    int i,j;
    PLF_type temp;
    
    for(i=0;i<MAX_BREAKPOINTS-1;i++)
      {
	for(j=0;j<MAX_BREAKPOINTS-1-i;j++)
	  if(S.B[j+1][1]<S.B[j][1]) 
	    {  
	      temp.B[1][1]=S.B[j][1];       
	      S.B[j][1]=S.B[j+1][1];
	      S.B[j+1][1]=temp.B[1][1];
	    }
      }

    for(i=0;i<MAX_BREAKPOINTS-1;i++)
      {
	for (j=0;j<MAX_BREAKPOINTS-1-i;j++)
	  if (S.B[j+1][2]<S.B[j][2]) 
	    {  
	      temp.B[1][2]=S.B[j][2];       
	      S.B[j][2]=S.B[j+1][2];
	      S.B[j+1][2]=temp.B[1][2];
	    }
      }
    S.nB=MAX_BREAKPOINTS;
    for (i=0;i<MAX_BREAKPOINTS-1;i++) if(S.B[i][1]<0) S.nB--; /* count number of breakpoints */
    
    return S;
  } 

   PLF_type GAalign(double A[MAX_DATAPOINTS][2],int nA)
   {
      PLF_type C[N_CANDIDATES], nextC[N_CANDIDATES], out_C;
      double F[N_CANDIDATES],temp; /* fitness vector */
      int F_index[N_CANDIDATES]; /* fitness vector index */
      int i,j,G,temp_index,parent;
     
      printf("aligning data.."); fflush(stdout);
      srand48(22069);
      for(i=0;i<N_CANDIDATES-1;i++)     // added -1 otherwise the program gets stuck in this loop ??? | 26-11-2008
	  {  
    	//if (!(i%50)) printf("GAalign = %d\n", i); fflush(stdout);    //used for debugging	  
	    C[i].B[MAX_BREAKPOINTS][1]=Xmax;
	    C[i].B[MAX_BREAKPOINTS][2]=Ymax;

	    for(j=0;j<MAX_BREAKPOINTS-1;j++)
	    {
	      C[i].B[j][1]=drand48()*Xmax-1;
	      C[i].B[j][2]=drand48()*Ymax-1;
	    }
	    C[i].nB=MAX_BREAKPOINTS; 
	    C[i]=sort_breakpoints(C[i]);
	  }

      for(G=0;G<=N_GENERATIONS;G++) /* loop through N_GENERATIONS generations */
	  {
	    if (!(G%100)) printf("."); fflush(stdout);
	    for(i=0;i<N_CANDIDATES;i++)
	    {
	      F[i]=fitness(A,nA,C[i]);
	      F_index[i]=i; 
	    }
	
	    /* sort candidates according to fitness */
	    for(i=0;i<N_CANDIDATES-1;i++)
	    {
	      for(j=0;j<N_CANDIDATES-1-i;j++)
		  {
		    if(F[j+1]>F[j]) 
		    {  
		      temp=F[j];       
		      F[j]=F[j+1];
		      F[j+1]=temp;
		      temp_index=F_index[j];       
		      F_index[j]=F_index[j+1];
		      F_index[j+1]=temp_index;
		    }
		  }
	    }
	  
	    if(G==N_GENERATIONS)
	    {
	      for(i=0;i<MAX_BREAKPOINTS;i++)
		  {
		    out_C.B[i][1]=C[F_index[0]].B[i][1];
		    out_C.B[i][2]=C[F_index[0]].B[i][2];
		    out_C.nB=C[F_index[0]].nB;
		  }
	      printf("done\n"); fflush(stdout);
	      return out_C; /* return after sorting and pulling out best candidate */
	    }
	
	    for(i=0;i<(int)floor((double)(FRACTION_KEPT*N_CANDIDATES+0.5));i++) /* keep FRACTION_KEPT best candidates */ //added typecasting and replaced round for floor and "+0.5" | 25-11-2008
	    {
	      nextC[i].nB=C[F_index[i]].nB;
	      for(j=0;j<MAX_BREAKPOINTS-1;j++)
		  {
		    nextC[i].B[j][1]=C[F_index[i]].B[j][1];
		    nextC[i].B[j][2]=C[F_index[i]].B[j][2];
		  }
	    }

	    for(i=(int)floor((double)(FRACTION_KEPT*N_CANDIDATES+0.5));i<N_CANDIDATES;i++)
	    { 
	      parent=(int)floor((double)(drand48()*FRACTION_KEPT*N_CANDIDATES+0.5)); /* choose fit parent */  //added typecasting and replaced round for floor and "+0.5" | 25-11-2008
	      nextC[i].nB=C[F_index[parent]].nB;
	      for(j=0;j<MAX_BREAKPOINTS-1;j++)
		  {
		    nextC[i].B[j][1]=C[F_index[parent]].B[j][1];
		    nextC[i].B[j][2]=C[F_index[parent]].B[j][2];
		  }
	    
	      for(j=0;j<MAX_BREAKPOINTS-1;j++) /* always leave top breakpoint? */
		  {
		    if(drand48()<0.1) /* major mutation/insertion */
		    {
		      nextC[i].B[j][1]=drand48()*Xmax-1; 
		      nextC[i].B[j][2]=drand48()*Ymax-1; 
		    }
		  
		    if(nextC[i].B[j][1]>0)
		    {
		      if(drand48()<0.5*(G/N_GENERATIONS))
			  {
			    nextC[i].B[j][1]+=10*(drand48()-0.5); /* nudge X */
			    nextC[i].B[j][2]+=3*(drand48()-0.5); /* nudge Y */
		      }

		      if(drand48()<0.1) /* delete breakpoint */
			  {
			    nextC[i].B[j][1]=-1;
			    nextC[i].B[j][2]=-1;
		      }
		    }
		  }
	      nextC[i]=sort_breakpoints(nextC[i]);
	    }
	 
	    for(i=0;i<N_CANDIDATES;i++) /* update C */
	    {
	      C[i].nB=nextC[i].nB;
	      for(j=0;j<MAX_BREAKPOINTS;j++)
		  {
		    C[i].B[j][1]=nextC[i].B[j][1];
		    C[i].B[j][2]=nextC[i].B[j][2];
	      }
	    }
 	  }      
    }


  /* parsing command line parameters */
  
  if( (argc==2) && ( (strcmp(argv[1],"--help")==0) || (strcmp(argv[1],"-help")==0) || (strcmp(argv[1],"-h")==0)) ) /* want help? */
    {
      printf("msalign - (c) Magnus Palmblad 2007-\n\nusage: msalign -1<LC-MS dataset 1 filename> -2<LC-MS dataset 2 filename> -e<max. mass error in MS-only data (in ppm)>  [-b<background> -l<typical standard deviation in LC retention time in LC-MS data> -X<Xmax> -Y<Ymax> -R<MS start scan>,<MS end scan> -o<output file> -nographics]\n\nfor more information, see http://www.ms-utils.org/msalign or e-mail magnus.palmblad@gmail.com\n");
      return 0;
    }
  
  if (argc<4 || argc>18) /* check for incorrect number of parameters */
    {
      printf("usage: msalign -1<LC-MS dataset 1 filename> -2<LC-MS dataset 2 filename> -e<max. mass error in MS-only data (in ppm)>  [-b<background> -l<typical standard deviation in LC retention time in LC-MS data> -X<Xmax> -Y<Ymax> -R<MS start scan>,<MS end scan> -o<output file> -nographics] (type msalign --help for more information)\n");
      return -1;
    }
  

  LC_sigma_set=0; Xmax_set=0; Ymax_set=0;
  for(i=1;i<argc;i++) {
    if( (argv[i][0]=='-') && (argv[i][1]=='1') ) 
      {
	strcpy(mzXML_filename_1,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); 
	strcpy(output_filename,mzXML_filename_1); strcat(output_filename,".alignment"); /* default output filename */
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='2') )
    {	
    	strcpy(mzXML_filename_2,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); 
	    
    }
    if( (argv[i][0]=='-') && (argv[i][1]=='e') ) max_mass_measurement_error=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
//    if( (argv[i][0]=='-') && (argv[i][1]=='b') ) background=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='c') ) costs=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='d') ) argFitness=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='f') ) SetFeatures=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='R') ) 
      {
	strcpy(temp,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); p=strtok(temp,","); 
	MS_start_scan=atof(p); p=strtok('\0',","); MS_end_scan=atof(p);
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='l') ) 
      {
	LC_sigma=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
	LC_sigma_set=1;
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='X') ) 
      {
	Xmax=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
	Xmax_set=1;
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='Y') ) 
      {
	Ymax=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
	Ymax_set=1;
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='o') ) strcpy(output_filename,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( strcmp(argv[i],"-nographics")==0 ) nographics=1;
  }
  
  
  /* initialize mzXML dataset 1 */
  
  printf("checking mzXML dataset 1..."); fflush(stdout);
  mzXML_file=rampOpenFile(mzXML_filename_1);

    /* Read the offset of the index */
    indexOffset = getIndexOffset (mzXML_file);
    
    /* Read the scan index into a vector, get LastScan */
    pScanIndex = readIndex(mzXML_file, indexOffset, &iLastScan);
    
    fflush(stdout);
    readRunHeader(mzXML_file, pScanIndex, &runHeader, iLastScan);

    n_unique_masses=MAX_DATAPOINTS;
    n_unique_masses1=MAX_DATAPOINTS;
    n_unique_masses2=MAX_DATAPOINTS;
    SIC_1_max=(double*) malloc(n_unique_masses*sizeof(double));
    SIC_1_argmax=(int*) malloc(n_unique_masses*sizeof(int));
    SIC_2_max=(double*) malloc(n_unique_masses*sizeof(double));
    SIC_2_argmax=(int*) malloc(n_unique_masses*sizeof(int));
    SIC_2a_max=(double*) malloc(n_unique_masses*sizeof(double));
    SIC_2a_argmax=(int*) malloc(n_unique_masses*sizeof(int));
    SIC_1_max_mz=(double*) malloc(n_unique_masses*sizeof(double));
    SIC_2_max_mz=(double*) malloc(n_unique_masses*sizeof(double));
    SIC_2a_max_mz=(double*) malloc(n_unique_masses*sizeof(double));
    n_masses = MAX_DATAPOINTS;                     // no value was asigned to n_masses so this was added!!!!

    printf("done\nreading MS scans between %d and %d...",MS_start_scan,MS_end_scan); fflush(stdout);
    go = 1;    //initalize the background checking
    CountFeatures = 0;
    background1 = 100000;

// ------------ determine background of first dataset
  CountAttempts = 0; // initialize the Attempts to find the features
  while (go == 1)
  {
    for(k=0;k<n_masses;k++) {SIC_1_max[k]=0.0; SIC_1_argmax[k]=-1; SIC_1_max_mz[k]=0.0;}
    for(k=0;k<n_masses;k++) {SIC_2_max[k]=0.0; SIC_2_argmax[k]=-1; SIC_2_max_mz[k]=0.0;}   // this comes back later, why??
    n_MS_spectra=0; n_unique_masses1=0;
    for(range[0]=MS_start_scan;range[0]<MS_end_scan;range[0]++)
    { 
      readHeader(mzXML_file, pScanIndex[range[0]], &scanHeader);
      if((scanHeader.msLevel==1) && (scanHeader.peaksCount>0)) /* non-empty MS spectrum? */
      {
        range[1]=range[0];
        if(getMsSpect(&mzXML_spectrum,mzXML_file,range)<0) continue; /* skip spectrum */
	for(i=0;i<mzXML_spectrum.size;i++) 
	{
	  if(mzXML_spectrum.yval[i]>background1)
	  {
	    new_mass=1;
	    for(k=0;k<n_unique_masses1;k++)
	    {
	      if(((fabs(mzXML_spectrum.xval[i]-SIC_1_max_mz[k]))/SIC_1_max_mz[k]*1e6)<max_mass_measurement_error) 
	      {
	        new_mass=0;
	        if(mzXML_spectrum.yval[i]>SIC_1_max[k]) {SIC_1_max[k]=mzXML_spectrum.yval[i]; SIC_1_argmax[k]=range[0];}
 	      }
	    }
	    if(new_mass==1) 
	    {
	      n_unique_masses1++;
	      SIC_1_max_mz[k]=mzXML_spectrum.xval[i]; SIC_1_max[k]=mzXML_spectrum.yval[i]; SIC_1_argmax[k]=range[0];
	    }
	  }
	}
	freeSpectStrct(mzXML_spectrum);
      }
    }
//    printf("%d - %#f\n", n_unique_masses1, background1); fflush(stdout);
    CountFeatures = n_unique_masses1;
    if (SetFeatures > CountFeatures)
    {
      oldbackground1 = background1;
      background1 = abs(background1 * 0.5);
      CountAttempts++;
    }
    else
    {
      oldbackground1 = background1;
      background1 = background1 + abs(background1 * 0.5);
      CountAttempts++;
    }
    if (CountAttempts > 100)
    {
      printf("\nCould not find background, aborted!"); fflush(stdout);
      exit(42); // it's not working so bye
    }
    if ((CountFeatures < (SetFeatures + abs(0.1*SetFeatures))) && (CountFeatures > (SetFeatures - abs(0.1*SetFeatures))))
    {
      go = 0;
    }
  }
// --------end determining background first dataset  

    printf("done (read %d MS spectra and %d unique masses above background)\n",n_MS_spectra,n_unique_masses1); fflush(stdout);

    /* initialize mzXML dataset 2 */
    
    printf("checking mzXML dataset 2..."); fflush(stdout);
    mzXML_file=rampOpenFile(mzXML_filename_2);

    /* Read the offset of the index */
    indexOffset = getIndexOffset (mzXML_file);
    
    /* Read the scan index into a vector, get LastScan */
    pScanIndex = readIndex(mzXML_file, indexOffset, &iLastScan);
    
    fflush(stdout);
    readRunHeader(mzXML_file, pScanIndex, &runHeader, iLastScan);

    printf("done\nreading MS scans between %d and %d...",MS_start_scan,MS_end_scan); fflush(stdout);
    go = 1;    //initalize the background checking
    CountFeatures = 0;
    background2 = 100000;

// ------------ determine background of second dataset
  CountAttempts = 0; // initialize the Attempts to find the features
  while (go == 1)
  {
    for(k=0;k<n_masses;k++) {SIC_2a_max[k]=0.0; SIC_2a_argmax[k]=-1; SIC_2a_max_mz[k]=0.0;}   // this comes back later, why??
    n_MS_spectra=0; n_unique_masses2=0;
    for(range[0]=MS_start_scan;range[0]<MS_end_scan;range[0]++)
    { 
      readHeader(mzXML_file, pScanIndex[range[0]], &scanHeader);
      if((scanHeader.msLevel==1) && (scanHeader.peaksCount>0)) /* non-empty MS spectrum? */
      {
        range[1]=range[0];
        if(getMsSpect(&mzXML_spectrum,mzXML_file,range)<0) continue; /* skip spectrum */
	for(i=0;i<mzXML_spectrum.size;i++) 
	{
	  if(mzXML_spectrum.yval[i]>background2)
	  {
	    new_mass=1;
	    for(k=0;k<n_unique_masses2;k++)
	    {
	      if(((fabs(mzXML_spectrum.xval[i]-SIC_2a_max_mz[k]))/SIC_2a_max_mz[k]*1e6)<max_mass_measurement_error) 
	      {
	        new_mass=0;
	        if(mzXML_spectrum.yval[i]>SIC_2a_max[k]) {SIC_2a_max[k]=mzXML_spectrum.yval[i]; SIC_2a_argmax[k]=range[0];}
 	      }
	    }
	    if(new_mass==1) 
	    {
	      n_unique_masses2++;
	      SIC_2a_max_mz[k]=mzXML_spectrum.xval[i]; SIC_2a_max[k]=mzXML_spectrum.yval[i]; SIC_2a_argmax[k]=range[0];
	    }
	  }
	}
	freeSpectStrct(mzXML_spectrum);
      }
    }
//    printf("%d - %#f\n", n_unique_masses2, background2); fflush(stdout);
    CountFeatures = n_unique_masses2;
    if (SetFeatures > CountFeatures)
    {
      oldbackground2 = background2;
      background2 = abs(background2 * 0.5);
      CountAttempts++;
    }
    else
    {
      oldbackground2 = background2;
      background2 = background2 + abs(background2 * 0.5);
      CountAttempts++;
    }
    if (CountAttempts > 100)
    {
      printf("Could not find background"); fflush(stdout);
      exit(42); // it's not working so bye
    }
    if ((CountFeatures < (SetFeatures + abs(0.1*SetFeatures))) && (CountFeatures > (SetFeatures - abs(0.1*SetFeatures))))
    {
      go = 0;
    }
  }
// --------end determining background second dataset  

    background2 = oldbackground2;    // get the determined background for the second dataset

// -------- find common unique masses
//  printf("%d - %d", n_unique_masses1, n_unique_masses2); fflush(stdout);
  n_unique_masses = n_unique_masses1;
  for (i=0;i<n_unique_masses1;i++)
  {
    for (k=0;k<n_unique_masses2;k++)
    {
      if(((fabs(SIC_1_max_mz[i]-SIC_2a_max_mz[k]))/SIC_2a_max_mz[k]*1e6)<max_mass_measurement_error) // is this m/z found back in the first file 
      {
        if(SIC_2a_max[k]>SIC_2_max[i]) {SIC_2_max[i]=SIC_2a_max[k]; SIC_2_max_mz[i]=SIC_2a_max_mz[k]; SIC_2_argmax[i]=SIC_2a_argmax[k];}
      }
    }
  }
   
    
    printf("done (read %d MS spectra and %d unique masses above background)\nsaving time pairs to .alignment file...",n_MS_spectra,n_unique_masses2); fflush(stdout);

    if ((outp=fopen(output_filename, "w"))==NULL) {printf("error opening output file %s",output_filename); return -1;}

    fprintf(outp,"m/z 1\tintensity 1\tscan 1\tm/z 2\tintensity 2\tscan 2\n");
    for(k=0;k<n_unique_masses;k++) 
    {		
	  fprintf(outp,"%1.4f %1.1f %d %1.4f %1.1f %d\n",SIC_1_max_mz[k], SIC_1_max[k], SIC_1_argmax[k], SIC_2_max_mz[k], SIC_2_max[k], SIC_2_argmax[k]); 
	  A[k][1]=SIC_1_argmax[k];
	  A[k][2]=SIC_2_argmax[k];
	  error[k]=1e6*(SIC_2_max_mz[k]-SIC_1_max_mz[k])/SIC_1_max_mz[k];
	  // printf("%f %f %f\n",unique_peptide[k].mz,SIC_max_mz[k],1e6*(SIC_max_mz[k]-unique_peptide[k].mz)/unique_peptide[k].mz);
    }
    fflush(outp);
    fclose(outp);
    printf("done\n");
    
    /* sort errors to calculate median error (for QC or crude recalibration) */

    for(i=0;i<n_unique_masses-1;i++)
    {
	  for(j=0;j<n_unique_masses-1-i;j++)
	  {
	    if(error[j+1]>error[j]) 
	    {  
		  swap=error[j];       
		  error[j]=error[j+1];
		  error[j+1]=swap;
	    }
	  }
    }
    printf("median mass measurement error %3.3f ppm\n",error[(int)floor(n_unique_masses/2)]);
    

    if(!LC_sigma_set)
    {

      /* estimate the 'standard deviation' of the retention time in the LC-MS dataset */
	
      /* sort peptides (A) according to LC retention time in LC-MS/MS */
	
	  for(i=0;i<n_unique_masses-1;i++)
	  {
	    for(j=0;j<n_unique_masses-1-i;j++)
	    {
		  if(A[j+1][1]<A[j][1]) 
		  {  
		    swap=A[j][2];       
		    A[j][2]=A[j+1][2];
		    A[j+1][2]=swap;

		    swap=A[j][1];       
		    A[j][1]=A[j+1][1];
		    A[j+1][1]=swap;
		  }
	    } 
	  }

	  /* find median and standard deviation of difference between two consecutive points */
	
	  non_matched=0;
      for(i=0;i<n_unique_masses-1;i++) if((A[i][2]>-1)&&(A[i+1][2]>-1)) scan_diff[i]=A[i+1][2]-A[i][2]; else {scan_diff[i]=99999999; non_matched+=1;}
	  
	  for(i=0;i<n_unique_masses-1;i++)
	  {
	    for(j=0;j<n_unique_masses-1-i;j++)
	    {
		  if(scan_diff[j+1]>scan_diff[j]) 
		  {  
		    swap=scan_diff[j];       
		    scan_diff[j]=scan_diff[j+1];
		    scan_diff[j+1]=swap;
		  }
	    }
	  }   

//	  for(i=0;i<n_unique_masses-1;i++) printf("%4.2f\n",scan_diff[i]);
	
	  scan_diff_median=scan_diff[(int)floor((n_unique_masses+non_matched)/2)];
      printf("median difference between consecutive peptides %3.3f (scans or time units)\n",scan_diff_median);

	  sum_squared_deviation=0;


	  for(i=non_matched+(int)floor((double)(n_unique_masses-non_matched)*0.2+0.5);i<n_unique_masses-1-floor((double)(n_unique_masses-non_matched)*0.2+0.5);i++) if(scan_diff[i]<99999998) sum_squared_deviation+=(scan_diff[i]-scan_diff_median)*(scan_diff[i]-scan_diff_median); /* this is a key line in the estimation of the residual standard deviation - keeping only the central 60% of the distribution of differences from the median */ // added (int) 25-11-2008
      
	  LC_sigma=sqrt(sum_squared_deviation/(n_unique_masses-2-non_matched));

	  printf("estimated standard deviation in LC retention time (in LC-MS data) %3.3f (scans or time units)\n",LC_sigma);

    }

    if(!Xmax_set) Xmax=floor((double)(Xmax_temp*1.2+0.5)); 
    if(!Ymax_set) Ymax=floor((double)(MS_end_scan*1.2+0.5));
    //printf("Xmax=%4.2f, Ymax=%4.2f\n",Xmax,Ymax); fflush(stdout);

    best_C=GAalign(A,n_unique_masses); /* align the LC-MS/MS and LC-MS datasets */

    strcat(output_filename,".plf");
    if ((outp=fopen(output_filename, "w"))==NULL) {printf("error opening output .plf file %s",output_filename); return -1;}
    
    printf("best fit = [0.0 0.0");
    fprintf(outp,"0.0 0.0\n");
    for(i=0;i<MAX_BREAKPOINTS;i++)
    {
      if(best_C.B[i][1]>0) 
	  {
	    printf("; %1.1f %1.1f",best_C.B[i][1],best_C.B[i][2]);
	    fprintf(outp,"%1.1f %1.1f\n",best_C.B[i][1],best_C.B[i][2]);
	  }
    }
    printf("]\nwriting best fit to .plf file...(done)\n"); 

    
    /* close files and free memory */

    fclose(outp);   
    rampCloseFile(mzXML_file);
    return 0;
}
