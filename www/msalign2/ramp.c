/***************************************************************************
                             RAMP


Non sequential parser for mzXML files
and mzData files, too!

                             -------------------
    begin                : Wed Oct 10
    copyright            : (C) 2003 by Pedrioli Patrick, ISB, Proteomics
    email                : ppatrick@student.ethz.ch
    additional work for C++, >2GB files in WIN32, and portability (C) 2004 by Brian Pratt, Insilicos LLC 
    additional work for mzData input (C) 2005 Brian Pratt Insilicos LLC
 ***************************************************************************/

/***************************************************************************
*    This program is free software; you can redistribute it and/or modify  *
*    it under the terms of the GNU Library or "Lesser" General Public      *
*    License (LGPL) as published by the Free Software Foundation;          *
*    either version 2 of the License, or (at your option) any later        *
*    version.                                                              *
***************************************************************************/

#define _O_BINARY O_BINARY // for gcc under Cygwin
#define _O_RDONLY O_RDONLY // for gcc under Cygwin

#include "ramp.h"

#include "base64.h"

#include "time.h"

#include "zlib.h"

#define SIZE_BUF 512

#if defined __LITTLE_ENDIAN
#define swapbytes(x) ntohl(x)  /* use system byteswap (ntohl is a noop on bigendian systems) */
#else
uint32_t swapbytes(uint32_t x) {
  return ((x & 0x000000ffU) << 24) | 
         ((x & 0x0000ff00U) <<  8) | 
         ((x & 0x00ff0000U) >>  8) | 
         ((x & 0xff000000U) >> 24);
}
#endif

uint64_t swapbytes64(uint64_t x) {
  return ((((uint64_t)swapbytes((uint32_t)(x & 0xffffffffU)) << 32) | 
            (uint64_t)swapbytes((uint32_t)(x >> 32))));
}

/****************************************************************
 * Utility functions					*
 ***************************************************************/

static char *findquot(const char *cp) { /* " and ' are both valid attribute delimiters */
   char *result = strchr(cp,'\"');
   if (!result) {
      result = strchr(cp,'\'');
   }
   return result;
}

static int isquot(const char c) { /* " and ' are both valid attribute delimiters */
      return ('\"'==c)||('\''==c);
}

static int setTagValue(const char* text,
      char* storage,
      int maxlen,
      const char* lead);

const char* skipspace(const char* pStr)
{
	while (isspace(*pStr))
		pStr++;
	if (*pStr == '\0')
		return NULL;
	return pStr;
}

static void getIsLittleEndian(const char *buf, int *result) {
   const char *p = strstr(buf,"byteOrder");
   if (p) {
      p = findquot(p);
      if (p++) {
         *result = (0!=strncmp(p,"network",7));
      }
   }
   // leave *result alone if we don't see anything!
}

/****************************************************************
 * mzML helper functions					*
 ***************************************************************/

/*
 * Reads a spectrum element starting at offset, up to stopTag or the
 * end of the spectrum, whichever comes first.
 * Returns a null terminated buffer that becomes property of the caller.
 */
static char *mzmlReadSpectrum(RAMPFILE *pFI, ramp_fileoffset_t offset, const char *stopTag)
{
   size_t alloc = 4096;
   size_t size = 0;
   size_t search = 0; // where to start looking for the stop tag
   size_t nread;
   char *result = (char *) malloc(alloc + 1);
   char *end;

   if (!result) {
      printf("Cannot allocate memory\n");
      return NULL;
   }
   ramp_fseek(pFI, offset, SEEK_SET);
   while (1) {
      if (size == alloc) {
         char *tmp = (char *) realloc(result, 2 * alloc + 1);
         if (!tmp) {
            printf("Cannot allocate memory\n");
            free(result);
            return NULL;
         }
         result = tmp;
         alloc *= 2;
      }
      nread = ramp_fread(result + size, alloc - size, pFI);
      if (nread == 0) {
         break; // end of file
      }
      size += nread;
      result[size] = 0;
      if (strstr(result + search, stopTag) || strstr(result + search, "</spectrum>")) {
         break;
      }
      search = size > 32 ? size - 32 : 0; // the tag may overhang the end of the buffer
   }
   result[size] = 0;
   // cut off everything after the stop tag, it belongs to the next element
   if ((end = strstr(result, stopTag))) {
      *end = 0;
   }
   if ((end = strstr(result, "</spectrum>"))) {
      *end = 0;
   }
   return result;
}

/*
 * Finds attribute attr in the tag starting at tag.
 * Returns a pointer to the (not null terminated) value and its length in *len,
 * or NULL if the tag doesn't have the attribute.
 */
static const char *mzmlAttr(const char *tag, const char *attr, int *len)
{
   const char *tagEnd = strchr(tag, '>');
   const char *p = tag;
   size_t attrLen = strlen(attr);

   while ((p = strstr(p + 1, attr)) && (!tagEnd || p < tagEnd)) {
      if (isspace((unsigned char) p[-1]) && (p[attrLen] == '=') && isquot(p[attrLen + 1])) {
         const char *value = p + attrLen + 2;
         const char *valueEnd = strchr(value, p[attrLen + 1]);
         if (!valueEnd) {
            return NULL;
         }
         *len = (int) (valueEnd - value);
         return value;
      }
   }
   return NULL;
}

/*
 * Finds the cvParam with the given accession between start and end (NULL = end of string).
 * Returns a pointer to the cvParam tag, or NULL if not found.
 */
static const char *mzmlFindCv(const char *start, const char *end, const char *accession)
{
   const char *p = start;
   int accLen = (int) strlen(accession);

   while ((p = strstr(p, "<cvParam")) && (!end || p < end)) {
      int len;
      const char *acc = mzmlAttr(p, "accession", &len);
      if (acc && (len == accLen) && !strncmp(acc, accession, accLen)) {
         return p;
      }
      p += 8;
   }
   return NULL;
}

/*
 * Reads the value of the cvParam with the given accession into *value.
 * Returns 1 if found, 0 otherwise.
 */
static int mzmlCvDouble(const char *start, const char *end, const char *accession, double *value)
{
   int len;
   const char *v;
   const char *cv = mzmlFindCv(start, end, accession);

   if (!cv) {
      return 0;
   }
   v = mzmlAttr(cv, "value", &len);
   *value = v ? atof(v) : 0.0;
   return 1;
}

/*
 * Reads a time cvParam, returns the time in seconds.
 */
static double mzmlTime(const char *cv)
{
   int len;
   const char *unit;
   const char *v = mzmlAttr(cv, "value", &len);
   double t = v ? atof(v) : 0.0;

   if (((unit = mzmlAttr(cv, "unitAccession", &len)) && !strncmp(unit, "UO:0000031", len)) ||
       ((unit = mzmlAttr(cv, "unitName", &len)) && !strncmp(unit, "minute", len))) {
      t *= 60.0;
   } else if (((unit = mzmlAttr(cv, "unitAccession", &len)) && !strncmp(unit, "UO:0000028", len)) ||
              ((unit = mzmlAttr(cv, "unitName", &len)) && !strncmp(unit, "millisecond", len))) {
      t /= 1000.0;
   }
   return t;
}

/*
 * Decodes one binaryDataArray (between arr and arrEnd) into count doubles.
 * Returns the number of values decoded.
 * Stops the program if the array uses a compression we can't read.
 */
static int mzmlDecodeArray(const char *arr, const char *arrEnd, int count, double *values)
{
   int endtest = 1;
   int weAreLittleEndian = *((char *)&endtest);
   int bytesPerValue = 4;
   int isFloat = 1;
   int isCompressed = 0;
   int n, len;
   size_t encodedLen = 0, decodedLen;
   const char *binary, *p;
   char *encoded, *decoded;
   unsigned char *data;

   // precision
   if (mzmlFindCv(arr, arrEnd, "MS:1000523")) { // 64-bit float
      bytesPerValue = 8;
   } else if (mzmlFindCv(arr, arrEnd, "MS:1000519")) { // 32-bit integer
      isFloat = 0;
   } else if (mzmlFindCv(arr, arrEnd, "MS:1000522")) { // 64-bit integer
      bytesPerValue = 8;
      isFloat = 0;
   }

   // compression
   if (mzmlFindCv(arr, arrEnd, "MS:1000574")) { // zlib
      isCompressed = 1;
   } else if (!mzmlFindCv(arr, arrEnd, "MS:1000576")) { // not "no compression" either
      p = strstr(arr, "Numpress");
      if (p && (p < arrEnd)) {
         printf("\nMS-Numpress compressed mzML files are not supported, convert the file without numpress, aborted!\n");
      } else {
         printf("\nUnsupported compression in mzML file, convert the file with zlib or no compression, aborted!\n");
      }
      fflush(stdout);
      exit(42);
   }

   // number of values in this array
   if ((p = mzmlAttr(arr, "arrayLength", &len)) && (atoi(p) < count)) {
      count = atoi(p);
   }
   if (count <= 0) {
      return 0;
   }

   // find the base64 text
   binary = strstr(arr, "<binary>");
   if (!binary || (binary > arrEnd)) {
      return 0;
   }
   binary += 8;
   p = strchr(binary, '<');
   if (!p) {
      return 0;
   }
   if ((encoded = (char *) malloc(p - binary + 1)) == NULL) {
      printf("Cannot allocate memory\n");
      return 0;
   }
   for (; binary < p; binary++) { // whitespace may be present in base64 char stream
      if (!isspace((unsigned char) *binary)) {
         encoded[encodedLen++] = *binary;
      }
   }
   encoded[encodedLen] = 0;

   decodedLen = (encodedLen / 4) * 3;
   if ((decoded = (char *) malloc(decodedLen + 1)) == NULL) {
      printf("Cannot allocate memory\n");
      free(encoded);
      return 0;
   }
   b64_decode_mio(decoded, encoded);
   // every '=' padding character at the end means one byte less
   for (len = (int) encodedLen - 1; (len >= 0) && (encoded[len] == '=') && (decodedLen > 0); len--) {
      decodedLen--;
   }
   free(encoded);

   if (isCompressed) {
      uLong uncomprLen = (uLong) count * bytesPerValue;
      data = (unsigned char *) calloc(uncomprLen + 1, 1);
      if (!data) {
         printf("Cannot allocate memory\n");
         free(decoded);
         return 0;
      }
      if (uncompress(data, &uncomprLen, (const Bytef *) decoded, decodedLen) != Z_OK) {
         printf("Error uncompressing mzML data\n");
      }
      free(decoded);
      if ((int) (uncomprLen / bytesPerValue) < count) {
         count = (int) (uncomprLen / bytesPerValue);
      }
   } else {
      data = (unsigned char *) decoded;
      if ((int) (decodedLen / bytesPerValue) < count) {
         count = (int) (decodedLen / bytesPerValue);
      }
   }

   // mzML binary data is always little endian
   for (n = 0; n < count; n++) {
      if (bytesPerValue == 4) {
         uint32_t u;
         memcpy(&u, data + 4 * n, 4);
         if (!weAreLittleEndian) {
            u = swapbytes(u);
         }
         if (isFloat) {
            float f;
            memcpy(&f, &u, 4);
            values[n] = (double) f;
         } else {
            values[n] = (double) (int32_t) u;
         }
      } else {
         uint64_t u;
         memcpy(&u, data + 8 * n, 8);
         if (!weAreLittleEndian) {
            u = swapbytes64(u);
         }
         if (isFloat) {
            double d;
            memcpy(&d, &u, 8);
            values[n] = d;
         } else {
            values[n] = (double) (int64_t) u;
         }
      }
   }
   free(data);
   return count;
}

/*
 * Reads the scan header of the mzML spectrum at lScanIndex.
 * scanHeader must be initialized with the defaults already.
 */
static void mzmlReadHeader(RAMPFILE *pFI,
                           ramp_fileoffset_t lScanIndex,
                           struct ScanHeaderStruct *scanHeader)
{
   char *spectrum = mzmlReadSpectrum(pFI, lScanIndex, "<binaryDataArrayList");
   const char *tag, *cv, *v;
   double value;
   int len;

   if (!spectrum) {
      return;
   }
   if ((tag = strstr(spectrum, "<spectrum ")) != NULL) {
      if ((v = mzmlAttr(tag, "index", &len))) {
         scanHeader->acquisitionNum = atoi(v) + 1; // mzML index is 0-based
      }
      if ((v = mzmlAttr(tag, "defaultArrayLength", &len))) {
         scanHeader->peaksCount = atoi(v);
      }
      if (mzmlCvDouble(tag, NULL, "MS:1000511", &value)) { // ms level
         scanHeader->msLevel = (int) value;
      } else if (mzmlFindCv(tag, NULL, "MS:1000579")) { // MS1 spectrum
         scanHeader->msLevel = 1;
      }
      if ((cv = mzmlFindCv(tag, NULL, "MS:1000016"))) { // scan start time
         scanHeader->retentionTime = mzmlTime(cv);
      }
      if (mzmlCvDouble(tag, NULL, "MS:1000501", &value) || // scan window lower limit
          mzmlCvDouble(tag, NULL, "MS:1000528", &value)) { // lowest observed m/z
         scanHeader->lowMZ = value;
      }
      if (mzmlCvDouble(tag, NULL, "MS:1000500", &value) || // scan window upper limit
          mzmlCvDouble(tag, NULL, "MS:1000527", &value)) { // highest observed m/z
         scanHeader->highMZ = value;
      }
      mzmlCvDouble(tag, NULL, "MS:1000285", &(scanHeader->totIonCurrent));
      mzmlCvDouble(tag, NULL, "MS:1000504", &(scanHeader->basePeakMZ));
      mzmlCvDouble(tag, NULL, "MS:1000505", &(scanHeader->basePeakIntensity));
      mzmlCvDouble(tag, NULL, "MS:1000744", &(scanHeader->precursorMZ)); // selected ion m/z
      if (mzmlCvDouble(tag, NULL, "MS:1000041", &value)) { // charge state
         scanHeader->precursorCharge = (int) value;
      }
      if (mzmlCvDouble(tag, NULL, "MS:1000042", &value)) { // peak intensity
         scanHeader->precursorIntensity = value;
      }
   }
   free(spectrum);
}

/*
 * Reads the peaks of the mzML spectrum at lScanIndex.
 * Returns a RAMPREAL* with m/z-int pairs, terminated by -1, which becomes
 * property of the caller, or NULL if the spectrum has no m/z or intensity array.
 */
static RAMPREAL *mzmlReadPeaks(RAMPFILE *pFI,
                               ramp_fileoffset_t lScanIndex)
{
   char *spectrum = mzmlReadSpectrum(pFI, lScanIndex, "</spectrum>");
   const char *tag, *arr, *arrEnd, *v;
   double *mz = NULL;
   double *inten = NULL;
   int nMz = 0, nInten = 0;
   int peaksCount = 0;
   int len, n;
   RAMPREAL *pPeaks = NULL;

   if (!spectrum) {
      return NULL;
   }
   if (!(tag = strstr(spectrum, "<spectrum ")) || !(v = mzmlAttr(tag, "defaultArrayLength", &len)) ||
       ((peaksCount = atoi(v)) <= 0)) {
      free(spectrum);
      return NULL; // No peaks in this scan!!
   }
   mz = (double *) malloc(peaksCount * sizeof(double));
   inten = (double *) malloc(peaksCount * sizeof(double));
   if (!mz || !inten) {
      printf("Cannot allocate memory\n");
      free(mz);
      free(inten);
      free(spectrum);
      return NULL;
   }

   // intensity and mz are written in two different arrays
   for (arr = tag; (arr = strstr(arr, "<binaryDataArray ")); arr = arrEnd) {
      if (!(arrEnd = strstr(arr, "</binaryDataArray>"))) {
         break;
      }
      if (mzmlFindCv(arr, arrEnd, "MS:1000514")) { // m/z array
         nMz = mzmlDecodeArray(arr, arrEnd, peaksCount, mz);
      } else if (mzmlFindCv(arr, arrEnd, "MS:1000515")) { // intensity array
         nInten = mzmlDecodeArray(arr, arrEnd, peaksCount, inten);
      }
   }
   free(spectrum);

   if ((nMz > 0) && (nInten > 0)) {
      peaksCount = nMz < nInten ? nMz : nInten;
      if ((pPeaks = (RAMPREAL *) malloc((peaksCount + 1) * 2 * sizeof(RAMPREAL))) == NULL) {
         printf("Cannot allocate memory\n");
      } else {
         for (n = 0; n < peaksCount; n++) {
            pPeaks[2 * n] = (RAMPREAL) mz[n];
            pPeaks[2 * n + 1] = (RAMPREAL) inten[n];
         }
         pPeaks[2 * peaksCount] = -1; // some callers want a terminator
      }
   }
   free(mz);
   free(inten);
   return pPeaks; // caller must free this pointer
}

/**************************************************
* open and close files *
**************************************************/
RAMPFILE *rampOpenFile(const char *filename) {
   RAMPFILE *result = (RAMPFILE *)calloc(1,sizeof(RAMPFILE));
   if (result) {
      int bOK;
#ifndef _LARGEFILE_SOURCE
     result->fileHandle = open(filename,_O_BINARY|_O_RDONLY);
     bOK = (result->fileHandle >= 0);
#else
     result->fileHandle = fopen(filename,"rb");
     bOK = (result->fileHandle != NULL);
#endif
     if (!bOK) {
        free(result);
        result = NULL;
     } else {
        char buf[1024];
        buf[sizeof(buf)-1] = 0;
        while (!ramp_feof(result)) {
           ramp_fgets(buf,sizeof(buf)-1,result);
           if (strstr(buf,"<mzML")) {
              result->bIsMzML = 1;
              break;
           } else if (strstr(buf,"<msRun")) {
              result->bIsMzData = 0;
              break;
           } else if (strstr(buf,"<mzData")) {
              result->bIsMzData = 1;
              break;
           }
        }
        ramp_fseek(result,0,SEEK_SET); // rewind
        // now set up the element tags we seek
     }
   }
   return result;
}

void rampCloseFile(RAMPFILE *pFI) {
   if (pFI) {
#ifndef _LARGEFILE_SOURCE
      close(pFI->fileHandle);
#else
      fclose(pFI->fileHandle);
#endif
      free(pFI);
   }
}

/**************************************************
* fgets() for win32 long files *
* TODO: this could be a LOT more efficient...
**************************************************/
#ifndef _LARGEFILE_SOURCE
char *ramp_fgets(char *buf,int len,RAMPFILE *handle) {
   int nread=0;
   int chunk;
   ramp_fileoffset_t pos = ramp_ftell(handle);
   buf[--len]=0; // nullterm for safety
   chunk = max(len/4,1); // usually all that's needed is a short read
   while (nread <= len) {
      char *newline;
      int nread_now = ramp_fread(buf+nread,chunk,handle);
      buf[nread+nread_now] = 0;
      if (!nread_now) {
         return nread?buf:NULL;
      }
      newline = strchr(buf+nread,'\n');
      if (newline) {
         *newline = 0;
         ramp_fseek(handle,pos+(newline-buf)+1,SEEK_SET); // so next read is at next line
         break;
      }
      nread+=nread_now;
      if (nread >= len) {
         break;
      }
      // apparently we need bigger reads
      chunk = len-nread;
   }
   return buf;
}
#endif

/****************************************************************
 * Find the Offset of the index					*
 ***************************************************************/
ramp_fileoffset_t getIndexOffset(RAMPFILE *pFI)
{
   int  i;
   ramp_fileoffset_t  indexOffset, indexOffsetOffset;
   char indexOffsetTemp[SIZE_BUF+1], buf;

   if (pFI->bIsMzData) {
      return -1; // no index in mzData
   }

   // in mzML the fileChecksum follows the indexListOffset, so look further back
   for (indexOffsetOffset = pFI->bIsMzML ? -400 : -120;  indexOffsetOffset++ < 0 ;)
   {
      char seekbuf[SIZE_BUF+1];
      char *target = pFI->bIsMzML ? "<indexListOffset>" : "<indexOffset>";
      int  nread;

      ramp_fseek(pFI, indexOffsetOffset, SEEK_END);
      nread = ramp_fread(seekbuf, strlen(target), pFI);
      seekbuf[nread] = '\0';
      
      if (!strcmp(seekbuf, target))
      {
         break;
      }
   }

   if (indexOffsetOffset >= 0) {
      return -1; // no answer
   }

   indexOffset = ramp_ftell(pFI);

   i = 0;
   while (ramp_fread(&buf, 1, pFI) && buf != '<')
   {
      indexOffsetTemp[i] = buf;
      i++;
   }
   indexOffsetTemp[i] = '\0';

   if (sizeof(ramp_fileoffset_t)==8)
      indexOffset = (atoll(indexOffsetTemp));
   else
      indexOffset = (atol(indexOffsetTemp));
   // now test this
   ramp_fseek(pFI, indexOffset, SEEK_SET);
   ramp_fread(indexOffsetTemp, sizeof(indexOffsetTemp), pFI);
   indexOffsetTemp[sizeof(indexOffsetTemp)-1] = 0;
   if (!strstr(indexOffsetTemp,"<index")) {
      indexOffset = -1; // broken somehow
   }
   return indexOffset;
}


/****************************************************************
 * Reads the Scan index in a list				*
 * Returns pScanIndex which becomes property of the caller	*
 * pScanIndex is -1 terminated					*
 ***************************************************************/
char buf[SIZE_BUF*16];

ramp_fileoffset_t *readIndex(RAMPFILE *pFI,
                ramp_fileoffset_t indexOffset,
                int *iLastScan)
{
   int  n;
   int  reallocSize = 8000;    /* initial # of scan indexes to expect */
   char *beginScanOffset;
   ramp_fileoffset_t *pScanIndex=NULL;
   int retryLoop;
   char* s;
   
   for (retryLoop = 2;retryLoop--;) {
     n = 1; // ramp is one based
     *iLastScan = 0;
     free(pScanIndex);
      if ((indexOffset < 0) || (retryLoop==0)) { // derive the index by inspection
         // no index found, derive it
         const char *scantag = pFI->bIsMzML?"<spectrum ":(pFI->bIsMzData?"<spectrum id":"<scan");
         int taglen = (int)strlen(scantag);
         ramp_fileoffset_t index = 0;
         pScanIndex = (ramp_fileoffset_t *)malloc( sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
         if (!pScanIndex) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         int nbytes;
         ramp_fseek(pFI,0,SEEK_SET);
         buf[sizeof(buf)-1] = 0;
         while ((nbytes = (int)ramp_fread(buf,sizeof(buf)-1,pFI))>taglen) {
            char *find;
            char *look=buf;
            int nread;
            buf[nbytes] = 0; // a short read at the end of the file leaves old data behind
            while (NULL != (find = strstr(look,scantag))) {
               pScanIndex[n++] = index+(find-buf); // ramp is 1-based
               (*iLastScan)++;
               look = find+taglen;
               if (reallocSize<=n) {
                  reallocSize+=500;
                  pScanIndex = (ramp_fileoffset_t *)realloc(pScanIndex, sizeof(ramp_fileoffset_t)*reallocSize); // allocate space for the scan index info
                  if (!pScanIndex) {
                     printf("Cannot allocate memory\n");
                     return NULL;
                  }
               }
            }
            nread = strlen(look)+(look-buf);
            if (*look && strchr(scantag,buf[nread-1])) { // check last char of buffer
               // possible that next scantag overhangs end of buffer
               ramp_fseek(pFI,-taglen,SEEK_CUR); // so next get includes it
            } 
            index = ramp_ftell(pFI);
         }
         break; // no need to retry
      } else {  // read the index out of the file (then check it!)
         struct ScanHeaderStruct scanHeader; // for test of index integrity
         int indexOK=1; // until we show otherwise

         if ((pScanIndex = (ramp_fileoffset_t *) malloc(reallocSize * sizeof(ramp_fileoffset_t))) == NULL) {
            printf("Cannot allocate memory\n");
            return NULL;
         }
         
         ramp_fseek(pFI, indexOffset, SEEK_SET);
         
         s = ramp_fgets(buf, SIZE_BUF, pFI);
         if (pFI->bIsMzML) { // skip to the spectrum index, the chromatogram index may come first
            while( s!=NULL && !strstr( buf , "name=\"spectrum\"" ) ) {
               s = ramp_fgets(buf, SIZE_BUF, pFI);
            }
            if (s != NULL) {
               s = ramp_fgets(buf, SIZE_BUF, pFI);
            }
         }
         while( s!=NULL && !strstr( buf , "<offset" ) ) {
            s = ramp_fgets(buf, SIZE_BUF, pFI);
         }
         
         if (s == NULL)
            continue;   // end of file reached, derive the index by inspection
         
         while (!strstr(buf, "/index")) {
            if ((beginScanOffset = (char *) (strstr(buf, ">"))) == NULL) {
               ramp_fgets(buf, SIZE_BUF, pFI);
               continue;
            }
            beginScanOffset++;
            
            if (sizeof(ramp_fileoffset_t)==8) {
               pScanIndex[n] = atoll(beginScanOffset);
            } else {
               pScanIndex[n] = atol(beginScanOffset);
            }
            n++;
            (*iLastScan)++;
            
            if (n == reallocSize) {
               ramp_fileoffset_t *pTmp;
               
               reallocSize = reallocSize + 500;
               
               pTmp = (ramp_fileoffset_t*) realloc(pScanIndex, reallocSize*sizeof(ramp_fileoffset_t));
               
               if (pTmp==NULL) {
                  printf("Cannot allocate memory\n");
                  return NULL;
               } else {
                  pScanIndex=pTmp;
               }
            }
            
            ramp_fgets(buf, SIZE_BUF, pFI);
         }
         if (n) {
            // OK, now test that to see if it's a good index or not
            readHeader(pFI,pScanIndex[1],&scanHeader);
            if (scanHeader.acquisitionNum == -1) { // first
               indexOK = 0; // bogus index
               free(pScanIndex);
               pScanIndex = NULL;
            } 
            if (indexOK && (n>1)) { // middle
               readHeader(pFI,pScanIndex[n/2],&scanHeader);
               if (scanHeader.acquisitionNum == -1) {
                  indexOK = 0; // bogus index
                  free(pScanIndex);
                  pScanIndex = NULL;
               }
            }
            if (indexOK && (n>1)) { // last
               readHeader(pFI,pScanIndex[n-1],&scanHeader);
               if (scanHeader.acquisitionNum == -1) {
                  indexOK = 0; // bogus index
                  free(pScanIndex);
                  pScanIndex = NULL;
               }
            }
         }
         if (indexOK) {
            break; // no retry
         }
      } // end if we claim to have an index
   } // end for retryloop
   pScanIndex[n] = -1;

   return (pScanIndex);
}

// helper func for reading mzData
const char *findMzDataTagValue(const char *pStr, const char *tag) {
   const char *find = strstr(pStr,tag);
   if (find) {
      find = strstr(find+1,"value=");
      if (find) {
         find = findquot(find);
         if (find) {
            find++; // pointing at value string
         }
      }
   }
   return find;
}

/*
 * Reads a time string, returns time in seconds.
 */
static double rampReadTime(RAMPFILE *pFI,const char *pStr) {
   double t=0;
   if (pFI->bIsMzData) {
      const char *tag = findMzDataTagValue(pStr, "TimeInMinutes");
      if (tag) {
         t = 60.0*atof(tag);
      } else if (NULL!=(tag = findMzDataTagValue(pStr, "TimeInSeconds"))) { // von Steffan Neumann
         t = atof(tag);
      }
   } else if (!sscanf(pStr, "PT%lfS", &t)) {  // usually this is elapsed run time
      /* but could be stored in for PxYxMxDTxHxMxS */
      struct tm fullTime; // apologies to those working after January 18, 19:14:07, 2038
      double secondsFrac=0;
      int bDate = 1; 
      while (!isquot(*++pStr)) {
         double val;
         if ('T'==*pStr) {
            pStr++;
            bDate = 0; // we're into the minutes:seconds portion
         }
         val = atof(pStr);
         while (('.'==*pStr)||isdigit(*pStr)) {
            pStr++;
         }
         switch (*pStr) {
         case 'Y':
            fullTime.tm_year = (int)val-1900; // years since 1900
            break;
         case 'M':
            if (bDate) {
               fullTime.tm_mon = (int)val-1; // range 0-11
            } else {
               fullTime.tm_min = (int)val; // range 0-59
            }            
            break;
         case 'D':
            fullTime.tm_mday = (int)val; // range 1-31
            break;
         case 'H':
            fullTime.tm_hour = (int)val; // range 0-23
            break;
         case 'S':
            fullTime.tm_sec = (int)val;
            secondsFrac = val-(double)fullTime.tm_sec;
            break;
         }
      }
      t = (double)mktime(&fullTime)+secondsFrac;
   }
   return t;
}

/*
 * helper func for faster parsing
 */
const char *matchAttr(const char *where,const char *attr,int len) {
   const char *look = where; // we assume this is pointed at '=', look back at attr
   while (len--) {
      if (*--look != attr[len]) {
         return NULL; // no match
      }
   }
   return where+2; // point past ="
}


//
// helper function to deal with mzdata with no newlines - breaks
// lines up at </...>
//
static char *ramp_nextTag(char *buf, int buflen, RAMPFILE *pFI) {
   char *result;
   buf[buflen-1]=1; // so we can detect full buffer read
   result = ramp_fgets(buf,buflen,pFI);
   if (result && !buf[buflen-1]) { // no newline found
      char *closer = strstr(buf+1,"</");
      if (closer) {
         *closer = 0; // temp. nullterm
         ramp_fseek(pFI,(1+closer-buf)-buflen,SEEK_CUR); // reposition for next read
      }
   }
   return result;
}


/*
 * Reads scan header information.
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE
 *    RETURNING !
 */
void readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex, // look here
                struct ScanHeaderStruct *scanHeader)
{
   char stringBuf[SIZE_BUF];
   char *pStr2;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   /*
    * initialize defaults
    */
   memset(scanHeader,0,sizeof(struct ScanHeaderStruct)); // mostly we want 0's
#define LOWMZ_UNINIT 1.111E6
   scanHeader->lowMZ =  LOWMZ_UNINIT;
   scanHeader->acquisitionNum = -1;
   scanHeader->seqNum = -1;


   if (pFI->bIsMzML) {
      mzmlReadHeader(pFI, lScanIndex, scanHeader);
   } else if (pFI->bIsMzData) {
      int bHasPrecursor = 0;
      while (ramp_nextTag(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         const char *closeTag = strstr(stringBuf, "</spectrumSettings>");
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if (closeTag && (closeTag < attrib)) {
               break; // into data territory now
            }
            if ((pStr = matchAttr(attrib, "spectrum id",11))) {
               sscanf(pStr, "%d", &(scanHeader->acquisitionNum));
            //} else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->basePeakMZ));      
            //} else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->totIonCurrent));      
            //} else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
            //   sscanf(pStr, "%lf", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d", &(scanHeader->msLevel));
            //} else if ((pStr = matchAttr(attrib, "length",6)))  { get this from array element
            //   sscanf(pStr, "%d", &(scanHeader->peaksCount));
            } else if ((pStr = findMzDataTagValue(attrib,"TimeInMinutes")))  {
               scanHeader->retentionTime = rampReadTime(pFI,stringBuf);
            } else if ((pStr = matchAttr(attrib, "mzRangeStart",12)))  {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "mzRangeStop",11)))  {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((pStr = findMzDataTagValue(attrib, "ScanMode"))) { 
               if ((pStr2 = (char *) findquot(pStr))) {
                  memcpy(&(scanHeader->scanType), pStr, pStr2-pStr);
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            }
         }
         if (closeTag) {
            break; // into data territory now
         }
      }
      do {
         
         /*
         * read precursor info
         */
         const char *pStr,*pStr2;
         if ((pStr = (char *) strstr(stringBuf, "<precursorList")))
         {
            bHasPrecursor = 1;
         } else if ((pStr = (char *) strstr(stringBuf, "</precursorList")))
         {
            bHasPrecursor = 0;
         }
         if (bHasPrecursor) { // in precursor section
            if (NULL!=(pStr2 = (char *) strstr(stringBuf, "spectrumRef="))) 
            {
               sscanf(pStr2 + 13, "%d", &(scanHeader->precursorScanNum));
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"ChargeState"))) 
            {
               scanHeader->precursorCharge = atoi(pStr2);
            }
            
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"MassToChargeRatio"))) 
            {
               scanHeader->precursorMZ = atoi(pStr2);
            }
            if (NULL!=(pStr2 = findMzDataTagValue(stringBuf,"mz"))) 
            {
               scanHeader->precursorMZ = atoi(pStr2);
            }
         }
         if (strstr(stringBuf, "</spectrumDesc>")) {
            break; // into data territory now
         }
         if (strstr(stringBuf, "</precursorList>")) {
            break; // into data territory now
         }
      } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
      // now read peaks count
      do {
         if (strstr(stringBuf, "ArrayBinary>")) {
            do {
               char *cp=strstr(stringBuf,"length=");
               if (cp) {
                  scanHeader->peaksCount = atoi(cp+8);
                  break;
               }
            } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
            break;
         }
      } while (ramp_nextTag(stringBuf, SIZE_BUF, pFI));
   } else { // mzXML
      while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
      {
         const char *attrib;
         const char *pStr;
         // find each attribute in stringBuf
         for (attrib=stringBuf-1;NULL!=(attrib=strchr(attrib+1,'='));) {
            if ((pStr = matchAttr(attrib, "num",3))) {
               sscanf(pStr, "%d", &(scanHeader->acquisitionNum));
            } else if ((pStr = matchAttr(attrib, "basePeakMz",10)))  {
               sscanf(pStr, "%lf", &(scanHeader->basePeakMZ));      
            } else if ((pStr = matchAttr(attrib, "totIonCurrent",13)))  {
               sscanf(pStr, "%lf", &(scanHeader->totIonCurrent));      
            } else if ((pStr = matchAttr(attrib, "basePeakIntensity",17)))  {
               sscanf(pStr, "%lf", &(scanHeader->basePeakIntensity));      
            } else if ((pStr = matchAttr(attrib, "msLevel",7)))  {
               sscanf(pStr, "%d", &(scanHeader->msLevel));
            } else if ((pStr = matchAttr(attrib, "peaksCount",10)))  {
               sscanf(pStr, "%d", &(scanHeader->peaksCount));
            } else if ((pStr = matchAttr(attrib, "retentionTime",13)))  {
               scanHeader->retentionTime = rampReadTime(pFI,pStr);
            } else if ((pStr = matchAttr(attrib, "lowMz",5)))  {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));
            } else if ((pStr = matchAttr(attrib, "highMz",6)))  {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((scanHeader->lowMZ==LOWMZ_UNINIT) &&  
               ((pStr = matchAttr(attrib, "startMz",7)))) {
               sscanf(pStr, "%lf", &(scanHeader->lowMZ));  
            } else if ((!scanHeader->highMZ) &&
               ((pStr = matchAttr(attrib, "endMz",5)))) {
               sscanf(pStr, "%lf", &(scanHeader->highMZ));
            } else if ((pStr = matchAttr(attrib, "scanType", 8))) { 
               if ((pStr2 = (char *) findquot(pStr))) {
                  memcpy(&(scanHeader->scanType), pStr, sizeof(char)*((pStr2-pStr)));
                  scanHeader->scanType[pStr2-pStr] = '\0';
               }
            } else if ((pStr = matchAttr(attrib, "collisionEnergy", 15)))  {
                 sscanf(pStr, "%lf", &(scanHeader->collisionEnergy));
            }
         }
         
         /*
         * read precursor mass
         */
         if ((pStr = (char *) strstr(stringBuf, "<precursorMz ")))
         {
            if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum="))) 
               sscanf(pStr2 + 18, "%d", &(scanHeader->precursorScanNum));
            if ((pStr2 = (char *) strstr(pStr, "precursorCharge="))) 
               sscanf(pStr2 + 17, "%d", &(scanHeader->precursorCharge));
            if ((pStr2 = (char *) strstr(pStr, "precursorIntensity="))) 
               sscanf(pStr2 + 20, "%lf", &(scanHeader->precursorIntensity));
            
            /*
            * Find end of tag.
            */
            while (!(pStr = strchr(pStr, '>')))
            {      
               ramp_fgets(stringBuf, SIZE_BUF, pFI);
               pStr = stringBuf;

               if ((pStr2 = (char *) strstr(stringBuf, "precursorScanNum="))) 
                  sscanf(pStr2 + 18, "%d", &(scanHeader->precursorScanNum));
               if ((pStr2 = (char *) strstr(pStr, "precursorCharge="))) 
                  sscanf(pStr2 + 17, "%d", &(scanHeader->precursorCharge));
               if ((pStr2 = (char *) strstr(pStr, "precursorIntensity="))) 
                  sscanf(pStr2 + 20, "%lf", &(scanHeader->precursorIntensity));
            }
            pStr++;	// Skip >
            
            /*
             * Skip past white space.
            */
            while (!(pStr = skipspace(pStr)))
            {
               ramp_fgets(stringBuf, SIZE_BUF, pFI);
               pStr = stringBuf;
            }
            
            sscanf(pStr, "%lf<", &(scanHeader->precursorMZ));
            //         printf("precursorMass = %lf\n", scanHeader->precursorMZ);
         }
         if (strstr(stringBuf, "<peaks")) {
            break; // into data territory now
         }
      }
   } // end else mzXML
   if (!scanHeader->retentionTime) { 
      scanHeader->retentionTime = scanHeader->acquisitionNum; // just some unique nonzero value
   }
   scanHeader->seqNum = scanHeader->acquisitionNum; //  default sequence number
   scanHeader->filePosition = lScanIndex;
}

/****************************************************************
 * Reads the MS level of the scan.				*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

int readMsLevel(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   int  msLevelLen;
   char stringBuf[SIZE_BUF+1];
   char szLevel[12];
   char *beginMsLevel, *endMsLevel;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   while (!(beginMsLevel = (char *) strstr(stringBuf, "msLevel=")))
   {
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

   beginMsLevel += 9;           // We need to move the length of msLevel="
   endMsLevel = (char *) findquot(beginMsLevel);
   msLevelLen = endMsLevel - beginMsLevel;

   strncpy(szLevel, beginMsLevel, msLevelLen);
   szLevel[msLevelLen] = '\0';

   return atoi(szLevel);
}


/****************************************************************
 * Reads startMz and endMz of the scan.				*
 * Returns 1.E6 if startMz was not set.*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readStartMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex)
{
  char stringBuf[SIZE_BUF+1];
  double startMz = 1.E6;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStart":"startMz";

  if (pFI->bIsMzML) {
     struct ScanHeaderStruct header;
     readHeader(pFI, lScanIndex, &header);
     return (header.lowMZ != LOWMZ_UNINIT) ? header.lowMZ : startMz;
  }

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf", &startMz);
        break;
      }
  }

  return startMz;
}


/****************************************************************
 * Reads startMz and endMz of the scan.				*
 * Returns 1.E6 if startMz was not set.	*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

double readEndMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex)
{
  char stringBuf[SIZE_BUF+1];
  double endMz = 0.0;
  char *pStr;
  const char *tag = pFI->bIsMzData?"mzRangeStop":"endMz";

  if (pFI->bIsMzML) {
     struct ScanHeaderStruct header;
     readHeader(pFI, lScanIndex, &header);
     return header.highMZ;
  }

  ramp_fseek(pFI, lScanIndex, SEEK_SET);
   
  while (ramp_fgets(stringBuf, SIZE_BUF, pFI))
  {
     if (strstr(stringBuf, pFI->bIsMzData?"</spectrumDesc":"<peaks"))
        break; // ran to end
      if ((pStr = strstr(stringBuf, tag))){
        sscanf(pStr + strlen(tag)+2, "%lf", &endMz);
        break;
      }
  }

  return endMz;
}

#define N_EXT_TYPES 4
static char *data_ext[N_EXT_TYPES]={".mzXML",".mzData",".mzxml",".mzdata"};


// construct a filename in inbuf from a basename, adding .mzXML or .mzData as exists
// returns inbuf, or NULL if neither .mzXML or .mzData file exists
char *rampConstructInputFileName(char *inbuf,int inbuflen,const char *basename_in) {
   return rampConstructInputPath(inbuf, inbuflen, "", basename_in);
}

char *rampConstructInputPath(char *inbuf,int inbuflen,const char *dir,const char *basename_in) {
   char *result = NULL;
   int i;
   char *basename = (char *)basename_in;
   char *tmpbuf = (char *)malloc(strlen(dir) + strlen(basename) + 20);
   char *append;
   if (basename_in==inbuf) { // same pointer
      basename = (char *)malloc(inbuflen);
      strncpy(basename,basename_in,inbuflen);
   }

   strcpy(tmpbuf, dir);
   strcat(tmpbuf, basename);

   append = tmpbuf + strlen(tmpbuf);

   for (i=0;i<N_EXT_TYPES;i++) {
      FILE *test;
      strcpy(append,data_ext[i]);
      test = fopen(tmpbuf,"r");
      if (test != NULL) {
         if (result) { // conflict! both mzXML and mzData are present
            if (strcasecmp(tmpbuf,result)) { // win32 isn't case sensitive
               printf("found both %s and %s, using %s\n",
                  tmpbuf,result,result);
            }
         } else { // found file, copy the constructed name
            result = strdup(tmpbuf);
         }
         fclose(test);
      }
   }
   if (!result) { // failed - caller can complain about lack of .mzXML
      strcpy(append,data_ext[0]);
      result = strdup(tmpbuf);
   }
   if (basename_in==inbuf) { // same pointer
      free(basename);
   }
   free(tmpbuf);

   if ((int) strlen(result) < inbuflen) {
      strcpy(inbuf, result);
      free(result);
      result = inbuf;
   } else {
      printf("buffer too small for file %s\n",
         result);
      free(result);
      result = NULL;
   }
   return result;
}

// return NULL if fname is not of extension type we handle,
// otherwise return pointer to .ext
char *rampValidFileType(const char *fname) {
   char *result;
   int i;
   for (i = N_EXT_TYPES;i--;) {
      result = strrchr(fname,'.');
      if (result && !strcasecmp(result,data_ext[i])) {
         break;
      }
      result = NULL; // try again
   }
   return result;
}

// remove the filename .ext, if found
// return NULL if no .ext found, else return fname
char *rampTrimBaseName(char *fname) {
   char *ext = rampValidFileType(fname);
   if (ext) {
      *ext = 0; // trim the extension
   }
   return ext?fname:NULL;
}


/****************************************************************
 * READS the number of peaks.			*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/
int readPeaksCount(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   char *stringBuf=(char *)malloc(SIZE_BUF+1);
   char *beginPeaksCount, *peaks;
   int result = 0;
   const char *tag = pFI->bIsMzData?"length=":"peaksCount=";
   ramp_fileoffset_t in_lScanIndex = lScanIndex;

   ramp_fseek(pFI, lScanIndex, SEEK_SET);

   // Get the num of peaks in the scan and allocate the space we need
   ramp_nextTag(stringBuf, SIZE_BUF, pFI);
   while (!(beginPeaksCount = (char *) strstr(stringBuf, tag)))
   {
      lScanIndex = ramp_ftell(pFI);
      ramp_nextTag(stringBuf, SIZE_BUF, pFI);
   }

   // We need to move forward the length of the tag
   beginPeaksCount += (strlen(tag)+1);
   result = atoi(beginPeaksCount);

   // mext call is probably to read the <peaks> section, position there if needed
   if (pFI->bIsMzData) {
      ramp_fseek(pFI, in_lScanIndex, SEEK_SET);
   } else {
      peaks = strstr(stringBuf,"<peaks");
      if (peaks) {
         ramp_fseek(pFI, lScanIndex+(peaks-stringBuf), SEEK_SET);
      }
   }
   free(stringBuf);
   return result;
}


/****************************************************************
 * READS the base64 encoded list of peaks.			*
 * Return a RAMPREAL* that becomes property of the caller!		*
 * The list is terminated by -1					*
 * !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
 *    RETURNING !!						*
 ***************************************************************/

RAMPREAL *readPeaks(RAMPFILE *pFI,
      ramp_fileoffset_t lScanIndex)
{
   int  n=0;
   int  peaksCount=0;
   int  peaksLen;       // The length of the base64 section
   int precision = 0;
   RAMPREAL *pPeaks = NULL;
   RAMPREAL *pPeaksDeRuled = NULL;
   
   int  endtest = 1;
   int  weAreLittleEndian = *((char *)&endtest);

   char *pData = NULL;
   char *pBeginData;
   char *pDecoded = NULL;

   char buf[1000];
   buf[sizeof(buf)-1] = 0;

   if (pFI->bIsMzML) {
      return mzmlReadPeaks(pFI, lScanIndex);
   }

   if (pFI->bIsMzData) {
      // intensity and mz are written in two different arrays
      int bGotInten = 0;
      int bGotMZ = 0;
      ramp_fseek(pFI,lScanIndex,SEEK_SET);
      while ((!(bGotInten && bGotMZ)) &&
           ramp_nextTag(buf,sizeof(buf)-1,pFI)) {
         int isArray = 0;
         int isInten = 0;
         int isLittleEndian = 0;
         int byteOrderOK;
         if (strstr(buf,"<mzArrayBinary")) {
            isArray = bGotMZ = 1;
         } else if (strstr(buf,"<intenArrayBinary")) {
            isArray = isInten = bGotInten = 1;
         }
         if (isArray) {
            int partial,triplets,bytes;
            const char *datastart;
            // now determine peaks count, precision
            while (!(datastart= (char *) strstr(buf, "<data")))
            {
               ramp_nextTag(buf, sizeof(buf)-1, pFI);
            }
            
            // find precision="xx"
            if( !(pBeginData = strstr( buf , "precision=" )))
            {
               precision = 32; // default value
            } else { // we found declaration
               precision = atoi(findquot(pBeginData)+1);
            }
            
            // find length="xx"
            if((!peaksCount) && (pBeginData = strstr( buf , "length=" )))
            {
               peaksCount = atoi(findquot(pBeginData)+1);
            }

            if (peaksCount <= 0)
            { // No peaks in this scan!!
               return NULL;
            }

            // find endian="xx"
            if((pBeginData = strstr( buf , "endian=" )))
            {
               isLittleEndian = !strncmp("little",findquot(pBeginData)+1,6);
            }
            
            // find close of <data>
            while( !(pBeginData = strstr( datastart , ">" )))
            {
               ramp_nextTag(buf, sizeof(buf)-1 , pFI);
               datastart = buf;
            }
            pBeginData++;	// skip the >
            
            // base64 has 4:3 bloat, precision/8 bytes per value
            bytes = (peaksCount*(precision/8));
            // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
            triplets = (bytes/3)+((bytes%3)!=0);
            peaksLen = (4*triplets)+1; // read the "<" from </data> too, to confirm lack of whitespace
            
            if ((pData = (char *) realloc(pData,1 + peaksLen)) == NULL)
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            
            // copy in any partial read of peak data, and complete the read
            strncpy(pData,pBeginData,peaksLen);
            pData[peaksLen] = 0;
            partial = strlen(pData);
            if (partial < peaksLen) {              
               ramp_fread(pData+partial,peaksLen-partial, pFI);
            }

            // whitespace may be present in base64 char stream
            while (pData[peaksLen-1]!='<') {
               char *cp;
               partial = 0;
               // didn't read all the peak info - must be whitespace
               for (cp=pData;*cp;) {
                  if (strchr("\t\n\r ",*cp)) {
                     memmove(cp,cp+1,peaksLen-(partial+cp-pData));
                     partial++;
                  } else {
                     cp++;
                  }
               }
               if (!ramp_fread(pData+peaksLen-partial,partial, pFI)) {
                  break;
               }
            }
            pData[peaksLen-1] = 0; // pure base64 now
                        
            if ((pDecoded = (char *) realloc(pDecoded,peaksCount * (precision/8) + 1)) == NULL)
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            // Base64 decoding
            b64_decode_mio(pDecoded, pData);
            
            if ((!pPeaks) && ((pPeaks = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL))
            {
               printf("Cannot allocate memory\n");
               return NULL;
            }
            
            // And byte order correction
            byteOrderOK = (isLittleEndian==weAreLittleEndian);
            if (32==precision) { // floats
               if (byteOrderOK) {
                  float *f = (float *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     pPeaks[isInten+(2*n)] = (RAMPREAL)*f++;
                  }
               } else {
                  uint32_t *u = (uint32_t *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     uint32_t tmp = swapbytes( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) (* (float *)(&tmp));
                  }
               }
            } else { // doubles
               if (byteOrderOK) {
                  double *d = (double *)pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     pPeaks[isInten+(2*n)] = (RAMPREAL)*d++;
                  }
               } else {
                  uint64_t *u = (uint64_t *) pDecoded;
                  for (n = 0; n < peaksCount; n++) {
                     uint64_t tmp = swapbytes64( *u++ );
                     pPeaks[isInten+(2*n)] = (RAMPREAL) (* (double *)(&tmp));
                  }
               }
            }
            if (bGotInten && bGotMZ) {
               break;
            }
         } // end if isArray
      } // end while we haven't got both inten and mz
      free(pData);
      free(pDecoded);
      pPeaks[peaksCount*2] = -1; // some callers want a terminator
   } else { // mzXML
       Byte *pUncompr;
       int isCompressed = 0;
      int partial,bytes,triplets;
      int isLittleEndian = 0; // default is network byte order (Big endian)
      int byteOrderOK;
      int       compressedLen = 0;
      int       decodedSize;
      char      *pToBeCorrected;
      e_contentType contType = mzInt; // default to m/z-int

      peaksCount = readPeaksCount(pFI, lScanIndex);
      if (peaksCount <= 0)
      { // No peaks in this scan!!
         return NULL;
      }
      
      // now determine peaks precision
      ramp_fgets(buf, sizeof(buf)-1, pFI);
      while (!(pBeginData = (char *) strstr(buf, "peaks")))
      {
         ramp_fgets(buf, sizeof(buf)-1, pFI);
      }
      getIsLittleEndian(buf,&isLittleEndian);

          // TODO ALL OF THE FOLLOWING CHECKS ASSUME THAT THE NAME AND THE VALUE OF THE
          // ATTRIBUTE ARE PRESENT AT THE SAME TIME IN THE BUFFER.
          // ADD A CHECK FOR THAT!
      while( 1 )
      { // Untill the end of the peaks element
          if( (pBeginData = strstr( buf , "precision=" )))
          { // read the precision attribute
              precision = atoi(strchr(pBeginData,'\"')+1);
          }
          if( (pBeginData = strstr( buf , "contentType=" )))
          { // read the contentType attribute
                  // we are only supporting m/z-int for the moment > return if it is something else
                  // TODO add support for the other content types
              if( (pBeginData = strstr( buf , "m/z-int" )))
              {
                  contType = mzInt;
              }
              else if( (pBeginData = strstr( buf , "m/z ruler" )))
              {
                  contType = mzRuler;
              }
              else
              {
                  char* pEndAttrValue;
                  pEndAttrValue = strchr( pBeginData + strlen( "contentType=\"") + 1 , '\"' );
                  pEndAttrValue  = '\0';
                  fprintf(stderr, "%s Unsupported content type\n" , pBeginData ); 
                  return NULL;
              }
          }
          if( (pBeginData = strstr( buf , "compressionType=" )))
          { // read the compressionType attribute.
              if( (pBeginData = strstr( buf , "zlib" )))
              {
                  isCompressed = 1;
              }
              else if( (pBeginData = strstr( buf , "none" )))
              {
                  isCompressed = 0;
              }
              else
              {
                  char* pEndAttrValue;
                  pEndAttrValue = strchr( pBeginData + strlen( "compressionType=\"") + 1 , '\"' );
                  pEndAttrValue = '\0';
                  fprintf(stderr, "%s Unsupported compression type\n" , pBeginData ); 
                  return NULL;
              }
          }
          if( (pBeginData = strstr( buf , "compressedLen=\"")))
          {
              compressedLen = atoi( pBeginData + strlen( "compressedLen=\"" ) );
          }
          if( !(pBeginData = strstr( buf , ">" )))
          { // There is more to read
              ramp_fgets(buf, sizeof(buf)-1 , pFI);
              getIsLittleEndian(buf,&isLittleEndian);
          }
          else
          {
              pBeginData++;	// skip the >
              break;
          }
      }
      if( !precision )
      { // precision attribute was not defined assume 32 by default
          precision = 32;
      }

      if( isCompressed )
      {
          bytes = compressedLen;
      }
      else
      {
          bytes = (peaksCount*(precision/4));
      }
      
      // base64 has 4:3 bloat, precision/8 bytes per value, 2 values per peak
      // for every 3 bytes base64 emits 4 characters - 1, 2 or 3 byte input emits 4 bytes
      triplets = (bytes/3)+((bytes%3)!=0);
      peaksLen = (4*triplets)+1; // read the "<" from </data> too, to confirm lack of whitespace
      
       if ((pData = (char *) malloc(1 + peaksLen)) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }
      pData[peaksLen] = 0;
      
      // copy in any partial read of peak data, and complete the read
      strncpy(pData,pBeginData,peaksLen);
      partial = strlen(pData);
      if (partial < peaksLen) {
         ramp_fread(pData+partial,peaksLen-partial, pFI);
      }
      // whitespace may be present in base64 char stream
      while (pData[peaksLen-1]!='<') {
         char *cp;
         partial = 0;
         // didn't read all the peak info - must be whitespace
         for (cp=pData;*cp;) {
            if (strchr("\t\n\r ",*cp)) {
               memmove(cp,cp+1,peaksLen+1-(partial+cp-pData));
               partial++;
            } else {
               cp++;
            }
         }
         if (!ramp_fread(pData+peaksLen-partial,partial, pFI)) {
            break;
         }
      }
      if( isCompressed )
      {
          decodedSize = compressedLen + 1;
      }
      else
      {
          // 2 values per peak, precision/8 bytes per value
          decodedSize = peaksCount * (precision/4) + 1;
      }
      pData[peaksLen-1] = 0; // pure base64 now
      
      if ((pDecoded = (char *) malloc( decodedSize )) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }
      // Base64 decoding
      b64_decode_mio(pDecoded, pData);
      free(pData);
      
      if ((pPeaks = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL)
      {
         printf("Cannot allocate memory\n");
         return NULL;
      }

          //Zlib decompression 
      if( isCompressed )
      {
          int err;
//        printf("Decompressing data\n");
          uLong uncomprLen = (peaksCount * precision/4 + 1);
			
          pUncompr = (Byte*)calloc((uInt) uncomprLen , 1);
			
          err = uncompress( pUncompr , &uncomprLen , (const Bytef*)pDecoded , decodedSize );
          free( pDecoded );
          pToBeCorrected = (char *)pUncompr;
      }
      else
      {
          pToBeCorrected = pDecoded;
      }
      
      // And byte order correction
      byteOrderOK = (isLittleEndian==weAreLittleEndian);
      if (32==precision) { // floats
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL) ((float *) pToBeCorrected)[n];
            } 
         } else {
            for (n = 0; n < (2 * peaksCount); n++) {
               uint32_t tmp = swapbytes(((uint32_t *) pToBeCorrected)[n]);
               pPeaks[n] = (RAMPREAL) * (float *)(&tmp);
            } 
         }
      } else { // doubles
         if (byteOrderOK) {
            for (n = 0; n < (2 * peaksCount); n++) {
               pPeaks[n] = (RAMPREAL)((double *) pToBeCorrected)[n];
            }
         } else {
            for (n = 0; n < (2 * peaksCount); n++) {
               uint64_t tmp = swapbytes64((uint64_t) ((uint64_t *) pToBeCorrected)[n]);
               pPeaks[n] = (RAMPREAL) * (double *)(&tmp);
            }
         }
      }

      if( contType == mzRuler )
      { // Convert back from m/z ruler contentType into m/z - int pairs
		RAMPREAL lastMass;
		RAMPREAL  deltaMass;
		  int multiplier;
          int j = 0;

          if ((pPeaksDeRuled = (RAMPREAL *) malloc((peaksCount+1) * 2 * sizeof(RAMPREAL) + 1)) == NULL)
          {
              printf("Cannot allocate memory\n");
              return NULL;
          }
         
          for (n = 0; n < (2 * peaksCount); )
          {
// printf("%f\n" , pPeaks[j] );
              if( (int) pPeaks[j] == -1 )
              { // Change in delta m/z
				++j;
				lastMass = (RAMPREAL) pPeaks[j++];
                deltaMass = pPeaks[j++];
				multiplier = 0;
//printf("%f %f\n" , lastMass , deltaMass );
              }
   		    pPeaksDeRuled[n++] = lastMass + (RAMPREAL) multiplier * deltaMass;
			++multiplier;
            pPeaksDeRuled[n++] = pPeaks[j++];
          }
          
          free(pToBeCorrected);
          pPeaksDeRuled[n] = -1;

          free( pPeaks );
          return (pPeaksDeRuled); // caller must free this pointer
      }
      
      free(pToBeCorrected);
      pPeaks[n] = -1;
   }

   return (pPeaks); // caller must free this pointer
}


/*
 * read just the info available in the msRun element
 */
void readMSRun(RAMPFILE *pFI,
                   struct RunHeaderStruct *runHeader)
{
   char stringBuf[SIZE_BUF+1];
   ramp_fseek(pFI, 0 , SEEK_SET); // rewind
   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   while((!strstr( stringBuf , pFI->bIsMzData?"<mzData":"<msRun" )) && !ramp_feof(pFI))  /* this should not be needed if index offset points to correct location */
   {
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }
   while(!ramp_feof(pFI))
   {
      const char *cp;
      if (NULL != (cp=strstr( stringBuf , pFI->bIsMzData?"spectrumList count":"scanCount" ))) {
         cp = findquot(cp);
         runHeader->scanCount = atoi(cp+1);
      }
      if (NULL != (cp=strstr( stringBuf , "startTime" ))) {
         cp = findquot(cp);
         runHeader->dStartTime = rampReadTime(pFI,cp+1);
      } 
      if (NULL != (cp=strstr( stringBuf , "endTime" ))) {
         cp = findquot(cp);
         runHeader->dEndTime = rampReadTime(pFI,cp+1);
      } 
      if (NULL != (cp=strstr( stringBuf , pFI->bIsMzData?"<spectrumDesc":"<scan" ))) {
         break; /* we're into data territory now */
      } 
      ramp_fgets(stringBuf, SIZE_BUF, pFI);
   }

}

/*
 * walk through each scan to find overall lowMZ, highMZ
 * sets overall start and end times also
 */
void readRunHeader(RAMPFILE *pFI,
                   ramp_fileoffset_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan)
{

   int i;
   struct ScanHeaderStruct scanHeader;
   double startMz = 0.0;
   double endMz = 0.0;

   readHeader(pFI, pScanIndex[1], &scanHeader);

   /*
    * initialize values to first scan
    */
   runHeader->lowMZ = scanHeader.lowMZ;
   runHeader->highMZ = scanHeader.highMZ;
   runHeader->dStartTime = scanHeader.retentionTime;
   runHeader->startMZ = readStartMz( pFI , pScanIndex[1] );
   runHeader->endMZ = readEndMz( pFI , pScanIndex[1] );
   for (i = 2; i <= iLastScan; i++)
   {
      readHeader(pFI, pScanIndex[i], &scanHeader);

      if (scanHeader.lowMZ < runHeader->lowMZ)
         runHeader->lowMZ = scanHeader.lowMZ;
      if (scanHeader.highMZ > runHeader->highMZ)
         runHeader->highMZ = scanHeader.highMZ;
      if( (startMz = readStartMz( pFI , pScanIndex[i] )) < runHeader->startMZ )
	runHeader->startMZ = startMz;
      if( (endMz = readEndMz( pFI , pScanIndex[i] )) > runHeader->endMZ )
	runHeader->endMZ = endMz;   
   }

   runHeader->dEndTime = scanHeader.retentionTime;
}



int setTagValue(const char* text,
   char* storage,
   int maxlen,
   const char* lead)
{
  char* result = NULL;
  char* term = NULL;
  int len = maxlen - 1;
  int leadlen = strlen(lead)+1; // include the opening quote

  result = strstr(text, lead);
  if(result != NULL)
  {
    char tail = *(result+leadlen-1); // nab the quote char (is it single or double quote?)
    term = strchr(result + leadlen, tail);
    if(term != NULL)
    {
      if((int)(strlen(result) - strlen(term) - leadlen) < len)
        len = strlen(result) - strlen(term) - leadlen;

      strncpy(storage, result + leadlen , len);
      storage[len] = 0;
      return 1;
    } // if term
  }
  return 0;
}


InstrumentStruct* getInstrumentStruct(RAMPFILE *pFI)
{
  InstrumentStruct* output = NULL;
  char* result = NULL;
  int found[] = {0, 0, 0, 0, 0};
  char stringBuf[SIZE_BUF+1];

   if ((output = (InstrumentStruct *) calloc(1,sizeof(InstrumentStruct))) == NULL)
   {
      printf("Cannot allocate memory\n");
      return NULL;
   } else {
      const char *cpUnknown="UNKNOWN";
      strncpy(output->analyzer,cpUnknown,sizeof(output->analyzer));
      strncpy(output->detector,cpUnknown,sizeof(output->analyzer));
      strncpy(output->ionisation,cpUnknown,sizeof(output->detector));
      strncpy(output->manufacturer,cpUnknown,sizeof(output->ionisation));
      strncpy(output->model,cpUnknown,sizeof(output->manufacturer));
   }

   ramp_fgets(stringBuf, SIZE_BUF, pFI);

   if (pFI->bIsMzData) {
   } else {
      int isAncient=0;
      while( !strstr( stringBuf , "<msInstrument" ) && 
            !(isAncient=(NULL!=strstr( stringBuf , "<instrument" ))) && 
            !strstr(stringBuf, "<dataProcessing") && 
            !ramp_feof(pFI))  /* this should not be needed if index offset points to correct location */
      {
         ramp_fgets(stringBuf, SIZE_BUF, pFI);
      }
      while(! strstr(stringBuf, isAncient?"</instrument":"</msInstrument") &&  ! strstr(stringBuf, "</dataProcessing") && !ramp_feof(pFI))
      {
        if(! found[0])
        {
          result = strstr(stringBuf, isAncient?"manufacturer=":"<msManufacturer");
          if(result != NULL && setTagValue(result, output->manufacturer, INSTRUMENT_LENGTH, isAncient?"manufacturer=":"value="))
	      found[0] = 1;
        }
        if(! found[1])
        {
           result = strstr(stringBuf, isAncient?"model=":"<msModel");
          if(result != NULL && setTagValue(result, output->model, INSTRUMENT_LENGTH, isAncient?"model=":"value="))
	      found[1] = 1;
        }
        if(! found[2])
        {
          result = strstr(stringBuf, isAncient?"ionisation=":"<msIonisation");
          if(result != NULL && setTagValue(result, output->ionisation, INSTRUMENT_LENGTH, isAncient?"ionisation=":"value="))
	      found[2] = 1;
        }
        if(! found[3])
        {
           result = strstr(stringBuf, isAncient?"msType=":"<msMassAnalyzer");
          if(result != NULL && setTagValue(result, output->analyzer, INSTRUMENT_LENGTH, isAncient?"msType=":"value="))
	      found[3] = 1;
        }
        if(! found[4])
        {
          result = strstr(stringBuf, "<msDetector");
          if(result != NULL && setTagValue(result, output->detector, INSTRUMENT_LENGTH, "value="))
	      found[4] = 1;
        }
        ramp_fgets(stringBuf, SIZE_BUF, pFI);

      } // while
   }

   if(found[0] || found[1] || found[2] || found[3] || found[4])
     return output;

   return NULL; // no data
}

// exercise at least some of the ramp interface - return non-0 on failure
int rampSelfTest(char *filename) { // if filename is non-null we'll exercise reader with it
#define N_TEST_NAME 5
   int result = 0; // assume success
   char buf[256];
   char buf2[256];
   int i;

   char *testname[N_TEST_NAME] = 
   {"foo.bar","foo.mzxml","foo.mzdata","foo.mzXML","foo.mzData"};

   // locate the .mzData or .mzXML extension in the buffer
   // return pointer to extension, or NULL if not found
   for (i=N_TEST_NAME;i--;) {
      result |= (!i) != !rampValidFileType(testname[i]); // 0th one in not a valid file type
   }

   // trim a filename of its .mzData or .mzXML extension
   // return trimmed buffer, or null if no proper .ext found
   for (i=N_TEST_NAME;i--;) {
      strncpy(buf,testname[i],sizeof(buf));
      result |= ((!i) != !rampTrimBaseName(buf));
      if (i) {
         result |= (strcmp(buf,"foo")!=0);
      }
   }

   if (filename && rampValidFileType(filename)) {
      // construct a filename in buf from a basename, adding .mzXML or .mzData as exists
      // returns buf, or NULL if neither .mzXML or .mzData file exists
      char *name;
      strncpy(buf,filename,sizeof(buf));
      rampTrimBaseName(buf);
      name = rampConstructInputFileName(buf,sizeof(buf),buf); // basename is in buf
      result |= (name==NULL);
      strncpy(buf,filename,sizeof(buf));
      rampTrimBaseName(buf);
      name = rampConstructInputFileName(buf2,sizeof(buf2),buf); // different buf, basename
      result |= (name==NULL);
   }
   return result;
}

// Cache support

// Get a new cache instance with a specified window size.  A larger window
// requires more memory, obviously.  Too small a window for the required
// function can lead perf comparable to no caching at all.

// Iterating over a range of 200 scans with a cache that contains 100 or
// fewer scans will yield no cache hits on each iteration.  Pick a cache
// size slightly larger than the window you expect to cover.

struct ScanCacheStruct *getScanCache(int size)
{
    struct ScanCacheStruct* cache = (struct ScanCacheStruct*) malloc(sizeof(struct ScanCacheStruct));
    cache->seqNumStart = 0;
    cache->size = size;
    cache->headers = (struct ScanHeaderStruct*) calloc(size, sizeof(struct ScanHeaderStruct));
    cache->peaks = (RAMPREAL**) calloc(size, sizeof(RAMPREAL*));
    return cache;
}

// Free all memory associated with a cache struct.
void freeScanCache(struct ScanCacheStruct* cache)
{
   if (cache) {
    int i;
    for (i = 0; i < cache->size; i++)
    {
        if (cache->peaks[i] != NULL)
            free(cache->peaks[i]);
    }
    free(cache->peaks);
    free(cache->headers);
    free(cache);
}
}

// Clear all cached values, freeing peaks, but not the cache arrays themselves.
void clearScanCache(struct ScanCacheStruct* cache)
{
    int i;
    for (i = 0; i < cache->size; i++)
    {
        if (cache->peaks[i] == NULL)
            continue;

        free(cache->peaks[i]);
        cache->peaks[i] = NULL;
    }
    memset(cache->headers, 0, cache->size * sizeof(struct ScanHeaderStruct));
}

// Shift the cache start index by a number of scans.  This moves the cache
// window left (negative) or right (positive).
void shiftScanCache(struct ScanCacheStruct* cache, int nScans)
{
    int i;
    cache->seqNumStart += nScans;
    if (abs(nScans) > cache->size)
    {
        // If the shift is larger than the size of the cache window,
        // just clear the whole cache.
        clearScanCache(cache);
    }
    else if (nScans > 0)
    {
        // Shifting window to the right.  Memory moves right, with new
        // empty scans on the end.

        // Free the peaks that memmove will overwrite.
        for (i = 0; i < nScans; i++)
        {
            if (cache->peaks[i] != NULL)
                free(cache->peaks[i]);
        }
        memmove(cache->peaks, cache->peaks + nScans,
            (cache->size - nScans) * sizeof(RAMPREAL*));
        memset(cache->peaks + cache->size - nScans, 0, nScans * sizeof(RAMPREAL*));
        memmove(cache->headers, cache->headers + nScans,
            (cache->size - nScans) * sizeof(struct ScanHeaderStruct));
        memset(cache->headers + cache->size - nScans, 0, nScans * sizeof(struct ScanHeaderStruct));
    }
    else if (nScans < 0)
    {
        // Shifting window to the left.  Memory moves right, with new
        // empty scans at the beginning.
        nScans = -nScans;

        // Free the peaks that memmove will overwrite.
        for (i = 0; i < nScans; i++)
        {
            if (cache->peaks[cache->size - 1 - i] != NULL)
                free(cache->peaks[cache->size - 1 - i]);
        }
        memmove(cache->peaks + nScans, cache->peaks,
            (cache->size - nScans) * sizeof(RAMPREAL*));
        memset(cache->peaks, 0, nScans * sizeof(RAMPREAL*));
        memmove(cache->headers  + nScans, cache->headers,
            (cache->size - nScans) * sizeof(struct ScanHeaderStruct));
        memset(cache->headers, 0, nScans * sizeof(struct ScanHeaderStruct));
    }
}

// Convert a scan index into a cache index, adjusting the cache window
// if necessary.
int getCacheIndex(struct ScanCacheStruct* cache, int seqNum)
{
    int seqNumStart = cache->seqNumStart;
    int size = cache->size;

    // First access, just set the start to seqNum.
    if (seqNumStart == 0)
        cache->seqNumStart = seqNum;
    // If requested scan is less than cache start, shift cache window
    // left to start at requested scan.
    else if (seqNum < seqNumStart)
        shiftScanCache(cache, (int) (seqNum - seqNumStart));
    // If requested scan is greater than cache end, shift cache window
    // right so last entry is requested scan.
    else if (seqNum >= seqNumStart + size)
        shiftScanCache(cache, (int) (seqNum - (seqNumStart + size - 1)));

    return (int) (seqNum - cache->seqNumStart);
}

struct ScanHeaderStruct* readHeaderCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    int i = getCacheIndex(cache, seqNum);
    if (cache->headers[i].msLevel == 0)
        readHeader(pFI, lScanIndex, cache->headers + i);
    return cache->headers + i;
}

int readMsLevelCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    struct ScanHeaderStruct* header = readHeaderCached(cache, seqNum, pFI, lScanIndex);
    return header->msLevel;
}

RAMPREAL *readPeaksCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex)
{
    int i = getCacheIndex(cache, seqNum);
    if (cache->peaks[i] == NULL)
        cache->peaks[i] = readPeaks(pFI, lScanIndex);
    return cache->peaks[i];
}

