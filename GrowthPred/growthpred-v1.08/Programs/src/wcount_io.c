/* ==================================================== */
/* Copyright (c) Atelier de BioInformatique		*/
/* @file: wCount_io.c					*/
/* util and IO functions				*/
/* history:						*/
/* April 97 <ed> first draft from Covgen		*/
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Gtypes.h"
#include "libwcount.h"

#if PROTO
typedef int (*QsortCompFunc)(const void *, const void *);
#else
typedef int (*QsortCompFunc)();
#endif

#define SKIP(buffer)	(*(buffer) == '/') || (*(buffer) == '#')

/* -------------------------------------------- */
/* open a file and assign it to stdin		*/
/* -------------------------------------------- */

static int AssignToStdin(filename)
	char *filename;
{
	char buffer[256];
	
	if (! freopen(filename, "r", stdin)) {
	   sprintf(buffer, "cant open file %s to read", filename);
	   return wCountError(buffer, 0);
	}

	return 1;
}
/* -------------------------------------------- */
/* clean a fgets buffer				*/
/* -------------------------------------------- */

static char *sCleanBuffer(char *buffer)
{
    int len;
    if (buffer) {
	len = strlen(buffer);
	if (buffer[len-1] == '\n')
	   buffer[len-1] = '\000';
    }
    return buffer;
}

/* -------------------------------------------- */
/* check if file can be accessed		*/
/* unix access() does not exist on Mac		*/
/* -------------------------------------------- */

int wCountCheckFile(char *filename, int toread)
{
	FILE *stream;
		
	if (! (stream = fopen(filename, (toread ? "r" : "w"))))
	  return 0;
	   
	fclose(stream);
	
	return 1;
}

/* -------------------------------------------- */
/* get # of lines starting with char in file	*/
/* -------------------------------------------- */

Int32 wCountGetNbLines(char *filename,int achar)
{
	Int32 count;
	FILE  *stream;
	char  buffer[1024];
	
	if (! (stream = fopen(filename, "r")))
	    return 0;
	
	count = 0;
	
	while (fgets(buffer, sizeof(buffer), stream)) {
	
	    if ((*buffer == '/') || (*buffer == '#'))
	    	continue;
	
	    if ((achar == 0) || (*buffer == achar))
	    	count++;
	}
	
	(void) fclose(stream);
		    	
	return count;
}


/* -------------------------------------------- */
/* uppercase sequence				*/
/* -------------------------------------------- */
char *wCountUpper(char *seq)
{
	char *s;
	
	for (s = seq; *s ; s++)
	   if (islower(*s))
	   	*s = toupper(*s);
	   	
	return seq;
}


/* -------------------------------------------- */
/* read restricted alphabet			*/
/* -------------------------------------------- */

int wCountReadAlphabet(char *filename, char *alphabet, int *count)
{
	int  icount, err;
	FILE *stream;
	char buffer[256], alpha[32];
	
	if (! (stream = fopen(filename, "r")))
	    return 0;
	err = 0;
	
	while (fgets(buffer, sizeof(buffer), stream)) {
	
	    if (SKIP(buffer))
	    	continue;
	
	    if (sscanf(buffer, "%1s%d", alpha, &icount) != 2)
		return err++;

	    if (islower(*alpha))
	        *alpha = toupper(*alpha);
	        
	    if (! isupper(*alpha))
	    	continue;
	    
	    *alphabet++ = *alpha;
	    
	    count[*alpha - 'A'] = ((icount<0) ? MAX_PAT_LEN : icount);
	}
	
	(void) fclose(stream);
	
	*alphabet = '\000';	    
		
	return (err ? 0 : 1);
}

/* -------------------------------------------- */
/* read error matrix				*/
/* -------------------------------------------- */

int wCountReadErrorMatrix(char *filename, SErrMat *errmat)
{
	int  i, j, in, err;
	FILE *stream;
	char *c, buffer[256];
	
	if (! (stream = fopen(filename, "r")))
	    return 0;
	    
	i = err = 0;
	
	while (fgets(buffer, sizeof(buffer), stream)) {
	
	    if (SKIP(buffer))
	    	continue;

	    for (c = buffer, in = j = 0 ; *c ; c++) {
	        
		if (*c == '/')
		    in = (in ? 0 : 1);    
		if (in)
		   continue;
		if (isdigit(*c)) {
		   errmat->allow[i][j] = ((*c == '0') ? 0 : 1);
		   j++;
		}
	    }	
	
	    if (j != ALPHA_LEN) 
	        return 0;
	    	
	    i++;
	}
	
	(void) fclose(stream);

	if (i != ALPHA_LEN)
		return 0;
		
	errmat->active = (err == 0);
	
	return (err ? 0 : 1);
}

/* -------------------------------------------- */
/* read generator				*/
/* -------------------------------------------- */

#define IS_RESTRICT(buf) (! strncmp(buf, "restrict:", 9))

int wCountReadGenerator(char *filename, SGenerator *gen)
{
	int  i, count, npos, err;
	FILE *stream;
	char buffer[256], alpha[32], skip[32];
	
	if (! (stream = fopen(filename, "r")))
	    return 0;
	   
	err = npos = 0;
	
	while (fgets(buffer, sizeof(buffer), stream)) {
	
	    (void) sCleanBuffer(buffer);
	    
	    if (SKIP(buffer))
	    	continue;
	    if (IS_RESTRICT(buffer))
	        break;
	    if (npos >= MAX_PAT_LEN)
		return 0;
	    	        
	    strcpy(gen->alpha[npos], buffer);
	
	    npos++;
	}

	gen->maxlen = npos;
	
	for (i = 0 ; i < ALPHA_LEN ; i++)
	    gen->count[i] = MAX_PAT_LEN;
	
	do {
	
	    if (! IS_RESTRICT(buffer))
	        continue;
	
	    if (sscanf(buffer, "%s%1s%d", skip, alpha, &count) != 3)
	    	return 0;

	    if (islower(*alpha))
	        *alpha = toupper(*alpha);
	        
	    if (! isupper(*alpha))
	    	return 0;
	    	    
	    gen->count[*alpha - 'A'] = (count < 0 ? MAX_PAT_LEN : count);
	
	} while (fgets(buffer, sizeof(buffer), stream));

	(void) fclose(stream);
	
	return 1;
}

#undef IS_RESTRICT

/* -------------------------------------------- */
/* clean user's oligo				*/
/* -------------------------------------------- */

static char *sCleanOligo(char *buf)
{
	char *c, *cc;
	
	for (c = cc = buf ; *c ; c++)
	    if (!isspace(*c))
	        *cc++ = *c;

	*cc = '\000';

	return buf;
}

/* -------------------------------------------- */
/* check user's oligo				*/
/* -------------------------------------------- */

static int sCheckOligo(char *buf)
{
	char *c;
	
	for (c = buf ; *c ; c++)
	    if (! isupper(*c))
		return 0;
		
	return 1;
}

/* -------------------------------------------- */
/* read user's oligo file			*/
/* -------------------------------------------- */

int wCountReadOligos(char *filename, SPattern *pat, Int32 npat,
    int nberr, Int32 overl)
{
	Int32    ipat;
	FILE     *stream;
	SPattern *p;
	char     buffer[256];

	if (! (stream = fopen(filename, "r")))
	    return -1;

	ipat = 0;
	p    = pat;

	while (fgets(buffer, sizeof(buffer), stream)) {
	
	    (void) sCleanBuffer(buffer);
	    
	    if (SKIP(buffer))
	    	continue;

	    (void) wCountUpper(buffer);
	    
	    p->patlen  = strlen(buffer);
	    p->maxerr  = nberr;
	    p->overlap = overl;

	    if (! sCheckOligo(sCleanOligo(buffer)))
	    	return -1;

	    if (p->patlen > MAX_PAT_LEN)
	    	return -2;
	    
	    strcpy(p->pat, buffer);

	    if (! CompileSMat(buffer, p))
	    	return -3;
	    p++;
	    ipat++;
	}
	
	(void) fclose(stream);

	if (ipat != npat)
	    return -1;
	
	return 1;
}

/* -------------------------------------------- */
/* get user oligo's from string			*/
/* -------------------------------------------- */


int wCountGetOligos(char *oligos, SPattern *pat, Int32 npat,
    int nberr, Int32 overl)
{
    Int32	    ipat, j, kmax;
    SPattern	    *p;
    char	    buffer[33], *excl;

    ipat = 0;
    p    = pat;
    kmax = 0;

    for (ipat=0; ipat<npat; ipat++, p++) {
    	    
		for (j=kmax; oligos[j]!= SEPARATOR; j++)
	    	;
	
		p->patlen  = (j - kmax );
		p->maxerr  = nberr;
		p->overlap = overl;
    
		strncpy(buffer, &(oligos[kmax]), p->patlen);
		buffer[p->patlen]='\000';

	    (void) wCountUpper(buffer);
	    if (! sCheckOligo(sCleanOligo(buffer)))
	    	return -1;
	    if (p->patlen > MAX_PAT_LEN)
	    	return -2;

	    strcpy(p->pat, buffer);

		kmax += 1 + (p->patlen);

		if (! CompileSMat(buffer, p))
	    	return -3;
    }

	if (ipat != npat)
	    return -1;
    
    return 1;
}

/* -------------------------------------------- */
/* read sequence file				*/
/* -------------------------------------------- */

int wCountReadSequences(char *filename, FastaSequencePtr *seq, 
	int nseq)
{
	Int32 i;

	if (! AssignToStdin(filename))
	    return 0;
	
	for (i = 0 ; i < nseq ; i++) {
	
	   if (! (seq[i] = NewFastaSequence()))
	   	return wCountError("not enough memory (NewFastaSequence)", 0);
	   
	   if (! ReadFastaSequence(stdin, seq[i]))
	   	return 0;
	
	    (void) wCountUpper(seq[i]->seq);

	}

	return 1;
}

/* -------------------------------------------- */
/* print out results for each oligo in heap		*/
/* -------------------------------------------- */

int PrintResPerOligo(CHeap *heap, int totlen, char *id, int count, int lab)
{
	int i, res;

				/* print out in order ------------------ */

	if (count==0){
		printf("/ gene\t");
	
		for (i = 0 ; i < heap->nbbest ; i++) 
	    	printf("#%s\t\%%s\t", heap->tree[i].pat, heap->tree[i].pat);
		
		printf("\n/ --------------------------------------\n");
	}
	
	printf("%s\t", id);
	
	for (i = 0 ; i < heap->nbbest ; i++) {
	
		res= lab * heap->tree[i].score;
		
	    printf("%ld\t%.2f\t", res, (float) 100 * res / totlen);
	}
	
	printf("\n");
	
	return 1;
}


/* -------------------------------------------- */
/* print out results for oligo in heap		*/
/* -------------------------------------------- */

int PrintResOligos(CHeap *heap, int totlen, int lab)
{
	int i, res;

				/* print out in order ------------------ */

	printf("/ oligo\t#\t%%\n");
	printf("/ --------------------------------------\n");
		
	for (i = 0 ; i < heap->nbbest ; i++) {
	
		res= lab * heap->tree[i].score;
		
	    printf("%s\t%ld\t%.2f\n", 
	   		heap->tree[i].pat,
	   		res, (float) 100 * res / totlen);
	}
	
	return 1;
}


/* ------------------------------------------------------------ */
/* function for scoring windows of cumulative counts		*/
/* ------------------------------------------------------------ */

int getWinOligo(FastaSequencePtr *seqdata, Int32 nseq, 
    Int32 winlen, Int32 step)
{
    Int32   i, j, k, count, nwin, res;
    char    *start, *end, *cover; 
    
    for (k=0; k<nseq; k++) {
    
		cover = (seqdata[k]->cover);
		nwin = floor((seqdata[k]->length - winlen)/step)+1;

		printf("/ #cov\t%%\n");
		printf("/ --------------------------------------\n");

				/* first window*/
		for (i=res=0; i<winlen;i++){
	    	if (cover[i]==1)
				res++;
		}

	   printf("%s\t%ld\t%ld\t%.2f\n", 
	   		res, (float) res / (float) winlen * 100.);
    
				/* setting markers*/
		start=&(cover[0]);
		end=&(cover[winlen-1]);
		
    			/* sliding window */
		for (i=1; i<nwin; i++){
	    	for (j=0, count=0; j<step; j++){
				if (*(start++) == 1)
		    		count--;
				if (*(end++) == 1)
		    		count++;
	    	}
	    	res = res + count;
	   		printf("%s\t%ld\t%ld\t%.2f\n", 
	   			res, (float) res / (float) winlen * 100.);
 	    }
    }	
    return 0;
}
