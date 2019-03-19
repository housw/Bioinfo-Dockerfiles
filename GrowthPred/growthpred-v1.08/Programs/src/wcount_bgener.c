/* ==================================================== */
/* Copyright (c) Atelier de BioInformatique		*/
/* @file: wCount_bgener.c				*/
/* oligo generator - brute force			*/
/* history:						*/
/* April 1997 <ed> modified from Covgen			*/
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _H_Gtypes
#include "Gtypes.h"
#endif

#ifndef _H_libwcount
#include "libwcount.h"
#endif

static Int32 sOligoNum = 0;
static Int32 sOligoTot = 0;

/* -------------------------------------------- */
/* eval oligo on sequences			*/
/* -------------------------------------------- */

static void sEvalOligo(char *oligo, CHeap *heap,
    FastaSequencePtr *seqdata, Int32 nseq, Int32 overl, int nerr, SErrMat *errmat,
    int coverflag, Int32 jump, Int32 startat, int overflag, int worstflag)
{
	Int32 count;

	count = wCountScanOneOligoAllSeq(oligo, seqdata, nseq, overl, nerr,
	    errmat, Vrai, coverflag, jump, overflag, startat);

	if (worstflag)
	    count=count*(-1);

	if (count > heap->tree->score)
	    wCountUpdateHeap(heap, sOligoNum, count, oligo);
}
	

/* -------------------------------------------- */
/* generation oligomeres			*/
/* -------------------------------------------- */

static void sRecursOligo(char *oligo, int size, int *count, SGenerator *gen,
	CHeap *heap, FastaSequencePtr *seqdata, Int32 nseq,  
	Int32 overl,int nerr, SErrMat *errmat, int coverflag, Int32 jump, 
	Int32 startat, int overflag, int worstflag) 
{
	int   i;
	char *c;
	
	if (size == gen->maxlen) {	/* end of recursion 	*/

	    oligo[size] = '\000';
	    sEvalOligo(oligo, heap, seqdata, nseq, overl, nerr,
		errmat, coverflag, jump, startat, overflag, worstflag);
	}

	else {				/* body of recursion	*/
	
	   for (c = gen->alpha[size] ; *c ; c++) {
	   
	       oligo[size] = *c;
	       
	       i = *c - 'A';
	       
	       count[i]++;
	       
	       if (count[i] <= gen->count[i])
	           sRecursOligo(oligo, size+1, count, gen, heap,
			seqdata, nseq, overl, nerr, errmat, coverflag,
			jump, startat, overflag, worstflag);
	           
	       count[i]--;
	   }
	}
}

/* -------------------------------------------- */

void wCountBruteGenereOligo(SGenerator *gen, CHeap *heap,
    FastaSequencePtr *seqdata, Int32 nseq, Int32 overl, int nerr, SErrMat *errmat,
    int coverflag, Int32 jump, int overflag, Int32 startat, int worstflag) 
{
	int	i, count[ALPHA_LEN];
	char	oligo[MAX_PAT_LEN+1];

	for (i = 0 ; i < ALPHA_LEN ; i++)
	    count[i] = 0;

	*oligo = '\000';
	
	sOligoNum = 0;
	sOligoTot = wCountNumberGen(gen);
	
	sRecursOligo(oligo, 0, count, gen, heap, seqdata, nseq,  
	    overl, nerr, errmat, coverflag, jump, startat, overflag, worstflag);
}

/* -------------------------------------------- */
/* default generator				*/
/* -------------------------------------------- */

int wCountDefaultGenerator(SGenerator *gen, char *alphabet, 
	int kmax, int *count)
{
	int  i;
	char *c;
	
	gen->maxlen = kmax;
	
	for (i = 0 ; i < kmax ; i++) {	
	    strcpy(gen->alpha[i], alphabet);
	    }

	for (i = 0 ; i < ALPHA_LEN ; i++) {	    
	    gen->count[i] = MAX_PAT_LEN;
	    }
	
	for (c = alphabet ; *c ; c++) {
	    i  = *c - 'A'; 
	    gen->count[i] = (count[i] < 0 ? MAX_PAT_LEN : count[i]);
	}
	
	return 1;
}

/* -------------------------------------------- */
/* compute # of oligo to generate		*/ 
/* ! don't take restrict: in account !		*/
/* -------------------------------------------- */

int wCountNumberGen(SGenerator *gen)
{
	int   i, n;
	
	for (i = 0 , n = 1 ; i < gen->maxlen ; i++) 
	   n *= strlen(gen->alpha[i]);
	    
	return n;
}

/* ----------------------------------------------------- */
/* Scan Stuff						 */
/* ----------------------------------------------------- */
/* -------------------------------------------- */
/* comput one sequence coverage			*/
/* -------------------------------------------- */

static Int32 SequenceCoverage(FastaSequencePtr *seqdata, long startat, Int32 nseq, 
    Int32 jump, int overflag, int patlen)
{
    Int32 i, j, count;
    char  *cov;

    for (j=count=0; j<nseq; j++){
	
	if (overflag){		    /* overlapping */
	    for (i = 0, cov = (seqdata[j]->cover) ; 
		i < (seqdata[j]->length)-startat ; i++)
		if (*cov++==1)
		    count++; 
	}
	else {			    /* not-overlapping */
	    for (i = 0, cov = (seqdata[j]->cover) ; 
		i < (seqdata[j]->length)-startat ; i++) {
		if (*cov++==1) {
		    count++;
		    i+=patlen;
		}
	    }
	}
    }
    return count;
}

/* -------------------------------------------- */
/* scan patterns				*/
/* and keep best in heap			*/
/* -------------------------------------------- */

void wCountScanPatterns(SPattern *pat, FastaSequencePtr *seqdata, Int32 nseq, 
    CHeap *heap, Int32 npat, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag, Int32 startat, int worstflag, Bool reset)
{
    Int32 i, j, coverscr, ok;
    char *seq, *cover;

    for (i = 0 ; i < npat ; i++) {
	
	for (j=ok=0; j<nseq; j++) {
	
	    seq = &(seqdata[j]->seq[startat]);
	    cover = (seqdata[j]->cover);
	    
	    if (reset||(i==0))
			memset(seqdata[j]->cover, 2, seqdata[j]->length);

	    ok += ManberSearch_S(seq, cover, pat + i, errmat,
		coverflag, jump, overflag);	
	}
	
	pat[i].hit = (ok ? 1 : 0);

	if (worstflag)	    
	    coverscr = -SequenceCoverage(seqdata, startat, nseq, jump, 
		overflag, pat->patlen);	
	else 
	    coverscr = SequenceCoverage(seqdata, startat, nseq, jump, 
		overflag, pat->patlen);

	if (coverscr > heap->tree->score)
	        wCountUpdateHeap(heap, i, coverscr, pat[i].pat);
    }
}


/* -------------------------------------------- */
/* scan one oligo all sequences			*/
/* -------------------------------------------- */

Int32 wCountScanOneOligoAllSeq(PatString cpat, FastaSequencePtr *seqdata,
    Int32 nseq, Int32 overl, Int32 nerr, SErrMat *errmat, Bool reset,
    int coverflag, Int32 jump, int overflag, Int32 startat)
{
    SPattern pat;
    Int32 i;
    char *seq, *cov;

    pat.maxerr  = nerr;
    pat.overlap = overl;
    pat.patlen  = strlen(cpat);

    strcpy(pat.pat, cpat);

    if (! CompileSMat(cpat, &pat))
       return wCountError("<internal> (Compile:ScanOneOligoAllSeq)", 0);

    for (i=0; i<nseq; i++) {
	
	cov = seqdata[i]->cover;
	seq = &(seqdata[i]->seq[startat]);
	
	if (reset)
	    memset(seqdata[i]->cover, 2, seqdata[i]->length);
	 
	(void) ManberSearch_S(seq, cov, &pat, errmat, coverflag, 
	    jump, overflag);
    }
	
    return SequenceCoverage(seqdata, startat, nseq, jump, overflag, pat.patlen);
}

/* ----------------------------------------------------- */
/* Search Stuff						 */
/* ----------------------------------------------------- */

#define IUPAC_DUAL
					
static long sDnaCode[]  =  {
#include "dna_code.h"
};

/* -------------------------------------------- */
/* test the pattern 				*/
/* -------------------------------------------- */

static Bool sCheckHit(char *frag, char *pat,int len, SErrMat * errmat)
{
    if (! errmat->active)
        return Vrai;
	    
    for (frag = frag - len + 1 ; *pat ; pat++, frag++) {
        if (errmat->allow[*frag-'A'][*pat-'A'] == 0)
            return Faux;
    }
	
    return Vrai;
}


/* -------------------------------------------- */
/* save the pattern     			*/
/* -------------------------------------------- */

static void sRecordHit(char *seq, char *sfnd, char *cover, Int32 overl,
    Int32 plen, int coverflag, Int32 jump, int overflag)
{
    Int32 pos, left, right, span, i, out=0;

    pos = sfnd - seq - plen + 1;	/* start pos of hit	*/

    if (pos < 0)			/* special case !	*/
	return;			/* first word		*/
	    
    left  = pos - overl;			/* leftmost	*/
    right = pos + overl + plen - 1;	 	/* rightmost	*/

    if (left < 0)
	left = 0;

    cover += left;
    seq   += left;
	
    span  = right - left + 1;
    
    if (!overflag)		/*tests overlap if !overflag*/
	{
	for(i=1; i<plen; i++)
	    if (*(cover-i)==1)
		out=1;
	if(out==1)
	    return;
	}
    
					/*tests phase */
    if (pos % jump != 0){
	return;	
	}
    					/* increments cover as coverage or count */
    if(coverflag) {  
	while ((span--) && (*seq++))
	    (*cover++) = 1;
	}
    else 
	*cover=1;

}	


/* -------------------------------------------- */
/* building S Matrix				*/
/* -------------------------------------------- */

int CompileSMat(char *cpat,SPattern *ppat)
{
	int 	i, j;
	long	code; 
	UInt32	pindx, amask, *smat;


	ppat->patlen = strlen(cpat);

	for (i = 0 ; i < ALPHA_LEN ; i++)
	   ppat->smat[i] = 0x0;

	for (i = ppat->patlen-1, amask = 0x1L ; i >= 0 ; i--, amask <<= 1) {

	    code = sDnaCode[cpat[i] - 'A'];

	    for (j = 0, pindx = 0x1L ; j < ALPHA_LEN ; j++, pindx <<= 1)
		if (code & pindx)
		    ppat-> smat[j] |= amask;
	}

	return 1;
}

	
/* -------------------------------------------- */
/* pattern search	 			*/
/* Substitutions only				*/
/* note: there should be a special case for a	*/
/*       k hit on the first k-1 letters of seq.	*/
/*       it is handled in sRecordHit	*/
/* -------------------------------------------- */

Int32 ManberSearch_S(char *seq, char *cover, SPattern *ppat,
    SErrMat *errmat, int coverflag, Int32 jump, int overflag)
{
    int		e, emax;
    Int32	count, pos;
    Ulong	smask, cmask, sindx;
    Ulong	*pr, r[2 * MAX_PAT_ERR + 2];
    char	*s, *ls;

					/* raz local r */
    emax = ppat->maxerr;
    r[0] = r[1] = r[2] = 0x0;
    cmask = smask = 0x1L << ppat->patlen;

    for (e = 0, pr = r + 3 ; e <= emax ; e++, pr += 2) {
	*pr = cmask;
	cmask = (cmask >> 1) | smask;
    }

				/* loop on text */
    count = 0;
    pos=0;
    
    for (ls = s = seq ; *s ; s++, pos++) {

	sindx  = ppat->smat[*s - 'A'];

        for (e = 0, pr = r ; e <= emax ; e++, pr += 2) {
	   	
	    pr[2]=  pr[3] | smask;
	    pr[3]=  (pr[0] >> 1)  		/* sub	 */
		    | ((pr[2] >> 1) & sindx);	/* ident */
	    
	    if ((pr[3] & 0x1L) && 
		sCheckHit(s, ppat->pat, ppat->patlen, errmat)) {
		
		sRecordHit(seq, s, cover, ppat->overlap, ppat->patlen, 
		    coverflag, jump, overflag);
		           
		if (ls != s)
	    	    count++;    	   
		ls = s;
	    }
	}
    }
	return count;
}

/* ------------------------------------------------------------ */
/* heap tree management						*/
/* ------------------------------------------------------------ */

/* -------------------------------------------- */
/* destroy tree					*/
/* -------------------------------------------- */

CHeap *wCountFreeHeap(CHeap *heap)
{
	if (heap) {
	    if (heap->tree)
		FREE(heap->tree);
	    FREE(heap);
	}
	
	return NULL;
}

/* -------------------------------------------- */
/* initialize a heap				*/
/* -------------------------------------------- */

CHeap *wCountInitHeap(CHeap *heap, long initheap)
{
	Int32 	i;
	Ctree 	*arbor;
	
	arbor = heap->tree;
	
	if (! arbor) 
	   return wCountFreeHeap(heap);
	   
	for (i = 0 ; i < heap->nbnodes ; i++, arbor++) {
	    arbor->score  = arbor->indx = initheap;
	    arbor->pat[0] = '\000';
	    arbor->hit    = 0;
	}
	    
	return heap;
}

/* -------------------------------------------- */
/* create new tree				*/
/* -------------------------------------------- */

CHeap *wCountNewHeap(Int32 nbnod, long initheap)
{
	Int32 	i;
	CHeap   *heap;
	Ctree 	*arbor;

	if (! (heap = NEW(CHeap)))
	    return NULL;

			/* compute actual # nodes */
			/* power of 2 >= nbnod	  */
				
	for (i = 1 ; (i - 1) < nbnod ; i *= 2)
		/* nop */ ;
	i--;
	
	heap->nbnodes = i;
		
	heap->depth = heap->nbnodes / 2;
	
	arbor = heap->tree = NEWN(Ctree, heap->nbnodes);

	if (! arbor) 
	   return wCountFreeHeap(heap);
	   
	for (i = 0 ; i < heap->nbnodes ; i++, arbor++) {
	    arbor->score  = arbor->indx = initheap;
	    arbor->pat[0] = '\000';
	    arbor->hit    = 0;
	}
	    
	return heap;
}

/* -------------------------------------------- */
/* update function of a minimiser		*/
/* (inversed version)				*/
/* -------------------------------------------- */

void wCountUpdateHeap(CHeap *heap,Int32 indx,Int32 score,PatString pat)
{
	Int32	i;
	Ctree 	*c_arbre, *l_arbre, *r_arbre, *e_arbre;

	for (c_arbre = heap->tree, i = 0 ; i < heap->depth ; ) {	
							
	    i = 2*i +1;				/* fils gauche 	*/
				
	    r_arbre = e_arbre = l_arbre = heap->tree + i;
	    r_arbre++;
				
	    if (r_arbre->score < l_arbre->score) {
		e_arbre++;
		i++;
	    }
				
	    if (score > e_arbre->score) {					
		*c_arbre = *e_arbre;		/* transfert cells	*/
		 c_arbre =  e_arbre;
	    }
	    else
		break;
	}					/* fin remise ˆ jour	*/

	c_arbre->score  = score;			
	c_arbre->indx   = indx;
	c_arbre->hit    = 1;
	strcpy(c_arbre->pat, pat);

}
