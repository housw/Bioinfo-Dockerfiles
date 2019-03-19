/* ==================================================== */
/* Copyright (c) Atelier de BioInformatique		*/
/* @file: wCount_hgener.c				*/
/* oligo generator - Hash array version			*/
/* history:						*/
/* April 1997 <ed> first version			*/
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _H_Gtypes
#include "Gtypes.h"
#endif

#ifndef _H_libwcount
#include "libwcount.h"
#endif

static Int32 sOligoNum = 0;
static Int32 sOligoTot = 0;

typedef struct {		/* Hcod array element  */
    Int32 score;
    Int32 right;
} wCountHCell;

#define TWO_BITS

static UInt32 sCode[] = {

#include "dna_code.h"

};
				/* hard wired codes	*/

#define NB_NUC	4	

static char sNuc[NB_NUC+1] = "ACGT";	

static UInt32 sNuCode[NB_NUC] = { 0, 1, 2, 3 };

/* -------------------------------------------- */
/* hascode a fragment				*/
/* -------------------------------------------- */

static int sHCode(char *pat, int patlen, UInt32 *h)
{
	int    i, bad;
	UInt32 cod, hh;

	bad = 0;
	hh  = 0;
		
	for (i = 0 ; i < patlen ; i++) {
	
	    if ((cod = sCode[*pat++ - 'A']) & 0x10)
	        bad = i+1;

	    hh = (hh <<= 2) | cod;
	}
	
	*h = hh;
	
	return bad;
}	

/* -------------------------------------------- */
/* update a hascode				*/
/* right sliding				*/
/* -------------------------------------------- */

static int sUpdateCode(int symb, int patlen, UInt32 *h)
{
	int    bad;
	UInt32 cod, msk;

	bad = 0;
		
	if ((cod = sCode[symb - 'A']) & 0x10)
	    bad = patlen;

	msk = 0xffffffff << (patlen << 1);

	*h = ((*h << 2) & ~msk) | cod;
	
	return bad;
}	

/* -------------------------------------------- */
/* un-hascode a number				*/
/* -------------------------------------------- */

static int sUnHCode(UInt32 h,int patlen,char *pat)
{
	int    i;
	char *pp;
	
	pp = pat + patlen;
	
	*(pp--) = '\000';
	
	for (i = 0 ; i < patlen ; i++) {
	    *pp-- = sNuc[h & 0x3L];
	    h >>= 2;
	}
	
	return 0;
}	

/* -------------------------------------------- */
/* update a hascode				*/
/* single mutation				*/
/* -------------------------------------------- */

static UInt32 sMutateCode(int isymb,int pos,UInt32 hcod)
{
	UInt32 msk, scod;
	
	pos <<= 1;			/* bit pos		*/

	msk = 0x00000003 << pos;

	scod = sNuCode[isymb] << pos;	/* symb is assumed OK   */

	return (hcod & ~msk) | scod;
}	

/* -------------------------------------------- */
/* get symbol at corresponding pos from code    */
/* -------------------------------------------- */

static int sDeMutate(int pos, UInt32 hcod)
{
	UInt32 msk;
	
	pos <<= 1;			/* bit pos		*/

	msk = 0x00000003 << pos;	/* bits			*/

	return (hcod & msk) >> pos;
}	

/* -------------------------------------------- */
/* record oligo in Hash array			*/
/* -------------------------------------------- */

static void sRecordOligo(wCountHCell *harr, UInt32 cod, Int32 pos, Int32 overl,
    Int32 seqlen,int patlen, int coverflag, int overflag)
{
    Int32 left, right, diff, lastright;

    harr += cod;			/* points to entry 	*/
	
		
    left = pos - overl;			/* left pos	*/
    if (left < 0)
	left = 0;
   
    right = pos + patlen + overl - 1;	/* right pos	*/
    if (right >= seqlen)
	right = seqlen - 1;
					/*test first time around */
    if (harr->right >= 0)
	lastright= harr->right;
    else 
	lastright=-1;
	
    					/* test on overlapping  */
    if (! (overflag==1 || (pos - lastright > 0)))
	return;
      
			/*  update according to mode COVER or COUNT*/
   if(coverflag){
	    				/* overlap with previous hit */
        if (left <= harr->right)
	    diff = harr->right - left + 1;
	else
	    diff = 0;

	    harr->score += (right - left + 1 - diff);
	}
    else 
	harr->score++;
    
    harr->right = right;
}

/* -------------------------------------------- */
/* oligo mutator				*/
/* recursive proc.				*/
/* -------------------------------------------- */

static void sMutateOligo(wCountHCell *harr, UInt32 cod, Int32 pos, 
    Int32 overl, Int32 seqlen, int patlen, int start, int nerr, int terr,
    SErrMat *errmat, int coverflag, Int32 jump, int overflag)
{
	int    ipos, isymb, osymb;
	UInt32 h;
	
				    /* test reading frame */ 
	if (pos % jump != 0)
	    return;
	   
	if (nerr == terr) {	    /* End of recursion	*/
	    sRecordOligo(harr, cod, pos, overl, seqlen, patlen, coverflag, 
		overflag);

	    return;
	}

	for (ipos = start ; ipos < patlen ; ipos++) {

	    osymb = sDeMutate(ipos, cod);
	    
	    for (isymb = 0 ; isymb < NB_NUC ; isymb++) {
	    
	        if (errmat->active && (! errmat->allow[isymb][osymb]))
	            continue;
	            
	    	h = sMutateCode(isymb, ipos, cod);

	        sMutateOligo(harr, h, pos, overl, seqlen, patlen, ipos+1, 
		    nerr+1, terr, errmat, coverflag, jump, overflag);
	    }
	}
}
		

/* -------------------------------------------- */
/* sequence driven generator			*/
/* internal gener for one sequence		*/
/* -------------------------------------------- */
static void sHashGenereSingle(wCountHCell *harr, char *seq, long seqlen, 
    Int32 overl, int patlen, int nerr, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag)
{
	int	   off, noff;
	Int32	   pos;
	UInt32     h;
	char       *sd, *sf, *ff;
	
	if (seqlen < patlen)
	    return;
	    
	sd = seq;
	sf = sd + patlen;
	ff = sd + seqlen;
					/* startup 	*/
	
	while (sf <= ff) {					    
	    if ((off = sHCode(sd, patlen, &h)) == 0)
	        break;
	    sd += off + 1;
	    sf += off + 1;
	}

	pos = sd - seq;
					/* elongation	*/
	do {
	    if (off == 0) {

		sMutateOligo(harr, h, pos, overl, seqlen, patlen, 0, 0, 
		    nerr, errmat, coverflag, jump, overflag);
		  
	          sOligoTot ++;
	    }

	    noff = sUpdateCode(*sf++, patlen, &h);
	    
	    off = (noff ? noff : (off ? off - 1 : 0));
	    
	    pos++;
	    
	} while (sf <= ff);

} 

/* -------------------------------------------- */
/* size of Hash array				*/
/* -------------------------------------------- */
static Int32 sSizeArray(char *alpha,int patlen)
{
	Int32  asiz, siz;
	
	asiz = strlen(alpha);
	
	siz = 1;
	
	while (patlen--) 
	    siz *= asiz;
	    
	return siz;
}

/* -------------------------------------------- */
/* reset right field of Hash array		*/
/* -------------------------------------------- */
static void sResetArray(wCountHCell *harr,Int32 siz)
{
	while (siz--) 
	    (harr++)->right = -kBigInt32;
}

/* -------------------------------------------- */
/* sequence driven generator			*/
/* -------------------------------------------- */

int wCountHashGenereOligo(CHeap *heap, FastaSequencePtr *seqdata, Int32 nseq,  
    Int32 overl, int patlen,int nerr, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag, Int32 startat, int worstflag) 
{
    Int32 	    i, siz, score;
    wCountHCell	    *harr;
    PatString	    pat;	
	
				/* allocate hash array */

    siz = sSizeArray(sNuc, patlen);
	
    if (! (harr = NEWN(wCountHCell, siz)))
	return 0;
	 
    memset(harr, 0 , siz * sizeof(wCountHCell));
    
				/* record actual oligos in hash array	*/

    sOligoTot = 0;

    for (i=0; i<nseq; i++){
	sResetArray(harr, siz);	
	sHashGenereSingle(harr, &(seqdata[i]->seq[startat]), 
	    (seqdata[i]->length) - startat, overl, 
	    patlen, nerr, errmat, coverflag, jump, overflag);
    }

				/* put best oligos in heap	  */

    sOligoNum = 0;

    for (i = 0 ; i < siz ; i++) {
	score = (worstflag ? -harr[i].score : harr[i].score);
	if (score > heap->tree->score) {
	    sUnHCode(i, patlen, pat);
	    wCountUpdateHeap(heap, sOligoNum++, score, pat);
	}
    }
				/* end	*/
    FREE(harr);

    return 1;
}
