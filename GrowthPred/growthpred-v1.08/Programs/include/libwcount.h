/* ==================================================== */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: libcount.h					*/
/* history:						*/
/* April 97 <ed> first draft (taken from covgen.c)	*/
/* ==================================================== */

#include "cov_fasta.h"

/* ----------------------------------------------- */
/* constants					   */
/* ----------------------------------------------- */

#define ALPHA_LEN	 26	/* length of alphabet		*/
#define MAX_PAT_LEN	 31	/* max pattern length		*/
#define MAX_PAT_ERR	 31	/* max # of errors		*/

#define OPEN_TO_WRITE	 0	/* for CheckFile()		*/
#define OPEN_TO_READ	 1	/*				*/

#define WANT_ROTATOR	 1	/* wants completion indicator	*/

#define SEPARATOR '.'

/* ----------------------------------------------- */
/* Errors					   */
/* ----------------------------------------------- */

#define ERR_NONE	      0	
#define ERR_FILOPEN	    -80
#define ERR_FILRESTRALPH    -81
#define ERR_FILERRMATR	    -82
#define ERR_FILPAT	    -83
#define ERR_FILUSROLIGO	    -84
#define ERR_FILSEQ	    -85
#define ERR_INCOPT	    -86
#define ERR_PATSIZE	    -87
#define ERR_PATPOS	    -88
#define ERR_NMEMSEQ	    -89
#define ERR_NMEMPAT	    -90
#define ERR_NMEMHEAP	    -91
#define ERR_NMEMHARR	    -92
#define ERR_MODE	    -93
#define ERR_NBEST	    -94
#define ERR_SMAT	    -95
#define ERR_JUMP	    -96
#define ERR_STARTAT	    -97
#define ERR_OVERL	    -98
#define ERR_USROLIGSIZE	    -99
#define ERR_WINLEN	    -100
#define ERR_WINSTEP	    -101
#define ERR_NSEQ	    -102

#define STR_FILOPEN	    "Error in open file"
#define STR_FILRESTRALPH    "Error in open file with restricted alphabet"
#define STR_FILERRMATR	    "Error in open file with error matrix"
#define STR_FILPAT	    "Error in open file with patterns"
#define STR_FILUSROLIGO	    "Error in open file with user oligos"
#define STR_FILSEQ	    "Error in open file with sequences"
#define STR_INCOPT	    "Error: Incompatible options"
#define STR_PATSIZE	    "Error: Pattern size is too long"
#define STR_PATPOS	    "Error: Too many positions in pattern (max: 31)"
#define STR_NMEMSEQ	    "Error: Not enough memory for sequences"
#define STR_NMEMPAT	    "Error: Not enough memory for patterns"
#define STR_NMEMHEAP	    "Error: Not enough memory for heap tree"
#define STR_NMEMHARR	    "Error: Not enough memory for hash array"
#define STR_MODE	    "Error: Bad mode option (expects COVER or COUNT)"
#define STR_NBEST	    "Zero oligo's to be kept, analysis terminated"
#define STR_SMAT	    "Error in S Matrix"
#define STR_JUMP	    "Error: jump cannot be smaller than 1"
#define STR_STARTAT	    "Error: reading cannot start at pos <0"
#define STR_OVERL	    "Error: overlap must be non-negative"
#define STR_USROLIGSIZE	    "Error: Wrong Oligo size"
#define STR_WINLEN	    "Error: sliding window must be positive"
#define STR_WINSTEP	    "Error: sliding window step must be positive"
#define STR_NSEQ	    "Error: number of sequences must be positive"

/* ----------------------------------------------- */
/* data structures				   */
/* ----------------------------------------------- */

typedef unsigned long 	Ulong;		/* to make sure		*/
typedef 	 long 	Slong;		/* across plateforms	*/

typedef char PatString[MAX_PAT_LEN+1];	/* pattern string	*/
typedef char AlphaString[ALPHA_LEN+1];	/* alphabet string	*/

typedef struct {			/* pattern structure 	*/
    Bool	hit;			/* hit			*/
    int		patlen;			/* length		*/
    int    	maxerr;			/* max # errors allowed */
    int		overlap;		/* overlap size		*/
    PatString	pat;			/* char pattern		*/
    Ulong	smat[ALPHA_LEN];	/* S matrix		*/
} SPattern;

typedef struct {			/* generator structure 	*/
    int		maxlen;			/* length		*/
    AlphaString	alpha[MAX_PAT_LEN];	/* alpha each position	*/
    int		count[ALPHA_LEN];	/* max count foreach a	*/    
} SGenerator;

typedef struct {			/* substitution error 	*/
    char  allow[ALPHA_LEN][ALPHA_LEN];	/* matrix		*/
    Bool  active;			/*			*/
} SErrMat;

typedef struct {			/* cellule pour tri des */
    Bool	hit;			/* has an hit		*/
    Int32	indx;			/* meilleurs oligos	*/
    Int32	score;			/* minimier		*/
    PatString	pat;			/* pattern string	*/
} Ctree;

typedef struct {			/* structure heap	*/
    Int32	nbnodes;		/* actual size		*/
    Int32	nbbest;			/* useful nodes		*/
    Int32	depth;			/* actual depth		*/
    Ctree	*tree;			/* actual heap		*/
} CHeap;

/* ----------------------------------------------- */
/* macros					   */
/* ----------------------------------------------- */

#ifndef NEW
#define NEW(typ)		(typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim) 		(typ*)malloc((unsigned long)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)	(typ*)realloc((void *) (ptr), (unsigned long)(dim) * sizeof(typ))
#define FREE(ptr)		free((void *) ptr)
#endif

#ifdef MACINTOSH
#define ROTATOR(s, t)	  fprintf(stderr, "%ld/%ld\n", s, t)
#else
#define ROTATOR(s, t)	  fprintf(stderr, "\r%ld/%ld", s, t)
#endif
	
#if WANT_ROTATOR
#define PLAY_ROTATOR(s,t,g) if (! (++(s) % (g))) ROTATOR(s, t);
#else
#define PLAY_ROTATOR(s,t,g) s++
#endif

/* ----------------------------------------------- */
/* prototypes					   */
/* ----------------------------------------------- */

				     /* wcount_io.c */
int wCountCheckFile(char *filename, int toread);
Int32 wCountGetNbLines(char *filename,int achar);
char *wCountUpper(char *seq);
int wCountReadAlphabet(char *filename, char *alphabet, int *count);
int wCountReadErrorMatrix(char *filename, SErrMat *errmat);
int wCountReadGenerator(char *filename, SGenerator *gen);
int wCountReadOligos(char *filename, SPattern *pat, Int32 npat,
    int nberr, Int32 overl);
int wCountGetOligos(char *oligos, SPattern *pat, Int32 npat,
    int nberr, Int32 overl);
int wCountReadSequences(char *filename,FastaSequencePtr *seq, 
	int nseq);
int PrintResOligos(CHeap *heap, int totlen, int lab);
int getWinOligo(FastaSequencePtr *seqdata, Int32 nseq, 
    Int32 winlen, Int32 step);
int PrintResPerOligo(CHeap *heap, int totlen, char *id, int count, int lab);



				     /* wcount_hgener.c */
int wCountHashGenereOligo(CHeap *heap, FastaSequencePtr *seqdata, Int32 nseq,  
    Int32 overl, int patlen,int nerr, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag, Int32 startat, int worstflag);


				     /* wcount_bgener.c */
void wCountBruteGenereOligo(SGenerator *gen, CHeap *heap,
    FastaSequencePtr *seqdata, Int32 nseq, Int32 overl, int nerr, SErrMat *errmat,
    int coverflag, Int32 jump, int overflag, Int32 startat, int worstflag); 
int wCountDefaultGenerator(SGenerator *gen, char *alphabet, 
	int kmax, int *count);
Int32 wCountNumberGen(SGenerator *gen);
void wCountScanPatterns(SPattern *pat, FastaSequencePtr *seqdata, Int32 nseq, 
    CHeap *heap, Int32 npat, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag, Int32 startat, int worstflag, Bool reset);
Int32 wCountScanOneOligoAllSeq(PatString cpat, FastaSequencePtr *seqdata,
    Int32 nseq, Int32 overl, Int32 nerr, SErrMat *errmat, Bool reset,
    int coverflag, Int32 jump, int overflag, Int32 startat);
int CompileSMat(char *cpat,SPattern *ppat);
Int32 ManberSearch_S(char *seq, char *cover, SPattern *ppat,
    SErrMat *errmat, int coverflag, Int32 jump, int overflag);
CHeap *wCountFreeHeap(CHeap *heap);
CHeap *wCountNewHeap(Int32 nbnod, long initheap);
void wCountUpdateHeap(CHeap *heap,Int32 indx,Int32 score,PatString pat);
CHeap *wCountInitHeap(CHeap *heap, long initheap);

				     /* wcount_hgenerSngl.c */
int wCountHashGenereEachOligo(CHeap *heap, FastaSequencePtr *seqdata, Int32 nseq,  
    Int32 overl, int patlen,int nerr, SErrMat *errmat, int coverflag, 
    Int32 jump, int overflag, Int32 startat, int worstflag);
