/* ==================================================== */
/* Copyright (c) Atelier de BioInformatique		*/
/* @file: wcount_main.c					*/
/* computes oligo <cover> or <counts> on sequences	*/
/* history:						*/
/* April 97 <ed> from covgen				*/
/* October 98 stand alone from Imagene version */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef SGI
#include <getopt.h>
#endif

#ifdef MACOSX
#include <unistd.h>
#endif
 
#include "Gtypes.h"

#include "libwcount.h"

#if PROTO
typedef int (*QsortCompFunc)(const void *, const void *);
#else
typedef int (*QsortCompFunc)();
#endif

#define WCOUNTSTANDALONE 1

#define nIncompOpt 7

#define QUARTET  "ACZCCZGCZGGZGTZCGZCTZTCZaczcczgczggzgtzcgzctztc"

#define LOCALPHA  "ACGT"

/* -------------------------------------------- */
/* local error notifier				*/
/* -------------------------------------------- */

int wCountError(char *msg, int errno)
{
    if (WCOUNTSTANDALONE)
	fprintf(stderr, "[libwCount] error : %s\n",  msg);
    
    return errno;
}

/* ------------------------------------------------------------ */
/* local variables						*/
/* ------------------------------------------------------------ */

static char sDftDnaAlpha[] = "ACGT"; 

static SErrMat sError={
    {
#include "error_mat.h"
    }, 
    Faux
};


/* ----------------------------------------------- */
/* gestion des arguments integer	   	   */
/* ----------------------------------------------- */

static int GetIArg(char *arg, int *value, int dft)
{
	int	varg;
	int	ok;

	ok = (sscanf(arg, "%d", &varg) == 1);

	*value = (ok ? varg : dft);
 
	return ok;
}

/* -------------------------------------------- */
/* printout header				*/
/* -------------------------------------------- */

#define PP fprintf(stdout, 

static void sPrintHeader(nseq, lseq, kmax, nerr, overlap, nbest, genflag, patflag,
			 oliflag, qckflag, alphabet, seqfile, patfile, olifile,
			 matfile, jump, startat, cumulflag)
	int	nseq, lseq, kmax, nerr, overlap, nbest;
	int	genflag, patflag, oliflag, qckflag, jump, startat, cumulflag;
	char	*alphabet, *seqfile, *patfile, *olifile, *matfile;
{
	PP	"/------------------------------------------			\n");
	PP	"/ WordCountQ Version 1 Mars 05	\n");
	PP	"/------------------------------------------			\n");
	PP	"/ sequence file          : %s\n", seqfile			   );
	PP	"/ nb sequences           : %ld\n", nseq			   );
	PP	"/ total length           : %ld\n", lseq			   );
	
	PP "/ internal generator     : k=%ld alpha=%s quick=%s\n", 
		      kmax, alphabet, 
		      qckflag ? "on" : "off"					   );
	PP	"/ # jump positions       : %ld\n", jump			   );
	PP	"/ start at position      : %ld\n", startat			   );
	PP	"/------------------------------------------			\n");
	
}

#undef PP
			
/* ----------------------------------------------- */
/* printout help				   */
/* ----------------------------------------------- */

#define PP fprintf(stdout, 

static void PrintHelp()
{
	PP	"------------------------------------------			\n");
	PP	" WordCountQ Version 0.1 Mars 05		\n");
	PP	"------------------------------------------			\n");
	PP	"synopsis :							\n");
	PP	"  compute nucleotides counts on sequences			\n");
	PP	"  at the 3rd codon position of quartets \n");
	PP	"								\n");
	PP	"usage: wcountq [options] seqfile				\n");
	PP	"------------------------------------------			\n");
	PP	"options:							\n");
	PP	"----- General options	 ------------------			\n");
	PP	"-h           : [H]elp - print <this> help			\n");
	PP	"------------------------------------------			\n");
	PP	"file formats							\n");
	PP	"								\n");
	PP	"seqfile:   fasta format sequences multple of 3 \n");
	PP	"--------------------------------				\n");
}

#undef PP

/* ----------------------------------------------- */
/* check options compatibility			   */
/* ----------------------------------------------- */

static int sCompatibleOption(int *vecopt)
{
	int   i, j, err=0;
			    /*a- alphabet		*/
			    /*A- restricted alphabet	*/
			    /*p- pattern file		*/
			    /*Q- quickmode 		*/
			    /*u- user oligos		*/
			    /*c- cover mode		*/
			    /*j- jump option		*/
			    /*l- cumulative counts	*/
			    
	static char sCompat[nIncompOpt][nIncompOpt] = { 	
		/*	a  A  p  u  c  j  l */
		/* a */ 1, 0, 0, 0, 1, 1, 0,
		/* A */ 0, 1, 0, 0, 1, 1, 0,
		/* p */ 0, 0, 1, 0, 1, 1, 0,
		/* u */ 0, 0, 0, 1, 1, 1, 1,
		/* c */ 1, 1, 1, 1, 1, 0, 1,
		/* j */ 1, 1, 1, 1, 0, 1, 1,
		/* l */ 0, 0, 0, 1, 1, 1, 1
	};

	for (i=0; i<nIncompOpt-1;i++){
	    if (*(vecopt+i)==1)
		for (j=i+1;j<nIncompOpt; j++)
		    if (sCompat[i][j]==0&& (*(vecopt+j)==1))
			err++;
	}
	return (err ? 0 : 1);	
}

/* ------------------------------------------------------------ */
/* total length							*/
/* ------------------------------------------------------------ */

static Int32 sTotalLength(FastaSequencePtr *seqdata, Int32 nseq)
{
    Int32 i, sum;
    
    for (i = sum = 0; i<nseq; i++)
		sum+=seqdata[i]->length;

    return sum;
}

/* -------------------------------------------- */
/* fonction de comparaison de cellules de 	*/
/* minimier pour quick_sort			*/
/* (pour tri en ordre decroissant par SCORE)		*/
/* -------------------------------------------- */

static int wCountCtreeComparS(Ctree *x, Ctree *y)
{
	if (y->score == x->score)
		return 0;
		
	return (y->score > x->score ? 1 : -1);
}
	
/* -------------------------------------------- */
/* fonction de comparaison de cellules de 	*/
/* minimier pour quick_sort			*/
/* (pour tri en ordre croissant par MOTS)		*/
/* -------------------------------------------- */

static int wCountCtreeComparM(Ctree *x, Ctree *y)
{
	if ((x->score < 0) && (y->score < 0))
		return 0;
	if (x->score < 0)
		return 1;
	if (y->score < 0)
		return -1;
	return strcmp(x->pat, y->pat);
}
	

/* ------------------------------------------------------------ */
/* main	call							*/
/* ------------------------------------------------------------ */

main( int argn, char *argv[])
{
	int		 i, j, errflag, patflag, coverflag, carg, value, oliflag,
			 qckflag, alpflag, nbest, jump, startat, overflag,
			 worstflag, genflag, overlap, nerr, kmax, nseq, ltot, allsortflag;
	int		 alcount[ALPHA_LEN], vecopt[nIncompOpt], cumulflag;
	long 	res[4], count;
	Int32		 npat;
	SPattern	 *pat;
	SGenerator	 generator;
	FastaSequencePtr *seqdata;
	CHeap		 *heap;
	AlphaString	 alphabet;
	char		 matfile[256], *c, lastopt[2], options[256], 
				seqfile[256],patfile[256], olifile[256], dinuc[3], *c1, *c2;

		    	    /* -------------------- */    
			    	/* INIT and defaults    */
			    	/* -------------------- */

	nbest       = 4;		/* default 50 best (worst)*/
	nerr   	    = 0;		/* default no errors allowed*/
	kmax        = 1;		/* default hexanucleotides */
	jump		= 3;		/* default jump 1 */
	startat		= 2;		/* default startat 0*/
	overflag	= 1;		/* default overlap counts */
	patflag	    = 0;		/* default no pattern file */
	alpflag	    = 0;		/* default std alphabet */
	cumulflag	= 0;		/* default not cumulative */
	worstflag	= 0;		/* default top */
	allsortflag = 1;		/* default not all*/
	
	errflag     = 0;		/* errors initialized at 0 */
	coverflag 	= 0;		/* just COUNT not COVER*/
	overlap 	= 0;		/* overlap for COVER off */
	genflag    = 1;		/* only generates oligos */
	oliflag 	= 0;		/* no User oligos */
	qckflag = 1; 		    /* setup auto-quick flag by default */
count = 0;
	strcpy(matfile, "<internal>");

    strcpy(alphabet, sDftDnaAlpha);
	
	for (i = 0 ; i < ALPHA_LEN ; i++)
	    alcount[i] = MAX_PAT_LEN;
	    
	for (i=0; i<nIncompOpt; i++)
	    vecopt[i]=0;
	
	lastopt[1] ='\000';
	*options   = '\000'; 

	/* ---------------------------- */    
	/* parse cmdline arguments	*/
	
	while ((carg = getopt(argn, argv, "hS")) != -1) {
	    
	    *lastopt = carg;
	    strcat(options, lastopt);
	    
	    switch(carg) {

		case 'h' :		/* Help		*/
		    PrintHelp();
		    exit(0);
		    break;

		case '?' :		/* misusage	*/
		    errflag++;
	    }
	}
	
	/* ---------------------------- */    
	/* should remain 1 arguments	*/

	if ((argn -= optind) != 1)
	    errflag++;

	if (errflag) 
	    return wCountError("INIT", 0);

	strcpy(seqfile, argv[optind]);
	
	    
			    /* ----------------------------------- */
			    /* update generator		*/
			    /* ----------------------------------- */	
			    
	wCountDefaultGenerator(&generator, alphabet, kmax, alcount);
			    
	npat = wCountNumberGen(&generator);

		/* ---------------------------- */
		/* open and read sequence files	*/

	if (! freopen(seqfile, "r", stdin))
	    (void) wCountError("cant open sequence file", 1);

	nseq = wCountGetNbLines(seqfile, '>');
	
	if (nseq<1)
	    return wCountError(STR_NSEQ, ERR_NSEQ);

	if (! (seqdata = NEWN(FastaSequencePtr, nseq)))
	    (void) wCountError("not enough memory (NEW Sequence)", 5);

	
	if (! wCountReadSequences(seqfile, seqdata, nseq))
	    (void) wCountError("reading sequence file", 15);

	ltot = sTotalLength(seqdata, nseq);
    
	sPrintHeader(nseq, ltot, kmax, nerr, overflag, nbest, 
		     genflag, patflag, oliflag, qckflag, 
		     alphabet, seqfile, patfile, olifile, matfile,
		     jump, startat, cumulflag);

		/* ---------------------------- */
		/* count				*/
	
	for (i=0;i<4; i++)
		res[i] = 0;
	
	for (i=0; i<nseq; i++){
		c1= seqdata[i]->seq;
		c2= &seqdata[i]->seq[startat];
		for (j=startat; j<seqdata[i]->length; j+=jump, c1+=jump, c2+=jump){
			strncpy(dinuc, c1, 2);
			dinuc[2]= '\0';
			if (strstr(QUARTET, dinuc)) {
/* 				printf("%s %c\n", dinuc, *c2); */
			    switch(*c2) {

				case 'A' :		/* A		*/
					res[0] ++;
				   	break;

				case 'C' :		/* A		*/
					res[1] ++;
				   	break;

				case 'G' :		/* A		*/
					res[2] ++;
				   	break;

				case 'T' :		/* A		*/
					res[3] ++;
				   	break;
			    }
			}	
		}
	}
	
	count = res[0]+res[1]+res[2] +res[3];
	
	for (i=0;i<4; i++){
		printf("%c\t%d\t%g\n", LOCALPHA[i],res[i], (float) res[i]/count );
	}

    exit(0);
}

