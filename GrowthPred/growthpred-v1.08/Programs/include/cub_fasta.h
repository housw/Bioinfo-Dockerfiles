/* ==================================================== */
/*	Copyright (c) GEM	*/
/* 	File: cub_fasta.h				*/
/*	History:					*/
/*	History:					*/
/*	1/09/09 : Alain Viari et Eduardo Rocha	*/
/* ==================================================== */

#ifndef _H_cub_fasta

#define _H_cub_fasta

/* ==================================================== */
/* Constantes						*/
/* ==================================================== */

//#define FASTA_NAMLEN  32  	/* max length of seq. name	 */
#define FASTA_NAMLEN  256      /* max length of seq. name       */
#define FASTA_COMLEN  550 	/* max length of seq. comment	 */

#define FASTA_CHAR_PER_LINE 50	/* # of chars per line in output */


#ifndef BUFSIZ
#define BUFSIZ 10000
#endif

#define ErrWarning	 0	    /* shell error codes		*/
#define ErrMemory	10
#define ErrSymbol	11
#define ErrReadSeq	12
#define ErrBadArg	13
#define ErrNotImp	14
#define ErrFilOpe	15
#define ErrFilAcc	16
#define ErrManySq	17
#define ErrSeqIo	18
#define ErrNoSeq	19

#define Seq_Rotator	100		/* be verbose every 100 seq.	*/

/* ==================================================== */
/* Macros						*/
/* ==================================================== */

#ifndef NEW
#define NEW(typ)		(typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim) 		(typ*)malloc((unsigned long)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)	(typ*)realloc((void *) (ptr), (unsigned long)(dim) * sizeof(typ))
#define FREE(ptr)		free((Ptr) ptr)
#endif

#define IS_UPPER(c) 	    (((c) >= 'A') && ((c) <= 'Z'))
#define IS_LOWER(c) 	    (((c) >= 'a') && ((c) <= 'z'))
#define IS_ALPHA(c) 	    (IS_UPPER(c) || IS_LOWER(c))
#define TO_UPPER(c) 	    ((c) - 'a' + 'A')
#define TO_LOWER(c) 	    ((c) - 'A' + 'a')

#define MIN(a, b) 	    ((a) < (b) ? (a) : (b))	
#define MAX(a, b) 	    ((a) > (b) ? (a) : (b))	
#define ABS(x) 		    ((x) >= 0  ? (x) : -(x))	

/* ==================================================== */
/* Structures de donnees				*/
/* ==================================================== */

typedef struct {			/* -- Sequence ---------------- */
	int	ok;			    /* error flag			*/
	int	length,			/* longueur			*/
			bufsize;		/* size of current seq buffer	*/
	char    name[FASTA_NAMLEN],	/* nom 				*/
		comment[FASTA_COMLEN],	/* commentaire			*/
		*seq,				/* sequence			*/
		*cover;				/* couverture <covgen added>	*/	
} FastaSequence, *FastaSequencePtr;

		    /* ------------------------ */
		    /* Program parameters	*/
		    /* ------------------------ */
typedef struct {
    int  nkeep;
    int  maxseq;
    int  seed;
    int  opt_verbose;
    char dbfile[FILENAME_MAX];
} Parameters;

typedef struct {
    int  ok;
    char id[BUFSIZ];
} idStructLG;



/* ==================================================== */
/*  Prototypes			*/
/* ==================================================== */

					/* cov_fasta.c 	*/

int 		 CountAlpha 	    ( char *buf );
char  		 *StrcpyAlpha 	    ( char *s1 , char *s2 );
char  		 *NextSpace 	    ( char *buffer );
char  		 *GetFastaName 	    ( char *buffer );
char 		 *GetFastaComment   ( char *buffer );

FastaSequencePtr FreeFastaSequence  ( FastaSequencePtr seq );
FastaSequencePtr NewFastaSequence   ( void );

int		 ReadFastaSequence  ( FILE *streamin, FastaSequencePtr seq );
int		 ReadFastaAlignment ( FILE *streamin, FastaSequencePtr seq );
void 		 WriteFastaSequence ( FILE *streamou, FastaSequencePtr seq , int char_per_line );

				/* fshuf_util.c    */
int	Erreur		    (char *msg , int stat);
int	MemoryErreur	    (char *where, int stat);
int	AssignToStdin	    (char *filename);
int	CheckAccess	    (char *filename, char *mode);
char	*StrTrimLeading	    (char *str);
char	*UpperStr	    (char *str);
char	*LowerStr	    (char *str);
char	*DeleteInStr	    (char *str, char *alphabet);
char	*ReplaceInStr	    (char *str, char *alphabet, int repchar);
float	UserCpuTime	    ();
float	SystemCpuTime	    ();
float	MemoryArena	    ();


#endif
