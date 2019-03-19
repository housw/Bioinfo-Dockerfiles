/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_translate.c					    */
/* @desc: util_translate a fasta nucleic sequence		    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* @+       <Gloup> : Feb 99 : -h -s -S options added	  	    */
/* @+       <Gloup> : Sep 99 : -m option added	(temp. patch)  	    */
/* @+       <Ed>    : June 01 : independent version   */
/* @+       <Ed>    : Oct 02 : -f added   */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "Gtypes.h"
#include "cub_fasta.h"

#define GENETIC_CODE_INSTANCE
#include "Genetic.h"
#undef GENETIC_CODE_INSTANCE

#define OPEN_TO_WRITE	 0	/* for CheckFile()		*/
#define OPEN_TO_READ	 1	/*				*/

static char sStarts[] = "ATG/GTG/TTG/AUG/GUG/UUG";

static char sNuc[]     = DNA_ALPHA;
static char sAnuc[]    = C_DNA_ALPHA;
static char sGenNuc[]  = GEN_NUC_ALPHA;
static char sGenPro[]  = GEN_PRO_ALPHA;
static int  sNucNum[5] = {0, 1, 2, 3, 3};
static char AA_DEPEP[] = "TVSAG";

/* -------------------------------------------- */
/* error reporting				*/
/* -------------------------------------------- */

static int Error(char *msg, int stat)
{
	fprintf(stderr, "*Error* [%d] %s\n", stat, msg);	
	
	if (stat)
	   exit(stat);
	   
	fflush(stderr);	
	return stat;
}


/* -------------------------------------------- */
/* uppercase sequence				*/
/* -------------------------------------------- */
static char *Upper(char *seq)
{
	char *s;
	
	for (s = seq; *s ; s++)
	   if (islower(*s))
	   	*s = toupper(*s);
	   	
	return seq;
}

/* -------------------------------------------- */
/* check if file can be accessed		*/
/* unix access() does not exist on MacOSNOTX	*/
/* -------------------------------------------- */

static int CheckFile(char *filename, int toread)
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

static int GetNbLines(char *filename, int achar)
{
	int count;
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
/* read sequence file				*/
/* -------------------------------------------- */

static int ReadSequences(char *filename, FastaSequencePtr *seq, 
	int nseq)
{
	int i;

	if (! AssignToStdin(filename))
	    return 0;
	
	for (i = 0 ; i < nseq ; i++) {
	
	   if (! (seq[i] = NewFastaSequence()))
	   		return Error("not enough memory (NewFastaSequence)", 0);
	   
	   if (! ReadFastaSequence(stdin, seq[i]))
	   		return 0;

	    (void) Upper(seq[i]->seq);

	}

	return 1;
}

/* ----------------------------------------------- */
/* copy s2 to s1 				   */

static void sStrCpy(char *s1, char *s2)
{
    while (*s2)
       *s1++ = *s2++;
    *s1 = '\000';
}

/* ----------------------------------------------- */
/* copy n chars from s2 to s1			   */

static void sStrNCpy(char *s1, char *s2, int n)
{
    while (*s2 && (n-- > 0))
       *s1++ = *s2++;
    *s1 = '\000';
}


/* ---------------------------------------------------- */
/* @Static : int * sTranslateCodon			*/
/* Purpose : translate codon -> aa			*/
/* see  pwg_codon_translate				*/
/* ---------------------------------------------------- */
static int sTranslateCodon(char *codon, int *code, int *frmshift_flag)
{
    int  i, base, hash;
    char *p;
    
    for (i = hash = 0 ; i < 3 ; i++) {
	if ((p = strchr(sGenNuc, *codon++)) != NULL) {
	    base = ((int) (p - sGenNuc)) / 2;
	    hash = (hash * 4) + sNucNum[base];
	}
	else {
	    hash = 64;	/* bad letter in codon	*/
	    break;	/* or incomplete codon	*/
	}
    }
    
    if ((frmshift_flag > 0) && (code[hash]== 20))
        (*frmshift_flag) ++ ;

    return (int) sGenPro[code[hash]];
}

/* ---------------------------------------------------- */
/* @Function : 	int * seq_translate			*/
/* Purpose : translate sequence to protein  		*/
/* ---------------------------------------------------- */
static char* seq_translate(char *seq, int codid, int *frmshift_flag)
{
    int  *code;
    char *ps, *ns;
    
    if ((codid < 0) || (codid >= GEN_MAX_CODES))
		return NULL;

    code = theGeneticCode[codid].code;

    for (ns = ps = seq ; ns[0] && ns[1] && ns[2] ; ns += 3)
		*ps++ = sTranslateCodon(ns, code, frmshift_flag);
    
    *ps = '\000';

    return seq;    
} 

/* ----------------------------------------------- */
/* printout help				   */
/* ----------------------------------------------- */

#define PP fprintf(stdout, 

static void PrintHelp()
{
	PP	"#------------------------------------------			\n");
	PP	"# translate Version 1.3 Sep 06			\n");
	PP	"#------------------------------------------			\n");
	PP	"#synopsis :							\n");
	PP	"#  gets gene sequences and translates to protein \n");
	PP	"#								\n");
	PP	"# usage: translate [-c (0-8)] [-s] [-S] seq.file				\n");
	PP	"#------------------------------------------			\n");
	PP	"#options:							\n");
	PP	"# -h          : [H]elp - print <this> help			\n");
	PP	"# -c code     :  \n");
	PP	"#             : 0 - universal \n");
	PP	"#             : 1 - mito yeast \n");
	PP	"#             : 2 - mito vertebrate \n");
	PP	"#             : 3 - filamentous fungi & Mycoplasma \n");
	PP	"#             : 4 - mito insects & platyhelminthes \n");
	PP	"#             : 5 - Candida cylindracea \n");
	PP	"#             : 6 - Ciliata \n");
	PP	"#             : 7 - Euplotes \n");
	PP	"#             : 8 - mito echinoderms \n");
	PP	"# -p          : remove Met as Map does \n");	
	PP	"# -f          : reject if there are inner stops (frameshifts) \n");	
	PP	"# -m          : force Met as start codon if [AGT]TG \n");
	PP	"#             : note: incompatible with -s & -M \n");
	PP	"# -M          : really force start to be Met \n");
	PP	"#             : note: incompatible with -s & -m\n");
	PP	"# -s          : ignore first (start) codon \n");
	PP	"#             : note: incompatible with -m & -M\n");
	PP	"# -S          : ignore last  (stop)  codon \n");	
	PP	"#------------------------------------------			\n");
	PP	"#file formats							\n");
	PP	"#								\n");
	PP	"# seqfile       :  fasta format				\n");
	PP	"#--------------------------------				\n");
}

#undef PP

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      minlen, opt, code, nseqlin, i, frmshift_flag, frmshift_cnt;
	int	      sa_flag, so_flag, m_flag, fm_flag, p_flag;
	FastaSequencePtr *seq;
	char	      start[4], fileseq[100];

	extern char *optarg;	/* externs for getopts (3C)	*/

	code = 0;
	sa_flag = so_flag = m_flag = fm_flag = p_flag  = frmshift_flag = frmshift_cnt = 0;
	
        while ((opt = getopt(argn, argv, "c:fhmMpsS")) != -1) {
	    
	    switch (opt) {

		case 'c':
		    if (   (sscanf(optarg, "%d", &code) != 1)
				|| (code < 0) || (code > 8)) 
		       return Error("bad code value: -c (0-8)", 0);
		    break;

		case 'f':
		    frmshift_flag = 1;
		    break;

		case 'h':
   		    PrintHelp();
		    exit(0);
		    break;

		case 'm':
		    m_flag  = 1;
		    break;

		case 'p':
		    p_flag  = 1;
		    break;

		case 'M':
		    fm_flag  = 1;
		    break;

		case 's':
		    sa_flag = 1;
		    break;

		case 'S':
		    so_flag = 1;
		    break;

		case '?':
		    (void) printf("usage: pwg_translate [-c (0-8)] [-s] [-S] < fastafile\n");
		    exit(6);
		    break;
	     }
	}

	if (sa_flag + m_flag + fm_flag > 1)
		return Error("m, s and p options are incompatible", 1);

	/* ---------------------------- */    
	/* should remain 1 arguments	*/

	if ((argn -= optind) != 1)
	    return Error("input file missing", 2);

	if (! CheckFile(argv[optind], OPEN_TO_READ))
	    return Error("input file unavailable", 2);
	else
            strcpy(fileseq, argv[optind]);
		
        nseqlin = GetNbLines(fileseq, '>');

	if (! (seq = NEWN(FastaSequencePtr, nseqlin)))
	    return Error("not enough memory (NEW Sequence)", 4);

	if (! ReadSequences(fileseq, seq, nseqlin))
	    return Error("reading sequence file", 5);

	for (i=0; i<nseqlin; i++) {

	    if (! seq[i]->ok) {
			(void) fprintf(stderr, "error at seq # %d\n", i+1);
			(void) fprintf(stderr,"bad length at seq # %d\n", i+1);
			continue;
	    }

	    minlen = 0;
	    if (sa_flag) minlen += 3;
	    if (so_flag) minlen += 3;

	    if (seq[i]->length <= minlen || (seq[i]->length % 3 != 0)) {
			(void) fprintf(stderr, "bad length at seq # %d\n", i+1);
			continue;
	    }

	    if (m_flag || fm_flag) {
	        sStrNCpy(start, seq[i]->seq, 3); /* keep start codon */
			(void) Upper(start);
	    }
	    
	    (void) seq_translate(seq[i]->seq, code, &frmshift_cnt); 

	    if (((m_flag) && (strstr(sStarts, start))) || fm_flag)
		 	seq[i]->seq[0] = 'M';

	    (void) strcat(seq[i]->comment, " (translation)");

	    seq[i]->length /= 3;

	    if (sa_flag) {
			sStrCpy(seq[i]->seq, seq[i]->seq + 1);
			seq[i]->length--;
	    }
	    
	    if (p_flag){
	     if (strchr(AA_DEPEP, seq[i]->seq[1])){
/*			printf("here is %c\n", seq[i]->seq[1]);*/
			sStrCpy(seq[i]->seq, seq[i]->seq + 1); 
			seq[i]->length--;				   

			}
	    }
	    
           /* printf("%c %d\n", seq[i]->seq[seq[i]->length-1], seq[i]->seq[seq[i]->length-1]); */
            if (seq[i]->seq[seq[i]->length-1] == 42) /* if the last is a stop */
                    frmshift_cnt -- ;
	    
        if (so_flag)
			seq[i]->length--;
            
        if (frmshift_flag > 0 && frmshift_cnt > 1) {
            (void) fprintf(stderr, "inner stops at seq # %d\n", i+1);
            frmshift_cnt=1;
			continue;
            }

	    WriteFastaSequence(stdout, seq[i], FASTA_CHAR_PER_LINE);

	}
	
	for (i=0; i<nseqlin; i++)
		FreeFastaSequence(seq[i]);

	exit(0);
}
