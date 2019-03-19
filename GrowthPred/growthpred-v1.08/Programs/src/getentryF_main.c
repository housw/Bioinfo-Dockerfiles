/* ==================================================== */
/* @File getentryF_main.c					*/
/* @history:						*/
/* @+	Aug. 02 <Ed> first draft			*/
/* @+	Sept. 02 <Ed> option -s added			*/
/* @+	Sept. 09 <Ed> option -f added			*/
/* ==================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>

#include "Gtypes.h"
#include "cub_fasta.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX    1024	/* max # of characters in a path name */
#endif

#ifndef BUFSIZ
#define BUFSIZ		1024
#endif

#define OPEN_TO_WRITE	 0	/* for CheckFile()		*/
#define OPEN_TO_READ	 1	/*				*/

#define MAX_SEQ		10	/* max # seq in db ONE DOES NOT NEED MORE THAN ONE	    */

#define DFT_KEEP	100	/* default # of seqs to keep NOT USED */

/* -------------------------------------------- */
/* get int argument				*/
/* -------------------------------------------- */

static int sGetIArg(char *opt, char *arg)
{
    int     iarg;
    char    buffer[BUFSIZ];

    if (   arg
	&& (sscanf(arg, "%d", &iarg) == 1)
	&& (iarg > 0))
	return iarg;

    (void) sprintf(buffer,  "[GetIArg] bad argument -%s %s",
			    opt, (arg ? arg : "<NULL>"));

    (void) Erreur(buffer, ErrBadArg);

    return -1; /* never reached */
}

/* -------------------------------------------- */
/* argument error				*/
/* -------------------------------------------- */

static int sBadArg(int carg, char *varg, int status)
{
    char    buffer[BUFSIZ];

    (void) sprintf(buffer, "option -%c Bad value: %s", carg, varg);

    (void) Erreur(buffer, status);

    return status;
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
/* check if file can be accessed		*/
/* unix access() does not exist on MacOSNOTX	*/
/* -------------------------------------------- */

static int word_len(char *s)
{
	int len;
	len = 0;
	
	while ((*s) && (*s != 32)){
		len ++;
		s++;
	}

	return len;
}


/* -------------------------------------------- */
/* read ids file				*/
/* -------------------------------------------- */

static int ReadData(char *filename, char *ids, int opt_case)
{
    int  i, err=0, counter=0, temp1, temp2;
    char *c, buffer[1000], oli1[1000], oliout[1000];
	char *oli2 = " ";    

    oli1[0] = '\000';
    if (! AssignToStdin(filename))
		return 0;
    
    while (fgets( buffer, sizeof(buffer), stdin )) {
		if (sscanf(buffer, "%s", &oli1)!=1) {
	    		(void) Erreur(filename, NULL);
	    		err++;
			}
			
		strcat(oli1, oli2);	
		
		if (opt_case ) /* has to be changed and applied to all ids*/
			UpperStr(oli1);

		strcat(ids, oli1);
		counter++;	    
    }
    return (err ? 0 : 1);
}

/* ----------------------------------------------- */
/* printout help				   */
/* ----------------------------------------------- */

#define PP fprintf(stdout, 

static void PrintGetEntryFHelp()
{
	PP	"------------------------------------------			\n");
	PP	" getentryF Version 0.1 Sept 09					\n");
	PP	"------------------------------------------			\n");
	PP	"synopsis :							\n");
	PP	"  gets entries in a database from a file with ids \n");
	PP	"								\n");
	PP	"usage: getentry [options] db_file					\n");
	PP	"------------------------------------------			\n");
	PP	"options:							\n");
	PP	"-h           : [H]elp - print <this> help			\n");
	PP	"-c string    : case insensitive	\n");
	PP	"-f string    : file with ids (one per line)	\n");
	PP	"------------------------------------------			\n");
	PP	"file formats							\n");
	PP	"db_file:        fasta format					\n");
	PP	" if n > MAX_SEQ it outputs an error \n");
	PP	"--------------------------------				\n");
	PP	"								\n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit			   */
/* ----------------------------------------------- */

#define PP fprintf(stderr, 

static int ExitUsageGE(int stat)
{
	PP	"usage: getentry [options] db_file		\n");
	PP	"type \"getentry -h\" for help		\n");

	if (stat)
	    exit(stat);
	    
	return stat;
}

#undef PP


/* ----------------------------------------------- */
/* main						   */
/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      carg, iarg, rarg, errflag, nseq, opt_number, opt_comment, opt_id;
	int 	      i, j, k, iseq, iran, nbseq, nbkept, opt_case, opt_file, inputok;
	FILE	      *dbfile, *idfile;
	Parameters    parms;
	FastaSequence **seq;
	char	      buffer[BUFSIZ], idstring[BUFSIZ], commstring[BUFSIZ], *s, *space, fileaddr[250], long_idstr[1000000];
	/*TODO memory allocation for long_idstr */
	
	/* ---------------------------- */    
	/* defaut options  		*/

	errflag		= 0;
	nseq        = 0;
	opt_number  = 0; /*no longer used*/
    opt_comment	= 0; /*no longer used*/
    opt_id 		= 1;
    opt_case	= 0;
	inputok 	= 0; /*i,nput file with names*/


	parms.nkeep       = DFT_KEEP;
	parms.maxseq      = MAX_SEQ;
	parms.seed        = 0;
	parms.opt_verbose = 0;
	parms.dbfile[0]   = '\000';

        /* ---------------------------- */    
        /* parse cmdline arguments      */
        
  	while ((carg = getopt(argn, argv, "chf:")) != -1) {
            
         switch(carg) {

		case 'c' :		/* case insensitive */
		    opt_case = 1;
		    break;

       	case 'h' :              /* help             */
		    PrintGetEntryFHelp();
		    exit(0);

		case 'f' :		/* file with ids */
		    strcpy(fileaddr, optarg);
		    if (! CheckFile(fileaddr, OPEN_TO_READ))
		        errflag++;
		    inputok = 1;
		    break;

		case '?' :		/* misusage	    */
		    errflag++;

	    }
	}

	/* ---------------------------- */    
	/* should remain 1 argument	*/

	rarg = argn - optind;

	if (rarg == 1) {
	    strcpy(parms.dbfile, argv[optind]);
	    if (! (dbfile = fopen(parms.dbfile, "r")))
			exit(ErrFilOpe);
	}
	else
	    errflag++;
	
	if (inputok == 0)
		errflag++;

	if (errflag) 
	    ExitUsageGE(1);

/*reading indexes*/

    if (!(ReadData(fileaddr, long_idstr, opt_case)))
		return Erreur("error reading index file", 2);

	/* ---------------------------- */    
	/* read DB			*/
	/* to count # of sequences	*/
	
	seq = NEWN(FastaSequence*, parms.maxseq);
	
	nbseq = 0;
		    
	while ( (seq[0] = NewFastaSequence())
	         && ReadFastaSequence(dbfile, seq[0])){

    	if (! seq[0]->ok) {
			sprintf(buffer, "IO error with sequence #%d (%s)\n", 
			        nbseq-1, seq[0]->name);
			Erreur(buffer, ErrSeqIo);
	    }
	    
		if (opt_case)
			UpperStr(seq[0]->name);

		if (s=strstr(long_idstr, seq[0]->name)) {
			if  (strlen(seq[0]->name)==word_len(s)) 
		   		WriteFastaSequence(stdout, seq[0], FASTA_CHAR_PER_LINE);
    	}
	    
	    nbseq++;
	    	    	
		FreeFastaSequence(seq[0]);
	}

	(void) fclose(dbfile);
		
	exit(0);
}
