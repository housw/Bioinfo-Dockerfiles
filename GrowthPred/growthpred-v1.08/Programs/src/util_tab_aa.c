/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_tab_aa.c						    */
/* @desc: tabulate aminoacid usage				    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* @+       <ed> : 2001 : standalone			    */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SGI
#include <getopt.h>
#endif

#include "Genetic.h"
/*#include "libbio.h"*/
#include "cub_fasta.h"

#ifndef Max
#define Max(i, j)  ((i) > (j) ? (i) : (j))
#endif

#define NB_AA 21
#define X_AA  (NB_AA-1)

static char sAA[] = "ACDEFGHIKLMNPQRSTVWY";

/* ----------------------------------------------- */
static int sCharIndex(char *s, int c)
{
    char *ss = strchr(s, c);
    return (ss ? (int) (ss - s) : -1);
}

/* ----------------------------------------------- */
static int sSum(int *t, char *s, int imax, char *ignore)
{
    int i, sum;
    
    for (i = sum = 0 ; i < imax ; i++)
	if ((ignore == 0) || (sCharIndex(ignore, s[i]) == -1))
	    sum += t[i];
	
    return Max(sum, 1);
}

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, k, imin, imax, nbseq, nbaa, opt;
	int	      f_flag, t_flag, r_flag, x_flag;
	int	      count[NB_AA], tot[NB_AA];
	FastaSequence *seq;
	char	      ignore[256];

	extern char *optarg;	/* externs for getopts (3C)	*/

	f_flag = 0;		/* consider first residue	*/
	t_flag = 0;		/* no total			*/
	r_flag = 0;		/* use counts not rel. freq.	*/
	x_flag = 0;		/* no unknown symbols		*/
	*ignore = '\000';
		
        while ((opt = getopt(argn, argv, "1hi:rtx")) != -1) {
	    
	    switch (opt) {

		case '1':
		    f_flag = 1;
		    break;

		case 'h':
		    (void) printf("tabulate amino-acid usage\n");
		    (void) printf("usage: util_tab_aa [-1] [-i alpha] [-x]\n");
		    (void) printf("   -1       : ignore first residue\n");
		    (void) printf("   -i alpha : ignore symbols in alpha\n");
		    (void) printf("   -r       : compute relative frequencies\n");
		    (void) printf("   -t       : print last total line\n");
		    (void) printf("   -x       : ignore unknown symbol\n");
		    exit(0);
		    break;

		case 'i':
		    (void) strcpy(ignore, optarg);
		    break;

		case 'r':
		    r_flag = 1;
		    break;

		case 't':
		    t_flag = 1;
		    break;

		case 'x':
		    x_flag = 1;
		    break;
		   
		case '?':
		    (void) printf("usage: tab_aa [-h] [-1] [-i alpha] [-r] [-t] [-x]\n");
		    exit(6);
		    break;
	     }
	}

	seq = NewFastaSequence();

	nbseq = 0;
	imin  = (f_flag ? 1 : 0);
	imax  = (x_flag ? NB_AA : NB_AA - 1);

	for (i = 0 ; i < NB_AA ; i++)
	    tot[i] = 0;

	printf("name");
	for (i = 0 ; i < imax ; i++) {
	    if (sCharIndex(ignore, sAA [i]) == -1)
			printf("\t%c", sAA[i]);
	}
	printf("\n");

	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok)
		(void) printf("error at seq # %d\n", nbseq);

	    for (i = 0 ; i < NB_AA ; i++)
		count[i] = 0;
	
	    for (i = imin ; i < seq->length ; i++) {
	    
		k = sCharIndex(sAA, seq->seq[i]);

		if (k >= 0)
		    count[k]++;
		else if (x_flag)
		    count[X_AA]++;
		else
		    fprintf(stderr, "Invalid symbol %c in sequence %s\n",
				     seq->seq[i], seq->name);
	    }

	    nbaa = (r_flag ? sSum(count, sAA, imax, ignore) : 0);

	    printf("%s", seq->name);

	    for (i = 0 ; i < imax ; i++) {
		if (sCharIndex(ignore, sAA[i]) == -1) {
		    if (r_flag)
			printf("\t%.1f", 100. * (float) count[i] / (float) nbaa);
		    else
			printf("\t%d", count[i]);
		}
	    }
	    printf("\n");

	    if (t_flag) {
		for (i = 0 ; i < NB_AA ; i++)
		    tot[i] += count[i];
	    }

	}

	if (t_flag) {
	    printf("total:");
	    nbaa = (r_flag ? sSum(tot, sAA, imax, ignore) : 0);
	    for (i = 0 ; i < imax ; i++) {
		if (sCharIndex(ignore, sAA [i]) == -1) {
		    if (r_flag)
			printf("\t%.1f", 100. * (float) tot[i] / (float) nbaa);
		    else
			printf("\t%d", tot[i]);
		}
	    }
	    printf("\n");
	}	    

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


