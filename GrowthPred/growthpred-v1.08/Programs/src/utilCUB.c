/* ---------------------------------------------------------------- */
/* Copyright (c) Eduardo Rocha @ Atelier de BioInformatique	Universite PIerre et Marie Curie */
/* @file: utilCUB.c					    */
/* @desc: tabulate codon usage strength ENC ENCp S					    */
/*								    */
/* @history:							    */
/* @+       <Ed>    : January 09 : first shot from util_tab_cod	*/
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define GENETIC_CODE_INSTANCE
#include "Genetic.h"
#undef GENETIC_CODE_INSTANCE


#include "cub_fasta.h"

#ifndef Max
#define Max(i, j)  ((i) > (j) ? (i) : (j))
#endif

#define NB_CODONS 65
#define X_CODON  (NB_CODONS-1)
#define NB_AA	  21
#define X_AA	 (NB_AA-1)

#define COUNT_MODE  0	    /* comptage absolu	    	*/
#define FREQU_MODE  1	    /* freq. relative	    	*/
#define SYNON_MODE  2	    /* freq. rel. sur synonymes	*/

//#define MIN_STORAGE 100
#define MIN_STORAGE 10000

#define USE_HASH    1	    /* quicker version 		*/
			    /* ! see listing !		*/

typedef char NameString[16];

typedef struct s_Storage {
    NameString  name;
    int 	counts[NB_CODONS];
} Storage;

static char sNuc[]     = DNA_ALPHA;
static char sAnuc[]    = C_DNA_ALPHA;
static char sGenNuc[]  = GEN_NUC_ALPHA;
static char sGenPro[]  = GEN_PRO_ALPHA;
static int  sNucNum[5] = {0, 1, 2, 3, 3};

static char sDna[] = "ACTG";
static char sAA[]  = "ACDEFGHIKLMNPQRSTVWY*";
static char mysDna[] = "ACGT";

/* -------------------------------------------- */
/* get float argument				*/
/* -------------------------------------------- */

static float sGetFArg(char *opt, char *arg)
{
    float     farg;
    char    buffer[BUFSIZ];

    if (   arg && (sscanf(arg, "%f", &farg) == 1))
	return farg;

    (void) sprintf(buffer,  "[GetFArg] bad argument -%s %s",
			    opt, (arg ? arg : "<NULL>"));

    return -1.0; /* never reached */
}



/* ----------------------------------------------- */
static char *sUpper(char *s)
{
    char *c;

    for (c = s ; *c ; c++) {
	if (islower(*c))
	    *c = toupper(*c);
    }

    return s;
}

/* ----------------------------------------------- */
static char *sUpperDna(char *s)
{
    char *c;

    for (c = sUpper(s) ; *c ; c++) {
	if (*c == 'U')
	    *c = 'T';	
    }

    return s;
}

/* ----------------------------------------------- */
static char *sUpperProt(char *s)
{
    char *c;

    for (c = sUpper(s) ; *c ; c++) {
	if (*c == '#')
	    *c = '*';	
    }

    return s;
}


/* ---------------------------------------------------- */
/* @Static : int * sTranslateCodon			*/
/* Purpose : translate codon -> aa			*/
/* see  bio_codon_translate				*/
/* ---------------------------------------------------- */
static int sTranslateCodon(char *codon, int *code)
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
    return (int) sGenPro[code[hash]];
}

/* ---------------------------------------------------- */
/* @Function : 	int * bio_codon_translate		*/
/* Purpose : return amino-acid				*/
/* input: codon char* 3 bases (in GEN_NUC_ALPHA)	*/
/*        codid int   (see Genetic.h)			*/
/* output: aa in one letter code (in GEN_PRO_ALPHA)	*/
/* ---------------------------------------------------- */
static int bio_codon_translate(char *codon, int codid)
{
    return sTranslateCodon(codon, theGeneticCode[codid].code);
}


/* ----------------------------------------------- */
static void sMakeStdCodons(char *buffer)
{
    int i, j, k;
    
    for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
    for (k = 0 ; k < 4 ; k++) {
	*buffer++ = sDna[i];
	*buffer++ = sDna[j];
	*buffer++ = sDna[k];
	*buffer++ = '/';
    }    
}

/* ----------------------------------------------- */
static void sMakeAaCodons(char *buffer, char *aa, char *codons, int code)
{
    int	    i, k;
    char    *c;

    while (*aa) {
	for (i = 0, c = codons ; i < X_CODON ; i++, c += 4) {
	    k = bio_codon_translate(c, code);
	    if (strchr(aa, k)) {
		strncpy(buffer, c, 3);
		buffer[3] = '/';
		buffer += 4;
	    }
	}
	aa++;
    } 
     
    *buffer = '\000';
}

/* ----------------------------------------------- */
static int sCodonIndex(char *s, char *c)
{
    int	i, ssa;
    char *sa, *sb;
    
    for (i = 0 ; i < 3 ; i++)
	if (! c[i]) return -1;
	
    sa  = c + 3;
    ssa = *sa;    
    *sa = '\000';
    sb = strstr(s, c);
    *sa = ssa;

    return (sb ? (int) (sb - s) / 4 : -1);
}

#if USE_HASH

/* ----------------------------------------------- */
/* this is a quicker version of sCodonIndex	   */
/* ! only works if alpha == sDna		   */

static int sCodonHash(char *alpha, char *c)
{
    int	i, h;
    char *p;
    
    for (i = h = 0 ; i < 3 ; i++) {
	if (! (p = strchr(alpha, c[i])))
	   return -1;
	h  = (h << 2) | (int) (p - alpha);
    }

    return h;
}

#endif

/* ----------------------------------------------- */
static int sAaIndex(char *saa, char *c, int code)
{
    int  aa;
    char *paa;
    
    aa = bio_codon_translate(c, code);
    if ((paa = strchr(saa, aa)) != 0) 
	return (int) (paa - saa);
    return X_AA;
}

/* ----------------------------------------------- */
static int sSumCounts(char *codons, int *counts, int zero)
{
    int  i, sum;
    char *trip;

    for (i = sum = 0, trip = codons ; i < X_CODON ; i++, trip += 4)
	sum += counts[i];

    return (sum ? sum : zero);
}

/* ----------------------------------------------- */
static void sSumAA(char *codons, char *aa, int code, int *counts, int *naa, int zero)
{
    int  i, k;
    char *trip;

    for (i = 0 ; i < NB_AA ; i++)
	naa[i] = 0;

    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
	k = sAaIndex(aa, trip, code);
	naa[k] += counts[i];
    }
    
    for (i = 0 ; i < NB_AA ; i++)
	if (naa[i] == 0)
	   naa[i] = zero;
}   


/* ----------------------------------------------- */
static void sPrintCodonUsage(int *tot, char *codons, int code)
{
    int	    k, n;
    char    *aa, *trip;
    char    bufaa[2],
	    bufco[4*NB_CODONS+1];

    bufaa[1] = '\000';
    
    for (aa = sAA ; *aa ; aa++) {
	*bufaa = *aa;
	sMakeAaCodons(bufco, bufaa, codons, code);
	printf ("%1s\n", bufaa);
	for (trip = bufco, n = 0 ; *trip ; trip += 4) {
	    k = sCodonIndex(codons, trip);
	    n += tot[k];
	}
	if (n == 0) n = 1;
	for (trip = bufco ; *trip ; trip += 4) {
	    k = sCodonIndex(codons, trip);
	    printf("  %3.3s\t%.2f %d\n", trip, (float) tot[k]/ (float) n, tot[k]);
	}
    }
}


/* ----------------------------------------------- */
static Storage *sIncreaseStorage(Storage *store, int *size)
{
    int      nsiz;
    Storage *new;

    nsiz = Max(*size * 2, MIN_STORAGE);
    

    if (store)
    	new = (Storage *) realloc(store, nsiz * sizeof(Storage));
    else
        new = (Storage *) malloc(nsiz * sizeof(Storage));

    if (new)
       *size = nsiz;
       
    return new;
}

/* ----------------------------------------------- */
static void sCopyStorage(Storage *store, int *counts, char *name)
{
    int i;
    
    for (i = 0 ; i < NB_CODONS ; i++)
        store->counts[i] = counts[i];
	
    (void) strncpy(store->name, name, sizeof(NameString));
    store->name[sizeof(NameString)-1] = '\000';
}    
   

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, j, k, l, imin, imax, nbseq, opt, code, sum, npart, total;
	int	      sa_flag, so_flag, r_flag, t_flag, p_flag, g_flag, a_flag, test_flag;
	int	      counts[NB_CODONS], totcd[NB_CODONS], corcd[NB_CODONS], frd_cnt[10], frd_cnt_out[10],
		      naa[NB_AA], totaa[NB_AA], tot_aa, aa_degeneracy[NB_AA], no_amin_error, sara_N, sara_O;
	double		faa[NB_AA],  fra[NB_AA], frd[10], faa_p[NB_AA], xaa_p[NB_AA], fra_p[NB_AA], frd_p[10], freq, ENC, ENCp;
	float 		nucfreq[4], expec_counts[NB_CODONS], expec_counts_aa[NB_AA], sara_S;
	Storage	      *partial;	       
	FastaSequence *seq;
	char	      *trip;
	char	      codons[4*NB_CODONS+1],
		      aas_ignore[256], 
		      usr_ignore[4*NB_CODONS+1]; 

	extern char *optarg;	/* externs for getopts (3C)	*/

	code    = 0;		/* universal genetic code	*/
	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	t_flag  = 0;		/* no total			*/
	r_flag  = 0;	/* no RSCU		*/
	p_flag  = 0;		/* no pretty print		*/
	g_flag  = 0;		/* no global correction		*/
	a_flag  = 0;		/* no aa names			*/
	test_flag = 0;		/* test if there is a stop, by default no */
	
	*aas_ignore = '\000';
	*usr_ignore = '\000';
	
	npart   = 0;
	partial = NULL;
	
	sMakeStdCodons(codons);	

	/* ---------------------------- */
	/* parse arguments 		*/
	/* ---------------------------- */
			
        while ((opt = getopt(argn, argv, "A:ac:C:G:hi:I:rsSt")) != -1) {
	    
	    switch (opt) {

		case 'a':
		    a_flag = 1;
		    break;

		case 'c':
		    if ((sscanf(optarg, "%d", &code) != 1)
				|| (code < 0) || (code > 8)) {
		       (void) printf("bad code value: -c (0-8)\n");
		       exit(5);
		    }
		    break;

		case 'A':
			nucfreq[0] = sGetFArg("A", optarg);
		    break;

		case 'C':
			nucfreq[1] = sGetFArg("C", optarg);
		    break;

		case 'G':
			nucfreq[2] = sGetFArg("G", optarg);
		    break;

		case 'h':
		    (void) printf("tabulate codon usage strength ENC ENC' P \n");
		    (void) printf("usage: utilCUB [-A freqA] [-C freqC] [-G freqG] [-c code] [-s] [-S] [-t]\n");
		    (void) printf("   -A   : frequency of A\n");
		    (void) printf("   -C   : frequency of C\n");
		    (void) printf("   -G   : frequency of G\n");
		    (void) printf("   -c code\n");
		    (void) printf("      0 : universal\n");
		    (void) printf("      1 : mito yeast\n");
		    (void) printf("      2 : mito vertebrate\n");
		    (void) printf("      3 : filamentous fungi\n");
		    (void) printf("      4 : mito insects & platyhelminthes\n");
		    (void) printf("      5 : Candida cylindracea\n");
		    (void) printf("      6 : Ciliata\n");
		    (void) printf("      7 : Euplotes\n");
		    (void) printf("      8 : mito echinoderms\n");
		    (void) printf("   -s       : ignore first (start) codon\n");
		    (void) printf("   -S       : ignore last  (stop)  codon\n");
		    (void) printf("   -t       : test if there is a stop and if so ignore \n");
		    exit(0);
		    break;

		case 'i':
		    (void) strcpy(usr_ignore, optarg);
		    (void) sUpperDna(usr_ignore);
		    break;

		case 'I':
		    (void) strcpy(aas_ignore, optarg);
		    (void) strcat(aas_ignore, "CW"); /*ignoring always Trp and Met*/
		    (void) sUpperProt(aas_ignore);
		    break;

		case 'r':
		    r_flag = 1;
		    break;

		case 's':
		    sa_flag = 1;
		    break;

		case 'S':
		    so_flag = 1;
		    break;

		case 't':
		    test_flag = 1;
		    break;
		
	     }
	}

	/* ---------------------------- */
	/* check usage			*/
	/* ---------------------------- */

	nucfreq[3] = 1 - nucfreq[0] -nucfreq[1] - nucfreq[2];
	if (nucfreq[3]<= 0 || nucfreq[0]<=0  || nucfreq[1]<= 0  || nucfreq[2]<=0){
		(void) printf("negative or absent frequencies: -A -C -G \n");
		exit(5);
	}

	if (a_flag && (! r_flag)){
		fprintf(stderr,"error: -a depends on -r\n");
		exit(5);
		}


	if (*aas_ignore && *usr_ignore) {
		fprintf(stderr,"tab_codon: -i and -I incompatible options\n");
		exit(5);
	}

	if (*aas_ignore)
	    sMakeAaCodons(usr_ignore, aas_ignore, codons, code);

	seq = NewFastaSequence();

	nbseq = 0;

	for (i = 0 ; i < NB_CODONS ; i++){
	    totcd[i] = 0;
	    expec_counts[i] = 1.0;
	    }

	for (i = 0 ; i < NB_AA ; i++){
	    aa_degeneracy[i] = 0;
	    expec_counts_aa[i] = 0.0;
	    }

    printf("name");
	    
	for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
		if (sCodonIndex(usr_ignore, trip) == -1) {
			if (r_flag)
			   	printf("\t%3.3s", trip);
		   	for (j=0; j<3; j++){
		   		for (k=0, l=0; k<4; k++)
		   			if (mysDna[k] == trip[j])
		   				l = k; 
		   		expec_counts[i] *= nucfreq[l];
		   	}
	    	k = sAaIndex(sAA, trip, code);
	    	aa_degeneracy[k] ++;
		   	expec_counts_aa[k] += expec_counts[i] ;
		 	if (a_flag) 
		 		printf("/%1c", sAA[k]);
		}
	}

	printf("\tENC\tENCp\tSi\n");

	for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
			if (sCodonIndex(usr_ignore, trip) == -1) {
	    	k = sAaIndex(sAA, trip, code);
			expec_counts[i] = expec_counts[i] / expec_counts_aa[k]; 
			if (expec_counts[i] <= 0){
				fprintf(stderr,"ERROR: expected value for codon %s3.3 of aa %c\n", trip, sAA[k]);
				exit(5);
				}
			}
		}
		

	/* ---------------------------- */
	/* loop on sequences		*/
	/* ---------------------------- */

while (ReadFastaSequence(stdin, seq)) {
            nbseq++;

	    if ((! seq->ok) || (seq->length%3 != 0)){
	    	(void) fprintf(stderr,"error at seq # %d\n", nbseq);
			continue;
			}

		seq->length -= (so_flag ? 3 : 0);
		
	    /* -------------------------------- */
	    /* compute counts			*/
	    /* -------------------------------- */

	    for (i = 0 ; i < NB_CODONS ; i++)
			counts[i] = 0;
	    
	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length;

  	    k = -1;
	    
	    for (i = imin ; i < imax ; i += 3) {

#if USE_HASH
			k = sCodonHash(sDna, seq->seq + i);
#else
			k = sCodonIndex(codons, seq->seq + i);
#endif

			if (k >= 0)
		    	counts[k]++;
			else
		    	fprintf(stderr, "invalid codon %3.3s at position %d in sequence %s\n",
				    seq->seq + i, i+1, seq->name);
	    }

	    /* remove ignored codons		*/

	    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
			if (sCodonIndex(usr_ignore, trip) >= 0)
		    	counts[i] = 0;
	    }

	    /* -------------------------------- */
	    /* compute totals			*/
	    /* -------------------------------- */
	    
	    for (i = 0 ; i < NB_CODONS ; i++)
			totcd[i] += counts[i];				

	    /* -------------------------------- */
	    /* store or print local values	*/
	    /* -------------------------------- */
		sum = sSumCounts(codons, counts, 1);

		for (i = 0 ; i < NB_AA ; i++)
		    faa[i] = xaa_p[i] = faa_p[i] = 0;

		sSumAA(codons, sAA, code, counts, naa, 0);

		printf("%s", seq->name);

		for (i = 0, trip = codons, sara_N = sara_O = 0 ; i < X_CODON ; i++, trip += 4) {
			if (sCodonIndex(usr_ignore, trip) == -1) {
					k = sAaIndex(sAA, trip, code);
					freq = -1;
				    if (naa[k] > 0) {
				    	freq = (float) counts[i] / (float) naa[k];
				    		/*computing S*/
				    	if (strchr("FYIN" , sAA[k])){
							if (trip[2] == 'T')
								sara_N += counts[i] ;
							else if (trip[2] == 'C')
								sara_O += counts[i] ;
							}
				    		/*computing for ENCp*/
				    	xaa_p[k] += (float) naa[k] * (freq - expec_counts[i]) * (freq - expec_counts[i])/expec_counts[i];
				    	faa[k] += freq*freq ;
			    		}
			    	else
						faa[k] = -1;
					if (r_flag)
						printf("\t%.1f", 100. * freq);
				
			}
		}

/* check for stops: if stops print NA NA NA and goes to next*/
		if (test_flag && naa[20] > 0){
			fprintf(stderr, "error %d STOPS in: %s\n", naa[20], seq->name);
			printf("\tNA\tNA\tNA\n"); 
			continue;
			}
			
		no_amin_error =  0; 
	
		for (i=0, ENC=0; i < NB_AA - 1; i++){
			fra[i] = 0;
			if (naa[i] > aa_degeneracy[i] ){
				fra[i] = (naa[i]*faa[i]-1)/ (naa[i] -1);
				faa_p[i] = (xaa_p[i] + naa[i] - aa_degeneracy[i] )/(aa_degeneracy[i] * (naa[i]-1));
				}
			else if ( aa_degeneracy[i] > 1){
				fra[i] = -1; /*aa does not exist*/
				faa_p[i] = -1;
				}
			else if (aa_degeneracy[i] == 1){
				fra[i] = 1;
				faa_p[i] = 1;
				}
			
			}
		
		for (i=0; i < 10; i++){
			frd[i] = frd_p[i] = 0;
			frd_cnt_out[i] = frd_cnt[i] = 0;
			}
/* frd is the F3, F4, etc  */
/* frd_cnt is to know for each how many aa contributed */

		for (i=0; i < NB_AA - 1; i++){
			if (fra[i] == -1) 
				frd_cnt_out[aa_degeneracy[i]] ++; 
			else {
				frd[aa_degeneracy[i]] += fra[i] ;
				frd_p[aa_degeneracy[i]] += faa_p[i];
				frd_cnt[aa_degeneracy[i]] ++; 
				}
			}

		for (i=0, ENC=0, ENCp=0; i < 10; i++){
			if (frd[i] > 0 && frd_p[i] > 0 && frd_cnt > 0 ){
				ENC += (frd_cnt[i]+ frd_cnt_out[i])/(frd[i]/frd_cnt[i]) ;
				ENCp += (frd_cnt[i]+ frd_cnt_out[i])/(frd_p[i]/frd_cnt[i]) ;

				}
			else if  (frd_cnt_out[i] > 0  && frd_cnt[i] == 0 )
				no_amin_error ++ ;			
			}
		
		sara_S = (float) sara_O / (float) (sara_O + sara_N);
		
		if (! no_amin_error)
			printf("\t%g\t%g\t%g\n", ENC, ENCp, sara_S);
		else 
			printf("\tNA\tNA\tNA\n");

	}
		 		

#if 1

	if (partial)
	   free(partial);

	FreeFastaSequence(seq);

#endif

	exit(0);
	
}
