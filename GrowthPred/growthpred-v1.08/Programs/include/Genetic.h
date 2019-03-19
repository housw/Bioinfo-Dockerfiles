/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: Genetic.h						    */
/* @desc: Genetic codes / include file				    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : first draft for PWG from acnuc	    */
/* ---------------------------------------------------------------- */

#ifndef _H_Genetic

#define _H_Genetic

/* ==================================================== */
/* Constants						*/
/* ==================================================== */

#define GEN_NUC_ALPHA	"AaCcGgTtUu"
#define GEN_PRO_ALPHA	"RLSTPAGVKNQHEDYCFIMW*X"
#define DNA_ALPHA       "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
#define C_DNA_ALPHA     "TVGHEFCDIJMLKNOPQYSAABWXRZtvghefcdijmlknopqysaabwxrz"

#define GEN_MAX_CODES	9

#define	GEN_CODE_UNIVL  0   /* Universal			*/
#define	GEN_CODE_MYEAS	1   /* mito yeast			*/
#define	GEN_CODE_MVERT  2   /* mito vertebrate			*/
#define	GEN_CODE_FUNGI  3   /* filamentous fungi		*/
#define	GEN_CODE_MINSE  4   /* mito insects & platyhelminthes	*/
#define	GEN_CODE_CANDI  5   /* Candida cylindracea		*/
#define	GEN_CODE_CILIA  6   /* Ciliata				*/
#define	GEN_CODE_EUPLO  7   /* Euplotes				*/
#define	GEN_CODE_MECHI  8   /* mito echinoderms			*/

/* ==================================================== */
/* Data Structures					*/
/* ==================================================== */

typedef struct GeneticCode {
    char   title[256];	/* nom du code		*/
    char   info[256];	/* informations		*/
    int	   code[65];	/* 64 codons + Error	*/
} GeneticCode;

/* ==================================================== */
/* Data Instanciation					*/
/* ==================================================== */

#ifdef GENETIC_CODE_INSTANCE

GeneticCode theGeneticCode[GEN_MAX_CODES] = {	/* Begin of codes */

{ /* 0: UNIVERSAL */
	{"Universal"},
	{""},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 20, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 1: MITOCHONDRIAL YEAST */
	{"Mitochondrial Yeast"}, 
	{"CUN=T  AUA=M  UGA=W"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 18, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  3,  3,  3,  3, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 19, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 2: MITOCHONDRIAL VERTEBRATES */
	{"Mitochondrial Vertebrates"}, 
	{"AGR=*  AUA=M  UGA=W"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	 20,  2, 20,  2, 18, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 19, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 3:	MITOCHONDRIAL FILAMENTOUS FUNGI */
	{"Mitochondrial Filamentous Fungi"}, 
	{"UGA=W"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 19, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 4:	MITOCHONDRIAL CODE OF INSECT AND PLATYHELMINTHES  */
	{"Mitochondrial Insects and Platyhelminthes"}, 
	{"AUA=M  UGA=W  AGR=S"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  2,  2,  2,  2, 18, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 19, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 5:	CANDIDA CYLINDRACEA (see nature 341:164) */
	{"Candida cylindracea"}, 
	{"CUG=S  CUA=?"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0, 21,  1,  2,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 20, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 6:	CILIATA	*/
	{"Ciliata"}, 
	{"UAR=Q"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 10, 14, 10, 14,  2,  2,  2,  2, 
	 20, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 7:	EUPLOTES */
	{"Euplotes"}, 
	{"UGA=C"},
	{ 8,  9,  8,  9,  3,  3,  3,  3, 
	  0,  2,  0,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 15, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}, 

{ /* 8:	MITOCHONDRIAL ECHINODERMS */
	{"Mitochondrial Echinoderms"}, 
	{"UGA=W  AGR=S  AAA=N"},
	{ 9,  9,  8,  9,  3,  3,  3,  3, 
	  2,  2,  2,  2, 17, 17, 18, 17, 
	 10, 11, 10, 11,  4,  4,  4,  4, 
	  0,  0,  0,  0,  1,  1,  1,  1, 
	 12, 13, 12, 13,  5,  5,  5,  5, 
	  6,  6,  6,  6,  7,  7,  7,  7, 
	 20, 14, 20, 14,  2,  2,  2,  2, 
	 19, 15, 19, 15,  1, 16,  1, 16, 
	 21}
}

/* end of codes */ };

#else

extern GeneticCode theGeneticCode[GEN_MAX_CODES];

#endif
#endif
