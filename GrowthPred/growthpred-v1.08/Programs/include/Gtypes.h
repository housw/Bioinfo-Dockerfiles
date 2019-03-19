/* ==================================================== */
/* Gloup Mar. 90 					*/
/* Copyright (c) Atelier de BioInformatique		*/
/* @file: Gtypes.h					*/
/* general (machine dependant types)			*/
/* ==================================================== */


#ifndef NULL
#include <stdio.h>
#endif

#define _H_Gtypes

/* ==================================================== */
/* constantes						*/
/* ==================================================== */

#ifndef PROTO
#define PROTO	1			/* prototypes flag		*/
#endif

#define Vrai	0x1			/* bool values	= TRUE		*/
#define Faux	0x0			/*              = FALSE		*/

#define Nil	NULL		 	/* nil pointer			*/

#define	kBigInt16	0x7fff	    	/* plus grand 16 bits signe	*/
#define	kBigInt32	0x7fffffff  	/* plus grand 32 bits signe	*/
#define	kBigUInt16	0xffff	    	/* plus grand 16 bits ~signe	*/
#define	kBigUInt32	0xffffffff  	/* plus grand 32 bits ~signe	*/

#define kBitsPerLong	32		/*  long = 32 bits		*/
#define kMaxShftLong	31		/*  BitsPerLong - 1 max shift	*/
#define kLog2BitLong	5		/*  =log2(BitsPerLong)		*/

#ifdef MACINTOSH
/* ==================================================== */
/*  Types (for Macintosh ThinK C)			*/
/* ==================================================== */

typedef long		Long;		/* plus grand mot signe		*/
typedef unsigned long	ULong;		/* plus grand mot signe		*/
typedef long		Int32;		/* Int32  = 32 bits signe	*/
typedef unsigned long	UInt32;		/* UInt32 = 32 bits ~signe	*/
typedef short		Int16;		/* Int16  = 16 bits signe	*/
typedef unsigned short	UInt16;		/* UInt32 = 16 bits ~signe	*/
typedef char		Int8;		/* Int8   = 8 bits signe	*/
typedef unsigned char	UInt8;		/* UInt8  = 8 bits ~signe	*/

typedef Boolean		Bool;		/* booleen  			*/

#else
/* ==================================================== */
/*  Types (for Sun & Iris)				*/
/* ==================================================== */

typedef long		Long;		/* plus grand mot signe		*/
typedef unsigned long	ULong;		/* plus grand mot signe		*/
typedef int		Int32;		/* Int32  = 32 bits signe	*/
typedef unsigned int	UInt32;		/* UInt32 = 32 bits ~signe	*/
typedef short		Int16;		/* Int16  = 16 bits signe	*/
typedef unsigned short	UInt16;		/* UInt32 = 16 bits ~signe	*/
typedef char		Int8;		/* Int8   = 8 bits signe	*/
typedef unsigned char	UInt8;		/* UInt8  = 8 bits ~signe	*/

typedef int		Bool;		/* booleen  (int for ANSI)	*/
		
typedef void 		*Ptr;		/* pointeur			*/
#endif


/* ==================================================== */
/*  special macro for prototypes			*/
/* ==================================================== */

#if PROTO
#define		P(s) 	s
#else
#define 	P(s) 	()
#endif
