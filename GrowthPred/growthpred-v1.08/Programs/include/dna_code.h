/* -------------------------------------------- */
/* @file: dna_code.h				*/
/* IUPAC encoding for DNA			*/
/* -------------------------------------------- */
/* -----------------------------------------	*/
/* direct encoding				*/
/* -----------------------------------------	*/
/* G/A/T/C					*/
/* U=T						*/
/* R=AG						*/
/* Y=CT						*/
/* M=AC						*/
/* K=GT						*/
/* S=CG						*/
/* W=AT						*/
/* H=ACT					*/
/* B=CGT					*/
/* V=ACG					*/
/* D=AGT					*/
/* N=ACGT					*/
/* X=ACGT					*/
/* EFIJLOPQZ  not recognized			*/
/* -----------------------------------------	*/
/* dual encoding				*/
/* -----------------------------------------	*/
/* A=ADHMNRVW					*/
/* B=BCDGHKMNRSTUVWY				*/
/* C=BCHMNSVY					*/
/* D=ABDGHKMNRSTUVWY				*/
/* G=BDGKNRSV					*/
/* H=ABCDHKMNRSTUVWY				*/
/* K=BDGHKNRSTUVWY				*/
/* M=ABCDHMNRSVWY				*/
/* N=ABCDGHKMNRSTUVWY				*/
/* R=ABDGHKMNRSVW				*/
/* S=BCDGHKMNRSVY				*/
/* T=BDHKNTUWY					*/
/* U=BDHKNTUWY					*/
/* V=ABCDGHKMNRSVWY				*/
/* W=ABDHKMNRTUVWY				*/
/* X=ABCDGHKMNRSTUVWY				*/
/* Y=BCDHKMNSTUVWY				*/
/* EFIJLOPQZ  not recognized			*/
/* -------------------------------------------- */

#ifdef IUPAC_DIRECT

	/* IUPAC Direct */

        0x00000001 /* A */, 0x00080044 /* B */, 0x00000004 /* C */, 
        0x00080041 /* D */, 0x00000000 /* E */, 0x00000000 /* F */, 
        0x00000040 /* G */, 0x00080005 /* H */, 0x00000000 /* I */, 
        0x00000000 /* J */, 0x00080040 /* K */, 0x00000000 /* L */, 
        0x00000005 /* M */, 0x00080045 /* N */, 0x00000000 /* O */, 
        0x00000000 /* P */, 0x00000000 /* Q */, 0x00000041 /* R */, 
        0x00000044 /* S */, 0x00080000 /* T */, 0x00080000 /* U */, 
        0x00000045 /* V */, 0x00080001 /* W */, 0x00080045 /* X */, 
        0x00080004 /* Y */, 0x00000000 /* Z */

#define IUPAC_DEFINED
#endif

#ifdef IUPAC_DUAL

	/* IUPAC Dual  */

	0x00623089 /* A */, 0x017e34ce /* B */, 0x01243086 /* C */, 
	0x017e34cb /* D */, 0x00000000 /* E */, 0x00000000 /* F */, 
	0x0026244a /* G */, 0x017e348f /* H */, 0x00000000 /* I */, 
	0x00000000 /* J */, 0x017e24ca /* K */, 0x00000000 /* L */, 
	0x0166308f /* M */, 0x017e34cf /* N */, 0x00000000 /* O */, 
	0x00000000 /* P */, 0x00000000 /* Q */, 0x006634cb /* R */, 
	0x012634ce /* S */, 0x0158248a /* T */, 0x0158248a /* U */, 
	0x016634cf /* V */, 0x017a348b /* W */, 0x017e34cf /* X */, 
	0x017c348e /* Y */, 0x00000000 /* Z */

#define IUPAC_DEFINED
#endif

#ifdef TWO_BITS

	/* just ACGT and 1 error bit (bit 5 ; msk=0x10) */

        0x00000000 /* A */, 0x00000010 /* B */, 0x00000001 /* C */, 
        0x00000010 /* D */, 0x00000010 /* E */, 0x00000010 /* F */, 
        0x00000002 /* G */, 0x00000010 /* H */, 0x00000010 /* I */, 
        0x00000010 /* J */, 0x00000010 /* K */, 0x00000010 /* L */, 
        0x00000010 /* M */, 0x00000010 /* N */, 0x00000010 /* O */, 
        0x00000010 /* P */, 0x00000010 /* Q */, 0x00000010 /* R */, 
        0x00000010 /* S */, 0x00000003 /* T */, 0x00000003 /* U */, 
        0x00000010 /* V */, 0x00000010 /* W */, 0x00000010 /* X */, 
        0x00000010 /* Y */, 0x00000010 /* Z */

#define IUPAC_DEFINED
#endif


#ifndef IUPAC_DEFINED
/*
 | #error  either IUPAC_DIRECT IUPAC_DUAL or TWO_BITS 
 | #error	should be defined before inclusion of this file
 */
#else
#undef  IUPAC_DEFINED
#endif



