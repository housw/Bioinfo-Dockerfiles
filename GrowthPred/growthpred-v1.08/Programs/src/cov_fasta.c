/* ==================================================== */
/*	Copyright (c) Atelier de BioInformatique	*/
/*	Dec. 92 					*/
/* 	File: cov_fasta.c				*/
/*	Purpose: sequence IO in fasta format		*/
/*     							*/
/*	History:					*/
/*	15/08/92 : <Gloup> first version		*/
/*	00/08/95 : <Gloup> modified version for covgen	*/
/* ====================================================	*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Gtypes.h"

#include "cov_fasta.h"

#define DEBUG	0

#define READ_NEXT Faux
#define PUSH_BACK Vrai

#ifdef MACINTOSH
#define LINE_FEED '\r'
#else
#define LINE_FEED '\n'
#endif

/* -------------------------------------------- */
/* @static: lecture bufferisee			*/
/* -------------------------------------------- */
static char * sNextIOBuffer(streamin, retain)
	FILE *streamin;
	Bool retain;
{
	char   *buf, *end;
	
	static char sBuffer[BUFSIZ];	/* in <stdio.h>	  */
	static Bool sRetained = Faux;
	
	buf = ((retain || sRetained)
			? sBuffer 
			: fgets(sBuffer, sizeof(sBuffer), streamin));	 	

	if (buf) {
	   end = buf + strlen(buf) - 1;
	   if (*end == LINE_FEED) *end = '\000';
	}

	sRetained = retain;
	
	return buf;
}

/* -------------------------------------------- */
/* compte le nombre de caracteres alpha	dans	*/
/* un buffer					*/
/* -------------------------------------------- */
Int32 CountAlpha(buf)
	char *buf;
{
	Int32 count;
	
	for (count = 0 ; *buf ; buf++)
	    if (isalpha(*buf))
		count++;
	
	return count;
}


/* -------------------------------------------- */
/* copy only alpha chars from s2 to s1		*/
/* -------------------------------------------- */
char * StrcpyAlpha(s1, s2)
	char *s1, *s2;
{
	for( ; *s2 ; s2++)
	    if (isalpha(*s2))
	    	*s1++ = *s2;

	*s1 = '\000';

	return s1;
}

/* -------------------------------------------- */
/* skip to next space in buffer			*/
/* -------------------------------------------- */
char * NextSpace(buffer)
	char *buffer;
{
	for (; *buffer ; buffer++)
	   if (isspace(*buffer))
	   	return buffer;
	
	return NULL;
}

/* -------------------------------------------- */
/* returns sequence name (FASTA)		*/
/* -------------------------------------------- */
char *GetFastaName(buffer)
	char *buffer;
{
	static char name[FASTA_NAMLEN];
	
	if (sscanf(buffer + 1, "%s", name) != 1)
	    strcpy(name, "<no Name>");

	return name;
}

/* -------------------------------------------- */
/* returns sequence comment (FASTA)		*/
/* -------------------------------------------- */
char *GetFastaComment(buffer)
	char *buffer;
{
	char   *space;
	static char comment[FASTA_COMLEN];
	
	space = NextSpace(buffer);
	
	if (space) {
	   if (strlen(space+1) >= FASTA_COMLEN-2) {
		strncpy(comment, space + 1, FASTA_COMLEN-2);
		comment[FASTA_COMLEN-1] = '\000';
	   }
	   else
		strcpy(comment, space + 1);
	}
	else
	   strcpy(comment, "<no comment>");

	return comment;
}

/* -------------------------------------------- */
/* liberation d'une sequence			*/
/* -------------------------------------------- */
FastaSequencePtr FreeFastaSequence(seq)
	FastaSequencePtr seq;
{
	if (seq) {
	    if (seq->seq)  FREE(seq->seq);
	    FREE(seq);
	}

	return NULL;
}
	
/* -------------------------------------------- */
/* allocation d'une sequence			*/
/* -------------------------------------------- */
FastaSequencePtr NewFastaSequence()
{
	FastaSequencePtr seq;
	
	if (! (seq = NEW(FastaSequence)))
	    return NULL;
	   
	seq->length   = 0;

	seq->seq = seq->cover = NULL;

	if (! (seq->seq = NEWN(char,  BUFSIZ)))
	    return FreeFastaSequence(seq);
	    
	if (! (seq->cover = NEWN(char,  BUFSIZ)))
	    return FreeFastaSequence(seq);

	seq->bufsize = BUFSIZ;

	*(seq->name)    = '\000';
	*(seq->comment) = '\000';

	seq->ok  = Vrai;
	
	return seq;
}

/* -------------------------------------------- */
/* lecture/redimensionnement d'une sequence au	*/
/* format Fasta					*/
/* returns : Faux -> last sequence		*/
/*	     Vrai -> more to read		*/
/*           <but> you must check seq->ok !	*/
/* -------------------------------------------- */
Bool ReadFastaSequence(streamin, seq)
	FILE *streamin;
	FastaSequencePtr seq;
{
	Int32	buflen;
	char 	*buffer, *tbuf;

	seq->ok = Faux;				/* assume error		*/

	buflen = seq->length = 0; 
	
	buffer = sNextIOBuffer(streamin, READ_NEXT);

	if (! (buffer && (*buffer == '>'))) 	/* sync error		*/
	    return Faux;			/* last sequence	*/
	
	strcpy(seq->name,    GetFastaName(buffer));
	
	strcpy(seq->comment, GetFastaComment(buffer));
	
	while (buffer = sNextIOBuffer(streamin, READ_NEXT)) {

	    if (*buffer == '>') {			   /* push buf	*/
	    	(void) sNextIOBuffer(streamin, PUSH_BACK); /* back	*/
		break;
	    }

	    buflen += CountAlpha(buffer);

	    if (buflen >= seq->bufsize) {
	    
	    	if (! (tbuf = REALLOC(char, seq->seq, buflen + 1)))
	    	   return Vrai;			/* but seq->ok is Faux	 */

		seq->seq = tbuf;
		
		seq->bufsize = buflen + 1;
		
	    }		

	    StrcpyAlpha(seq->seq + seq->length, buffer);

	    seq->length = buflen;
	
	}

	if (! (tbuf = REALLOC(char, seq->cover, seq->length + 1)))
	    return Vrai;			/* but seq->ok is Faux	*/

	seq->cover = tbuf;
	
	seq->seq[seq->length] = '\000';

	return (seq->ok = Vrai);
}

/* -------------------------------------------- */
/* ecriture d'une sequence au format Fasta	*/
/* -------------------------------------------- */
void WriteFastaSequence(streamou, seq, char_per_line)
	FILE  *streamou;
	FastaSequencePtr seq;
	Int32 char_per_line;
{
	Int32 i, nlines, rest;
	char *buf, *end, tempo;

	fputc('>', streamou);
	fputs((*(seq->name)    ? seq->name    : "<no name>")   , streamou);
	fputc(' ', streamou);
	fputs((*(seq->comment) ? seq->comment : "<no comment>"), streamou);
	fputc(LINE_FEED, streamou);

	nlines = seq->length / char_per_line;

	buf = seq->seq;

	for (i = 0 ; i < nlines ; i++) {
	    end = buf + char_per_line;
	    tempo = *end;
	    *end = '\000';
	    fputs(buf, streamou);
	    fputc(LINE_FEED , streamou);
	    *end = tempo;
	    buf += char_per_line;
	}

	if (rest = (seq->length % char_per_line)) {
	   end = buf + rest;
	   tempo = *end;
	   *end = '\000';
	   fputs(buf, streamou);
	   fputc(LINE_FEED , streamou);
	   *end = tempo;
	}
}

