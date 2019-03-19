#
# (getentryF:translate:utilCUB:util_tab_aa) Makefile
# 
#

#
# Machine dependant flags
# defaults
#

MACHINE = MACOSX 

MALLOC = 

PROTO = PROTO=1

CCKR = 

#
# Machine independant flags
#

INCDIR   = ../include
BINDIR   = ..
	
PROG	 = 	utilCUB \
	translate \
	util_tab_aa \
	getentryF \
	wcountq 

OBJS     =	cub_fasta.o

OBJ     =	wcount_bgener.o		\
		wcount_hgener.o		\
		wcount_io.o		\
		cov_fasta.o
		
OPTIM = -O -w

CFLAGS 	 = 	$(OPTIM) -D$(MACHINE) -D$(PROTO) -I$(INCDIR)  $(CCKR)

LDFLAGS  = 	$(MALLOC) -lm 

#
# Rules
#

.KEEP_STATE :

all	   : $(PROG) 

utilCUB :  utilCUB.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(@).c $(OBJS) $(LDFLAGS)	
	
util_tab_aa :  util_tab_aa.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(@).c $(OBJS) $(LDFLAGS)	

translate :  translate.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(@).c $(OBJS) $(LDFLAGS)
			
getentryF:	getentryF_main.c $(OBJS) 
	$(CC) $(CFLAGS) -o $@  $(@)_main.c $(OBJS) 

wcountq :  wcountq.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(@).c $(OBJ) $(LDFLAGS)	


install :   $(PROG)
	-\mv utilCUB $(BINDIR)/utilCUB
	-\mv translate $(BINDIR)/translate
	-\mv util_tab_aa $(BINDIR)/util_tab_aa
	-\mv getentryF $(BINDIR)/getentryF
	-\mv wcountq $(BINDIR)/wcountq
	
clean :
	-\rm -f $(PROG) test_? cvstatic* *% *.o
