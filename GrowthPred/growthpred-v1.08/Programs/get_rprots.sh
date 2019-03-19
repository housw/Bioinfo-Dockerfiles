#!/bin/bash 

# $1 is genome from which to find ribosomal-protein genes / fasta format (*.gen)
# $2 is the ribosomal database with full path
# $3 is the outfile with full path
# $4 blastall + full path
# $5 getentryF + full path
# needs blastall, rprot database, getentryF

# EDITING TO ADD LENGTHS

awk 'BEGIN{name="";seq=""} {if (NR == 1) {name=$1} else {if (substr($1,1,1) == ">") {print name"__"length(seq)/3; print seq; name=$1; seq="";long=0} else {seq=seq$0}}} END{print name"__"length(seq); print seq}' $1 > $3.plen.tmp

# BLASTING

$4 -p blastx -d $2 -i $3.plen.tmp -m 8 -e 0.00001 > $3.blast.all

sed 's/__/ /g' $3.blast.all | awk '{if ($2 < $4*1.10 && $2 > $4*0.90 && $5 > 40) {print $1}}' - | sort -u > $3.ids.lst

# RETRIEVING GENES FROM IDS 

$5 -f $3.ids.lst $1 > $3.ribs

rm -f $3.plen.tmp $3.ids.lst $3.blast.all

exit

