#!/bin/bash

#based on frequency of IVYWREL (Zeldovich, 2007, Plos Computational Biology)
# $1 = genetic code
# $2 = dirfile/file
# $3 = dirbin

$3translate -c $1 $2 2>&1| grep -v "bad length" | $3util_tab_aa -r -t  2>&1| grep -v "Invalid" | tail -1 | awk '{var= $9+$19+$21+$20+$16+$5+$11; print (937*(var/100))-335}'

exit



