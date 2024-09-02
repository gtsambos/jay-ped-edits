##Generate input for the tool out of hapIBD output
##Input: Argument 1: name of input HapIBD output, Argument 2: prefix of output file name
##Output: two files, ${2}.first_order.seg and ${2}.last_order.seg, two files ordered by base coordinate of first and last markers of segments, respectively.
##NOTE: Currently expects single-chrom input

sort -k 5 -n $1 > $1.first_order
sort -k 6 -n $1 > $1.last_order
