gcc -g -Wall  -O2 -o samtools_qual -I ../htslib qual_image.c libbam.a -Lbcftools -lbcf  ../htslib/libhts.a -lcurses  -lm -lz -lpthread 
