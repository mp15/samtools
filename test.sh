gcc -g -Wall  -O2 -o samtools_qual -I /opt/local/include/libpng15 -I ../htslib qual_image.c libbam.a -Lbcftools -lbcf -L/opt/local/lib  ../htslib/libhts.a -lcurses  -lm -lz -lpthread -lpng
