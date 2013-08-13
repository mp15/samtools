// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of samtools.
//
// samtools is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"

typedef struct parsed_readname {
	char *machine_name;
	int run_id;
	int lane;
	int surface; // 0 based
	int swathe; // 0 based
	int tile; // 0 based
	int x;
	int y;
	int index;
} parsed_readname_t;

static parsed_readname_t* parse_readname(const char* readname)
{
	// very crude parser to parse read names of form
	// <machine_name>_<run_id>:<lane>:<tile>:<x>:<y>#<index>
	parsed_readname_t* retval = (parsed_readname_t*)malloc(sizeof(parsed_readname_t));
	char* tofree;
	char* toparse;
	char* token;
	tofree = toparse = strdup(readname);
	
	token = strsep(&toparse,":");
	retval->machine_name = strdup(strsep(&token,"_"));
	retval->run_id = atoi(token);
	retval->lane = atoi(strsep(&toparse,":"));
	token = strsep(&toparse,":");
	retval->surface = token[0] - '1';
	retval->swathe = token[1] - '1';
	retval->tile = atoi(token+2);
	retval->x = atoi(strsep(&toparse,":"));
	retval->y = atoi(strsep(&toparse,"#"));
	retval->index = atoi(toparse);
	free(tofree);
	
	return retval;
}

static void parsed_readname_destroy(parsed_readname_t* destroy)
{
	free(destroy->machine_name);
	free(destroy);
}

static void bam_qualview_core(samFile* in, FILE* out)
{
	bam1_t* b = bam_init1();
	bam_hdr_t* hdr = sam_hdr_read(in);
	int surface_1 = 0, surface_2 = 0;
	while (sam_read1(in, hdr, b) != 0) {
		char* readname = bam_get_qname(b);
		parsed_readname_t* parse = parse_readname(readname);
		switch (parse->surface) {
			case 1:
				++surface_1;
				break;
			case 2:
				++surface_2;
				break;
		}
		parsed_readname_destroy(parse);
	}
	printf("Surface 1: %d\n"
		   "Surface 2: %d\n", surface_1, surface_2);

	bam_destroy1(b);
}

static void usage()
{
	fprintf(stderr,"Usage information\n");
}

#if 0
int bam_qualview(int argc, char *argv[])
#else
int main(int argc, char *argv[])
#endif
{
	int c;
	samFile* in;
	FILE* out;
	while ((c = getopt(argc, argv, "")) >= 0) {
		switch (c) {
		}
	}
	if (optind+1 >= argc) usage();
	in = sam_open(argv[optind], "rb", NULL);
	out = (strcmp(argv[optind+1], "-") == 0)? stdout : fopen(argv[optind+1], "w");
	bam_qualview_core(in, out);
	sam_close(in);
	fclose(out);
	return 0;
}

