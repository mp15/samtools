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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"

typedef struct welford_stat {
	int samples;
	double current_mean;
	double current_sd;
} welford_stat_t;

void welford_add(welford_stat_t* stat, double value) {
	stat->samples++;
	if ( stat->samples == 1 ) {
		stat->current_mean = value;
		stat->current_sd = 0.0;
	} else {
		double old_mean = stat->current_mean;
		stat->current_mean = stat->current_mean + (value - stat->current_mean)/stat->samples;
		stat->current_sd = stat->current_sd + (value - old_mean)*(value - stat->current_mean);
	}
}

double welford_sd(welford_stat_t* stat) {
	if (stat->samples == 0) return 0.0;
	double variance = stat->current_sd/(stat->samples-1);
	return sqrt(variance);
}

typedef struct parsed_readname {
	char *machine_name;
	int run_id;
	int lane;
	int surface; // 0 based
	int swath; // 0 based
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
	retval->swath = token[1] - '1';
	retval->tile = atoi(token+2)-1;
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

static void dump_tile( int tile_grid[2][3][16] )
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid[surface][swath][0])
			{
				printf( "%d", tile_grid[surface][swath][0] );
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid[surface][swath][tile] == 0) break;
					printf( "\t%d", tile_grid[surface][swath][tile] );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	printf( "\n" );
}

static void dump_tile_avg( welford_stat_t tile_grid_qual[2][3][16][100], int qual )
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_qual[surface][swath][0]->samples != 0)
			{
				printf( "%f", tile_grid_qual[surface][swath][0][qual].current_mean );
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_qual[surface][swath][tile]->samples == 0) break;
					printf( "\t%f", tile_grid_qual[surface][swath][tile][qual].current_mean );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_qual[surface][swath][0]->samples != 0)
			{
				printf( "%f", welford_sd(&tile_grid_qual[surface][swath][0][qual]));
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_qual[surface][swath][tile]->samples == 0) break;
					printf( "\t%f", welford_sd(&tile_grid_qual[surface][swath][tile][qual]) );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	printf( "\n" );
}

static void dump_tile_mq_avg( int tile_grid[2][3][16], int tile_grid_mq[2][3][16] )
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid[surface][swath][0])
			{
				printf( "%d", tile_grid_mq[surface][swath][0]/tile_grid[surface][swath][0] );
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid[surface][swath][tile] == 0) break;
					printf( "\t%d", tile_grid_mq[surface][swath][tile]/tile_grid[surface][swath][tile] );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	printf( "\n" );
}

static void bam_glomp(bam1_t *b, welford_stat_t tile_grid_qual[100]) {
	int qual;
	for (qual = 0; qual < b->core.l_qseq; ++qual) {
		 welford_add(&tile_grid_qual[qual], bam_get_qual(b)[qual]);
	}
}

static void clear_tq(welford_stat_t tile_grid_qual[2][3][16][100])
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			for (tile = 0; tile < 16; ++tile) {
				memset(tile_grid_qual[surface][swath][tile], 0, sizeof(welford_stat_t)*100);
			}
		}
	}
}

static void bam_qualview_core(samFile* in, FILE* out)
{
	bam1_t* b = bam_init1();
	bam_hdr_t* hdr = sam_hdr_read(in);
	int count = 0;
	int surface[2] = {0, 0};
	int swath[2][3] = {
		{ 0, 0, 0, },
		{ 0, 0, 0, },
	};
	int tile_grid [2][3][16] = {
		{ {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, },
		{ {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, },
	};
	int tile_grid_mq [2][3][16] = {
		{ {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, },
		{ {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, },
	};

	welford_stat_t tile_grid_qual[2][3][16][100];
	clear_tq(tile_grid_qual);
				
	while (sam_read1(in, hdr, b) >= 0) {
		++count;
		char* readname = bam_get_qname(b);
		parsed_readname_t* parse = parse_readname(readname);
		++surface[parse->surface];
		++swath[parse->surface][parse->swath];
		++tile_grid[parse->surface][parse->swath][parse->tile];
		bam_glomp(b, tile_grid_qual[parse->surface][parse->swath][parse->tile]);
		if ((b->core.flag&BAM_FUNMAP) == 0) {
			tile_grid_mq[parse->surface][parse->swath][parse->tile] += b->core.qual;
		}
		
		parsed_readname_destroy(parse);
	}
	
	// summary
	printf("Count: %d\n"
		   "Surface 1: %d\n"
		   "\tSwath 1: %d\n"
		   "\tSwath 2: %d\n"
		   "\tSwath 3: %d\n"
		   "Surface 2: %d\n"
		   "\tSwath 1: %d\n"
		   "\tSwath 2: %d\n"
		   "\tSwath 3: %d\n",
		   count,
		   surface[0],
		   swath[0][0],
		   swath[0][1],
		   swath[0][2],
		   surface[1],
		   swath[1][0],
		   swath[1][1],
		   swath[1][2]
		   );
	dump_tile(tile_grid);
	dump_tile_avg(tile_grid_qual,99);
	dump_tile_mq_avg(tile_grid,tile_grid_mq);
	
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

