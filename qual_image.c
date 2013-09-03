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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"
#include <png.h>

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
	int full_tile;
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
	retval->full_tile = atoi(token);
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

static void dump_tile_isize_avg( welford_stat_t tile_grid_isize[2][3][16] )
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_isize[surface][swath][0].samples != 0)
			{
				printf( "%f", tile_grid_isize[surface][swath][0].current_mean );
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_isize[surface][swath][tile].samples == 0) break;
					printf( "\t%f", tile_grid_isize[surface][swath][tile].current_mean );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_isize[surface][swath][0].samples != 0)
			{
				printf( "%f", welford_sd(&tile_grid_isize[surface][swath][0]));
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_isize[surface][swath][tile].samples == 0) break;
					printf( "\t%f", welford_sd(&tile_grid_isize[surface][swath][tile]) );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}

}

static void dump_tile_qual_avg( welford_stat_t tile_grid_qual[2][3][16][100], int qual )
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_qual[surface][swath][0][qual].samples != 0)
			{
				printf( "%f", tile_grid_qual[surface][swath][0][qual].current_mean );
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_qual[surface][swath][tile][qual].samples == 0) break;
					printf( "\t%f", tile_grid_qual[surface][swath][tile][qual].current_mean );
				}
				printf( "\n" );
			}
		}
		printf( "\n" );
	}
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			if (tile_grid_qual[surface][swath][0][qual].samples != 0)
			{
				printf( "%f", welford_sd(&tile_grid_qual[surface][swath][0][qual]));
				for (tile = 1; tile < 16; ++tile) {
					if (tile_grid_qual[surface][swath][tile][qual].samples == 0) break;
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

static void clear_isize(welford_stat_t tile_grid_isize[2][3][16])
{
	int surface, swath;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			memset(tile_grid_isize[surface][swath], 0, sizeof(welford_stat_t)*16);
		}
	}
}


static void bam_qualview_core(samFile* in, FILE* output[2][3][16], png_structp png_ptr[2][3][16], png_infop info_ptr[2][3][16] )
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
	welford_stat_t tile_grid_isize[2][3][16];

	clear_tq(tile_grid_qual);
	clear_isize(tile_grid_isize);
	int full_tile;
	png_structp current_png;
	png_infop current_png_info;
	png_bytepp current_bitmap;
	
	if (sam_read1(in, hdr, b) >= 0) {
		parsed_readname_t* parse = parse_readname(bam_get_qname(b));
		full_tile = parse->full_tile;
		current_png = png_ptr[parse->surface][parse->swath][parse->tile];
		current_png_info = info_ptr[parse->surface][parse->swath][parse->tile];
		parsed_readname_destroy(parse);
	
		current_bitmap = (png_bytepp)calloc(2048, sizeof(png_bytep));
		int i;
		for (i = 0; i < 2048; ++i) {
			current_bitmap[i] = (png_bytep)calloc(10000,png_get_rowbytes(current_png,current_png_info));
		}
	}
	else {
		bam_destroy1(b);
		return;
	}
	do {
		if (b->core.flag&BAM_FSECONDARY)
			continue;
		++count;
		parsed_readname_t* parse = parse_readname(bam_get_qname(b));
		int read = 0;
		if ((b->core.flag&(BAM_FREAD1|BAM_FREAD2)) == (BAM_FREAD1|BAM_FREAD2)) {
		} else if (b->core.flag&BAM_FREAD1) {
			read = 1;
		} else if (b->core.flag&BAM_FREAD2) {
			read = 2;
		}
		if (full_tile != parse->full_tile)
		{
			png_write_image(current_png, current_bitmap);
			png_write_end(current_png, current_png_info);
			int i;
			for (i = 0; i < 2048; ++i) {
				memset(current_bitmap[i], 0, png_get_rowbytes(current_png,current_png_info)*10000);
			}
			full_tile = parse->full_tile;
			current_png = png_ptr[parse->surface][parse->swath][parse->tile];
			current_png_info = info_ptr[parse->surface][parse->swath][parse->tile];
		}

		++surface[parse->surface];
		++swath[parse->surface][parse->swath];
		++tile_grid[parse->surface][parse->swath][parse->tile];
		bam_glomp(b, tile_grid_qual[parse->surface][parse->swath][parse->tile]);
		if ((b->core.flag&BAM_FUNMAP) == 0) {
			if ((b->core.flag&BAM_FPROPER_PAIR) != 0 && (b->core.isize > 0)) {
				welford_add(&tile_grid_isize[parse->surface][parse->swath][parse->tile], (double)b->core.isize);
				fprintf(output[parse->surface][parse->swath][parse->tile], "%d\t%d\t%d\t%d\n", bam_get_qual(b)[99], parse->x, parse->y, read);
				
				current_bitmap[parse->x][parse->y] = bam_get_qual(b)[99];
			}

			tile_grid_mq[parse->surface][parse->swath][parse->tile] += b->core.qual;
		}
		
		parsed_readname_destroy(parse);
	} while (sam_read1(in, hdr, b) >= 0);

	int i;
	for (i = 0; i < 2048; ++i) {
		free(current_bitmap[i]);
	}
	free(current_bitmap);
	
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
	dump_tile_qual_avg(tile_grid_qual,99);
	dump_tile_isize_avg(tile_grid_isize);
	dump_tile_mq_avg(tile_grid,tile_grid_mq);
	
	bam_destroy1(b);
}

static void usage()
{
	fprintf(stderr,"Usage information\n");
}

static bool open_files(const char* prefix, FILE* out[2][3][16], FILE* png[2][3][16], png_structp png_ptr[2][3][16], png_infop info_ptr[2][3][16])
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			for (tile = 0; tile < 16; ++tile) {
				char buf[255];
				char buf_png[255];
				sprintf(buf, "%s_%d_%d_%d.tsv",prefix, surface, swath, tile);
				sprintf(buf_png, "%s_%d_%d_%d.png",prefix, surface, swath, tile);
				out[surface][swath][tile] = fopen(buf, "w");
				png[surface][swath][tile] = fopen(buf_png, "wb");
				
				png_ptr[surface][swath][tile] = png_create_write_struct
				(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
				
				if (!png_ptr[surface][swath][tile])
					return false;
				
				info_ptr[surface][swath][tile] = png_create_info_struct(png_ptr[surface][swath][tile]);
				if (!info_ptr[surface][swath][tile])
				{
					png_destroy_write_struct(&png_ptr[surface][swath][tile],
											 (png_infopp)NULL);
					return false;
				}
				
				png_init_io(png_ptr[surface][swath][tile], png[surface][swath][tile]);
				png_set_IHDR(png_ptr[surface][swath][tile], info_ptr[surface][swath][tile], 2048,10000, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
			}
		}
	}
	return true;
}

static void close_files(FILE* out[2][3][16], FILE* png[2][3][16], png_structp png_ptr[2][3][16])
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			for (tile = 0; tile < 16; ++tile) {
				fclose(out[surface][swath][tile]);
				fclose(png[surface][swath][tile]);
			}
		}
	}
}


#if 0
int bam_qualview(int argc, char *argv[])
#else
int main(int argc, char *argv[])
#endif
{
	int c;
	samFile* in;
	FILE* out[2][3][16];
	FILE* png[2][3][16];
	png_structp png_ptr[2][3][16];
	png_infop info_ptr[2][3][16];
	while ((c = getopt(argc, argv, "")) >= 0) {
		switch (c) {
		}
	}
	if (optind+1 >= argc) usage();
	in = sam_open(argv[optind], "rb", NULL);
	open_files(argv[optind+1],out, png, png_ptr, info_ptr);
	bam_qualview_core(in, out, png_ptr, info_ptr);
	sam_close(in);
	close_files(out, png, png_ptr);
	return 0;
}

