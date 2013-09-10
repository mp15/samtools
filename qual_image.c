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

const int X_LEN = 2048+100;
const int Y_LEN = 10000+100;

#ifdef FULL_PALETTE
#define QUAL_PALETTE_LEN 51
static const png_color QUAL_PALETTE[QUAL_PALETTE_LEN] = {
	{ .red = 0x00, .green = 0x00, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x00, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x07, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x0E, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x15, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x1C, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x22, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x29, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x30, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x37, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x3E, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x45, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x4C, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x53, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x5A, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x60, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x67, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x6E, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x75, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x7C, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x83, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x8A, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x91, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x98, .blue = 0x00 },
	{ .red = 0xFF, .green = 0x9F, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xA5, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xAC, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xB3, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xBA, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xC1, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xC8, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xCF, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xD6, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xDD, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xE3, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xEA, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xF1, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xF8, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x00 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x0B },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x20 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x35 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x4A },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x60 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x75 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x8A },
	{ .red = 0xFF, .green = 0xFF, .blue = 0x9F },
	{ .red = 0xFF, .green = 0xFF, .blue = 0xB5 },
	{ .red = 0xFF, .green = 0xFF, .blue = 0xCA },
	{ .red = 0xFF, .green = 0xFF, .blue = 0xDF },
	{ .red = 0xFF, .green = 0xFF, .blue = 0xF4 },
};
#else
#define QUAL_PALETTE_LEN 5
static const png_color QUAL_PALETTE[QUAL_PALETTE_LEN] = {
	{ .red = 0xFF, .green = 0xFF, .blue = 0xFF }, // None
	{ .red = 0xFF, .green = 0x00, .blue = 0x00 }, // 0-10
	{ .red = 0x00, .green = 0xFF, .blue = 0x00 }, // 10-20
	{ .red = 0x00, .green = 0x00, .blue = 0x7F }, // 20-30
	{ .red = 0x00, .green = 0x00, .blue = 0xFF }, // 30+
};
#endif



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


static bool bam_qualview_core(samFile* in, FILE* output[2][3][16], const char* prefix )
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
	png_structp current_png[2];
	png_infop current_png_info[2];
	png_bytepp current_bitmap[2];
	FILE* png[2];
	
	if (sam_read1(in, hdr, b) >= 0) {
		parsed_readname_t* parse = parse_readname(bam_get_qname(b));
		full_tile = parse->full_tile;
		
		char buf_png[2][255];
		sprintf(buf_png[0], "%s_%d_%d_%d_fwd.png", prefix, parse->surface, parse->swath, parse->tile);
		sprintf(buf_png[1], "%s_%d_%d_%d_rev.png", prefix, parse->surface, parse->swath, parse->tile);
		parsed_readname_destroy(parse);
		// Create images
		png[0] = fopen(buf_png[0], "wb");
		png[1] = fopen(buf_png[1], "wb");
		
		current_png[0] = png_create_write_struct
		(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		current_png[1] = png_create_write_struct
		(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		
		if (!current_png[0] || !current_png[1])
			return false;
		
		current_png_info[0] = png_create_info_struct(current_png[0]);
		if (!current_png_info[0])
		{
			png_destroy_write_struct(&current_png[0],
									 (png_infopp)NULL);
			return false;
		}
		current_png_info[1] = png_create_info_struct(current_png[1]);
		if (!current_png_info[1]) {
			png_destroy_write_struct(&current_png[0],
									 &current_png_info[0]);
			png_destroy_write_struct(&current_png[1],
									 (png_infopp)NULL);
			return false;
		}
		
		png_init_io(current_png[0], png[0]);
		png_init_io(current_png[1], png[1]);
		png_set_IHDR(current_png[0], current_png_info[0], X_LEN, Y_LEN, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
		png_set_IHDR(current_png[1], current_png_info[1], X_LEN, Y_LEN, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
		png_set_PLTE(current_png[0], current_png_info[0], QUAL_PALETTE, 51);
		png_set_PLTE(current_png[1], current_png_info[1], QUAL_PALETTE, 51);
		png_write_info(current_png[0], current_png_info[0]);
		png_write_info(current_png[1], current_png_info[1]);

		current_bitmap[0] = (png_bytepp)png_malloc(current_png[0], Y_LEN * sizeof(png_bytep));
		current_bitmap[1] = (png_bytepp)png_malloc(current_png[1], Y_LEN * sizeof(png_bytep));
		int i;
		for (i = 0; i < Y_LEN; ++i) {
			current_bitmap[0][i] = (png_bytep)png_malloc(current_png[0], png_get_rowbytes(current_png[0],current_png_info[0]));
			current_bitmap[1][i] = (png_bytep)png_malloc(current_png[1], png_get_rowbytes(current_png[1],current_png_info[1]));
			memset(current_bitmap[0][i], 0, png_get_rowbytes(current_png[0],current_png_info[0]));
			memset(current_bitmap[1][i], 0, png_get_rowbytes(current_png[1],current_png_info[1]));
		}
	}
	else {
		bam_destroy1(b);
		return true;
	}
	do {
		if (b->core.flag&BAM_FSECONDARY)
			continue;
		parsed_readname_t* parse = parse_readname(bam_get_qname(b));
		int read = -1;
		if ((b->core.flag&(BAM_FREAD1|BAM_FREAD2)) == (BAM_FREAD1|BAM_FREAD2)) {
		} else if (b->core.flag&BAM_FREAD1) {
			read = 0;
		} else if (b->core.flag&BAM_FREAD2) {
			read = 1;
		}
		if (full_tile != parse->full_tile)
		{
			png_write_image(current_png[0], current_bitmap[0]);
			png_write_image(current_png[1], current_bitmap[1]);
			png_write_end(current_png[0], current_png_info[0]);
			png_write_end(current_png[1], current_png_info[1]);

			int i;
			for (i = 0; i < Y_LEN; ++i) {
				png_free(current_png[0], current_bitmap[0][i]);
				png_free(current_png[1], current_bitmap[1][i]);
			}
			png_free(current_png[0], current_bitmap[0]);
			png_free(current_png[1], current_bitmap[1]);

			png_destroy_write_struct(&current_png[0], &current_png_info[0]);
			png_destroy_write_struct(&current_png[1], &current_png_info[1]);
			fclose(png[0]);
			fclose(png[1]);

			current_bitmap[0] = (png_bytepp)png_malloc(current_png[0], Y_LEN * sizeof(png_bytep));
			current_bitmap[1] = (png_bytepp)png_malloc(current_png[1], Y_LEN * sizeof(png_bytep));
			for (i = 0; i < Y_LEN; ++i) {
				current_bitmap[0][i] = (png_bytep)png_malloc(current_png[0], png_get_rowbytes(current_png[0],current_png_info[0]));
				current_bitmap[1][i] = (png_bytep)png_malloc(current_png[1], png_get_rowbytes(current_png[1],current_png_info[1]));
				memset(current_bitmap[0][i], 0, png_get_rowbytes(current_png[0],current_png_info[0]));
				memset(current_bitmap[1][i], 0, png_get_rowbytes(current_png[1],current_png_info[1]));
			}
			full_tile = parse->full_tile;
			char buf_png[2][255];
			sprintf(buf_png[0], "%s_%d_%d_%d_fwd.png", prefix, parse->surface, parse->swath, parse->tile);
			sprintf(buf_png[1], "%s_%d_%d_%d_rev.png", prefix, parse->surface, parse->swath, parse->tile);

			// Create images
			png[0] = fopen(buf_png[0], "wb");
			png[1] = fopen(buf_png[1], "wb");
			
			current_png[0] = png_create_write_struct
			(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
			current_png[1] = png_create_write_struct
			(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
			
			if (!current_png[0] || !current_png[1])
				return false;
			
			current_png_info[0] = png_create_info_struct(current_png[0]);
			current_png_info[1] = png_create_info_struct(current_png[1]);
			if (!current_png_info[0] || !current_png_info[1])
			{
				png_destroy_write_struct(&current_png[0],
										 (png_infopp)NULL);
				png_destroy_write_struct(&current_png[1],
										 (png_infopp)NULL);
				return false;
			}
			
			png_init_io(current_png[0], png[0]);
			png_init_io(current_png[1], png[1]);
			png_set_IHDR(current_png[0], current_png_info[0], X_LEN, Y_LEN, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
			png_set_IHDR(current_png[1], current_png_info[1], X_LEN, Y_LEN, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
			png_set_PLTE(current_png[0], current_png_info[0], QUAL_PALETTE, QUAL_PALETTE_LEN);
			png_set_PLTE(current_png[1], current_png_info[1], QUAL_PALETTE, QUAL_PALETTE_LEN);
			png_write_info(current_png[0], current_png_info[0]);
			png_write_info(current_png[1], current_png_info[1]);
		}
		if ((b->core.flag&BAM_FUNMAP) == 0) {
			if ((b->core.flag&BAM_FPROPER_PAIR) != 0 && (b->core.isize > 0)) {
				++count;
				bam_glomp(b, tile_grid_qual[parse->surface][parse->swath][parse->tile]);
				
				++surface[parse->surface];
				++swath[parse->surface][parse->swath];
				++tile_grid[parse->surface][parse->swath][parse->tile];

				welford_add(&tile_grid_isize[parse->surface][parse->swath][parse->tile], (double)b->core.isize);
				fprintf(output[parse->surface][parse->swath][parse->tile], "%d\t%d\t%d\t%d\n", bam_get_qual(b)[99], parse->x, parse->y, read);
#ifdef FULL_PALETTE
				current_bitmap[read][parse->y/10][parse->x/10] = bam_get_qual(b)[99];
#else
				if (bam_get_qual(b)[99] > 30) {
					current_bitmap[read][parse->y/10][parse->x/10] = 4;
				} else if (bam_get_qual(b)[99] > 20) {
					current_bitmap[read][parse->y/10][parse->x/10] = 3;
				} else if (bam_get_qual(b)[99] > 10) {
					current_bitmap[read][parse->y/10][parse->x/10] = 2;
				} else {
					current_bitmap[read][parse->y/10][parse->x/10] = 1;
				}
#endif
			}

			tile_grid_mq[parse->surface][parse->swath][parse->tile] += b->core.qual;
		}
		
		parsed_readname_destroy(parse);
	} while (sam_read1(in, hdr, b) >= 0);
	png_write_image(current_png[0], current_bitmap[0]);
	png_write_image(current_png[1], current_bitmap[1]);
	png_write_end(current_png[0], current_png_info[0]);
	png_write_end(current_png[1], current_png_info[1]);

	int i;
	for (i = 0; i < Y_LEN; ++i) {
		png_free(current_png[0], current_bitmap[0][i]);
		png_free(current_png[1], current_bitmap[1][i]);
	}
	png_free(current_png[0], current_bitmap[0]);
	png_free(current_png[1], current_bitmap[1]);

	png_destroy_write_struct(&current_png[0], &current_png_info[0]);
	png_destroy_write_struct(&current_png[1], &current_png_info[1]);
	fclose(png[0]);
	fclose(png[1]);

	
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
	return true;
}

static void usage()
{
	fprintf(stderr,"Usage information\n");
}

static bool open_files(const char* prefix, FILE* out[2][3][16])
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			for (tile = 0; tile < 16; ++tile) {
				char buf[255];
				sprintf(buf, "%s_%d_%d_%d.tsv",prefix, surface, swath, tile);
				out[surface][swath][tile] = fopen(buf, "w");
			}
		}
	}
	return true;
}

static void close_files(FILE* out[2][3][16])
{
	int surface, swath, tile;
	for (surface = 0; surface < 2; ++surface) {
		for (swath = 0; swath < 3; ++swath) {
			for (tile = 0; tile < 16; ++tile) {
				fclose(out[surface][swath][tile]);
			}
		}
	}
}


int bam_qualview(int argc, char *argv[])
{
	int c;
	samFile* in;
	FILE* out[2][3][16];
	while ((c = getopt(argc, argv, "")) >= 0) {
		switch (c) {
		}
	}
	if (optind+1 >= argc) {
		usage();
		return -1;
	}
	in = sam_open(argv[optind], "rb", NULL);
	open_files(argv[optind+1],out);
	bam_qualview_core(in, out, argv[optind+1]);
	sam_close(in);
	close_files(out);
	return 0;
}

