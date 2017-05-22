/* bam_extract.c -- samtools command to extract specified regions for assembly.
 
 Copyright (c) 2013, 2015, 2016 Genome Research Limited.
 
 Author: Martin O. Pollard <mp15@sanger.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include "sam_opts.h"
#include "samtools.h"

struct parsed_opts {
    char* input_name;
    char* filter_region;
    char* output_name;
    sam_global_args ga;
    htsThreadPool p;
};

struct state;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

struct state {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
    uint16_t flag_off;
    int beg0, end0, tid0;
};

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    free(opts->output_name);
    free(opts->input_name);
    if (opts->p.pool) hts_tpool_destroy(opts->p.pool);
    sam_global_args_free(&opts->ga);
    free(opts);
}

static void cleanup_state(state_t* state)
{
    if (!state) return;
    if (state->output_file) sam_close(state->output_file);
    bam_hdr_destroy(state->output_header);
    if (state->input_file) sam_close(state->input_file);
    bam_hdr_destroy(state->input_header);
    free(state);
}

static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools extract [options] [-o <output.bam>] <input.bam> <region>\n"
            "\n"
            "Options:\n"
            "  -o FILE   Where to write output to [stdout]\n"
            );
    sam_global_opt_help(fp, "..O..@");
}

static int parse_args(int argc, char** argv, parsed_opts_t** opts)
{
    *opts = NULL;
    int n;

    if (argc == 1) { usage(stdout); return true; }

    parsed_opts_t* retval = calloc(1, sizeof(parsed_opts_t));
    if (! retval ) {
        fprintf(stderr, "[%s] Out of memory allocating parsed_opts_t\n", __func__);
        return false;
    }
    // Set defaults
    sam_global_args_init(&retval->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, 0, 'O', 0, 0, '@'),
        { NULL, 0, NULL, 0 }
    };

    while ((n = getopt_long(argc, argv, "o:O:h@:", lopts, NULL)) >= 0) {
        switch (n) {
//            case '':
//                break;
            case 'o':
                retval->output_name = strdup(optarg);
                break;
            case 'h':
                usage(stdout);
                free(retval);
                return 0;
            case '?':
                usage(stderr);
                free(retval);
                return 1;
            default:
                if (parse_sam_global_opt(n, optarg, lopts, &retval->ga) == 0) break;
                usage(stderr);
                free(retval);
                return 1;
        }
    }
    if (argc-optind < 2) {
        fprintf(stderr, "You must specify an input file and region.\n");
        usage(stderr);
        cleanup_opts(retval);
        return 1;
    }
    retval->input_name = strdup(argv[optind+0]);
    retval->filter_region = strdup(argv[optind+1]);

    *opts = retval;
    return 0;
}

static int init(const parsed_opts_t* opts, state_t** state_out) {
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));
    if (retval == NULL) {
        fprintf(stderr, "[init] Out of memory allocating state struct.\n");
        return false;
    }
    *state_out = retval;

    // Open files
    retval->input_file = sam_open_format(opts->input_name, "r", &opts->ga.in);
    if (retval->input_file == NULL) {
        print_error_errno("addreplacerg", "could not open \"%s\"", opts->input_name);
        return 1;
    }
    retval->input_header = sam_hdr_read(retval->input_file);
    hts_idx_t *idx = sam_index_load(retval->input_file, opts->input_name);
    if (idx == NULL) {
        fprintf(stderr, "[%s] fail to load index for %s\n", __func__, opts->input_name);
        exit(EXIT_FAILURE);
    }
    hts_itr_t *iter;
    if ( (iter=sam_itr_querys(idx, retval->input_header, opts->filter_region)) == 0) {
        fprintf(stderr, "[E::%s] fail to parse region '%s' with %s\n", __func__, opts->filter_region, opts->input_name);
        exit(EXIT_FAILURE);
    }
    retval->beg0 = iter->beg, retval->end0 = iter->end, retval->tid0 = iter->tid;
    hts_idx_destroy(idx);

    retval->output_header = bam_hdr_dup(retval->input_header);
    retval->output_file = sam_open_format(opts->output_name == NULL?"-":opts->output_name, "w", &opts->ga.out);
    retval->flag_off = 2;

    if (retval->output_file == NULL) {
        print_error_errno("extract", "could not create \"%s\"", opts->output_name);
        return 1;
    }

    if (opts->p.pool) {
        hts_set_opt(retval->input_file,  HTS_OPT_THREAD_POOL, &opts->p);
        hts_set_opt(retval->output_file, HTS_OPT_THREAD_POOL, &opts->p);
    }
    return 0;
}

static int extract(state_t* state)
{
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("extract", "[%s] Could not write header to output file", __func__);
        return 1;
    }

    bam1_t* file_read = bam_init1();
    int ret;
    while ((ret = sam_read1(state->input_file, state->input_header, file_read)) >= 0) {
        if (file_read->core.flag&(state->flag_off) &&
            !(state->tid0 == file_read->core.tid &&
             state->beg0 <= file_read->core.pos &&
             state->end0 >= bam_endpos(file_read))) continue;
        if (sam_write1(state->output_file, state->output_header, file_read) < 0) {
            print_error_errno("extract", "[%s] Could not write read to output file", __func__);
            bam_destroy1(file_read);
            return 1;
        }
    }
    bam_destroy1(file_read);
    if (ret != -1) {
        print_error_errno("extract", "[%s] Error reading from input file", __func__);
        return 1;
    } else {
        return 0;
    }
}

int main_extract(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;

    if (parse_args(argc, argv, &opts)) goto error;
    if (opts == NULL) return EXIT_SUCCESS; // Not an error but user doesn't want us to proceed
    if (!opts || init(opts, &state)) goto error;

    if (extract(state)) goto error;

    cleanup_state(state);
    cleanup_opts(opts);

    return EXIT_SUCCESS;
error:
    cleanup_state(state);
    cleanup_opts(opts);

    return EXIT_FAILURE;
}
