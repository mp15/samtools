/* bam_fastxam.c -- samtools command to convert FASTA/FASTQ to SAM/BAM/CRAM.

   Copyright (c) 2016 Genome Research Limited.

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
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>

#include "sam_opts.h"
#include "samtools.h"

KSEQ_INIT(gzFile, gzread)

struct parsed_opts {
    char* input_name_1;
    char* input_name_2;
    char* output_name;
    sam_global_args ga;
};

struct state;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

struct state {
    gzFile input_file_1;
    gzFile input_file_2;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
    void (*mode_func)(const state_t*, bam1_t*);
};

static void cleanup_opts(parsed_opts_t* opts)
{
}

static void cleanup_state(state_t* state)
{
}

static void usage(FILE *fp, const char* command)
{
    fprintf(fp,
            "Usage: samtools %s [options] [-D] [-o <output.bam>] <input_fwd.fq> <input_rev.fq>\n"
            "\n"
            "Options:\n"
            "  -D STRING dummy option\n",
            command
            );
    sam_global_opt_help(fp, "..O..");
}

static bool parse_args(int argc, char** argv, parsed_opts_t** opts)
{
    *opts = NULL;
    int n;

    if (argc == 1) { usage(stdout, argv[0]); return true; }

    parsed_opts_t* retval = calloc(1, sizeof(parsed_opts_t));
    if (! retval ) {
        fprintf(stderr, "[%s] Out of memory allocating parsed_opts_t\n", __func__);
        return false;
    }
    // Set defaults
    sam_global_args_init(&retval->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, 0, 'O', 0, 0),
        { NULL, 0, NULL, 0 }
    };
    if (strcmp(argv[0],"cram") == 0) {
        retval->ga.out.format = cram;
    } else if (strcmp(argv[0],"bam") == 0) {
        retval->ga.out.format = bam;
    } else if (strcmp(argv[0],"sam") == 0) {
        retval->ga.out.format = sam;
    } else {
        fprintf(stderr, "[%s] Unknown conversion type (We should never get here).\n", __func__);
        free(retval);
        return true;
    }

    while ((n = getopt_long(argc, argv, "D:", lopts, NULL)) >= 0) {
        switch (n) {
            case 'h':
                usage(stdout, argv[0]);
                free(retval);
                return true;
            case '?':
                usage(stderr, argv[0]);
                free(retval);
                return false;
            case 'O':
            default:
                if (parse_sam_global_opt(n, optarg, lopts, &retval->ga) == 0) break;
                usage(stderr, argv[0]);
                free(retval);
                return false;
        }
    }
    if (argc-optind < 2) {
        fprintf(stderr, "You must specify an input file.\n");
        usage(stderr, argv[0]);
        cleanup_opts(retval);
        return false;
    }
    retval->input_name_1 = strdup(argv[optind+0]);
    retval->input_name_2 = strdup(argv[optind+1]);

    *opts = retval;
    return true;

}

static bool init(const parsed_opts_t* opts, state_t** state_out)
{
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));
    if (retval == NULL) {
        fprintf(stderr, "[init] Out of memory allocating state struct.\n");
        return false;
    }
    *state_out = retval;
    return true;
}

static bool fastxam(state_t* state)
{
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("fastxam", "[%s] Could not write header to output file", __func__);
        return false;
    }
    kseq_t *seq_1 = kseq_init(state->input_file_1);
    kseq_t *seq_2 = kseq_init(state->input_file_2);

    while (feof(state->input_file_1) && feof(state->input_file_2)) {
        ssize_t l = 0;
    
        if ((l = kseq_read(seq_1)) >= 0) {
            bam1_t b;
            b.core.tid = -1;
            b.core.pos = -1;
            b.core.bin = 0;
            b.core.qual = 0;
            b.core.l_qname = ks_len(&seq_1->name);
            b.core.flag = BAM_FPAIRED|BAM_FUNMAP|BAM_FMUNMAP|BAM_FREAD1;
            b.core.n_cigar = 0;
            b.core.l_qseq = ks_len(&seq_1->seq);//read_one_seq_l;
            b.core.mtid = -1;
            b.core.mpos = -1;
            b.data = calloc(1,(b.core.l_qname + ((b.core.l_qseq + 1)>>1) + b.core.l_qseq));
            int i;
            for (i = 0; i < b.core.l_qseq; ++i) {
                bam_seqi(bam_get_seq(&b),i) = seq_nt16_table[(int)ks_str(&seq_1->seq)[i]];
                bam_get_qual(&b)[i] = ks_str(&seq_1->qual)[i] - 33;
            }
            if(sam_write1(state->output_file, state->output_header, &b) < 0) return false;
        } else {
            return false;
        }

        if ((l = kseq_read(seq_2)) >= 0) {
        } else {
            return true;
        }
    }

    return true;
}

int main_fastxam(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;

    if (!parse_args(argc, argv, &opts)) goto error;
    if (opts == NULL) return EXIT_SUCCESS; // Not an error but user doesn't want us to proceed
    if (!opts || !init(opts, &state)) goto error;

    if (!fastxam(state)) goto error;

    cleanup_opts(opts);
    cleanup_state(state);

    return EXIT_SUCCESS;
error:
    cleanup_opts(opts);
    cleanup_state(state);

    return EXIT_FAILURE;
}

