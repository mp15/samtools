/* bam_hicdigest.c -- samtools command to filter hic reads.

   Copyright (c) 2013, 2015-2017, 2019-2021, 2025 Genome Research Limited.

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

#include <htslib/sam.h>
#include "kvec.h"
#include "samtools.h"
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

typedef enum {
    sam_output,
    tsv_output,
} output_mode_enum;

/// @brief A pair of BAM records (read1 and read2)
typedef struct bam1_t_pair {
    bam1_t* read1;
    bam1_t* read2;
} bam1_t_pair_t;

// Forward declarations
struct state;
struct parsed_opts;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

static bool sam_output_pair(const bam1_t_pair_t pair, state_t *state);
static bool tsv_output_pair(const bam1_t_pair_t pair, state_t *state);

// Utility functions for bam1_t
#define bam_is_mapped(b) (((b)->core.flag&BAM_FUNMAP) == 0)
#define bam_get_mapq(b) ((b)->core.qual)
#define bam_get_flag(b) ((b)->core.flag)
#define bam_get_tid(b) ((b)->core.tid)
#define bam_get_pos(b) ((b)->core.pos)
#define bam_get_chr(b, hdr) sam_hdr_tid2name(hdr, b->core.tid)
#define bam_get_ncigar(b) ((b)->core.n_cigar)

typedef kvec_t(bam1_t*) bam_vec_t;

struct parsed_opts {
    char* input_name;
    char* output_name;
    int no_pg;
    sam_global_args ga;
    htsThreadPool p;
    int uncompressed;
    output_mode_enum output_mode;
};

struct state {
    const parsed_opts_t* opts;
    samFile* input_file;
    sam_hdr_t* input_header;
    union {
        samFile* output_file;
        FILE* output_tsv;
    };
    sam_hdr_t* output_header;
    bool (*output_func)(const bam1_t_pair_t pair, state_t *state);
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

static bool cleanup_state(state_t* state)
{
    if (!state) return false;
    bool clean_exit = true;
    switch(state->opts->output_mode) {
        case sam_output:
            if (state->output_file){
                if (sam_close(state->output_file) != 0) {
                    print_error_errno("hicdigest", "[cleanup_state] Error closing SAM output file");
                    clean_exit = false;
                }
            }
            sam_hdr_destroy(state->output_header);
            break;
        case tsv_output:
            if (fclose(state->output_tsv) != 0) {
                print_error_errno("hicdigest", "[cleanup_state] Error closing TSV output file");
                clean_exit = false;
            }
        default:
            abort(); // This should never happen
    }
    if (state->input_file) {
        if (sam_close(state->input_file) != 0) {
            print_error_errno("hicdigest", "[cleanup_state] Error closing SAM input file");
            clean_exit = false;
        }
    }
    sam_hdr_destroy(state->input_header);
    free(state);
    return clean_exit;
}

static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools hicdigest [options] [-o <output.bam>] <input.bam>\n"
            "\n"
            "Options:\n"
            "  -m MODE   Set the type of output from one of sam, tsv output modes [sam]\n"
            "  -o FILE   Where to write output to [stdout]\n"
            "  -u        Output uncompressed data\n"
            "  --no-PG   Do not add a PG line\n"
            );
    sam_global_opt_help(fp, "..O..@..");
}

static bool parse_args(int argc, char** argv, parsed_opts_t** opts)
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
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };
    retval->output_mode = sam_output;

    while ((n = getopt_long(argc, argv, "m:o:O:h@:u", lopts, NULL)) >= 0) {
        switch (n) {
            case 'o':
                retval->output_name = strdup(optarg);
                break;
            case 'h':
                usage(stdout);
                free(retval);
                return true;
            case 1:
                retval->no_pg = 1;
                break;
            case 'u':
                retval->uncompressed = 1;
                break;
            case 'm': {
                if (strcmp(optarg, "sam") == 0) {
                    retval->output_mode = sam_output;
                } else if (strcmp(optarg, "tsv") == 0) {
                    retval->output_mode = tsv_output;
                } else {
                    usage(stderr);
                    return false;
                }
                break;
            }
            case '?':
                usage(stderr);
                free(retval);
                return false;
            case 'O':
            default:
                if (parse_sam_global_opt(n, optarg, lopts, &retval->ga) == 0) break;
                usage(stderr);
                free(retval);
                return false;
        }
    }

    if (argc-optind < 1) {
        fprintf(stderr, "You must specify an input file.\n");
        usage(stderr);
        cleanup_opts(retval);
        return false;
    }

    retval->input_name = strdup(argv[optind+0]);

    if (retval->ga.nthreads > 0) {
        if (!(retval->p.pool = hts_tpool_init(retval->ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            return false;
        }
    }

    *opts = retval;
    return true;
}

static bool init(const parsed_opts_t* opts, state_t** state_out) {
    char output_mode[9] = "w";
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));

    if (retval == NULL) {
        fprintf(stderr, "[init] Out of memory allocating state struct.\n");
        return false;
    }
    *state_out = retval;

    // Open files
    retval->opts = opts;
    retval->input_file = sam_open_format(opts->input_name, "r", &opts->ga.in);
    if (retval->input_file == NULL) {
        print_error_errno("hicdigest", "could not open \"%s\"", opts->input_name);
        return false;
    }
    retval->input_header = sam_hdr_read(retval->input_file);

    if (opts->p.pool) {
        hts_set_opt(retval->input_file,  HTS_OPT_THREAD_POOL, &opts->p);
    }

    switch (opts->output_mode) {
        case sam_output:
        {
            // Only do sam opens if we are in SAM mode
            retval->output_header = sam_hdr_dup(retval->input_header);

            if (opts->uncompressed) {
                strcat(output_mode, "0");
            }
            if (opts->output_name) { // File format auto-detection
                sam_open_mode(output_mode + strlen(output_mode),
                            opts->output_name, NULL);
            }
            retval->output_file = sam_open_format(opts->output_name == NULL?"-":opts->output_name, output_mode, &opts->ga.out);

            if (retval->output_file == NULL) {
                print_error_errno("hicdigest", "could not create \"%s\"", opts->output_name);
                return false;
            }
            if (opts->p.pool) {hts_set_opt(retval->output_file, HTS_OPT_THREAD_POOL, &opts->p);}

            retval->output_func = &sam_output_pair;
            break;
        }
        case tsv_output:
            retval->output_tsv = opts->output_name ? fopen(opts->output_name, "w") : stdout;
            retval->output_func = &tsv_output_pair;
            break;
    }

    return true;
}

#define bam_get_rlen(rec) bam_cigar2rlen(rec->core.n_cigar, bam_get_cigar(rec))

int sort_bam_mq_rlen(const void *a, const void *b) {
    const bam1_t *bam_a = *(const bam1_t **)a;
    const bam1_t *bam_b = *(const bam1_t **)b;
    if (bam_get_mapq(bam_b)  != bam_get_mapq(bam_a) ) return bam_get_mapq(bam_b) - bam_get_mapq(bam_a);
    if (bam_get_rlen(bam_a) != bam_get_rlen(bam_b)) return (bam_get_rlen(bam_b) > bam_get_rlen(bam_a)) ? 1 : -1;

    return 0;
}

static bool bam_is_r1(const bam1_t *b) {
    return bam_get_flag(b) & BAM_FREAD1;
}
static bool bam_is_r2(const bam1_t *b) {
    return bam_get_flag(b) & BAM_FREAD2;
}

#define SAME_LOCUS_MAX 1000     // 1 kb
#define NEAR3_MAX     20000     // 20 kb

/// @brief get read start position
/// @param b read to get start pos for
/// @return 
static inline hts_pos_t bam_get_rstart_pos(const bam1_t *b) {
    if (!b) return -1;
    return (!bam_is_rev(b)) ? b->core.pos : (b->core.pos + (bam_get_rlen(b) > 0 ? bam_get_rlen(b) - 1 : 0));
}

// This relies on the fact that the reads grouped by name coming out of BWA go PrimaryR1, SecondaryR1, PrimaryR2, SecondaryR2 
static bam1_t_pair_t is_valid_hic_cluster(const bam_vec_t* accum) {
    // return true if the reads pass the filter
    if (kv_size(*accum) < 2 || kv_size(*accum) > 4) return (bam1_t_pair_t){0};

    // if 2 reads just take them
    if (kv_size(*accum) == 2) {
        return (bam1_t_pair_t){ .read1 = kv_A(*accum, 0), .read2 = kv_A(*accum, 1) };
    }

    // select R1
    bam_vec_t r1_reads;
	kv_init(r1_reads);
    kv_select(*accum, r1_reads, bam_is_r1);
    // select R2
    bam_vec_t r2_reads;
	kv_init(r2_reads);
    kv_select(*accum, r2_reads, bam_is_r2);

    if (kv_size(*accum) == 3) {
        bam1_t *r1 = kv_A(*accum, 1);
        bam1_t *r2 = kv_A(*accum, 2);
        bam1_t *r3 = kv_A(*accum, 3);

        // Check the distance they are apart
        hts_pos_t dist_12 = abs(bam_get_tid(r1)-bam_get_tid(r2))*10000000 + abs(bam_get_rstart_pos(r1)-bam_get_rstart_pos(r2));
        hts_pos_t dist_23 = abs(bam_get_tid(r2)-bam_get_tid(r3))*10000000 + abs(bam_get_rstart_pos(r2)-bam_get_rstart_pos(r3));
        hts_pos_t dist_13 = abs(bam_get_tid(r1)-bam_get_tid(r3))*10000000 + abs(bam_get_rstart_pos(r1)-bam_get_rstart_pos(r3));

        // if none of them are < 1000 drop the read
        if (dist_12 > SAME_LOCUS_MAX && dist_23 > SAME_LOCUS_MAX && dist_13 > SAME_LOCUS_MAX) {
            return (bam1_t_pair_t){0};
        }
        if ((r1->core.flag & (BAM_FREAD1|BAM_FREAD2)) ==  (r2->core.flag & (BAM_FREAD1|BAM_FREAD2))) {
            return (bam1_t_pair_t){.read1 = dist_13 > dist_23 ? r1: r2, .read2 = r3};
        } else if ((r1->core.flag & (BAM_FREAD1|BAM_FREAD2)) ==  (r3->core.flag & (BAM_FREAD1|BAM_FREAD2))) {
            return (bam1_t_pair_t){.read1 = dist_12 > dist_23 ? r1: r3, .read2 = r2};
        } else if ((r2->core.flag & (BAM_FREAD1|BAM_FREAD2)) ==  (r3->core.flag & (BAM_FREAD1|BAM_FREAD2))) {
            return (bam1_t_pair_t){.read1 = dist_12 > dist_13 ? r2: r3, .read2 = r1};
        } else {
            fprintf(stderr, "[%s] Logic error with 3 reads\n", __func__);
            abort();
        }
    } else if (kv_size(*accum) == 4) {
        bam1_t *r1 = kv_A(r1_reads, 1);
        bam1_t *r2 = kv_A(r1_reads, 2);
        bam1_t *r3 = kv_A(r2_reads, 1);
        bam1_t *r4 = kv_A(r2_reads, 2);

        hts_pos_t dist_13 = abs(bam_get_tid(r1)-bam_get_tid(r3))*10000000 + abs(bam_get_rstart_pos(r1)-bam_get_rstart_pos(r3));
	    hts_pos_t dist_24 = abs(bam_get_tid(r2)-bam_get_tid(r4))*10000000 + abs(bam_get_rstart_pos(r2)-bam_get_rstart_pos(r4));
	    hts_pos_t dist_14 = abs(bam_get_tid(r1)-bam_get_tid(r4))*10000000 + abs(bam_get_rstart_pos(r1)-bam_get_rstart_pos(r4));
	    hts_pos_t dist_23 = abs(bam_get_tid(r2)-bam_get_tid(r3))*10000000 + abs(bam_get_rstart_pos(r2)-bam_get_rstart_pos(r3));

	    if ((dist_13 < SAME_LOCUS_MAX && dist_24 < SAME_LOCUS_MAX) || (dist_14 < SAME_LOCUS_MAX && dist_23 < SAME_LOCUS_MAX)) {
            return (bam1_t_pair_t){.read1 = r1, .read2 = r3};
        } else {
            return (bam1_t_pair_t){0};
        }
    }
    return (bam1_t_pair_t){0};
}

static bool sam_output_pair(const bam1_t_pair_t pair, state_t *state)
{
    if (sam_write1(state->output_file, state->output_header, pair.read1) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write read1 to output file", __func__);
        return false;
    }
    if (sam_write1(state->output_file, state->output_header, pair.read2) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write read2 to output file", __func__);
        return false;
    }
    return true;
}

static int fprint_cigar(FILE* fp, const uint32_t* cigar, const int n_cigar)
{
    int retval = 0;
    for (int i = 0; i < n_cigar; i++) {
        int chars;
        if ((chars = fprintf(fp, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]))) < 0) {
            return chars;
        }
        retval += chars;
    }
    return retval;
}

static bool tsv_output_bam1(const bam1_t* b1, const bam1_t* b2, state_t *state)
{
    //fprintf(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1]"$"fname"$"(prevlinenum+1),name[read2]"$"fname"$"linenum)
    if (fprintf(state->output_tsv, "%i\t%s\t%ld\t%i\t%s\t%ld\t%d\t",
            bam_is_rev(b1),
            bam_get_chr(b1, state->input_header),
            bam_get_pos(b1),
            bam_is_rev(b2),
            bam_get_chr(b2, state->input_header),
            bam_get_pos(b2),
            bam_get_mapq(b1)) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write to output file", __func__);
        return false;
    }

    if (fprint_cigar(state->output_tsv, bam_get_cigar(b1), bam_get_ncigar(b1)) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write to output file", __func__);
        return false;
    }

    if (fprintf(state->output_tsv, "\t%s\t%d\t",
            bam_get_seq(b1),
            bam_get_mapq(b2)) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write to output file", __func__);
        return false;
    }

    if (fprint_cigar(state->output_tsv, bam_get_cigar(b2), bam_get_ncigar(b2)) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write to output file", __func__);
        return false;
    }

    if (fprintf(state->output_tsv, "\t%s\t%s\t%s\n",
            bam_get_seq(b2),
            bam_get_qname(b1),
            bam_get_qname(b2)) < 0) {
        print_error_errno("hicdigest", "[%s] Could not write to output file", __func__);
        return false;
    }
    return true;
}

// Compares two bam records and see which comes first
static bool bam_less_than(const bam1_t* read1, const bam1_t* read2) {
    if (bam_get_tid(read1) < bam_get_tid(read2)) { return true; }
        else if (bam_get_tid(read1) > bam_get_tid(read2)) { return false; }
    // chr1 == chr2
    if (bam_is_rev(read1) < bam_is_rev(read2)) { return true; }
        else if (bam_is_rev(read1) > bam_is_rev(read2)) { return false; }
    // s1 == s2 && c1 == c2
    if (bam_get_pos(read1) < bam_get_pos(read2)) { return true; }
        else if (bam_get_pos(read1) > bam_get_pos(read2)) { return false; }
    // all are equal, doesn't matter
    return true;
}

static bool tsv_output_pair(const bam1_t_pair_t pair, state_t *state)
{
    if (bam_less_than(pair.read1, pair.read2)) {
        return tsv_output_bam1(pair.read1, pair.read2, state);
    }
    else {
        return tsv_output_bam1(pair.read2, pair.read1, state);
    }
}

// Process a batch of reads with same readname accumulated
// Note this destroys the list passed in
static bool process_batch(bam_vec_t* accum, state_t *state)
{
    // Do the accum reads pass the filter?
    bam1_t_pair_t pair = is_valid_hic_cluster(accum);
    if (pair.read1 && pair.read2) {
        // if so write them out
        if (!(state->output_func(pair, state))) {
            kv_destroy(*accum);
            return false;
        }
    }
    // get rid of items either way
    for (size_t item = 0; item != kv_size(*accum); item++) {
        bam_destroy1(kv_A(*accum, item));
    }
    kv_destroy(*accum);

    return true;
}

static bool hicdigest(parsed_opts_t *opts, state_t* state, char *arg_list)
{
    // Write PG line to header
    if (!opts->no_pg && sam_hdr_add_pg(state->output_header, "samtools",
                                       "VN", samtools_version(),
                                       arg_list ? "CL": NULL,
                                       arg_list ? arg_list : NULL,
                                       NULL))
        return false;

    // Write header to output
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("hicdigest", "[%s] Could not write header to output file", __func__);
        return false;
    }
    // Write index if requested
    char *idx_fn = NULL;
    if (opts->ga.write_index) {
        if (!(idx_fn = auto_index(state->output_file, opts->output_name, state->output_header))) {
            return false;
        }
    }

    // Requires namesorted reads for this to work
    bam_vec_t accum;
	kv_init(accum);

    bam1_t* file_read = bam_init1();
    // Preload the first read
    int ret = sam_read1(state->input_file, state->input_header, file_read);
    if (ret < 0) {
        print_error_errno("hicdigest", "[%s] Error reading from input file", __func__);
        free(idx_fn);
        return false;
    }
    char *prev_qname = bam_get_qname(file_read);

    if (ret > 0) {
        do {
            // do we have the same read name as last read?
            if (strcmp(bam_get_qname(file_read), prev_qname))
            {
                // process the batch
                if (!process_batch(&accum, state)) {
                    bam_destroy1(file_read);
                    free(idx_fn);
                    return false;
                }
                kv_init(accum);
            }
            // keep only well mapped reads
            if (bam_is_mapped(file_read) && bam_get_mapq(file_read) > 0) {
                // Add the read to accum
                kv_push(bam1_t*, accum, bam_dup1(file_read));
            }
        } while ((ret = sam_read1(state->input_file, state->input_header, file_read)) >= 0);
        // process last batch
        process_batch(&accum, state);
    }

    bam_destroy1(file_read);
    if (ret != -1) {
        print_error_errno("hicdigest", "[%s] Error reading from input file", __func__);
        free(idx_fn);
        return false;
    } else {
        if (opts->ga.write_index) {
            if (sam_idx_save(state->output_file) < 0) {
                print_error_errno("hicdigest", "[%s] Writing index failed", __func__);
                free(idx_fn);
                return false;
            }
        }
        free(idx_fn);
        return true;
    }
}

int main_hicdigest(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;
    char *arg_list = stringify_argv(argc+1, argv-1);
    if (!arg_list)
        return EXIT_FAILURE;

    if (!parse_args(argc, argv, &opts)) goto error;
    if (opts) { // Not an error but user doesn't want us to proceed
        if (!init(opts, &state) || !hicdigest(opts, state, arg_list))
            goto error;
    }

    int final_state = EXIT_SUCCESS;
    if (!cleanup_state(state))
        final_state = EXIT_FAILURE;
    cleanup_opts(opts);
    free(arg_list);

    return final_state;
error:
    cleanup_state(state);
    cleanup_opts(opts);
    free(arg_list);

    return EXIT_FAILURE;
}
