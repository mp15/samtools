/*  bam_rmrg.c -- remove readgroup subcommand.

    Copyright (C) 2013, 2014, 2016 Genome Research Ltd.

    Author: Martin Pollard <mp15@sanger.ac.uk>

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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include "sam_opts.h"


KHASH_MAP_INIT_STR(c2i, int)

struct parsed_opts {
    sam_global_args ga;
    bool verbose;
    char* input_name; // from argv, don't free
    char* output_name;
    kh_c2i_t* rg_hash;
    kh_c2i_t* pu_hash;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
    kh_c2i_t* rg_hash;
};

typedef struct state state_t;

static int cleanup_state(state_t* status);
static void cleanup_opts(parsed_opts_t* opts);

static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools rmrg [-r <RGID>] [-R <rgids.fofn>] [-p <PUID>] [-P <puids.fofn>] [-o <output.bam>] <input.bam>\n"
"Options:\n"
"  -o FILE         output filename\n"
"  -r STRING       RG ID of readgroup to remove\n"
"  -R FILE1        file of RG ID's to be removed from input file\n"
"  -p STRING       RG PU of readgroup to remove\n"
"  -P FILE2        file of RG ID's to be removed from input file\n"
            );
    sam_global_opt_help(write_to, "-....");
}

int read_list(kh_c2i_t * hash, const char *restrict filename)
{
    FILE* list_file = fopen(filename, "r");
    if (list_file == NULL) return 0;
    char* buf = NULL;
    size_t linecap = 0;
    ssize_t linelen = 0;
    while( (linelen = getline(&buf, &linecap, list_file)) > 0) {
        if (buf[linelen-1] == '\n') buf[linelen-1]= '\0';
        int ret = 0;
        kh_put(c2i, hash, strdup(buf), &ret);
    }
    fclose(list_file);

    return 1;
}

// Takes the command line options and turns them into something we can understand
static parsed_opts_t* parse_args(int argc, char** argv)
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "vr:R:p:P:";
    char* delim;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    parsed_opts_t* retval = calloc(sizeof(parsed_opts_t), 1);
    if (! retval ) { perror("cannot allocate option parsing memory"); return NULL; }
    retval->rg_hash = kh_init_c2i();
    retval->pu_hash = kh_init_c2i();

    sam_global_args_init(&retval->ga);

    int opt;
    while ((opt = getopt_long(argc, argv, optstring, lopts, NULL)) != -1) {
        switch (opt) {
            case 'v':
                retval->verbose = true;
                break;
            case 'o':
                retval->output_name = optarg;
                break;
            case 'r':
            {
                int ret;
                kh_put(c2i, retval->rg_hash, optarg, &ret);
                break;
            }
            case 'R':
            {
                if ( read_list(retval->rg_hash, optarg) == 0) {
                    free(retval);
                    return NULL;
                }
                break;
            }
            case 'p':
            {
                int ret;
                kh_put(c2i, retval->pu_hash, optarg, &ret);
                break;
            }
            case 'P':
            {
                if ( read_list(retval->pu_hash, optarg) == 0) {
                    free(retval);
                    return NULL;
                }
                break;
            }
            default:
                if (parse_sam_global_opt(opt, optarg, lopts, &retval->ga) == 0) break;
                /* else fall-through */
            case '?':
                usage(stdout);
                free(retval);
                return NULL;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1) {
        fprintf(stderr, "Invalid number of arguments: %d\n", argc);
        usage(stderr);
        free(retval);
        return NULL;
    }

    retval->input_name = argv[0];

    return retval;
}


// Filters a header of @RG lines where ID != id_keep
// TODO: strip @PG's descended from other RGs and their descendants
static bool filter_header_rg(bam_hdr_t* hdr, const char* id_keep)
{
    kstring_t str = {0, 0, NULL};

    regex_t rg_finder;

    if (regcomp(&rg_finder, "^@RG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE) != 0) {
        return false;
    }

    // regex vars
    char* header = hdr->text;
    regmatch_t* matches = (regmatch_t*)calloc(sizeof(regmatch_t),2);
    kstring_t found_id = { 0, 0, NULL };
    int error;

    while ((error = regexec(&rg_finder, header, 2, matches, 0)) == 0) {
        kputsn(header, matches[0].rm_so, &str); // copy header up until the found RG line

        found_id.l = 0;
        kputsn(header+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &found_id); // extract ID
        // if it matches keep keep it, else we can just ignore it
        if (strcmp(ks_str(&found_id), id_keep) == 0) {
            kputsn(header+matches[0].rm_so, (matches[0].rm_eo+1)-matches[0].rm_so, &str);
        }
        // move pointer forward
        header += matches[0].rm_eo+1;
    }
    // cleanup
    free(found_id.s);
    free(matches);
    regfree(&rg_finder);
    // Did we leave loop because of an error?
    if (error != REG_NOMATCH) {
        return false;
    }

    // Write remainder of string
    kputs(header, &str);

    // Modify header
    hdr->l_text = ks_len(&str);
    free(hdr->text);
    hdr->text = ks_release(&str);

    return true;
}

// Set the initial state
static state_t* init(parsed_opts_t* opts)
{
    state_t* retval = calloc(sizeof(state_t), 1);
    if (!retval) {
        fprintf(stderr, "Out of memory");
        return NULL;
    }

    retval->input_file = sam_open_format(opts->input_name, "rb", &opts->ga.in);
    if (!retval->input_file) {
        fprintf(stderr, "Could not open input file (%s)\n", opts->input_name);
        free(retval);
        return NULL;
    }
    retval->input_header = sam_hdr_read(retval->input_file);
    if (retval->input_header == NULL) {
        fprintf(stderr, "Could not read header for file '%s'\n",
                opts->input_name);
        cleanup_state(retval);
        return NULL;
    }

    // Open output file
    retval->output_file = sam_open_format(opts->output_name?opts->output_name:"-", "wb", &opts->ga.out);
    retval->output_header = bam_hdr_dup(retval->input_header);
    if (!retval->output_file || !retval->output_header) {
        fprintf(stderr, "Could not open output file (%s)\n", opts->output_name);
        cleanup_state(retval);
        return NULL;
    }
    
    retval->rg_hash = opts->rg_hash;
    opts->rg_hash = NULL;
    // TODO: Look up all PUs and transfer them to the RG kill hash

    return retval;
}

static bool rm_rg(state_t* state)
{
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        fprintf(stderr, "Could not write output file header\n");
        return false;
    }

    bam1_t* file_read = bam_init1();
    // Read the first record
    int r;
    if ((r=sam_read1(state->input_file, state->input_header, file_read)) < 0) {
        // Nothing more to read?  Ignore this file
        bam_destroy1(file_read);
        file_read = NULL;
        if (r < -1) {
            fprintf(stderr, "Could not write read sequence\n");
            return false;
        }
    }

    while (file_read != NULL) {
        // Get RG tag from read and look it up in hash to find if we write it out
        uint8_t* tag = bam_aux_get(file_read, "RG");
        khiter_t iter;
        if ( tag != NULL ) {
            char* rg = bam_aux2Z(tag);
            iter = kh_get_c2i(state->rg_hash, rg);
        } else {
            iter = kh_end(state->rg_hash);
        }

        // Write the read out if not in hash
        if (iter == kh_end(state->rg_hash)) {
            if (sam_write1(state->output_file, state->output_header, file_read) < 0) {
                fprintf(stderr, "Could not write sequence\n");
                return false;
            }
        }
        // Replace written read with the next one to process
        if ((r=sam_read1(state->input_file, state->input_header, file_read)) < 0) {
            // Nothing more to read?  Ignore this file in future
            bam_destroy1(file_read);
            file_read = NULL;
            if (r < -1) {
                fprintf(stderr, "Could not write read sequence\n");
                return false;
            }
        }
    }

    return true;
}

static int cleanup_state(state_t* status)
{
    int ret = 0;

    if (!status) return 0;
    sam_close(status->input_file);
    bam_hdr_destroy(status->output_header);
    ret |= sam_close(status->output_file);
    bam_hdr_destroy(status->input_header);
    kh_destroy_c2i(status->rg_hash);
    free(status);

    return ret;
}

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    sam_global_args_free(&opts->ga);
    free(opts);
}

int main_rmrg(int argc, char** argv)
{
    int ret = 1;
    parsed_opts_t* opts = parse_args(argc, argv);
    if (!opts ) goto cleanup_opts;
    state_t* status = init(opts);
    if (!status) goto cleanup_opts;

    if (rm_rg(status)) ret = 0;

    ret |= (cleanup_state(status) != 0);
cleanup_opts:
    cleanup_opts(opts);

    return ret;
}
