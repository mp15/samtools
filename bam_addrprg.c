/* bam_addrprg.c -- samtools command to add or replace readgroups.

   Copyright (c) 2013, 2015 Genome Research Limited.

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

#include <htslib/sam.h>
#include <htslib/kstring.h>
#include "samtools.h"
#include "sam_opts.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

typedef enum {
    overwrite_all,
    orphan_only,
} rg_mode;

struct parsed_opts {
    char* input_name;
    char* output_name;
    char* rg_id;
    char* rg_line;
    rg_mode mode;
    sam_global_args ga;
};

struct state;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

struct state {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
    char* rg_id;
    void (*mode_func)(const state_t*, bam1_t*);
};

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    free(opts->rg_id);
    free(opts->output_name);
    free(opts->input_name);
    sam_global_args_free(&opts->ga);
    free(opts);
}

static void cleanup_state(state_t* state)
{
    if (!state) return;
    free(state->rg_id);
    if (state->output_file) sam_close(state->output_file);
    bam_hdr_destroy(state->output_header);
    if (state->input_file) sam_close(state->input_file);
    bam_hdr_destroy(state->input_header);
    free(state);
}

// Converts \t and \n into real tabs and newlines
static char* basic_unescape(const char* in)
{
    assert(in);
    char *ptr, *out;
    out = ptr = malloc(strlen(in)+1);
    size_t size = 0;
    while (*in) {
        if (*in == '\\') {
            ++in;
            if (*in == '\0') {
                fprintf(stderr, "[%s] Unterminated escape sequence.\n", __func__);
                free(out);
                return NULL;
            }
            switch (*in) {
            case '\\':
                *ptr = '\\';
                break;
            case 't':
                *ptr = '\t';
                break;
            case 'n':
                fprintf(stderr, "[%s] \\n in escape sequence is not supported.\n", __func__);
                free(out);
                return NULL;
            default:
                fprintf(stderr, "[%s] Unsupported escape sequence.\n", __func__);
                free(out);
                return NULL;
            }
        } else {
            *ptr = *in;
        }
        ++in;
        ++ptr;
        ++size;
    }
    *ptr = '\0';
    ++size;
    char* tmp = (char*)realloc(out, size);
    if (!tmp) {
        free(out);
    }
    return tmp;
}

// These are to be replaced by samtools header parser
// Extracts the first @RG line from a string.
static char* get_rg_line(const char* text, size_t* last)
{
    const char* rg = text;
    if (rg[0] != '@' || rg[1] != 'R' || rg[2] != 'G' ) {
        if ((rg = (const char*)strstr(text,"\n@RG")) == NULL) {
            return NULL;
        }
        rg++;//skip initial \n
    }
    // duplicate the line for return
    char* line;
    char* end = strchr(rg, '\n');
    if (end) {
        line = strndup(rg,(end-rg));
        *last = end - rg;
    } else {
        line = strdup(rg);
        *last = strlen(rg);
    }
    return line;
}

// Given a @RG line return the id
static char* get_rg_id(const char* input)
{
    assert(input!=NULL);
    char* line = strdup(input);
    char *next = line;
    char* token = strsep(&next, "\t");
    token = strsep(&next,"\t"); // skip first token it should always be "@RG"
    while (next != NULL) {
        char* key = strsep(&token,":");
        if (!strcmp(key,"ID")) {
            char* retval = strdup(token);
            free(line);
            return retval;
        }
        token = strsep(&next,"\t");
    }
    free(line);
    return NULL;
}

// Confirms the existance of an RG line with a given ID in a bam header
static bool confirm_rg( const bam_hdr_t *hdr, const char* rgid )
{
    assert( hdr != NULL && rgid != NULL );

    char *ptr, *start;
    bool found = false;
    start = ptr = strndup(hdr->text, hdr->l_text);
    while (ptr != NULL && *ptr != '\0' && found == false ) {
        size_t end = 0;
        char* line = get_rg_line(ptr, &end);
        if (line == NULL) break; // No more @RG
        char* id;
        if (((id = get_rg_id(line)) != NULL) && !strcmp(id, rgid)) {
            found = true;
        }
        free(id);
        free(line);
        ptr += end;
    }
    free(start);
    return found;
}

static char* get_first_rgid( const bam_hdr_t *hdr )
{
    assert( hdr != NULL );
    char *ptr, *start;
    char* found = NULL;
    start = ptr = strndup(hdr->text, hdr->l_text);
    while (ptr != NULL && *ptr != '\0' && found == NULL ) {
        size_t end = 0;
        char* line = get_rg_line(ptr, &end);
        if ( line ) {
            found = get_rg_id(line);
        } else break;
        free(line);
        ptr += end;
    }
    free(start);
    return found;
}

static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools addreplacerg [options] [-r <@RG line> | -R <existing id>] <input.bam> <output.bam>\n"
            "\n"
            "Options:\n"
            "  -r STRING @RG line text\n"
            "  -R STRING ID of @RG line in existing header to use\n"
            "  -m MODE   Set the mode of operation from one of overwrite_all, orphan_only [overwrite_all]\n");
    sam_global_opt_help(fp, "..O..");
}

static bool parse_args(int argc, char** argv, parsed_opts_t** opts)
{
    *opts = NULL;
    int n;

    parsed_opts_t* retval = calloc(1, sizeof(parsed_opts_t));
    if (! retval ) {
        fprintf(stderr, "[%s] Out of memory allocating parsed_opts_t\n", __func__);
        return true;
    }
    // Set defaults
    retval->mode = overwrite_all;
    sam_global_args_init(&retval->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, 0, 'O', 0, 0),
        { NULL, 0, NULL, 0 }
    };
    kstring_t rg_line = {0,0,NULL};

    while ((n = getopt_long(argc, argv, "r:R:m:o:O:l:h", lopts, NULL)) >= 0) {
        switch (n) {
            case 'r':
                // Are we adding to existing rg line?
                if (ks_len(&rg_line) == 0) {
                    if (strlen(optarg)<3 || (optarg[0] != '@' && optarg[1] != 'R' && optarg[2] != 'G')) {
                        kputs("@RG\t", &rg_line);
                    }
                } else {
                    kputs("\t", &rg_line);
                }
                kputs(optarg, &rg_line);
                break;
            case 'R':
                retval->rg_id = strdup(optarg);
                break;
            case 'm': {
                if (strcmp(optarg, "overwrite_all") == 0) {
                    retval->mode = overwrite_all;
                } else if (strcmp(optarg, "orphan_only") == 0) {
                    retval->mode = orphan_only;
                } else {
                    usage(stderr);
                    return false;
                }
                break;
            }
            case 'o':
                retval->input_name = strdup(optarg);
                break;
            case 'h':
                usage(stdout);
                return false;
            case '?':
                usage(stderr);
                return true;
            case 'O':
            default:
                if (parse_sam_global_opt(n, optarg, lopts, &retval->ga) == 0) break;
                usage(stderr);
                return true;
        }
    }
    retval->rg_line = ks_release(&rg_line);

    if (argc-optind < 1) {
        usage(stdout);
        cleanup_opts(retval);
        return true;
    }
    if (retval->rg_id && retval->rg_line) {
        fprintf(stderr, "The options -r and -R are mutually exclusive.\n");
        cleanup_opts(retval);
        return true;
    }

    if (retval->rg_line)
    {
        char* tmp = basic_unescape(retval->rg_line);

        if ((retval->rg_id = get_rg_id(tmp)) == NULL) {
            fprintf(stderr, "[%s] The supplied RG line lacks an ID tag.\n", __func__);
            free(tmp);
            cleanup_opts(retval);
            return true;
        }
        retval->rg_line = tmp;
    }
    retval->input_name = strdup(argv[optind+0]);

    *opts = retval;
    return false;
}

static void overwrite_all_func(const state_t* state, bam1_t* file_read)
{
    uint8_t* data = (uint8_t*)strdup(state->rg_id);
    int len = strlen(state->rg_id)+1;
    // If the old exists delete it
    uint8_t* old = bam_aux_get(file_read, "RG");
    if (old != NULL) {
        bam_aux_del(file_read, old);
    }

    bam_aux_append(file_read, "RG", 'Z', len, data);
    free(data);
}

static void orphan_only_func(const state_t* state, bam1_t* file_read)
{
    uint8_t* data = (uint8_t*)strdup(state->rg_id);
    int len = strlen(state->rg_id)+1;
    // If the old exists don't do anything
    uint8_t* old = bam_aux_get(file_read, "RG");
    if (old == NULL) {
        bam_aux_append(file_read, "RG",'Z',len,data);
    }
    free(data);
}

static bool init(const parsed_opts_t* opts, state_t** state_out) {
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));
    if (retval == NULL) {
        fprintf(stderr, "[init] Out of memory allocating state struct.\n");
        return true;
    }
    *state_out = retval;

    // Open files
    retval->input_file = sam_open_format(opts->input_name, "r", &opts->ga.in);
    if (retval->input_file == NULL) {
        fprintf(stderr, "[init] Could not open input file: %s\n", opts->input_name);
        return true;
    }
    retval->input_header = sam_hdr_read(retval->input_file);

    retval->output_header = bam_hdr_dup(retval->input_header);
    retval->output_file = sam_open_format(opts->output_name == NULL?"-":opts->output_name, "w", &opts->ga.out);
    
    if (retval->output_file == NULL) {
        print_error_errno("Could not open output file: %s\n", opts->output_name);
        return true;
    }

    if (opts->rg_line) {
        // Append new RG line to header.
        // Check does not already exist
        if ( confirm_rg(retval->output_header, opts->rg_id) ) {
            fprintf(stderr, "[init] ID of new RG line specified conflicts with that of an existing header RG line. Overwrite not yet implemented.\n");
            return true;
        }
        retval->rg_id = strdup(opts->rg_id);
        size_t new_len = strlen( retval->output_header->text ) + strlen( opts->rg_line ) + 2;
        char* new_header = malloc(new_len);
        if (!new_header) {
            fprintf(stderr, "[init] Out of memory whilst writing new header.\n");
            return true;
        }
        sprintf(new_header,"%s%s\n", retval->output_header->text, opts->rg_line);
        free(retval->output_header->text);
        retval->output_header->text = new_header;
        retval->output_header->l_text = (int)new_len;
    } else {
        if (opts->rg_id) {
            // Confirm what has been supplied exists
            if ( !confirm_rg(retval->output_header, opts->rg_id) ) {
                fprintf(stderr, "RG ID supplied does not exist in header. Supply full @RG line with -r instead?\n");
                return true;
            }
            retval->rg_id = strdup(opts->rg_id);
        } else {
            if ((retval->rg_id = get_first_rgid(retval->output_header)) == NULL ) {
                fprintf(stderr, "No RG specified on command line or in existing header.\n");
                return true;
            }
        }
    }

    switch (opts->mode) {
        case overwrite_all:
            retval->mode_func = &overwrite_all_func;
            break;
        case orphan_only:
            retval->mode_func = &orphan_only_func;
            break;
    }

    return false;
}

static bool readgroupise(state_t* state)
{
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("[%s] Could not write header to output file", __func__);
        return true;
    }

    bam1_t* file_read = bam_init1();
    int ret;
    while ((ret = sam_read1(state->input_file, state->input_header, file_read)) >= 0) {
        state->mode_func(state, file_read);

        if (sam_write1(state->output_file, state->output_header, file_read) < 0) {
            print_error_errno("[%s] Could not write read to output file", __func__);
            bam_destroy1(file_read);
            return true;
        }
    }
    bam_destroy1(file_read);
    if (ret != -1) {
        print_error_errno("[%s] Error reading from input file", __func__);
        return true;
    } else {
        return false;
    }
}

int main_addreplacerg(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;

    if (parse_args(argc, argv, &opts)) goto error;
    if (!opts) return 0;
    if (!opts || init(opts, &state)) goto error;

    if (readgroupise(state)) goto error;

    cleanup_opts(opts);
    cleanup_state(state);

    return 0;
error:
    cleanup_opts(opts);
    cleanup_state(state);

    return 1;
}
