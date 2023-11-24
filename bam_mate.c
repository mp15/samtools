/*  bam_mate.c -- fix mate pairing information and clean up flags.

    Copyright (C) 2009, 2011-2017, 2019, 2022 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.
    Portions copyright (C) 2012 Peter Cock, The James Hutton Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "samtools.h"


#define MD_MIN_QUALITY 15

/*
 * This function calculates ct tag for two bams, it assumes they are from the same template and
 * writes the tag to the first read in position terms.
 */
static void bam_template_cigar(bam1_t *b1, bam1_t *b2, kstring_t *str)
{
    bam1_t *swap;
    int i;
    hts_pos_t end;
    uint32_t *cigar;
    str->l = 0;
    if (b1->core.tid != b2->core.tid || b1->core.tid < 0 || b1->core.pos < 0 || b2->core.pos < 0 || b1->core.flag&BAM_FUNMAP || b2->core.flag&BAM_FUNMAP) return; // coordinateless or not on the same chr; skip
    if (b1->core.pos > b2->core.pos) swap = b1, b1 = b2, b2 = swap; // make sure b1 has a smaller coordinate
    kputc((b1->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
    kputc((b1->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
    for (i = 0, cigar = bam_get_cigar(b1); i < b1->core.n_cigar; ++i) {
        kputw(bam_cigar_oplen(cigar[i]), str);
        kputc(bam_cigar_opchr(cigar[i]), str);
    }
    end = bam_endpos(b1);
    kputw(b2->core.pos - end, str);
    kputc('T', str);
    kputc((b2->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
    kputc((b2->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
    for (i = 0, cigar = bam_get_cigar(b2); i < b2->core.n_cigar; ++i) {
        kputw(bam_cigar_oplen(cigar[i]), str);
        kputc(bam_cigar_opchr(cigar[i]), str);
    }

    uint8_t* data;
    if ((data = bam_aux_get(b1,"ct")) != NULL) bam_aux_del(b1, data);
    if ((data = bam_aux_get(b2,"ct")) != NULL) bam_aux_del(b2, data);

    bam_aux_append(b1, "ct", 'Z', str->l+1, (uint8_t*)str->s);
}

/*
 * What This Program is Supposed To Do:
 * Fill in mate coordinates, ISIZE and mate related flags from a name-sorted
 * alignment.
 *
 * How We Handle Input
 *
 * Secondary and supplementary Reads:
 * -write to output unchanged
 * All Reads:
 * -if pos == 0 (1 based), tid == -1 set UNMAPPED flag
 * single Reads:
 * -if pos == 0 (1 based), tid == -1, or UNMAPPED then set UNMAPPED, pos = 0,
 *  tid = -1
 * -clear bad flags (PAIRED, MREVERSE, PROPER_PAIR)
 * -set mpos = 0 (1 based), mtid = -1 and isize = 0
 * -write to output
 * Paired Reads:
 * -if read is unmapped and mate is not, set pos and tid to equal that of mate
 * -sync mate flags (MREVERSE, MUNMAPPED), mpos, mtid
 * -recalculate ISIZE if possible, otherwise set it to 0
 * -optionally clear PROPER_PAIR flag from reads where mapping or orientation
 *  indicate this is not possible (Illumina orientation only)
 * -calculate ct and apply to lowest positioned read
 * -write to output
 * Limitations
 * -Does not handle tandem reads
 * -Should mark supplementary reads the same as primary.
 * Notes
 * -CT definition appears to be something else in spec, this was in here before
 *  I started tampering with it, anyone know what is going on here? To work
 *  around this I have demoted the CT this tool generates to ct.
 */

static void sync_unmapped_pos_inner(bam1_t* src, bam1_t* dest) {
    if ((dest->core.flag & BAM_FUNMAP) && !(src->core.flag & BAM_FUNMAP)) {
        // Set unmapped read's RNAME and POS to those of its mapped mate
        // (recommended best practice, ensures if coord sort will be together)
        dest->core.tid = src->core.tid;
        dest->core.pos = src->core.pos;
    }
}

static void sync_mate_inner(bam1_t* src, bam1_t* dest)
{
    // sync mate pos information
    dest->core.mtid = src->core.tid; dest->core.mpos = src->core.pos;
    // sync flag info
    if (src->core.flag&BAM_FREVERSE)
        dest->core.flag |= BAM_FMREVERSE;
    else
        dest->core.flag &= ~BAM_FMREVERSE;
    if (src->core.flag & BAM_FUNMAP) {
        dest->core.flag |= BAM_FMUNMAP;
    }
}

// Is it plausible that these reads are properly paired?
// Can't really give definitive answer without checking isize
static bool plausibly_properly_paired(bam1_t* a, bam1_t* b)
{
    if ((a->core.flag & BAM_FUNMAP) || (b->core.flag & BAM_FUNMAP)) return false;
    assert(a->core.tid >= 0); // This should never happen if FUNMAP is set correctly

    if (a->core.tid != b->core.tid) return false;

    bam1_t* first = a;
    bam1_t* second = b;
    hts_pos_t a_pos = a->core.flag&BAM_FREVERSE ? bam_endpos(a) : a->core.pos;
    hts_pos_t  b_pos = b->core.flag&BAM_FREVERSE ? bam_endpos(b) : b->core.pos;
    if (a_pos > b_pos) {
        first = b;
        second = a;
    } else {
        first = a;
        second = b;
    }

    if (!(first->core.flag&BAM_FREVERSE) && (second->core.flag&BAM_FREVERSE))
        return true;
    else
        return false;
}

// Returns 0 on success, -1 on failure.
static int bam_format_cigar(const bam1_t* b, kstring_t* str)
{
    // An empty cigar is a special case return "*" rather than ""
    if (b->core.n_cigar == 0) {
        return (kputc('*', str) == EOF) ? -1 : 0;
    }

    const uint32_t *cigar = bam_get_cigar(b);
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; ++i) {
        if (kputw(bam_cigar_oplen(cigar[i]), str) == EOF) return -1;
        if (kputc(bam_cigar_opchr(cigar[i]), str) == EOF) return -1;
    }

    return 0;
}

// Returns 0 on success, -1 on failure.
static int sync_mq_mc(bam1_t* src, bam1_t* dest)
{
    if ( (src->core.flag & BAM_FUNMAP) == 0 ) { // If mapped
        // Copy Mate Mapping Quality
        uint32_t mq = src->core.qual;
        uint8_t* data;
        if ((data = bam_aux_get(dest,"MQ")) != NULL) {
            bam_aux_del(dest, data);
        }

        bam_aux_append(dest, "MQ", 'i', sizeof(uint32_t), (uint8_t*)&mq);
    }
    // Copy mate cigar if either read is mapped
    if ( (src->core.flag & BAM_FUNMAP) == 0 || (dest->core.flag & BAM_FUNMAP) == 0 ) {
        uint8_t* data_mc;
        if ((data_mc = bam_aux_get(dest,"MC")) != NULL) {
            bam_aux_del(dest, data_mc);
        }

        // Convert cigar to string
        kstring_t mc = { 0, 0, NULL };
        if (bam_format_cigar(src, &mc) < 0) return -1;

        bam_aux_append(dest, "MC", 'Z', ks_len(&mc)+1, (uint8_t*)ks_str(&mc));
        free(mc.s);
    }
    return 0;
}

// Copy flags.
// Returns 0 on success, -1 on failure.
static int sync_mate(bam1_t* a, bam1_t* b)
{
    sync_unmapped_pos_inner(a,b);
    sync_unmapped_pos_inner(b,a);
    sync_mate_inner(a,b);
    sync_mate_inner(b,a);
    if (sync_mq_mc(a,b) < 0) return -1;
    if (sync_mq_mc(b,a) < 0) return -1;
    return 0;
}


static uint32_t calc_mate_score(bam1_t *b)
{
    uint32_t score = 0;
    uint8_t  *qual = bam_get_qual(b);
    int i;

    for (i = 0; i < b->core.l_qseq; i++) {
        if (qual[i] >= MD_MIN_QUALITY) score += qual[i];
    }

    return score;
}


static int add_mate_score(bam1_t *src, bam1_t *dest)
{
    uint8_t *data_ms;
    uint32_t mate_score = calc_mate_score(src);

    if ((data_ms = bam_aux_get(dest, "ms")) != NULL) {
        bam_aux_del(dest, data_ms);
    }

    if (bam_aux_append(dest, "ms", 'i', sizeof(uint32_t), (uint8_t*)&mate_score) == -1) {
        return -1;
    }

    return 0;
}

// Completely delete the CIGAR field
static void clear_cigar(bam1_t *b) {
    memmove(bam_get_cigar(b), bam_get_seq(b),
            b->data + b->l_data - bam_get_seq(b));
    b->l_data -= 4*b->core.n_cigar;
    b->core.n_cigar = 0;
}

// Trim a CIGAR field to end on reference position "end".  Remaining bases
// are turned to soft clips.
static int bam_trim(bam1_t *b, hts_pos_t end) {
    hts_pos_t pos = b->core.pos;
    int n_cigar = b->core.n_cigar, i;
    uint32_t new_cigar_a[1024];
    uint32_t *new_cigar = new_cigar_a;
    uint32_t *cigar = bam_get_cigar(b);

    // Find end of alignment or end of ref
    int op = 0, oplen = 0;
    for (i = 0; i < n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        if (!(bam_cigar_type(op) & 2))
            continue;
        pos += oplen;
        if (pos > end)
            break;
    }

    if (i == n_cigar)
        // looks fine already
        return 0;

    int old_i = i, j = 0;
    // At worst we grow by 1 element (eg 100M -> 70M30S)
    if (n_cigar-i >= 1024-1) {
        new_cigar = malloc(4*(n_cigar-i+1));
        if (!new_cigar)
            return -1;
    }

    // We fill out to new_cigar from here on.
    if (pos-oplen < end) {
        // Partial CIGAR op?  Split existing tag.
        cigar[old_i++] = bam_cigar_gen(end - (pos-oplen), op);
        new_cigar[j++] = bam_cigar_gen(pos-end, BAM_CSOFT_CLIP);
    } else if (pos > end) {
        // entirely off the chromosome; this will trigger CIGAR *, MQUAL 0
        b->core.flag |= BAM_FUNMAP;
        b->core.flag &= ~BAM_FPROPER_PAIR;
    } else {
        // CIGAR op started on the trim junction
        new_cigar[j++] = bam_cigar_gen(oplen, BAM_CSOFT_CLIP);
    }

    // Replace trailing elements.
    for (i++; i < n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CHARD_CLIP) {
            new_cigar[j++] = cigar[i];
        } else {
            new_cigar[j-1] =
                bam_cigar_gen(bam_cigar_oplen(new_cigar[j-1]) + oplen,
                              BAM_CSOFT_CLIP);
        }
    }

    // We now have cigar[0..old_i-1] for existing CIGAR
    // and new_cigar[0..j-1] for new CIGAR trailing component.

    if (old_i+j == n_cigar) {
        // Fits and no data move needed
        memcpy(&cigar[old_i], new_cigar, j*4);
    } else {
        uint8_t *seq_old = bam_get_seq(b);
        uint8_t *aux_end = b->data + b->l_data;
        int nshift;
        if (old_i+j < n_cigar) {
            // Smaller, and can move data down
            nshift = -4*(n_cigar - (old_i+j));
        } else {
            // Bigger, so grow BAM and move data up
            nshift = 4*(old_i+j - n_cigar);
            // FIXME: make htslib's sam_realloc_bam_data public
            if (b->l_data + nshift > b->m_data) {
                uint8_t *new_data = realloc(b->data, b->l_data + nshift);
                if (!new_data) {
                    if (new_cigar != new_cigar_a)
                        free(new_cigar);
                    return -1;
                }
                b->m_data = b->l_data + nshift;
                if (b->data != new_data) {
                    b->data = new_data;
                    seq_old = bam_get_seq(b);
                    aux_end = b->data + b->l_data;
                    cigar = bam_get_cigar(b);
                }
            }
        }
        memmove(seq_old+nshift, seq_old, aux_end - seq_old);
        b->l_data += nshift;
        memcpy(&cigar[old_i], new_cigar, j*4);
        b->core.n_cigar = old_i+j;
    }

    if (new_cigar != new_cigar_a)
        free(new_cigar);

    return 0;
}

// Parses a comma-separated list of "pos", "mqual", "unmap", "cigar", and "aux"
// keywords for the bam sanitizer.
int bam_sanitize_options(const char *str) {
    int opt = 0;

    while (str && *str) {
        const char *str_start;
        while(*str && *str == ',')
            str++;

        for (str_start = str; *str && *str != ','; str++);
        int len = str - str_start;
        if (strncmp(str_start, "all", 3) == 0 || *str_start == '*')
            opt = FIX_ALL;
        else if (strncmp(str_start, "none", 4) == 0 ||
                 strncmp(str_start, "off", 3) == 0)
            opt = 0;
        else if (strncmp(str_start, "on", 2) == 0)
            // default for position sorted data
            opt = FIX_MQUAL | FIX_UNMAP | FIX_CIGAR | FIX_AUX;
        else if (strncmp(str_start, "pos", 3) == 0)
            opt |= FIX_POS;
        else if (strncmp(str_start, "mqual", 5) == 0)
            opt |= FIX_MQUAL;
        else if (strncmp(str_start, "unmap", 5) == 0)
            opt |= FIX_UNMAP;
        else if (strncmp(str_start, "cigar", 5) == 0)
            opt |= FIX_CIGAR;
        else if (strncmp(str_start, "aux", 3) == 0)
            opt |= FIX_AUX;
        else {
            print_error("sanitize", "Unrecognised keyword %.*s\n",
                        len, str_start);
            return -1;
        }
    }

    return opt;
}

int bam_sanitize(sam_hdr_t *h, bam1_t *b, int flags) {
    if ((flags & FIX_POS) && b->core.tid < 0) {
        // RNAME * => pos 0. NB can break alignment chr/pos sort order
        b->core.pos = -1;
        if (flags & FIX_UNMAP)
            b->core.flag |= BAM_FUNMAP;
    }

    if ((flags & FIX_CIGAR) && !(b->core.flag & BAM_FUNMAP)) {
        // Mapped => unmapped correction
        if (b->core.pos < 0 && (flags & FIX_UNMAP)) {
            b->core.flag |= BAM_FUNMAP;
        } else {
            hts_pos_t cur_end, rlen = sam_hdr_tid2len(h, b->core.tid);
            if (b->core.pos >= rlen && (flags & FIX_UNMAP)) {
                b->core.flag |= BAM_FUNMAP;
                if (flags & FIX_POS)
                    b->core.tid = b->core.pos = -1;
            } else if ((cur_end = bam_endpos(b)) > rlen) {
                if (bam_trim(b, rlen) < 0)
                    return -1;
            }
        }
    }

    if (b->core.flag & BAM_FUNMAP) {
        // Unmapped -> cigar/qual correctoins
        if ((flags & FIX_CIGAR) && b->core.n_cigar > 0)
            clear_cigar(b);

        if (flags & FIX_MQUAL)
            b->core.qual = 0;

        // Remove NM, MD, CG, SM tags.
        if (flags & FIX_AUX) {
            uint8_t *from = bam_aux_first(b);
            uint8_t *end = b->data + b->l_data;
            uint8_t *to = from ? from-2 : end;

#define XTAG(a) (((a)[0]<<8) + (a)[1])
            while (from) {
                uint8_t *next = bam_aux_next(b, from);
                if (!next && errno != ENOENT)
                    return -1;

                // Keep tag unless one of a specific set.
                // NB "to" always points to an aux tag start, while
                // "from" is after key.
                from -= 2;
                int key = (int)from[0]<<8 | from[1];
                if (key != XTAG("NM") && key != XTAG("MD") &&
                    key != XTAG("CG") && key != XTAG("SM")) {
                    ptrdiff_t len = (next ? next-2 : end) - from;
                    if (from != to)
                        memmove(to, from, len);
                    to += len;
                }
                from = next;
            }
            b->l_data = to - b->data;
        }
    }

    return 0;
}

// Look for 3 tags in one pass, for efficiencies sake
static inline void find_tags(bam1_t *b,
                             char *t1, uint8_t **t1p,
                             char *t2, uint8_t **t2p,
                             char *t3, uint8_t **t3p) {
    *t1p = *t2p = *t3p = NULL;
    uint8_t *aux = bam_aux_first(b);

    while (aux) {
        if (aux[-2] == t1[0] && toupper(aux[-1]) == t1[1])
            *t1p = aux;
        else if (aux[-2] == t2[0] && toupper(aux[-1]) == t2[1])
            *t2p = aux;
        else if (aux[-2] == t3[0] && toupper(aux[-1]) == t3[1])
            *t3p = aux;
        aux = bam_aux_next(b, aux);
    }
}

// Return 5' and 3' CIGAR hard-clip counts
static inline void hard_clips(bam1_t *b, int *end5, int *end3) {
    uint32_t *cigar = bam_get_cigar(b);
    int ncigar = b->core.n_cigar;
    int endL = 0, endR = 0, nh = 0;

    if (ncigar && bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)
        endL = bam_cigar_oplen(cigar[0]), nh=1;
    if (ncigar > nh && bam_cigar_op(cigar[ncigar-1]) == BAM_CHARD_CLIP)
        endR = bam_cigar_oplen(cigar[ncigar-1]);

    if (b->core.flag & BAM_FREVERSE) {
        *end5 = endR;
        *end3 = endL;
    } else {
        *end5 = endL;
        *end3 = endR;
    }
}

// Get MM, ML and MN tags, and 5' and 3' hard-clip lengths.
// MNi is integer copy of MN, or -1 if absent/invalid
void get_mod_info(bam1_t *b, uint8_t **MM, uint8_t **ML, uint8_t **MN,
                  int *MNi, int *end5, int *end3) {
    find_tags(b, "MM", MM, "ML", ML, "MN", MN);
    if (*MN) {
        *MNi = bam_aux2i(*MN);
        if (errno == EINVAL)
            *MNi = -1;
    } else {
        *MNi = -1;
    }

    if (*MM)
        hard_clips(b, end5, end3);
    else
        *end5 = *end3 = 0; // don't need if MM not found
}

typedef struct MM_state {
    // tags found on "pre" BAM
    uint8_t *MM, *ML, *MN;
} MM_state;

uint8_t *MN_enc(uint8_t *tag, uint32_t n) {
    if (n > UINT16_MAX) {
        tag[0] = 'I';
        i32_to_le(n, tag+1);
        tag += 5;
    } else if (n > UINT8_MAX) {
        tag[0] = 'S';
        i16_to_le(n, tag+1);
        tag += 3;
    } else {
        *tag++ = 'C';
        *tag++ = n;
    }

    return tag;
}

// Trim 5'/3' bases off MM and ML tags, using a previous sequence as a guide.
int trim_MM(bam1_t *pre, bam1_t *cur, int end5, int end3,
            uint8_t *MM, uint8_t *ML, uint8_t *MN) {
    // seqi index / seq_nt16_table[]
    int counts5[16] = {0}, counts3[16] = {0};

    uint8_t *seq = bam_get_seq(pre);
    int i, j;
    for (i = 0; i < end5; i++)
        counts5[bam_seqi(seq, i)]++;
    memcpy(counts3, counts5, 16 * sizeof(*counts3));
    for (; i < pre->core.l_qseq - end3; i++)
        counts3[bam_seqi(seq, i)]++;

    // "p" is position in pre.
    // "q" is position in cur.
    // Hence move up "p" to start and copy from there to "q".
    uint8_t *MMp, *MLp, *MLend, *MMq = NULL, *MLq = NULL;
    if (ML && ML[0] == 'B' && ML[1] == 'C') {
        MLp = ML+6;
        MLend = MLp + le_to_u32(ML+2);
    } else {
        ML = MLp = NULL;
    }
    MMq = MM+1;
    MLq = MLp;
    for (MMp = MM+1; *MMp && (!ML || *MLp); ) {
        int fundamental = seq_nt16_table[*MMp];
        while (*MMp && *MMp != ',')
            *MMq++ = *MMp++;
        if (*MMp)
            *MMq++ = *MMp++;

        // Now on comma separated list for MM and BC array for ML. Skip
        int n = 0;
        while (n < counts5[fundamental]) {
            char *endptr;
            long delta = strtol(MMp, &endptr, 10);
            if (counts5[fundamental] - n > delta) {
                // Skip entire delta in MM and ML.
                // Eg counts[]=10, MM=3,10 ML=<10><20> => MM=10 ML=<20>
                n += delta+1;
                MLp++;
            } else {
                // Shrink delta, writing MM and ML is unchanged.
                // Eg counts[]=3, MM=10,4 ML=<10><20> => MM=7,4 ML=<10><20>
                //int l = sprintf(MMq, "%ld,", delta - counts5[fundamental]);
                int l = sprintf(MMq, "%ld,", delta - (counts5[fundamental]-n));
                //n = counts5[fundamental];
                n += delta+1;
                MMq += l;
                *MLq++ = *MLp++;
            }

            MMp = endptr;
            if (*MMp != ',')
                // error?  if not ; also?
                break;
            MMp++;
        }

        // Copy
        while (n < counts3[fundamental]) {
            char *endptr;
            long delta = strtol(MMp, &endptr, 10);
            if (counts3[fundamental] - n > delta) {
                // Copy entire delta in MM and ML including [,;]
                memmove(MMq, MMp, (uint8_t *)endptr - MMp + 1);
                MMq += (uint8_t *)endptr - MMp + 1;
                n += delta+1;
                if (ML)
                    *MLq++ = *MLp++;
            } else {
                // Next mod is into 3' cutoff, so can terminate MM/ML now
                n = counts3[fundamental];
                MLp++;
            }

            MMp = endptr;
            if (*MMp != ',')
                break;
            MMp++;
        }

        // Skip
        while (*MMp && *MMp != ';') {
            while (*MMp && *MMp != ',' && *MMp != ';')
                MMp++;
            if (*MMp == ',')
                MMp++;

            if (ML)
                MLp++;
        }
        MMq[-1] = ';'; // replaces , with ; if clipping right
        MMp++;
    }

    MMp++; // skip nul
    *MMq++ = 0;

    // Adjust ML B array length
    u32_to_le(MLq-(ML+6), ML+2);

    // Move MM and ML down to include their MM:Z and ML:B bits
    if (MM) MM-=2;
    if (ML) ML-=2;

    // Now MM/ML are start of tags, MMq/MLq are ends of edited tags,
    // and MMp/MLp are ends of original tags.  Walk through tags taking up
    // any gaps
    //
    // Eg XXXXXXmmmmm--YYYlllll-ZZ (m and l are edited MM and ML tags)
    // => XXXXXXmmmmmYYYlllllZZ

    uint8_t *tag = bam_get_aux(cur), *tag_end = cur->data + cur->l_data;
    uint8_t *to = tag;
    //int mn_done = 0;
    while (tag < tag_end) {
        if (tag[0] == 'M' && (tag[1] == 'M' || tag[1] == 'm')) {
            // Slow but easy
            
            memmove(to, MM, MMq-MM); // length of new tag
            to += MMq-MM;
            tag = MMp; // size of old tag
        } else if (tag[0] == 'M' && (tag[1] == 'L' || tag[1] == 'l')) {
            memmove(to, ML, MLq-ML);
            to += MLq-ML;
            tag = MLp;
//        } else if (tag[0] == 'M' && tag[1] == 'N') {
//            tag = MN_enc(tag, cur->core.l_qseq);
//            mn_done = 1;
        } else {
            // Want aux_skip, but it's private.
            // So we use bam_aux_next with work-arounds. :(
            uint8_t *from = tag;
            tag = bam_aux_next(cur, tag+2);
            tag = tag ? tag-2 : tag_end;
            memmove(to, from, tag-from);
        }
    }
    cur->l_data = to - cur->data;

//    // Set MN tag
//    if (!mn_done) {
//        // inefficient, but minimal compared to everything else
//        if (bam_aux_update_int(cur, "MN", cur->core.l_qseq) < 0)
//            return -1;
//    }

    return 0;
}

// Fix base modification tags MM, ML and MN.
// For supplementary-style alignments we may have hard-clipped the sequence
// and just duplicated the MM/ML tags.  Use the primary alignment to get the
// clipped sequence so we can trim MM/ML accordingly.
int fix_MM(bam1_t *pre, bam1_t *cur, int name_match, MM_state *state) {
    int end5, end3;
    int MNi = 0; // MN of -1 is used as indicator for no valid mods

    if (!name_match) {
        // First time we've see this name.
        // Look for base modification tags and sanity check.
        fprintf(stderr, "New read %s, len %d\n", bam_get_qname(cur), cur->core.l_qseq);
        get_mod_info(cur, &state->MM, &state->ML, &state->MN, &MNi, &end5, &end3);
        if (!state->MM)
            return 0;

        if (MNi > 0 && cur->core.l_qseq - MNi == end5) {
            // MN and SEQ lengths disagree by 5' hard-clip size.
            // Let's assume that's the cause.
            fprintf(stderr, "TODO: Trim MM, ML, MN by %d\n", end5);
        } else if (end5 == 0 && MNi <= 0) {
            // No MN tag, but no clipping. Assume MM is valid
            MNi = cur->core.l_qseq;
            fprintf(stderr, "TODO: create MN tag\n");
        } else {
            // No MN tag and hard clips present - fail
            fprintf(stderr, "TODO: drop MN / ML tags\n");
            MNi = -1;
        }
    } else {
        uint8_t *cur_MM = NULL, *cur_ML = NULL, *cur_MN = NULL;
        MNi = -1;
        get_mod_info(cur, &cur_MM, &cur_ML, &cur_MN, &MNi, &end5, &end3);

        if (!cur_MM)
            return 0;

        fprintf(stderr, "MNi=%d seq=%d L %d R %d\n",
                MNi, cur->core.l_qseq, end5, end3);

        // Does MN match seq length?  If so, we believe it's already valid
        if (MNi == cur->core.l_qseq)
            return 0;

        // Length mismatch and/or no known length, so check vs full seq.
        if (pre->core.l_qseq != cur->core.l_qseq + end3 + end5) {
            fprintf(stderr, "    Mismatch: CLEAR MM, ML, MN\n");
            uint8_t *tag = bam_aux_first(cur);
            while (tag) {
                if (tag[-2] == 'M' &&
                    (tag[-1] == 'M' || tag[-1] == 'm' ||
                     tag[-1] == 'L' || tag[-1] == 'l' ||
                     tag[-1] == 'N')) {
                    // NB: not as efficient as it could be as we move tags
                    // multiple times, but simple.
                    fprintf(stderr, "        Del %.2s\n", tag-2);
                    tag = bam_aux_remove(cur, tag);
                } else {
                    tag = bam_aux_next(cur, tag);
                }
            }
            return 0;
        } else if (end5 != 0) {
            fprintf(stderr, "    Match: FIX MM\n");
            trim_MM(pre, cur, end5, end3, cur_MM, cur_ML, cur_MN);
        } // else no 5' hard clips so MM is already valid

        // Set MN so we've validated it.
        fprintf(stderr, "    Add MN\n");
        // inefficient, but minimal compared to everything else
        if (bam_aux_update_int(cur, "MN", cur->core.l_qseq) < 0)
            return -1;
    }

    return 0;
}

// currently, this function ONLY works if each read has one hit
static int bam_mating_core(samFile *in, samFile *out, int remove_reads,
                           int proper_pair_check, int add_ct,
                           int do_mate_scoring, char *arg_list, int no_pg,
                           int sanitize_flags, int base_mods)
{
    sam_hdr_t *header;
    bam1_t *b[2] = { NULL, NULL };
    int curr, has_prev, result;
    hts_pos_t pre_end = 0, cur_end = 0;
    kstring_t str = KS_INITIALIZE;
    struct MM_state mm_state;

    header = sam_hdr_read(in);
    if (header == NULL) {
        fprintf(stderr, "[bam_mating_core] ERROR: Couldn't read header\n");
        return 1;
    }

    // Accept unknown, unsorted, or queryname sort order, but error on coordinate sorted.
    if (!sam_hdr_find_tag_hd(header, "SO", &str) && str.s && !strcmp(str.s, "coordinate")) {
        fprintf(stderr, "[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname.\n");
        goto fail;
    }
    ks_free(&str);

    if (!no_pg && sam_hdr_add_pg(header, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto fail;

    if (sam_hdr_write(out, header) < 0) goto write_fail;

    b[0] = bam_init1();
    b[1] = bam_init1();
    curr = 0; has_prev = 0;
    while ((result = sam_read1(in, header, b[curr])) >= 0) {
        bam1_t *cur = b[curr], *pre = b[1-curr];
        if (bam_sanitize(header, cur, sanitize_flags) < 0)
            goto fail;

        int name_match = has_prev
            ? strcmp(bam_get_qname(cur), bam_get_qname(pre)) == 0
            : 0;

        if (name_match)
            fprintf(stderr, "%s %5d %s %5d\n",
                    bam_get_qname(cur), cur->core.flag,
                    bam_get_qname(pre), pre->core.flag);

        if (base_mods)
            fix_MM(pre, cur, name_match, &mm_state);

        if (cur->core.flag & BAM_FSECONDARY)
        {
            if ( !remove_reads ) {
                if (sam_write1(out, header, cur) < 0) goto write_fail;
            }
            continue; // skip secondary alignments
        }
        if (cur->core.flag & BAM_FSUPPLEMENTARY)
        {
            if (sam_write1(out, header, cur) < 0) goto write_fail;
            continue; // pass supplementary alignments through unchanged (TODO:make them match read they came from)
        }
        if ((cur->core.flag&BAM_FUNMAP) == 0) // If mapped calculate end
        {
            cur_end = bam_endpos(cur);
        }

        if (has_prev) { // do we have a pair of reads to examine?
            if (name_match) { // identical pair name
                pre->core.flag |= BAM_FPAIRED;
                cur->core.flag |= BAM_FPAIRED;
                if (sync_mate(pre, cur)) goto fail;

                if (pre->core.tid == cur->core.tid && !(cur->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))
                    && !(pre->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))) // if safe set TLEN/ISIZE
                {
                    hts_pos_t cur5, pre5;
                    cur5 = (cur->core.flag&BAM_FREVERSE)? cur_end : cur->core.pos;
                    pre5 = (pre->core.flag&BAM_FREVERSE)? pre_end : pre->core.pos;
                    cur->core.isize = pre5 - cur5; pre->core.isize = cur5 - pre5;
                } else cur->core.isize = pre->core.isize = 0;
                if (add_ct) bam_template_cigar(pre, cur, &str);
                // TODO: Add code to properly check if read is in a proper pair based on ISIZE distribution
                if (proper_pair_check && !plausibly_properly_paired(pre,cur)) {
                    pre->core.flag &= ~BAM_FPROPER_PAIR;
                    cur->core.flag &= ~BAM_FPROPER_PAIR;
                }

                if (do_mate_scoring) {
                    if ((add_mate_score(pre, cur) == -1) || (add_mate_score(cur, pre) == -1)) {
                        fprintf(stderr, "[bam_mating_core] ERROR: unable to add mate score.\n");
                        goto fail;
                    }
                }

                // Write out result
                if ( !remove_reads ) {
                    if (sam_write1(out, header, pre) < 0) goto write_fail;
                    if (sam_write1(out, header, cur) < 0) goto write_fail;
                } else {
                    // If we have to remove reads make sure we do it in a way that doesn't create orphans with bad flags
                    if(pre->core.flag&BAM_FUNMAP) cur->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(cur->core.flag&BAM_FUNMAP) pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(!(pre->core.flag&BAM_FUNMAP)) {
                        if (sam_write1(out, header, pre) < 0) goto write_fail;
                    }
                    if(!(cur->core.flag&BAM_FUNMAP)) {
                        if (sam_write1(out, header, cur) < 0) goto write_fail;
                    }
                }
                has_prev = 0;
            } else { // unpaired?  clear bad info and write it out
                pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
                pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                if ( !remove_reads || !(pre->core.flag&BAM_FUNMAP) ) {
                    if (sam_write1(out, header, pre) < 0) goto write_fail;
                }
            }
        } else has_prev = 1;
        curr = 1 - curr;
        pre_end = cur_end;
    }
    if (result < -1) goto read_fail;
    if (has_prev && !remove_reads) { // If we still have a BAM in the buffer it must be unpaired
        bam1_t *pre = b[1-curr];
        if (pre->core.tid < 0 || pre->core.pos < 0 || pre->core.flag&BAM_FUNMAP) { // If unmapped
            pre->core.flag |= BAM_FUNMAP;
            pre->core.tid = -1;
            pre->core.pos = -1;
        }
        pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
        pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);

        if (sam_write1(out, header, pre) < 0) goto write_fail;
    }
    sam_hdr_destroy(header);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    ks_free(&str);
    return 0;

 read_fail:
    print_error("fixmate", "Couldn't read from input file");
    goto fail;

 write_fail:
    print_error_errno("fixmate", "Couldn't write to output file");
 fail:
    sam_hdr_destroy(header);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    ks_free(&str);
    return 1;
}

void usage(FILE* where)
{
    fprintf(where,
"Usage: samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n"
"Options:\n"
"  -r           Remove unmapped reads and secondary alignments\n"
"  -p           Disable FR proper pair check\n"
"  -c           Add template cigar ct tag\n"
"  -m           Add mate score tag\n"
"  -u           Uncompressed output\n"
"  -z, --sanitize FLAG[,FLAG]\n"
"               Sanitize alignment fields [defaults to all types]\n"
"  --no-PG      do not add a PG line\n");

    sam_global_opt_help(where, "-.O..@-.");

    fprintf(where,
"\n"
"As elsewhere in samtools, use '-' as the filename for stdin/stdout. The input\n"
"file must be grouped by read name (e.g. sorted by name). Coordinated sorted\n"
"input is not accepted.\n");
}

int bam_mating(int argc, char *argv[])
{
    htsThreadPool p = {NULL, 0};
    samFile *in = NULL, *out = NULL;
    int c, remove_reads = 0, proper_pair_check = 1, add_ct = 0, res = 1,
        mate_score = 0, no_pg = 0, sanitize_flags = FIX_ALL, base_mods = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    char wmode[4] = {'w', 'b', 0, 0};
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };
    char *arg_list = NULL;

    // parse args
    if (argc == 1) { usage(stdout); return 0; }
    while ((c = getopt_long(argc, argv, "rpcmMO:@:uz:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'r': remove_reads = 1; break;
        case 'p': proper_pair_check = 0; break;
        case 'c': add_ct = 1; break;
        case 'm': mate_score = 1; break;
        case 'M': base_mods = 1; break;
        case 'u': wmode[2] = '0'; break;
        case 1: no_pg = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': usage(stderr); goto fail;
        case 'z':
            if ((sanitize_flags = bam_sanitize_options(optarg)) < 0)
                exit(1);
            break;
        }
    }
    if (optind+1 >= argc) { usage(stderr); goto fail; }

    if (!no_pg && !(arg_list =  stringify_argv(argc+1, argv-1)))
        goto fail;

    // init
    if ((in = sam_open_format(argv[optind], "rb", &ga.in)) == NULL) {
        print_error_errno("fixmate", "cannot open input file");
        goto fail;
    }
    sam_open_mode(wmode+1, argv[optind+1], NULL);
    if ((out = sam_open_format(argv[optind+1], wmode, &ga.out)) == NULL) {
        print_error_errno("fixmate", "cannot open output file");
        goto fail;
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            goto fail;
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }

    // run
    res = bam_mating_core(in, out, remove_reads, proper_pair_check, add_ct,
                          mate_score, arg_list, no_pg, sanitize_flags,
                          base_mods);

    // cleanup
    sam_close(in);
    if (sam_close(out) < 0) {
        fprintf(stderr, "[bam_mating] error while closing output file\n");
        res = 1;
    }

    if (p.pool) hts_tpool_destroy(p.pool);
    free(arg_list);
    sam_global_args_free(&ga);
    return res;

 fail:
    if (in) sam_close(in);
    if (out) sam_close(out);
    if (p.pool) hts_tpool_destroy(p.pool);
    free(arg_list);
    sam_global_args_free(&ga);
    return 1;
}


