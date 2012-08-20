#ifndef ERRMOD_H
#define ERRMOD_H

#include <stdint.h>

struct __errmod_coef_t;

typedef struct {
	double depcorr;
	struct __errmod_coef_t *coef;
} errmod_t;

/* parameters for the model used to fit heterozygotes */
typedef struct {
    int model_sel;
    union {
        double p;
        struct {
            double alpha;
            double beta;
        } ab;
    } param;
} call_model_t;

#define MODEL_SEL_BINOM 0
#define MODEL_SEL_BETABINOM 1

errmod_t *errmod_init(const double depcorr, const call_model_t* model);
void errmod_destroy(errmod_t *em);

/*
	n: number of bases
	m: maximum base
	bases[i]: qual:6, strand:1, base:4
	q[i*m+j]: phred-scaled likelihood of (i,j)
 */
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q);

#endif
