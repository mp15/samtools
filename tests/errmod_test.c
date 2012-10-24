#include <stdio.h>
#include "errmod.h"

#define MAKEBASE(base, strand, qual) base|(strand<<4)|(qual<<5)

void printq(float *q, const int m);
void printlhet(errmod_t* error_model);
/* table of constants generated for given depcorr and eta */
typedef struct __errmod_coef_t {
	double *fk, *beta, *lhet;
} errmod_coef_t;

int main(int argc, char** argv)
{
    /* error model */
    errmod_t* error_model_a;
    errmod_t* error_model_b;
    errmod_t* error_model_c;
    /* observed bases */
    uint16_t bases[] = { MAKEBASE(0,0,20), MAKEBASE(0,0,20), MAKEBASE(0,0,20), MAKEBASE(0,0,20),
        MAKEBASE(0,0,20), MAKEBASE(0,0,20), MAKEBASE(1,0,20), MAKEBASE(1,0,20) };
    float q_a[25], q_b[25], q_c[25];
    /* number of observations */
    const int n = 8;
    const int m = 2;

    /* setup with depcorr of 1-0.83 taken from BAM2BCF */
    call_model_t model;
    model.model_sel = MODEL_SEL_BINOM;
    model.param.p = 0.5;
    error_model_a = errmod_init(1.0 - 0.83, &model);
    model.param.p = 0.25;
    error_model_b = errmod_init(1.0 - 0.83, &model);
    model.model_sel = MODEL_SEL_BETABINOM;
    model.param.ab.alpha = 0.3;
    model.param.ab.beta = 0.7;
    error_model_c = errmod_init(1.0 - 0.83, &model);

    /* run test */
    errmod_cal(error_model_a, n, m, bases, q_a);
    errmod_cal(error_model_b, n, m, bases, q_b);
    errmod_cal(error_model_c, n, m, bases, q_c);

    /* write out results */
    printq(q_a,m);
    printq(q_b,m);
    printq(q_c,m);
return 0;
#if 0
    uint16_t bases_t[] = { MAKEBASE(1,0,63), MAKEBASE(1,0,63), MAKEBASE(1,0,63), MAKEBASE(1,0,63) };

    /* run test */
    errmod_cal(error_model, n, m, bases_t, q);

    /* write out results */
    printq(q);

    uint16_t bases_g[] = { MAKEBASE(2,0,63), MAKEBASE(2,0,63), MAKEBASE(2,0,63), MAKEBASE(2,0,63) };
    /* run test */
    errmod_cal(error_model, n, m, bases_g, q);

    /* write out results */
    printq(q);

    uint16_t bases_c[] = { MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63) };
    /* run test */
    errmod_cal(error_model, n, m, bases_c, q);

    /* write out results */
    printq(q);

    /* observed bases */
    printf("50-50 het test\r\n");
    uint16_t bases_het[] = { MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63),
        MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63),
        MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63),
        MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63) };
    /* run test */
    errmod_cal(error_model, n*4, m, bases_het, q);

    /* write out results */
    printq(q);

    /* observed bases */
    printf("slightly off 50-50 het test\r\n");
    uint16_t bases_hetoff[] = { MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63),
        MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(0,0,63), MAKEBASE(3,0,63),
        MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63),
        MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63), MAKEBASE(3,0,63) };
    /* run test */
    errmod_cal(error_model, n*4, m, bases_hetoff, q);

    /* write out results */
    printq(q);
    
    printlhet(error_model);

    /* clean up */
    errmod_destroy(error_model);

    return 0;
#endif
}

void printq(float *q, const int m)
{
    /* write out results */
    printf("q is:\r\n");
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%7.3f ", q[i*m+j]);
        }
        printf("\r\n");
    }
}

void printlhet(errmod_t* error_model)
{
    int i, j;
    for ( i = 0; i < 256; i++)
    {
        for ( j = 0; j < 256; j++)
        {
            printf("%d\t%d\t%03.3f\r\n", i, j, error_model->coef->lhet[i|j<<8]);
        }
    }

}
