
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "array.h"
#include "random.h"

static unsigned int seeds[] = {1, 0xC75E464, 0xC0EDA55,
                               -1 /* setinel */};

static int sizes[] = {100, 2000,
                      -1 /* setinel */};

static double sigmas[] = {1.0, 1e-4,
                          -1 /*setinel*/};

static inline double
ftime_to_double(struct timeval t)
{
        return (double)t.tv_sec + t.tv_usec*1e-6;
}

static inline void
rstat(array_t *arr, double *meanp, double *sigmap)
{
        REQ_VARR(arr);
        double sum_x=0.0, sum_x2=0.0;
        double * CEX_RESTRICT data = ARR_DATA_AS(double, arr);
        unsigned int N=3*ARR_LENGTH(arr);
        for (int i=0; i<N; i++) {
                sum_x += data[i];
                sum_x2 += data[i] * data[i];
        }
        *meanp = sum_x / (double)N;
        *sigmap = sqrt(sum_x2 / (double)(N) - *meanp * *meanp);
}

int
main(int argc, char **argv)
{
        if (sizeof(vec_t)!=3*sizeof(double)) {
                Fatal("assume no padding of vector stuct");
        }
        for (unsigned int *seedp=seeds; *seedp!=-1; seedp++) {
                for (int *sizep=sizes; *sizep!=-1; sizep++) {
                        for (double *sigmap=sigmas; *sigmap>0; sigmap++) {
                                CEX_seed_random(*seedp);
                                array_t *arr = CEX_make_vec_array(0);
                                CEX_align_array(arr, sizeof(double));
                                CEX_prealloc_array(arr, *sizep);
                                ARR_LENGTH(arr) = *sizep;
                                struct timeval start, end;
                                gettimeofday(&start, NULL);
                                CEX_generate_gauss(arr, *sigmap);
                                gettimeofday(&end, NULL);
                                double rmean, rsigma;
                                rstat(arr, &rmean, &rsigma);
                                printf("seed=%08X size=%d sigma=%.4e => tm = %.3f ms mean=%.4f sigma=%.4e\n",
                                       *seedp, *sizep, *sigmap,
                                       1e3*(ftime_to_double(end) - ftime_to_double(start)),
                                       rmean, rsigma);
                                fflush(stdout);
                                CEX_free_array(arr);

                        }
                }
        }
        return 0;
}
