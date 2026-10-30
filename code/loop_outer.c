// The flattened-array loop from the case study, Example 2.
// Sequential:   gcc -O3 -fopenmp -DSERIAL loop_outer.c -o serial
// GCC automatic: gcc -O3 -fopenmp -DSERIAL -ftree-parallelize-loops=4 loop_outer.c -o autopar
// Outer loop:   gcc -O3 -fopenmp loop_outer.c -o outer
// GCC's report: gcc -O3 -fopenmp -DSERIAL -ftree-parallelize-loops=4 -c loop_outer.c \
//               -fdump-tree-parloops-details -o /dev/null
//               grep FAILED loop_outer.c.*parloops*
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define ROWS 4000
#define REPEAT 300
int main(void) {
    double *A = malloc(sizeof(double) * 1000 * ROWS);
    for (long i = 0; i < 1000L * ROWS; i++) A[i] = i % 7;
    double start = omp_get_wtime();
    for (int r = 0; r < REPEAT; r++) {
#ifndef SERIAL
#pragma omp parallel for num_threads(4)
#endif
        for (int i = 0; i < ROWS; i++)
            for (int j = 0; j < 500; j++)
                A[1000 * i + j] = A[1000 * i + j + 500] * 0.5 + 1.0;
    }
    double time = omp_get_wtime() - start, sum = 0;
    for (long i = 0; i < 1000L * ROWS; i++) sum += A[i];
    printf("%.3f s  checksum %.6f\n", time, sum);
}
