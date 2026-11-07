// Nussinov CPU baseline: single thread, -O3, inner loop vectorizable.
// Keeps a transposed copy of the table so the column read in the split term is
// contiguous, which lets the compiler emit SIMD add+max.  This is a fair baseline,
// not a strawman: a naive column-strided loop is several times slower.
//
// usage: cpu_baseline <seqfile> <L> [reps]     (seqfile: one line of ACGU)
// prints: n, N(0,n-1), seconds per fold, Gops/s (ops = n(n^2-1)/6 add-max steps)
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int pairs(char a, char b) {
    return (a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'G' && b == 'C') ||
           (a == 'C' && b == 'G') || (a == 'G' && b == 'U') || (a == 'U' && b == 'G');
}

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static int32_t fold(const char *s, int n, int L, int32_t *N, int32_t *T) {
    // N[i*n + j] = N(i,j); T[j*n + i] = N(i,j)
    memset(N, 0, sizeof(int32_t) * n * n);
    memset(T, 0, sizeof(int32_t) * n * n);
    for (int span = 1; span < n; span++) {
        for (int i = 0; i + span < n; i++) {
            int j = i + span;
            int32_t inner = (i + 1 <= j - 1) ? N[(i + 1) * n + (j - 1)] : 0;
            int32_t best = inner + ((j - i > L) && pairs(s[i], s[j]));
            const int32_t *row = &N[i * n];          // N(i, k),   k = i .. j-1
            const int32_t *col = &T[j * n + 1];      // N(k+1, j), k = i .. j-1
            for (int k = i; k < j; k++) {
                int32_t v = row[k] + col[k];
                best = v > best ? v : best;
            }
            N[i * n + j] = best;
            T[j * n + i] = best;
        }
    }
    return N[n - 1];
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <seqfile> <L> [reps]\n", argv[0]);
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    static char s[1 << 16];
    if (!fgets(s, sizeof s, f)) { fprintf(stderr, "empty file\n"); return 1; }
    fclose(f);
    int n = (int)strcspn(s, "\r\n");
    int L = atoi(argv[2]);
    int reps = argc > 3 ? atoi(argv[3]) : 3;
    int32_t *N = malloc(sizeof(int32_t) * n * n), *T = malloc(sizeof(int32_t) * n * n);
    int32_t r = fold(s, n, L, N, T);                 // warm-up
    double best = 1e30;
    for (int k = 0; k < reps; k++) {
        double t0 = now();
        r = fold(s, n, L, N, T);
        double dt = now() - t0;
        if (dt < best) best = dt;
    }
    double ops = (double)n * ((double)n * n - 1) / 6.0;
    printf("%d %d %.6f %.3f\n", n, r, best, ops / best / 1e9);
    free(N); free(T);
    return 0;
}
