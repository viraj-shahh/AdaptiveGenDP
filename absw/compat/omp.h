/* stub: build ksw-test without OpenMP (pragmas ignored, single thread) */
#pragma once
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_get_num_threads(void) { return 1; }
static inline void omp_set_num_threads(int n) { (void)n; }
