#ifndef _H_TEST_
#define _H_TEST_


int jacobian_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int la_test(void);
int newton_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int residual_test(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int sw2_wd_test(int npx, int npy, int nt);
int timeloop_test(int npx, int npy, int nt);
int run_tests(void);

#endif
