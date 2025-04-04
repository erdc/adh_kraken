#ifndef _H_TEST_
#define _H_TEST_


int test_jacobian(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int test_la(void);
int test_newton(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int test_residual(int npx, int npy, double xmin, double xmax, double ymin, double ymax);
int test_sw2_nb(int nt);
int test_sw2_wd(int npx, int npy, int nt);
int test_timeloop(int npx, int npy, int nt);
int test_comm_update(int comm_type);
int run_tests(void);

#endif
