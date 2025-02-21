// An AdH SuperModel
#ifndef H_SDT_
#define H_SDT_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    double dt, dt_old, dt_err, dt_prev;
    double t_init, t_prev, t_final;
    double tau_temporal;
    double time;
    int nt; // current time-step
    int t_adpt_flag;
    int ientry_out;
    bool green_ampt;
} SDT;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void sdt_init(SDT *t);
void sdt_check_units(int uin, int uout);
int sdt_get_conversion_factor(int units,int direct);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
