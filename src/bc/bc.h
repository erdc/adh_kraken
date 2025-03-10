#ifndef H_BC_
#define H_BC_

// int get_id(SIO *io, char **token, int nmax, char *error_descript);
// void read_bc(SMODEL_SUPER *mod);
// void read_bc_AWRITE(SMODEL_SUPER *mod, char **token);
// void read_bc_CN(SMODEL_SUPER *mod, char **token);
// void read_bc_DB(SMODEL_SUPER *mod, char **token);
// void read_bc_DEBUG(SMODEL_SUPER *mod, char **token);
// void read_bc_FILE_OUTPUT(SMODEL_SUPER *mod, char **token);
// void read_bc_FR(SMODEL_SUPER *mod, char **token);
// void read_bc_IP(SMODEL_SUPER *mod, char **token);
// void read_bc_MP(SMODEL_SUPER *mod, char **token);
// void read_bc_NB(SMODEL_SUPER *mod, char **token);
// void read_bc_NOTERM(SMODEL_SUPER *mod, char **token);
// void read_bc_OFF(SMODEL_SUPER *mod, char **token);
// void read_bc_OP(SMODEL_SUPER *mod, char **token);
// void read_bc_PC(SMODEL_SUPER *mod, char **token);
// void read_bc_SCREEN_OUTPUT(SMODEL_SUPER *mod, char **token);
// void read_bc_SERIES(SSERIES *series_out, char **token);
// void read_bc_STRINGS(SMODEL_SUPER *mod, char **token);
// void read_bc_TC(SMODEL_SUPER *mod, char **token);
// void read_bc_TEST(SMODEL_SUPER *mod, char **token);
// void read_bc_prep(SMODEL_SUPER *mod);
void read_bc_MODEL(SMODEL_SUPER *sm, FILE *fp, char codes[][10]);


#endif
