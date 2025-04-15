#include "adh.h"


double get_node_uda(SMODEL_SUPER *sm, int inode) {
    return sm->ivars[sm->ivar_pos.var[adh_def._UDA]][inode];
}
void put_node_uda(SMODEL_SUPER *sm, int inode, double value) {
    sm->ivars[sm->ivar_pos.var[adh_def._UDA]][inode] = value;;
}


double get_node_vda(SMODEL_SUPER *sm, int inode) {
    return sm->ivars[sm->ivar_pos.var[adh_def._VDA]][inode];
}
void put_node_vda(SMODEL_SUPER *sm, int inode, double value) {
    sm->ivars[sm->ivar_pos.var[adh_def._VDA]][inode] = value;;
}

double get_node_h(SMODEL_SUPER *sm, int inode) {
    return sm->ivars[sm->ivar_pos.var[adh_def._H]][inode];
}
void put_node_h(SMODEL_SUPER *sm, int inode, double value) {
    sm->ivars[sm->ivar_pos.var[adh_def._H]][inode] = value;;
}
