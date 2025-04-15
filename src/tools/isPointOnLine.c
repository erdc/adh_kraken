#include "adh.h"
bool isPointOnLine(SVECT p1, SVECT p2, SVECT pointCheck) {
    double d1 = sqrt(pow(p1.x - pointCheck.x,2) + pow(p1.y - pointCheck.y,2));
    double d2 = sqrt(pow(p2.x - pointCheck.x,2) + pow(p2.y - pointCheck.y,2));
    double d3 = sqrt(pow(p1.x - p2.x,2) + pow(p1.y - p2.y,2));
    if (fabs((d1 + d2) - d3) < SMALL6) {return TRUE;}
    return FALSE;
}
